"""
能带对齐（Band Offset）分析 — 平均势法 (Van de Walle & Martin)

输出:
  output/band_alignment.csv       能带对齐参数表 (Origin/Excel 用)
  output/locpot_planar_*.csv      静电势曲线数据
  output/band_alignment.json      完整数值存档
  output/band_alignment.png       能带对齐图
  output/locpot_bulk_regions.png  体相区标识图
"""

import argparse
import gzip
import json
import logging
import os
import shutil
import tempfile
import warnings
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from monty.json import MontyDecoder
from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Locpot

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 11, "axes.linewidth": 1.2, "figure.dpi": 100,
    "mathtext.default": "regular",
})

DATA_DIR = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S/data")
OUTPUT_DIR = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S/postprocessing/output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 异质结 -> (LDH slab, CMCH slab) — Level 2 应变 slab（共同晶格）
HETERO_MAP = {
    "hetero_intrinsic": ("LDH_strained", "CMCH_strained"),
    "hetero_s_doped": ("LDH_S_strained", "CMCH_strained"),
    "hetero_s_exposed": ("LDH_S_flip_strained", "CMCH_strained"),
}

SYSTEM_LABELS = {
    "CMCH_strained": "CMCH (Strained)",
    "LDH_strained": "CoNiOH2 (Strained)",
    "LDH_S_strained": "CoNiOH2S (Strained)",
    "LDH_S_flip_strained": "CoNiOH2S (Strained-flip)",
    "hetero_intrinsic": "Intrinsic Het.",
    "hetero_s_doped": "S-Doped Het.",
    "hetero_s_exposed": "S-Exposed Het.",
}


# ═══════════════════════════════════════════════
# Utility: read files
# ═══════════════════════════════════════════════

def read_structure(name: str) -> Structure | None:
    static = DATA_DIR / name / "3-static"
    contcar = static / "CONTCAR.gz"
    if not contcar.exists():
        return None
    try:
        with gzip.open(contcar, "rt") as f:
            return Structure.from_str(f.read(), fmt="poscar")
    except Exception:
        return None


def read_locpot_planar(name: str) -> tuple[np.ndarray, int, float | None]:
    """返回 (planar_avg, ngz, efermi)."""
    static = DATA_DIR / name / "3-static"
    locpot = static / "LOCPOT_dipole.gz"
    if not locpot.exists():
        locpot = static / "LOCPOT.gz"
    if not locpot.exists():
        return None, 0, None

    # 解压
    try:
        if str(locpot).endswith(".gz"):
            tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".LOCPOT")
            with gzip.open(locpot, "rb") as f_in:
                shutil.copyfileobj(f_in, tmp)
            tmp_path = tmp.name
            tmp.close()
        else:
            tmp_path = str(locpot)

        lp = Locpot.from_file(tmp_path)
        avg = lp.get_average_along_axis(2)
        ngx, ngy, ngz = lp.dim
        if str(locpot).endswith(".gz"):
            os.unlink(tmp_path)

        # E_Fermi: 从 OUTCAR_dipole 或 OUTCAR
        outcar = static / "OUTCAR_dipole.gz"
        if not outcar.exists():
            outcar = static / "OUTCAR.gz"
        if outcar.exists():
            from pymatgen.io.vasp.outputs import Outcar
            o = Outcar(outcar)
            efermi = o.efermi
        else:
            efermi = None

        return avg, ngz, efermi
    except Exception as e:
        return None, 0, None


def read_vbm_cbm(name: str) -> tuple[float | None, float | None]:
    """从 band_out.json 读取 VBM/CBM."""
    band_path = DATA_DIR / name / "band_out.json"
    if not band_path.exists():
        log.warning("  No band_out.json for %s", name)
        return None, None

    try:
        with open(band_path) as f:
            raw = json.load(f)
        bd_dict = raw.get("vasp_objects", {}).get("bandstructure")
        if bd_dict is None:
            log.warning("  No bandstructure for %s", name)
            return None, None

        bd = MontyDecoder().process_decoded(bd_dict)
        vbm = bd.get_vbm()
        cbm = bd.get_cbm()

        v_e = vbm.get("energy") if vbm else None
        c_e = cbm.get("energy") if cbm else None
        log.info("  %s: VBM=%.4f, CBM=%s", name, v_e, f"{c_e:.4f}" if c_e else "N/A")
        return v_e, c_e
    except Exception as e:
        log.warning("  Failed to read band for %s: %s", name, e)
        return None, None


# ═══════════════════════════════════════════════
# Core: find bulk-like regions
# ═══════════════════════════════════════════════

def find_slab_bulk_region(avg: np.ndarray, ngz: int, struct: Structure) -> tuple[int, int]:
    """slab 体相区域: z 方向中间 40%."""
    c = struct.lattice.c
    grid = np.linspace(0, c, ngz)
    z_all = [s.coords[2] for s in struct]
    z_min, z_max = min(z_all), max(z_all)
    margin = (z_max - z_min) * 0.3  # 两端各去掉 30%
    start_z = z_min + margin
    end_z = z_max - margin
    start_idx = max(0, int(np.searchsorted(grid, start_z)))
    end_idx = min(ngz, int(np.searchsorted(grid, end_z)))
    # 确保至少 5 个格点
    if end_idx - start_idx < 5:
        mid = (z_min + z_max) / 2
        half = (z_max - z_min) * 0.15
        start_idx = max(0, int(np.searchsorted(grid, mid - half)))
        end_idx = min(ngz, int(np.searchsorted(grid, mid + half)))
    return start_idx, end_idx


def find_hetero_bulk_regions(avg: np.ndarray, ngz: int, struct: Structure) -> dict:
    """在 hetero 中定位 LDH 和 CMCH 的体相区域.

    用 marker 元素自动定位界面:
      CMCH marker = Mn（或 Co，若 Mn 不存在）
      LDH  marker = Ni（或 S，若 Ni 不存在）

    Returns: {"cmch": (start_idx, end_idx), "ldh": (start_idx, end_idx)}
    """
    grid = np.linspace(0, struct.lattice.c, ngz)

    # 从结构确定各材料的 z 范围
    cmch_markers = [s.coords[2] for s in struct if s.species_string in ("Mn", "C")]
    if not cmch_markers:
        cmch_markers = [s.coords[2] for s in struct if s.species_string == "Co"]
    ldh_markers = [s.coords[2] for s in struct if s.species_string in ("Ni", "S")]
    if not ldh_markers:
        ldh_markers = [s.coords[2] for s in struct if s.species_string == "O"]

    if not cmch_markers or not ldh_markers:
        # 备选: 用 z 坐标中位数分割
        z_all = sorted([s.coords[2] for s in struct])
        mid_z = z_all[len(z_all) // 2]
        cmch_z = [z for z in z_all if z <= mid_z]
        ldh_z = [z for z in z_all if z > mid_z]
    else:
        cmch_z = cmch_markers
        ldh_z = ldh_markers

    cmch_z_min, cmch_z_max = min(cmch_z), max(cmch_z)
    ldh_z_min, ldh_z_max = min(ldh_z), max(ldh_z)

    # 界面中点
    interface_mid = (cmch_z_max + ldh_z_min) / 2

    # 缓冲区: 从界面中点各退 1.5 Å
    buffer = 1.5
    cmch_top = interface_mid - buffer
    ldh_bottom = interface_mid + buffer

    # CMCH 体相: [cmch_z_min, cmch_top] 的中间 40%
    cmch_low = cmch_z_min + (cmch_top - cmch_z_min) * 0.3
    cmch_high = cmch_z_min + (cmch_top - cmch_z_min) * 0.7
    cmch_start = max(0, int(np.searchsorted(grid, cmch_low)))
    cmch_end = min(ngz, int(np.searchsorted(grid, cmch_high)))

    # LDH 体相: [ldh_bottom, ldh_z_max] 的中间 40%
    ldh_low = ldh_bottom + (ldh_z_max - ldh_bottom) * 0.3
    ldh_high = ldh_bottom + (ldh_z_max - ldh_bottom) * 0.7
    ldh_start = max(0, int(np.searchsorted(grid, ldh_low)))
    ldh_end = min(ngz, int(np.searchsorted(grid, ldh_high)))

    # 保护: 确保至少 5 个格点
    if cmch_end - cmch_start < 5:
        mid_c = (cmch_z_min + cmch_top) / 2
        span = (cmch_top - cmch_z_min) * 0.25
        cmch_start = max(0, int(np.searchsorted(grid, mid_c - span)))
        cmch_end = min(ngz, int(np.searchsorted(grid, mid_c + span)))
    if ldh_end - ldh_start < 5:
        mid_l = (ldh_bottom + ldh_z_max) / 2
        span = (ldh_z_max - ldh_bottom) * 0.25
        ldh_start = max(0, int(np.searchsorted(grid, mid_l - span)))
        ldh_end = min(ngz, int(np.searchsorted(grid, mid_l + span)))

    log.info("    Regions: CMCH z=[%.2f,%.2f] idx=[%d,%d], LDH z=[%.2f,%.2f] idx=[%d,%d], interface=%.2f",
             grid[cmch_start], grid[min(cmch_end, ngz-1)], cmch_start, cmch_end,
             grid[ldh_start], grid[min(ldh_end, ngz-1)], ldh_start, ldh_end, interface_mid)

    return {
        "cmch": (cmch_start, cmch_end),
        "ldh": (ldh_start, ldh_end),
    }


# ═══════════════════════════════════════════════
# Band alignment for one hetero
# ═══════════════════════════════════════════════

def compute_band_offset(het_name: str) -> dict:
    """对一个异质结计算能带对齐."""
    if het_name not in HETERO_MAP:
        return {"name": het_name, "error": "Unknown hetero"}
    ldh_name, cmch_name = HETERO_MAP[het_name]

    # 1. 读取各体系 LOCPOT
    avg_het, ngz_het, ef_het = read_locpot_planar(het_name)
    avg_ldh, ngz_ldh, _ = read_locpot_planar(ldh_name)
    avg_cmch, ngz_cmch, _ = read_locpot_planar(cmch_name)
    if avg_het is None or avg_ldh is None or avg_cmch is None:
        return {"name": het_name, "error": "Missing LOCPOT data"}

    struct_het = read_structure(het_name)
    struct_ldh = read_structure(ldh_name)
    struct_cmch = read_structure(cmch_name)
    if not all([struct_het, struct_ldh, struct_cmch]):
        return {"name": het_name, "error": "Missing structure data"}

    # 2. VBM/CBM from DOS
    vbm_ldh, cbm_ldh = read_vbm_cbm(ldh_name)
    vbm_cmch, cbm_cmch = read_vbm_cbm(cmch_name)
    gap_ldh = (cbm_ldh - vbm_ldh) if (vbm_ldh and cbm_ldh) else None
    gap_cmch = (cbm_cmch - vbm_cmch) if (vbm_cmch and cbm_cmch) else None

    if vbm_ldh is None or vbm_cmch is None:
        return {"name": het_name, "error": "Missing VBM/CBM data"}

    # 3. slab 体相平均势
    s_ldh, e_ldh = find_slab_bulk_region(avg_ldh, ngz_ldh, struct_ldh)
    s_cmch, e_cmch = find_slab_bulk_region(avg_cmch, ngz_cmch, struct_cmch)
    v_ref_ldh_slab = float(avg_ldh[s_ldh:e_ldh].mean())
    v_ref_cmch_slab = float(avg_cmch[s_cmch:e_cmch].mean())

    # 4. hetero 体相平均势
    regions = find_hetero_bulk_regions(avg_het, ngz_het, struct_het)
    sc, ec = regions["cmch"]
    sl, el = regions["ldh"]
    v_ref_cmch_het = float(avg_het[sc:ec].mean())
    v_ref_ldh_het = float(avg_het[sl:el].mean())

    # 5. 静电势偏移
    dv_cmch = v_ref_cmch_het - v_ref_cmch_slab
    dv_ldh = v_ref_ldh_het - v_ref_ldh_slab

    # 6. 对齐后的 VBM/CBM
    het_vbm_cmch = vbm_cmch + dv_cmch
    het_cbm_cmch = cbm_cmch + dv_cmch if cbm_cmch else None
    het_vbm_ldh = vbm_ldh + dv_ldh
    het_cbm_ldh = cbm_ldh + dv_ldh if cbm_ldh else None

    # 7. Band offsets
    vb_offset = het_vbm_ldh - het_vbm_cmch
    cb_offset = (het_cbm_ldh - het_cbm_cmch) if (het_cbm_ldh and het_cbm_cmch) else None

    # Type
    # vb_offset > 0: LDH VBM 高于 CMCH VBM (straddling gap → Type I)
    # vb_offset < 0: staggered → Type II
    band_type = "Type I (straddling)"
    if cb_offset is not None:
        if vb_offset > 0 and cb_offset > 0:
            band_type = "Type I (straddling)"
        elif vb_offset > 0 and cb_offset < 0:
            band_type = "Type II (staggered)"
        elif vb_offset < 0 and cb_offset > 0:
            band_type = "Type II (staggered)"
        else:
            band_type = "Type III (broken gap)"

    cbm_cmch_str = f"{het_cbm_cmch:.4f}" if het_cbm_cmch else "N/A"
    cbm_ldh_str  = f"{het_cbm_ldh:.4f}"  if het_cbm_ldh  else "N/A"
    cb_off_str   = f"{cb_offset:+.4f}"   if cb_offset    else "N/A"
    gap_ldh_str  = f"{gap_ldh:.4f}"      if gap_ldh      else "N/A"

    print(f"\n  {'':─^60}")
    print(f"  {het_name}  ({SYSTEM_LABELS.get(het_name, '')})")
    print(f"  {'':─^60}")
    print(f"  {'':<15s}  {'CMCH':>12s}  {'LDH':>12s}")
    print(f"  {'VBM (eV)':<15s}  {het_vbm_cmch:>12.4f}  {het_vbm_ldh:>12.4f}")
    print(f"  {'CBM (eV)':<15s}  {cbm_cmch_str:>12s}  {cbm_ldh_str:>12s}")
    print(f"  {'Gap (eV)':<15s}  {(cbm_cmch - vbm_cmch) if cbm_cmch else 0:>12.4f}  {gap_ldh_str:>12s}")
    print(f"  {'VB offset (eV)':<15s}  {vb_offset:+12.4f}")
    cb_text = f"{cb_offset:+.4f}" if cb_offset is not None else "N/A"
    print(f"  {'CB offset (eV)':<15s}  {cb_text:>12s}")
    print(f"  {'Type':<15s}  {band_type:>24s}")
    print(f"  {'DV_CMCH (eV)':<15s}  {dv_cmch:>+12.4f}")
    print(f"  {'DV_LDH (eV)':<15s}  {dv_ldh:>+12.4f}")

    return {
        "name": het_name,
        "label": SYSTEM_LABELS.get(het_name, het_name),
        "ldh_name": ldh_name,
        "cmch_name": cmch_name,
        "vbm_cmch": float(het_vbm_cmch),
        "cbm_cmch": float(het_cbm_cmch) if het_cbm_cmch else None,
        "vbm_ldh": float(het_vbm_ldh),
        "cbm_ldh": float(het_cbm_ldh) if het_cbm_ldh else None,
        "vb_offset": float(vb_offset),
        "cb_offset": float(cb_offset) if cb_offset else None,
        "band_type": band_type,
        "dv_cmch": float(dv_cmch),
        "dv_ldh": float(dv_ldh),
        "v_ref_cmch_slab": v_ref_cmch_slab,
        "v_ref_ldh_slab": v_ref_ldh_slab,
        "v_ref_cmch_het": v_ref_cmch_het,
        "v_ref_ldh_het": v_ref_ldh_het,
        # For plotting
        "avg_het": avg_het.tolist(),
        "ngz_het": ngz_het,
        "c_het": struct_het.lattice.c,
        "cmch_region": [int(sc), int(ec)],
        "ldh_region": [int(sl), int(el)],
    }


# ═══════════════════════════════════════════════
# Plot
# ═══════════════════════════════════════════════

def plot_band_alignment(results: dict):
    """绘制所有异质结的能带对齐图."""
    hets = ["hetero_intrinsic", "hetero_s_doped", "hetero_s_exposed"]
    data = []
    for h in hets:
        r = results.get(h)
        if r and "error" not in r:
            data.append(r)

    if not data:
        log.warning("No band alignment data to plot")
        return

    fig, ax = plt.subplots(figsize=(9, 5.5))

    n = len(data)
    bar_width = 0.25
    x_base = np.arange(n)

    labels = [d["label"] for d in data]
    colors = ["#2171b5", "#2ca02c"]  # CMCH=blue, LDH=green

    for i, d in enumerate(data):
        x = x_base[i]
        # CMCH band
        vbm_c = d["vbm_cmch"]
        cbm_c = d.get("cbm_cmch")
        gap_c = (cbm_c - vbm_c) if cbm_c else 0.3
        if not cbm_c:
            gap_c = 0.3
        ax.bar(x - bar_width / 2, gap_c, bar_width, bottom=vbm_c,
               color=colors[0], alpha=0.85, edgecolor="k", lw=0.5)

        # LDH band
        vbm_l = d["vbm_ldh"]
        cbm_l = d.get("cbm_ldh")
        gap_l = (cbm_l - vbm_l) if cbm_l else 0.3
        if not cbm_l:
            gap_l = 0.3
        ax.bar(x + bar_width / 2, gap_l, bar_width, bottom=vbm_l,
               color=colors[1], alpha=0.85, edgecolor="k", lw=0.5)

        # Label VBM inside bottom of bar, CBM above top of bar
        ax.text(x - bar_width / 2, vbm_c - 0.12,
                f"VBM {vbm_c:.2f}", ha="center", va="top", fontsize=7.5)
        if cbm_c:
            ax.text(x - bar_width / 2, vbm_c + gap_c + 0.08,
                    f"CBM {cbm_c:.2f}", ha="center", va="bottom", fontsize=7.5)
        else:
            ax.text(x - bar_width / 2, vbm_c + gap_c + 0.08,
                    "CBM N/A", ha="center", va="bottom", fontsize=7.5)

        ax.text(x + bar_width / 2, vbm_l - 0.12,
                f"VBM {vbm_l:.2f}", ha="center", va="top", fontsize=7.5)
        if cbm_l:
            ax.text(x + bar_width / 2, vbm_l + gap_l + 0.08,
                    f"CBM {cbm_l:.2f}", ha="center", va="bottom", fontsize=7.5)
        else:
            ax.text(x + bar_width / 2, vbm_l + gap_l + 0.08,
                    "CBM N/A", ha="center", va="bottom", fontsize=7.5)

        # Gap annotation inside each bar
        if gap_c > 0:
            ax.text(x - bar_width / 2, vbm_c + gap_c / 2,
                    f"{gap_c:.2f}", ha="center", va="center", fontsize=6.5,
                    color="white", fontweight="bold")
        if gap_l > 0:
            ax.text(x + bar_width / 2, vbm_l + gap_l / 2,
                    f"{gap_l:.2f}", ha="center", va="center", fontsize=6.5,
                    color="white", fontweight="bold")

        # VB offset: horizontal double-headed arrow between VBM levels
        vo = d["vb_offset"]
        ax.annotate("", xy=(x + bar_width / 2, vbm_l), xytext=(x - bar_width / 2, vbm_c),
                    arrowprops=dict(arrowstyle="<->", color="#555555", lw=1.0,
                                    shrinkA=0, shrinkB=0))
        mid_y = (vbm_c + vbm_l) / 2
        ax.text(x, mid_y + 0.18, f"VB off: {vo:+.3f} eV",
                ha="center", va="bottom", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.2", fc="lightyellow", alpha=0.85, ec="none"))

        # Type annotation
        ax.text(x, max(vbm_c, vbm_l) + max(gap_c, gap_l) + 0.35,
                d["band_type"], ha="center", fontsize=8, fontstyle="italic",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, ec="gray"))

    ax.set_xticks(x_base)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("Energy (eV) vs. CMCH bulk reference", fontsize=11)
    ax.set_title("Band Alignment — Heterojunctions", fontweight="bold", fontsize=12)
    ax.legend(["CMCH (substrate)", "LDH (film)"], fontsize=9, loc="upper center",
              ncol=2, framealpha=0.8)
    ax.axhline(0, color="k", ls="--", lw=0.6, alpha=0.4)
    ax.grid(True, alpha=0.2, axis="y")

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "band_alignment.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("Plot → output/band_alignment.png")

    # Second plot: electrostatic potential showing bulk regions
    plot_bulk_regions(results)


def plot_bulk_regions(results: dict):
    """静电势曲线 + 体相区域标注（含界面位置与原始网格点）."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    hets = ["hetero_intrinsic", "hetero_s_doped", "hetero_s_exposed"]

    for ax_idx, het_name in enumerate(hets):
        ax = axes[ax_idx]
        r = results.get(het_name)
        if not r or "error" in r:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(het_name)
            continue

        grid = np.linspace(0, r["c_het"], r["ngz_het"])
        avg = np.array(r["avg_het"])

        # Raw grid points (faint) + smooth line
        ax.scatter(grid, avg, s=1.5, c="#cccccc", zorder=1, label="Raw grid")
        ax.plot(grid, avg, "k-", lw=0.8, zorder=2)

        # Shade bulk regions
        sc, ec = r["cmch_region"]
        sl, el = r["ldh_region"]
        ax.axvspan(grid[sc], grid[ec], color="blue", alpha=0.08, zorder=0)
        ax.axvspan(grid[sl], grid[el], color="green", alpha=0.08, zorder=0)

        # Bulk region edge markers
        for idx, c, lab in [(sc, "blue", "CMCH bulk"), (ec, "blue", None),
                            (sl, "green", "LDH bulk"), (el, "green", None)]:
            ax.axvline(grid[idx], color=c, ls=":", lw=0.7, alpha=0.5, zorder=0)

        # Interface midpoint vertical dashed line
        cmch_top_z = grid[ec] if ec < len(grid) else grid[-1]
        ldh_bot_z  = grid[sl] if sl < len(grid) else grid[0]
        interface_z = (cmch_top_z + ldh_bot_z) / 2
        ax.axvline(interface_z, color="red", ls="--", lw=1.0, alpha=0.7, zorder=3)
        ax.text(interface_z, avg.max(), "Interface",
                ha="center", va="bottom", fontsize=7.5,
                color="red", rotation=90, zorder=4)

        # Legend entries
        legend_handles = [
            Patch(facecolor="blue", alpha=0.15, edgecolor="blue", label="CMCH bulk region"),
            Patch(facecolor="green", alpha=0.15, edgecolor="green", label="LDH bulk region"),
            plt.Line2D([0], [0], color="red", ls="--", lw=1.0, label="Interface"),
        ]

        ax.set_xlabel("z / Å", fontsize=11)
        ax.set_ylabel("V(z) / eV", fontsize=11)
        ax.set_title(r["label"], fontsize=10, fontweight="bold")
        ax.tick_params(labelsize=9)
        ax.grid(True, alpha=0.2)
        if ax_idx == 2:
            ax.legend(handles=legend_handles, fontsize=8, loc="upper right")

    # Shared footnote
    fig.text(0.5, -0.02,
             "Bulk regions: central 40% of each slab z-range, excluding 3 Å interfacial buffer.",
             ha="center", fontsize=8, fontstyle="italic", color="gray")

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "locpot_bulk_regions.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("Plot → output/locpot_bulk_regions.png")


# ═══════════════════════════════════════════════
# Data export (CSV / JSON)
# ═══════════════════════════════════════════════

def export_band_alignment_csv(results: dict):
    """导出能带对齐汇总表 → output/band_alignment.csv"""
    rows = []
    for h in ["hetero_intrinsic", "hetero_s_doped", "hetero_s_exposed"]:
        r = results.get(h)
        if not r or "error" in r:
            continue
        rows.append({
            "heterojunction": r["label"],
            "VBM_CMCH (eV)": r["vbm_cmch"],
            "CBM_CMCH (eV)": r.get("cbm_cmch"),
            "gap_CMCH (eV)": (r["cbm_cmch"] - r["vbm_cmch"]) if r.get("cbm_cmch") else None,
            "VBM_LDH (eV)": r["vbm_ldh"],
            "CBM_LDH (eV)": r.get("cbm_ldh"),
            "gap_LDH (eV)": (r["cbm_ldh"] - r["vbm_ldh"]) if r.get("cbm_ldh") else None,
            "VB_offset (eV)": r["vb_offset"],
            "CB_offset (eV)": r.get("cb_offset"),
            "dV_CMCH (eV)": r["dv_cmch"],
            "dV_LDH (eV)": r["dv_ldh"],
            "band_type": r["band_type"],
        })

    header = list(rows[0].keys())
    with open(OUTPUT_DIR / "band_alignment.csv", "w") as f:
        f.write(",".join(header) + "\n")
        for row in rows:
            vals = [f"{row[k]}" if row[k] is not None else "N/A" for k in header]
            f.write(",".join(vals) + "\n")
    log.info("Data  → output/band_alignment.csv")


def export_locpot_csv(results: dict):
    """每个异质结导出 z / V(z) / region_mask → output/locpot_planar_{het}.csv"""
    hets = ["hetero_intrinsic", "hetero_s_doped", "hetero_s_exposed"]
    for het_name in hets:
        r = results.get(het_name)
        if not r or "error" in r:
            continue

        ngz = r["ngz_het"]
        c   = r["c_het"]
        avg = np.array(r["avg_het"])
        sc, ec = r["cmch_region"]
        sl, el = r["ldh_region"]

        grid = np.linspace(0, c, ngz)
        # region_mask: 0=none, 1=CMCH, 2=LDH
        mask = np.zeros(ngz, dtype=int)
        mask[sc:ec] = 1
        mask[sl:el] = 2

        with open(OUTPUT_DIR / f"locpot_planar_{het_name}.csv", "w") as f:
            f.write("z_Ang,V_eV,region\n")
            region_labels = {0: "none", 1: "CMCH_bulk", 2: "LDH_bulk"}
            for iz in range(ngz):
                f.write(f"{grid[iz]:.4f},{avg[iz]:.4f},{region_labels[mask[iz]]}\n")
    log.info("Data  → output/locpot_planar_*.csv")


def _to_serializable(obj):
    """递归转换 numpy 类型为 Python 原生类型，供 json.dump 使用."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj) if isinstance(obj, np.floating) else int(obj)
    if isinstance(obj, dict):
        return {k: _to_serializable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_serializable(v) for v in obj]
    return obj


def export_results_json(results: dict):
    """完整计算结果 → output/band_alignment.json（含 avg_het 曲线数据）"""
    with open(OUTPUT_DIR / "band_alignment.json", "w") as f:
        json.dump(_to_serializable(results), f, indent=2, ensure_ascii=False)
    log.info("Data  → output/band_alignment.json")


def load_results() -> dict | None:
    """从 JSON 加载上次的结果（含 avg_het 曲线数据）."""
    path = OUTPUT_DIR / "band_alignment.json"
    if not path.exists():
        return None
    try:
        with open(path) as f:
            data = json.load(f)
        # avg_het 从 list 转为 np.array（JSON 反序列化为 list）
        for h in data.values():
            if "avg_het" in h and isinstance(h["avg_het"], list):
                h["avg_het"] = np.array(h["avg_het"])
        log.info("Loaded from %s", path.name)
        return data
    except Exception as e:
        log.warning("Cannot load %s: %s", path.name, e)
        return None


# ═══════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Band alignment analysis")
    parser.add_argument("--het", type=str, help="指定异质结名称（默认全部）")
    parser.add_argument("--no-plot", action="store_true", help="仅导出数据，不生成图片")
    parser.add_argument("--force", action="store_true", help="强制重新计算（忽略已有 JSON）")
    args = parser.parse_args()

    hets = [args.het] if args.het else list(HETERO_MAP.keys())

    # 优先从 JSON 恢复（avg_het 曲线数据已包含在 JSON 中）
    results = None if args.force else load_results()

    if results is None:
        log.info("Band alignment for %d heterojunctions", len(hets))
        results = {}
        for h in hets:
            log.info("Processing %s ...", h)
            r = compute_band_offset(h)
            results[h] = r
            if "error" in r:
                log.warning("  ✗ %s: %s", h, r["error"])
            else:
                log.info("  ✓ %s: VB offset = %+.4f eV (%s)", h, r["vb_offset"], r["band_type"])
    else:
        log.info("Loaded %d heterojunction results from JSON", len(results))

    # Final summary
    print("\n" + "=" * 80)
    print("BAND ALIGNMENT SUMMARY")
    print("=" * 80)
    print(f"  {'Het.':<22s}  {'VB off.':>10s}  {'CB off.':>10s}  {'DV_CMCH':>10s}  {'DV_LDH':>10s}  {'Type':<22s}")
    print(f"  {'─'*22}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*22}")
    for h in hets:
        r = results.get(h, {})
        if "error" in r:
            print(f"  {SYSTEM_LABELS.get(h, h):<22s}  ERROR: {r['error']}")
        else:
            cb = f"{r.get('cb_offset', 0):+.4f}" if r.get('cb_offset') is not None else "N/A"
            print(f"  {r['label']:<22s}  {r['vb_offset']:>+10.4f}  {cb:>10s}  "
                  f"{r['dv_cmch']:>+10.4f}  {r['dv_ldh']:>+10.4f}  {r['band_type']:<22s}")

    # Data export
    log.info("Exporting data...")
    export_band_alignment_csv(results)
    export_locpot_csv(results)
    export_results_json(results)

    # Plots
    if not args.no_plot:
        log.info("Generating plots...")
        plot_band_alignment(results)
    log.info("Done. Output → %s", OUTPUT_DIR)


if __name__ == "__main__":
    main()
