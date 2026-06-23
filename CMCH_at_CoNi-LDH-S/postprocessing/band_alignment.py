"""
能带对齐（Band Offset）分析 — 平均势法

原理 (Van de Walle & Martin):
  孤立 slab 的 VBM/CBM 以自身 E_F 为参考。异质结中两个材料的 E_F 对齐后，
  需要找到与 E_F 无关的参考势来比较二者的带边偏移。

  参考势 = LOCPOT planar average 在体相区域的平均值

  对于 slab: V_ref = V(z) 在 slab 中心体相区的平均值
  对于 hetero: 分别提取两材料体相区的 V_ref_het
  ΔV = V_ref_het - V_ref_slab  → 该材料在异质结中的静电势偏移
  VBM_het = VBM_slab + ΔV,  CBM_het = CBM_slab + ΔV

用法:
  python band_alignment.py                          # 分析3个异质结
  python band_alignment.py --het hetero_intrinsic    # 单个异质结

输出:
  终端: 能带对齐参数表
  output/band_alignment.png / .pdf (能带对齐图)
  output/locpot_bulk_regions.png / .pdf (体相区标识图)
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
from matplotlib.patches import FancyBboxPatch
from monty.json import MontyDecoder
from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Locpot

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.linewidth": 1.2, "figure.dpi": 100})

DATA_DIR = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S/data")
OUTPUT_DIR = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S/postprocessing/output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 异质结 -> (LDH slab, CMCH slab)
HETERO_MAP = {
    "hetero_intrinsic": ("CoNiOH2-slab", "CoMnH2CO5-slab"),
    "hetero_s_doped": ("CoNiOH2S-noH-slab", "CoMnH2CO5-slab"),
    "hetero_s_exposed": ("CoNiOH2S-noH-slab-flip", "CoMnH2CO5-slab"),
}

SYSTEM_LABELS = {
    "CoNiOH2-slab": "CoNiOH2 (Slab)",
    "CoNiOH2S-noH-slab": "CoNiOH2S (Slab)",
    "CoNiOH2S-noH-slab-flip": "CoNiOH2S (Slab-flip)",
    "CoMnH2CO5-slab": "CoMnH2CO5 (Slab)",
    "hetero_intrinsic": "Intrinsic Het.",
    "hetero_s_doped": "S-Doped Het.",
    "hetero_s_exposed": "S-Exposed Het.",
}

SLAB_CMCH = "CoMnH2CO5-slab"


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

def _get_element_summary(struct: Structure) -> dict:
    """返回每种元素在 z 方向的分布范围."""
    elem_info = {}
    for i, site in enumerate(struct):
        sym = site.species_string
        z = site.coords[2]
        if sym not in elem_info:
            elem_info[sym] = {"z_min": z, "z_max": z, "indices": []}
        elem_info[sym]["z_min"] = min(elem_info[sym]["z_min"], z)
        elem_info[sym]["z_max"] = max(elem_info[sym]["z_max"], z)
        elem_info[sym]["indices"].append(i)
    return elem_info


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

    cbm_str = f"{het_cbm_ldh:.4f}" if het_cbm_ldh else "N/A"
    cb_off_str = f"{cb_offset:+.4f}" if cb_offset else "N/A"
    gap_str = f"{gap_ldh:.4f}" if gap_ldh else "N/A"

    print(f"\n  {'':─^60}")
    print(f"  {het_name}  ({SYSTEM_LABELS.get(het_name, '')})")
    print(f"  {'':─^60}")
    print(f"  {'':<15s}  {'CMCH':>12s}  {'LDH':>12s}")
    print(f"  {'VBM (eV)':<15s}  {het_vbm_cmch:>12.4f}  {het_vbm_ldh:>12.4f}")
    print(f"  {'CBM (eV)':<15s}  {cbm_str:>12s}  {cbm_str if het_cbm_ldh else 'N/A':>12s}")
    print(f"  {'Gap (eV)':<15s}  {gap_cmch if gap_cmch else 0:>12.4f}  {gap_str:>12s}")
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

    fig, ax = plt.subplots(figsize=(9, 5))

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

        # Label values
        ax.text(x - bar_width / 2, vbm_c - 0.15, f"{vbm_c:.2f}", ha="center", va="top", fontsize=7)
        ax.text(x - bar_width / 2, vbm_c + gap_c + 0.05, f"{cbm_c:.2f}" if cbm_c else "N/A",
                ha="center", va="bottom", fontsize=7)
        ax.text(x + bar_width / 2, vbm_l - 0.15, f"{vbm_l:.2f}", ha="center", va="top", fontsize=7)
        ax.text(x + bar_width / 2, vbm_l + gap_l + 0.05, f"{cbm_l:.2f}" if cbm_l else "N/A",
                ha="center", va="bottom", fontsize=7)

        # VB offset annotation
        mid = x
        vo = d["vb_offset"]
        ax.annotate(f"VB off: {vo:+.3f} eV",
                    xy=(mid, (vbm_c + vbm_l) / 2),
                    xytext=(mid + 0.5, (vbm_c + vbm_l) / 2 + 0.5),
                    fontsize=8, ha="center",
                    arrowprops=dict(arrowstyle="->", color="gray", lw=0.8),
                    bbox=dict(boxstyle="round,pad=0.2", fc="lightyellow", alpha=0.8, ec="none"))

        # Type annotation
        ax.text(x, max(vbm_c, vbm_l) + max(gap_c, gap_l) + 0.3,
                d["band_type"], ha="center", fontsize=8, fontstyle="italic",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, ec="gray"))

    ax.set_xticks(x_base)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("Energy (eV)")
    ax.set_title("Band Alignment — Heterojunctions", fontweight="bold")
    ax.legend(["CMCH (substrate)", "LDH (film)"], fontsize=9)
    ax.axhline(0, color="k", ls="--", lw=0.6, alpha=0.4)
    ax.grid(True, alpha=0.2, axis="y")

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "band_alignment.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "band_alignment.pdf", bbox_inches="tight")
    plt.close(fig)
    log.info("Plot → output/band_alignment.{png,pdf}")

    # Second plot: electrostatic potential showing bulk regions
    plot_bulk_regions(results)


def plot_bulk_regions(results: dict):
    """静电势曲线 + 体相区域标注."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    hets = ["hetero_intrinsic", "hetero_s_doped", "hetero_s_exposed"]

    for ax_idx, het_name in enumerate(hets):
        ax = axes[ax_idx]
        r = results.get(het_name)
        if not r or "error" in r:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(het_name)
            continue

        grid = np.linspace(0, r["c_het"], r["ngz_het"])
        ax.plot(grid, r["avg_het"], "k-", lw=0.8)

        # Shade bulk regions
        sc, ec = r["cmch_region"]
        sl, el = r["ldh_region"]
        ax.axvspan(grid[sc], grid[ec], color="blue", alpha=0.1, label="CMCH bulk")
        ax.axvspan(grid[sl], grid[el], color="green", alpha=0.1, label="LDH bulk")

        ax.set_xlabel("z (Å)")
        ax.set_ylabel("V(z) (eV)")
        ax.set_title(r["label"], fontsize=10, fontweight="bold")
        ax.grid(True, alpha=0.2)
        if ax_idx == 2:
            ax.legend(fontsize=7)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "locpot_bulk_regions.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "locpot_bulk_regions.pdf", bbox_inches="tight")
    plt.close(fig)
    log.info("Plot → output/locpot_bulk_regions.{png,pdf}")


# ═══════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Band alignment analysis")
    parser.add_argument("--het", type=str, help="指定异质结名称（默认全部）")
    args = parser.parse_args()

    hets = [args.het] if args.het else list(HETERO_MAP.keys())

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

    # Plots
    log.info("Generating plots...")
    plot_band_alignment(results)
    log.info("Done. Output → %s", OUTPUT_DIR)


if __name__ == "__main__":
    main()
