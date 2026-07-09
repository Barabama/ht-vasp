"""
界面差分电荷分析 — LOCPOT 静电势 Poisson 方程法

原理:
  CHGCAR 直接减法不可行 → 核心电子密度在原子核附近极高 (>1000 e⁻/Å³),
  异质结与孤立 slab 之间的微小原子位移会在此区域产生 4-5 个数量级的伪影.

  改用 LOCPOT 静电势法:
  1. 从 LOCPOT 提取平面平均静电势 V(z) (光滑函数, 无核心尖峰)
  2. ΔV(z) = V_hetero(z) − V_CMCH(z) − V_LDH_shifted(z)
  3. 泊松方程 Δρ(z) = −ε₀ · d²(ΔV)/dz² → 电荷重分布
  4. 界面区域积分 → 净电荷转移

参考: 所有 7 个 Level 2 体系共享同一 FFT 网格 (160,140,420),
      无需插值对齐.

用法:
  python postprocessing/charge_difference_3d.py
  python postprocessing/charge_difference_3d.py --force
  python postprocessing/charge_difference_3d.py --no-plot
"""

import argparse
import json
import gzip
import logging
import tempfile
import shutil
import os
import warnings
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from pymatgen.io.vasp.outputs import Locpot
from pymatgen.core import Structure

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

plt.rcParams.update({"font.family": "sans-serif", "font.size": 11, "axes.linewidth": 1.2})

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "data"
OUTPUT_DIR = SCRIPT_DIR / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

HETERO_MAP = {
    "hetero_intrinsic": ("LDH_strained", "CMCH_strained"),
    "hetero_s_doped": ("LDH_S_strained", "CMCH_strained"),
    "hetero_s_exposed": ("LDH_S_flip_strained", "CMCH_strained"),
}

SYSTEM_LABELS = {
    "CMCH_strained": "CMCH_slab (Strained)",
    "LDH_strained": "LDH_slab (Strained)",
    "LDH_S_strained": "LDH_S_slab (Strained)",
    "LDH_S_flip_strained": "LDH_S_slab (Strained-flip)",
    "hetero_intrinsic": "Intrinsic Het.",
    "hetero_s_doped": "S-Doped Het.",
    "hetero_s_exposed": "S-Exposed Het.",
}

COLORS = {"hetero_intrinsic": "#1f77b4", "hetero_s_doped": "#2ca02c", "hetero_s_exposed": "#ff7f0e"}

# ε₀ in e²/(eV·Å) — for Poisson: Δρ = −ε₀ × d²V/dz²
EPS0 = 0.00552635


# ═══════════════════════════════════════════════
# LOCPOT 加载
# ═══════════════════════════════════════════════

def load_locpot(name: str):
    """加载 LOCPOT, 返回 (planar_avg, ngz, c)."""
    static = DATA_DIR / name / "3-static"
    locpot = static / "LOCPOT_dipole.gz"
    if not locpot.exists():
        locpot = static / "LOCPOT.gz"
    if not locpot.exists():
        return None, 0, 0

    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".LOCPOT")
    try:
        with gzip.open(locpot, "rb") as f_in:
            shutil.copyfileobj(f_in, tmp)
        tmp_path = tmp.name
        tmp.close()
        lp = Locpot.from_file(tmp_path)
        avg = lp.get_average_along_axis(2)
        ngz = lp.dim[2]
        with gzip.open(static / "CONTCAR.gz", "rt") as f:
            struct = Structure.from_str(f.read(), fmt="poscar")
        c = struct.lattice.c
        return avg, ngz, c
    finally:
        try:
            os.unlink(tmp_path)
        except Exception:
            pass


def get_z_offset(s_ldh: Structure, s_het: Structure) -> float:
    """计算 LDH slab 在孤立 slab 与异质结之间的 z 偏移 (Å)."""
    z_ldh_list = [s.coords[2] for s in s_ldh if s.species_string in ("Ni", "S")]
    z_het_list = [s.coords[2] for s in s_het if s.species_string in ("Ni", "S")]
    if not z_ldh_list or not z_het_list:
        return 0.0
    z_ldh = np.mean(z_ldh_list)
    z_het = np.mean(z_het_list)
    return float(z_het - z_ldh)


# ═══════════════════════════════════════════════
# 核心计算 — LOCPOT Poisson 法
# ═══════════════════════════════════════════════

def compute_interface_analysis(het_name: str) -> dict:
    """界面静电势差分 + Poisson 反推电荷重分布."""
    ldh_name, cmch_name = HETERO_MAP[het_name]

    avg_het, ngz_het, c_het = load_locpot(het_name)
    avg_cmch, ngz_cmch, c_cmch = load_locpot(cmch_name)
    avg_ldh, ngz_ldh, c_ldh = load_locpot(ldh_name)
    if avg_het is None:
        return {"name": het_name, "error": "Missing LOCPOT"}

    # 结构加载 (用于 z 偏移和界面标记)
    s_het = Structure.from_file(str(DATA_DIR / het_name / "3-static" / "CONTCAR.gz"))
    s_cmch = Structure.from_file(str(DATA_DIR / cmch_name / "3-static" / "CONTCAR.gz"))
    s_ldh = Structure.from_file(str(DATA_DIR / ldh_name / "3-static" / "CONTCAR.gz"))

    grid_het = np.linspace(0, c_het, ngz_het)

    # 1. CMCH: 底部对齐, 无需偏移
    f_cmch = interp1d(np.linspace(0, c_cmch, ngz_cmch), avg_cmch,
                      kind="linear", bounds_error=False, fill_value=np.nan)
    v_cmch = f_cmch(grid_het)

    # 2. LDH: 需要 z 方向偏移
    z_offset = get_z_offset(s_ldh, s_het)
    f_ldh = interp1d(np.linspace(0, c_ldh, ngz_ldh), avg_ldh,
                     kind="linear", bounds_error=False, fill_value=np.nan)
    v_ldh = f_ldh(grid_het - z_offset)

    # 3. ΔV(z)
    if np.any(np.isnan(v_ldh)):
        log.warning("NaNs found in interpolated LDH potential — z_offset may be out of range")
        v_ldh = np.nan_to_num(v_ldh, nan=0.0)
    delta_v = avg_het - v_cmch - v_ldh

    # 4. 界面标记
    cmch_coords = [s.coords[2] for i, s in enumerate(s_het)
                   if s.species_string in ("Mn", "C")]
    ldh_coords = [s.coords[2] for i, s in enumerate(s_het)
                  if s.species_string in ("Ni", "S")]
    cmch_top = max(cmch_coords) if cmch_coords else 0.0
    ldh_bot = min(ldh_coords) if ldh_coords else 0.0
    interface_z = float((cmch_top + ldh_bot) / 2)
    area = s_het.lattice.a * s_het.lattice.b

    # 5. Poisson: Δρ = −ε₀ · d²(ΔV)/dz²
    dv_smooth = gaussian_filter1d(delta_v, sigma=2.0)
    d2v = np.gradient(np.gradient(dv_smooth, grid_het), grid_het)
    delta_rho = -EPS0 * d2v

    # 6. 累计电荷 Q(z) = ∫₀ᶻ Δρ(z') dz' · area
    q_cum = np.array([np.trapezoid(delta_rho[:i + 1], grid_het[:i + 1]) * area
                      for i in range(len(grid_het))])
    q_total = float(q_cum[-1])

    # 7. 界面区域净电荷
    imask = (grid_het >= cmch_top) & (grid_het <= ldh_bot)
    q_intf = float(np.trapezoid(delta_rho[imask], grid_het[imask]) * area) if imask.sum() > 2 else 0.0

    # 8. 界面偶极矩
    dv_intf = delta_v[imask]
    dipole_v = float(np.trapezoid(dv_intf, grid_het[imask])) if imask.sum() > 2 and np.isfinite(dv_intf).all() else 0.0

    log.info("  z_offset=%.2f Å | ΔV=[%.2f, %.2f] eV | Q_intf=%+.6f e⁻ | Dipole=%.4f eV·Å",
             z_offset, delta_v.min(), delta_v.max(), q_intf, dipole_v)

    return {
        "name": het_name,
        "label": SYSTEM_LABELS.get(het_name, het_name),
        "c": float(c_het), "ngz": ngz_het,
        "grid_z": grid_het.tolist(),
        "V_planar_avg": delta_v.tolist(),
        "rho_planar_avg": delta_rho.tolist(),
        "q_cumulative": q_cum.tolist(),
        "z_offset": float(z_offset),
        "interface_z": interface_z,
        "cmch_top_z": float(cmch_top),
        "ldh_bot_z": float(ldh_bot),
        "charge_transfer_total": q_total,
        "charge_transfer_interface": q_intf,
        "dipole_V": float(dipole_v),
        "delta_V_range": [float(delta_v.min()), float(delta_v.max())],
        "area": float(area),
    }


# ═══════════════════════════════════════════════
# 数据导出
# ═══════════════════════════════════════════════

def _serialize(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj) if isinstance(obj, np.floating) else int(obj)
    if isinstance(obj, dict):
        return {k: _serialize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_serialize(v) for v in obj]
    return obj


def export_summary_csv(results: dict):
    """汇总表 → output/charge_diff_summary.csv"""
    header = ["system", "label", "z_offset_A", "Q_total_e", "Q_interface_e",
              "dipole_V_eV_A", "delta_V_min_eV", "delta_V_max_eV", "area_A2", "interface_z_A"]
    with open(OUTPUT_DIR / "charge_diff_summary.csv", "w") as f:
        f.write(",".join(header) + "\n")
        for het_name in HETERO_MAP:
            r = results.get(het_name, {})
            if "error" in r:
                continue
            f.write(f"{het_name},{r['label']},{r['z_offset']:.4f},{r['charge_transfer_total']:.6f},"
                    f"{r['charge_transfer_interface']:.6f},{r['dipole_V']:.4f},"
                    f"{r['delta_V_range'][0]:.4f},{r['delta_V_range'][1]:.4f},"
                    f"{r['area']:.4f},{r['interface_z']:.4f}\n")
    log.info("CSV → charge_diff_summary.csv")


def export_profiles_csv(results: dict):
    """曲线数据 (长格式) → output/charge_diff_profiles.csv, Origin/Excel 用"""
    rows = []
    for het_name in HETERO_MAP:
        r = results.get(het_name, {})
        if "error" in r or "grid_z" not in r:
            continue
        g, dv, dr, qc = r["grid_z"], r["V_planar_avg"], r["rho_planar_avg"], r["q_cumulative"]
        for i in range(len(g)):
            rows.append([het_name, g[i], dv[i], dr[i], qc[i]])
    header = ["system", "z_A", "delta_V_eV", "delta_rho_e_per_A3", "Q_cumulative_e"]
    with open(OUTPUT_DIR / "charge_diff_profiles.csv", "w") as f:
        f.write(",".join(header) + "\n")
        for row in rows:
            f.write(",".join(str(v) for v in row) + "\n")
    log.info("CSV → charge_diff_profiles.csv")


def export_json(results: dict):
    """完整数值 → output/charge_diff.json"""
    path = OUTPUT_DIR / "charge_diff.json"
    path.write_text(json.dumps(_serialize(results), indent=2, ensure_ascii=False))
    log.info("JSON → charge_diff.json")


def load_results() -> dict | None:
    """从 JSON 恢复上次运行结果."""
    path = OUTPUT_DIR / "charge_diff.json"
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


# ═══════════════════════════════════════════════
# 绘图
# ═══════════════════════════════════════════════

def plot_interface_profile(name: str, r: dict):
    """单体系三面板: ΔV(z), Δρ(z), Q(z)."""
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 7), sharex=True)
    grid = np.array(r["grid_z"])
    dv = np.array(r["V_planar_avg"])
    dr = np.array(r["rho_planar_avg"])
    qc = np.array(r["q_cumulative"])
    cmch_top, ldh_bot, iface = r["cmch_top_z"], r["ldh_bot_z"], r["interface_z"]

    def shade(ax):
        ax.axvspan(0, cmch_top, color="blue", alpha=0.05)
        ax.axvspan(ldh_bot, r["c"], color="green", alpha=0.05)
        ax.axvline(iface, color="purple", ls="--", lw=0.8, alpha=0.6)

    ax1.plot(grid, dv, "k-", lw=1.2); shade(ax1)
    ax1.axhline(0, color="gray", ls="--", lw=0.5)
    ax1.set_ylabel("ΔV(z) (eV)")
    ax1.set_title(f"Interface Analysis — {r['label']}", fontweight="bold")
    ax1.grid(True, alpha=0.2)

    ax2.plot(grid, dr, "k-", lw=1.2); shade(ax2)
    ax2.fill_between(grid, 0, dr, where=(dr > 0), color="red", alpha=0.15, label="Accumulation")
    ax2.fill_between(grid, 0, dr, where=(dr < 0), color="blue", alpha=0.15, label="Depletion")
    ax2.axhline(0, color="gray", ls="--", lw=0.5)
    ax2.set_ylabel("Δρ(z) (e⁻/Å³)"); ax2.legend(fontsize=8); ax2.grid(True, alpha=0.2)

    ax3.plot(grid, qc, "k-", lw=1.2); shade(ax3)
    ax3.axhline(0, color="gray", ls="--", lw=0.5)
    ax3.set_xlabel("z (Å)"); ax3.set_ylabel("Cumulative Q(z) (e⁻)"); ax3.grid(True, alpha=0.2)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / f"charge_diff_{r['name']}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("PNG → charge_diff_%s.png", r['name'])


def plot_comparison(all_results: dict):
    """三异质结 ΔV(z) 叠加."""
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for het_name in HETERO_MAP:
        r = all_results.get(het_name)
        if not r or "error" in r:
            continue
        ax.plot(np.array(r["grid_z"]), np.array(r["V_planar_avg"]),
                color=COLORS.get(het_name, "#333"), lw=1.2,
                label=SYSTEM_LABELS.get(het_name, het_name))
    ax.axhline(0, color="gray", ls="--", lw=0.5)
    ax.set_xlabel("z / Å"); ax.set_ylabel("ΔV(z) / eV")
    ax.set_title("Interface Potential Difference — Comparison", fontweight="bold")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.2)
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "charge_diff_comparison.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("PNG → charge_diff_comparison.png")


# ═══════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Interface charge difference analysis (LOCPOT Poisson)")
    parser.add_argument("--no-plot", action="store_true", help="仅导出数据不绘图")
    parser.add_argument("--force", action="store_true", help="强制重算（忽略缓存的 JSON）")
    args = parser.parse_args()

    # 从 JSON 恢复或重新计算
    all_results = None if args.force else load_results()

    if all_results is None:
        log.info("LOCPOT-based interface analysis for %d heterojunctions", len(HETERO_MAP))
        all_results = {}
        for het_name in HETERO_MAP:
            log.info("Processing %s ...", het_name)
            r = compute_interface_analysis(het_name)
            all_results[het_name] = r
            if "error" in r:
                log.warning("  ✗ %s: %s", het_name, r["error"])
            else:
                log.info("  ✓ %s", het_name)
    else:
        log.info("Loaded %d results from cached JSON", len(all_results))

    # 打印结果表
    print("\n" + "=" * 85)
    print("INTERFACE ANALYSIS SUMMARY (LOCPOT Poisson method)")
    print("=" * 85)
    print(f"  {'Het.':<22s}  {'Q_total (e⁻)':>15s}  {'Q_intf (e⁻)':>15s}  {'Dipole_V (eV·Å)':>17s}")
    print(f"  {'─'*22}  {'─'*15}  {'─'*15}  {'─'*17}")
    for het_name in HETERO_MAP:
        r = all_results.get(het_name, {})
        if "error" in r:
            print(f"  {SYSTEM_LABELS.get(het_name, het_name):<22s}  ERROR: {r['error']}")
        else:
            print(f"  {SYSTEM_LABELS.get(het_name, het_name):<22s}  "
                  f"{r['charge_transfer_total']:>+15.6f}  "
                  f"{r['charge_transfer_interface']:>+15.6f}  "
                  f"{r['dipole_V']:>+17.4f}")

    # 导出数据
    log.info("Exporting data...")
    export_summary_csv(all_results)
    export_profiles_csv(all_results)
    export_json(all_results)

    # 绘图
    if not args.no_plot:
        log.info("Generating plots...")
        for het_name in HETERO_MAP:
            r = all_results.get(het_name)
            if not r or "error" in r:
                continue
            plot_interface_profile(het_name, r)
        plot_comparison(all_results)

    log.info("Done. Output → %s", OUTPUT_DIR)


if __name__ == "__main__":
    main()
