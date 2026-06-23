"""
功函数分析 — 7 体系 (4 slab + 3 hetero)

从 LOCPOT_dipole.gz 提取 planar average 静电势，计算 Φ = V_vacuum - E_Fermi

用法:
  python work_function_analysis.py                         # 全部分析
  python work_function_analysis.py --systems CoNiOH2-slab  # 指定体系

输出:
  终端: 功函数对比表
  output/wf_profile.png / .pdf  (静电势曲线叠加图)
  output/wf_comparison.png / .pdf  (柱状图)
"""

import argparse
import gzip
import logging
import os
import shutil
import tempfile
import warnings
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pymatgen.io.vasp.outputs import Locpot, Outcar
from pymatgen.analysis.surface_analysis import WorkFunctionAnalyzer
from pymatgen.core import Structure

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.linewidth": 1.2, "figure.dpi": 100})

DATA_DIR = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S/data")
OUTPUT_DIR = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S/postprocessing/output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SYSTEMS_WF = [
    "CoNiOH2-slab", "CoNiOH2S-noH-slab", "CoNiOH2S-noH-slab-flip", "CoMnH2CO5-slab",
    "hetero_intrinsic", "hetero_s_doped", "hetero_s_exposed",
]

SYSTEM_LABELS = {
    "CoNiOH2-slab": "CoNiOH2 (Slab)",
    "CoNiOH2S-noH-slab": "CoNiOH2S (Slab)",
    "CoNiOH2S-noH-slab-flip": "CoNiOH2S (Slab-flip)",
    "CoMnH2CO5-slab": "CoMnH2CO5 (Slab)",
    "hetero_intrinsic": "Intrinsic Het.",
    "hetero_s_doped": "S-Doped Het.",
    "hetero_s_exposed": "S-Exposed Het.",
}

SYSTEM_COLORS = {
    "CoNiOH2-slab": "#1f77b4",
    "CoNiOH2S-noH-slab": "#2ca02c",
    "CoNiOH2S-noH-slab-flip": "#ff7f0e",
    "CoMnH2CO5-slab": "#d62728",
    "hetero_intrinsic": "#1f77b4",
    "hetero_s_doped": "#2ca02c",
    "hetero_s_exposed": "#9467bd",
}


# ═══════════════════════════════════════════════
# Work Function Extraction
# ═══════════════════════════════════════════════

def extract_wf(name: str) -> dict:
    """从 LOCPOT_dipole 提取功函数."""
    static_dir = DATA_DIR / name / "3-static"

    # LOCPOT: 优先偶极修正版
    locpot = static_dir / "LOCPOT_dipole.gz"
    if not locpot.exists():
        locpot = static_dir / "LOCPOT.gz"
    if not locpot.exists():
        return {"name": name, "error": "No LOCPOT found"}

    # OUTCAR: 对应版本
    outcar_dipole = static_dir / "OUTCAR_dipole.gz"
    outcar = outcar_dipole if outcar_dipole.exists() else static_dir / "OUTCAR.gz"
    if not outcar.exists():
        return {"name": name, "error": "No OUTCAR found"}

    # Structure
    contcar = static_dir / "CONTCAR.gz"
    if not contcar.exists():
        return {"name": name, "error": "No CONTCAR"}
    with gzip.open(contcar, "rt") as f:
        struct = Structure.from_str(f.read(), fmt="poscar")

    # Decompress LOCPOT
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

    except Exception as e:
        return {"name": name, "error": f"LOCPOT read: {e}"}

    # E_Fermi
    try:
        o = Outcar(outcar)
        efermi = o.efermi
    except Exception as e:
        return {"name": name, "error": f"OUTCAR read: {e}"}

    if efermi is None:
        return {"name": name, "error": "No E_Fermi"}

    # WorkFunctionAnalyzer
    try:
        wfa = WorkFunctionAnalyzer(
            structure=struct,
            locpot_along_c=avg,
            efermi=efermi,
        )
        wf = wfa.work_function
        v_vac = wfa.vacuum_locpot
    except Exception as e:
        # fallback: smart plateau detection
        c = struct.lattice.c
        grid = np.linspace(0, c, ngz)
        mat_top = max(s.frac_coords[2] for s in struct) * c
        vac_start = np.searchsorted(grid, mat_top + 3.0)
        vac_window = min(40, ngz - vac_start)
        if vac_window > 5:
            min_std = float("inf")
            best = 0
            for i in range(vac_start, ngz - vac_window):
                s = avg[i:i+vac_window].std()
                if s < min_std:
                    min_std = s
                    best = i
            v_vac = float(avg[best:best+vac_window].mean())
            wf = v_vac - efermi
        else:
            v_vac = float(avg[-int(ngz*0.15):].mean())
            wf = v_vac - efermi

    return {
        "name": name,
        "label": SYSTEM_LABELS.get(name, name),
        "work_function": float(wf),
        "vacuum_level": float(v_vac),
        "efermi": float(efermi),
        "locpot_avg": avg.tolist(),
        "ngz": ngz,
        "has_dipole_correction": "dipole" in str(locpot),
    }


# ═══════════════════════════════════════════════
# Plot: electrostatic potential profiles
# ═══════════════════════════════════════════════

def plot_wf_profiles(results: dict):
    """静电势曲线叠加 (含真空水平标注)."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    groups = [("Slab", ["CoNiOH2-slab", "CoNiOH2S-noH-slab", "CoNiOH2S-noH-slab-flip", "CoMnH2CO5-slab"]),
              ("Hetero", ["hetero_intrinsic", "hetero_s_doped", "hetero_s_exposed"])]

    for ax, (gname, systems) in zip(axes, groups):
        for name in systems:
            r = results.get(name)
            if not r or "error" in r or "locpot_avg" not in r:
                continue
            struct = _get_structure(name)
            if struct is None:
                continue
            c = struct.lattice.c
            grid = np.linspace(0, c, r["ngz"])
            ax.plot(grid, r["locpot_avg"], lw=1.0, color=SYSTEM_COLORS[name], label=SYSTEM_LABELS[name])
            # vacuum level line
            ax.axhline(r["vacuum_level"], color=SYSTEM_COLORS[name], ls=":", lw=0.7, alpha=0.5)

        ax.axhline(0, color="gray", ls="--", lw=0.5, alpha=0.3)
        ax.set_xlabel("z (Å)")
        ax.set_ylabel("V(z) (eV)")
        ax.set_title(f"Electrostatic Potential — {gname}", fontweight="bold")
        ax.legend(fontsize=7, framealpha=0.8)
        ax.grid(True, alpha=0.2)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "wf_profile.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "wf_profile.pdf", bbox_inches="tight")
    plt.close(fig)
    log.info("Plot → output/wf_profile.{png,pdf}")


def _get_structure(name: str) -> Structure | None:
    static = DATA_DIR / name / "3-static"
    contcar = static / "CONTCAR.gz"
    if not contcar.exists():
        return None
    try:
        with gzip.open(contcar, "rt") as f:
            return Structure.from_str(f.read(), fmt="poscar")
    except Exception:
        return None


# ═══════════════════════════════════════════════
# Plot: bar chart comparison
# ═══════════════════════════════════════════════

def plot_wf_comparison(results: dict):
    """功函数柱状图对比."""
    names = [n for n in SYSTEMS_WF if n in results and "error" not in results[n]]
    vals = [results[n]["work_function"] for n in names]
    labels = [SYSTEM_LABELS[n] for n in names]
    colors = [SYSTEM_COLORS[n] for n in names]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    bars = ax.bar(range(len(names)), vals, color=colors, alpha=0.85, edgecolor="k", linewidth=0.5)
    for i, (bar, val) in enumerate(zip(bars, vals)):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.08,
                f"{val:.3f}", ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=9)
    ax.set_ylabel("Work Function Φ (eV)")
    ax.set_title("Work Function Comparison", fontweight="bold")
    ax.grid(True, alpha=0.2, axis="y")
    ax.set_ylim(0, max(vals) * 1.25)

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "wf_comparison.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "wf_comparison.pdf", bbox_inches="tight")
    plt.close(fig)
    log.info("Plot → output/wf_comparison.{png,pdf}")


# ═══════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Work function analysis")
    parser.add_argument("--systems", type=str, nargs="+", help="体系名")
    args = parser.parse_args()

    systems = args.systems or SYSTEMS_WF

    log.info("Work function analysis for %d systems", len(systems))
    results = {}
    for name in systems:
        r = extract_wf(name)
        results[name] = r
        if "error" in r:
            log.warning("  ✗ %s: %s", name, r["error"])
        else:
            dc = " (dipole)" if r.get("has_dipole_correction") else ""
            log.info("  ✓ %s: Φ = %.4f eV%s", name, r["work_function"], dc)

    # Print table
    print("\n" + "=" * 80)
    print("WORK FUNCTION SUMMARY")
    print("=" * 80)
    print(f"  {'System':<28s}  {'Φ (eV)':>10s}  {'V_vac':>10s}  {'E_Fermi':>10s}  {'Correction':>12s}")
    print(f"  {'─'*28}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*12}")
    for name in systems:
        r = results.get(name, {})
        if "error" in r:
            print(f"  {SYSTEM_LABELS.get(name, name):<28s}  ERROR: {r['error']}")
        else:
            dc = "dipole" if r.get("has_dipole_correction") else "standard"
            print(f"  {SYSTEM_LABELS.get(name, name):<28s}  {r['work_function']:>10.4f}  "
                  f"{r['vacuum_level']:>10.4f}  {r['efermi']:>10.4f}  {dc:>12s}")

    # Key comparisons
    print("\n── Key Contrasts ──")
    names_slab = ["CoNiOH2-slab", "CoNiOH2S-noH-slab", "CoNiOH2S-noH-slab-flip"]
    vals = [results[n].get("work_function") for n in names_slab]
    labels = ["Pristine LDH", "S-doped LDH", "S-doped LDH (flip)"]
    if all(v is not None for v in vals):
        print(f"  LDH: {labels[0]}={vals[0]:.3f}, {labels[1]}={vals[1]:.3f} (Δ={vals[1]-vals[0]:+.3f})")
        print(f"  LDH surface S effect: flip={vals[2]:.3f} vs internal={vals[1]:.3f}")

    names_het = ["hetero_intrinsic", "hetero_s_doped", "hetero_s_exposed"]
    vh = [results[n].get("work_function") for n in names_het]
    if all(v is not None for v in vh):
        print(f"  Hetero: intrinsic={vh[0]:.3f}, s_doped={vh[1]:.3f}, s_exposed={vh[2]:.3f}")
        print(f"  S doping effect: Δ_intr_doped={vh[1]-vh[0]:+.3f}, Δ_intr_exposed={vh[2]-vh[0]:+.3f}")

    # Plots
    log.info("Generating plots...")
    plot_wf_profiles(results)
    plot_wf_comparison(results)
    log.info("Done. Output → %s", OUTPUT_DIR)


if __name__ == "__main__":
    main()
