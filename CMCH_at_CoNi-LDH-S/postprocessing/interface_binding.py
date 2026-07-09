"""
界面结合能分析 — Level 2 应变 slab (共同晶格)

E_bind = E_hetero − E_cmch_slab − E_ldh_slab
E_bind_per_area = E_bind / (a × b)

用法:
  python postprocessing/interface_binding.py
  python postprocessing/interface_binding.py --force
  python postprocessing/interface_binding.py --no-plot

输出:
  output/interface_binding.png        柱状图
  output/interface_binding.csv        汇总表 (Origin/Excel 用)
  output/interface_binding.json       完整数值存档
"""

import argparse
import gzip
import json
import logging
import warnings
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
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
    "hetero_intrinsic": "Intrinsic Het.",
    "hetero_s_doped": "S-Doped Het.",
    "hetero_s_exposed": "S-Exposed Het.",
}

COLORS = {"hetero_intrinsic": "#1f77b4", "hetero_s_doped": "#2ca02c", "hetero_s_exposed": "#ff7f0e"}


# ═══════════════════════════════════════════════
# 数据读取
# ═══════════════════════════════════════════════

def read_energy(name: str) -> float | None:
    path = DATA_DIR / name / "static_out.json"
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())["output"]["energy"]
    except Exception as e:
        log.warning("Failed to read energy for %s: %s", name, e)
        return None


def get_area(name: str) -> float | None:
    path = DATA_DIR / name / "3-static" / "CONTCAR.gz"
    if not path.exists():
        return None
    try:
        with gzip.open(path, "rt") as f:
            struct = Structure.from_str(f.read(), fmt="poscar")
        return np.linalg.norm(np.cross(struct.lattice.matrix[0], struct.lattice.matrix[1]))
    except Exception as e:
        log.warning("Failed to get area for %s: %s", name, e)
        return None


def compute_binding() -> dict:
    results = {}
    for het_name, (ldh_name, cmch_name) in HETERO_MAP.items():
        e_het = read_energy(het_name)
        e_cmch = read_energy(cmch_name)
        e_ldh = read_energy(ldh_name)

        if e_het is None or e_cmch is None or e_ldh is None:
            log.warning("Skip %s: missing energy data", het_name)
            continue

        e_bind = e_het - e_cmch - e_ldh
        area = get_area(het_name)
        e_bind_area = e_bind / area if area else None

        results[het_name] = {
            "label": SYSTEM_LABELS.get(het_name, het_name),
            "E_hetero": float(e_het),
            "E_CMCH": float(e_cmch),
            "E_LDH": float(e_ldh),
            "E_bind": float(e_bind),
            "E_bind_per_area": float(e_bind_area) if e_bind_area else None,
            "area": float(area) if area else None,
        }
        area_str = f" ({e_bind_area:.6f} eV/Å²)" if e_bind_area else ""
        log.info("  ✓ %s: E_bind = %.4f eV%s", het_name, e_bind, area_str)
    return results


# ═══════════════════════════════════════════════
# 数据导出 (CSV / JSON)
# ═══════════════════════════════════════════════

def export_csv(results: dict):
    """汇总表 → output/interface_binding.csv"""
    rows = []
    for het_name in HETERO_MAP:
        r = results.get(het_name, {})
        if not r:
            continue
        rows.append({
            "system": het_name,
            "label": r.get("label", het_name),
            "E_hetero (eV)": r["E_hetero"],
            "E_CMCH (eV)": r["E_CMCH"],
            "E_LDH (eV)": r["E_LDH"],
            "E_bind (eV)": r["E_bind"],
            "E_bind_per_area (eV/A^2)": r["E_bind_per_area"] if r["E_bind_per_area"] is not None else "",
            "area (A^2)": r["area"] if r["area"] is not None else "",
        })

    if not rows:
        return
    header = list(rows[0].keys())
    with open(OUTPUT_DIR / "interface_binding.csv", "w") as f:
        f.write(",".join(header) + "\n")
        for row in rows:
            f.write(",".join(str(row[k]) for k in header) + "\n")
    log.info("Data  → output/interface_binding.csv")


def export_json(results: dict):
    with open(OUTPUT_DIR / "interface_binding.json", "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    log.info("Data  → output/interface_binding.json")


def load_results() -> dict | None:
    path = OUTPUT_DIR / "interface_binding.json"
    if not path.exists():
        return None
    try:
        with open(path) as f:
            data = json.load(f)
        log.info("Loaded from %s", path.name)
        return data
    except Exception as e:
        log.warning("Cannot load %s: %s", path.name, e)
        return None


# ═══════════════════════════════════════════════
# Plot
# ═══════════════════════════════════════════════

def plot_bar(results: dict):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    hets = [h for h in HETERO_MAP if h in results]
    labels = [SYSTEM_LABELS.get(h, h) for h in hets]
    bind_vals = [results[h].get("E_bind") for h in hets]
    area_vals = [results[h].get("E_bind_per_area") for h in hets]
    colors = [COLORS.get(h, "#333") for h in hets]
    x = np.arange(len(hets))
    width = 0.5

    ax1.bar(x, bind_vals, width, color=colors, alpha=0.85, edgecolor="k", lw=0.5)
    for i, v in enumerate(bind_vals):
        ax1.text(i, v + (0.05 * abs(max(bind_vals))),
                 f"{v:.4f}", ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, fontsize=9, rotation=15, ha="right")
    ax1.set_ylabel("E_bind / eV")
    ax1.set_title("Interface Binding Energy", fontweight="bold")
    ax1.axhline(0, color="k", ls="--", lw=0.5, alpha=0.4)
    ax1.grid(True, alpha=0.2, axis="y")

    ax2.bar(x, area_vals, width, color=colors, alpha=0.85, edgecolor="k", lw=0.5)
    for i, v in enumerate(area_vals):
        if v is not None:
            ax2.text(i, v + (0.05 * abs(max(area_vals))),
                     f"{v:.6f}", ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, fontsize=9, rotation=15, ha="right")
    ax2.set_ylabel("E_bind / area / (eV/Å²)")
    ax2.set_title("Binding Energy per Area", fontweight="bold")
    ax2.axhline(0, color="k", ls="--", lw=0.5, alpha=0.4)
    ax2.grid(True, alpha=0.2, axis="y")

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "interface_binding.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("Plot → output/interface_binding.png")


# ═══════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Interface binding energy analysis")
    parser.add_argument("--no-plot", action="store_true", help="仅导出数据，不生成图片")
    parser.add_argument("--force", action="store_true", help="强制重新计算（忽略已有 JSON）")
    args = parser.parse_args()

    results = None if args.force else load_results()

    if results is None:
        log.info("Computing interface binding energies ...")
        results = compute_binding()
    else:
        log.info("Loaded %d heterojunction results from JSON", len(results))

    # Print table
    print("\n" + "=" * 90)
    print("INTERFACE BINDING ENERGY (Level 2, common lattice)")
    print("=" * 90)
    print(f"  {'Het.':<22s}  {'E_het':>12s}  {'E_CMCH':>12s}  {'E_LDH':>12s}  "
          f"{'E_bind':>12s}  {'E_bind/area':>12s}")
    print(f"  {'─'*22}  {'─'*12}  {'─'*12}  {'─'*12}  {'─'*12}  {'─'*12}")
    for het_name in HETERO_MAP:
        r = results.get(het_name, {})
        if not r:
            continue
        eb = r.get("E_bind")
        eba = r.get("E_bind_per_area")
        print(f"  {SYSTEM_LABELS.get(het_name, het_name):<22s}  "
              f"{r['E_hetero']:>12.4f}  {r['E_CMCH']:>12.4f}  {r['E_LDH']:>12.4f}  "
              f"{f'{eb:.4f}' if eb else '—':>12s}  {f'{eba:.4f}' if eba else '—':>12s}")

    # Data export
    log.info("Exporting data...")
    export_csv(results)
    export_json(results)

    # Plots
    if not args.no_plot and results:
        plot_bar(results)

    log.info("Done. Output → %s", OUTPUT_DIR)


if __name__ == "__main__":
    main()
