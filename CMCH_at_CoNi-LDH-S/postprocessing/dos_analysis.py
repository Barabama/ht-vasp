"""
DOS/PDOS 分析 — 10 体系电子结构对比

功能:
  1. TDOS 对比图 (Bulk / Slab / Hetero 三组)
  2. PDOS 元素分辨 (Co-3d, Ni-3d, O-2p, S-3p)
  3. d-band center 汇总
  4. PDOS at E_F (费米能级处轨道贡献)
  5. Key contrast: 本征 vs S 掺杂带隙崩塌

用法:
  python dos_analysis.py                          # 全量分析
  python dos_analysis.py --systems CoNiOH2 CoNiOH2S-noH  # 指定体系
"""

import argparse
import json
import logging
import warnings
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from monty.json import MontyDecoder
from pymatgen.electronic_structure.dos import CompleteDos, Spin, OrbitalType

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.linewidth": 1.2})

DATA_DIR = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S/data")
OUTPUT_DIR = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S/postprocessing/output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

SYSTEMS_BULK = ["CoNiOH2", "CoNiOH2S-noH", "CoMnH2CO5"]
SYSTEMS_SLAB = ["CoNiOH2-slab", "CoNiOH2S-noH-slab", "CoNiOH2S-noH-slab-flip", "CoMnH2CO5-slab"]
SYSTEMS_HETERO = ["hetero_intrinsic", "hetero_s_doped", "hetero_s_exposed"]
SYSTEMS_ALL = SYSTEMS_BULK + SYSTEMS_SLAB + SYSTEMS_HETERO

SYSTEM_LABELS = {
    "CoNiOH2": "CoNiOH2 (Bulk)", "CoNiOH2S-noH": "CoNiOH2S (Bulk)", "CoMnH2CO5": "CoMnH2CO5 (Bulk)",
    "CoNiOH2-slab": "CoNiOH2 (Slab)", "CoNiOH2S-noH-slab": "CoNiOH2S (Slab)",
    "CoNiOH2S-noH-slab-flip": "CoNiOH2S (Slab-flip)", "CoMnH2CO5-slab": "CoMnH2CO5 (Slab)",
    "hetero_intrinsic": "Intrinsic Het.", "hetero_s_doped": "S-Doped Het.", "hetero_s_exposed": "S-Exposed Het.",
}

SYSTEM_COLORS = {
    "CoNiOH2": "#1f77b4", "CoNiOH2S-noH": "#2ca02c", "CoMnH2CO5": "#d62728",
    "CoNiOH2-slab": "#1f77b4", "CoNiOH2S-noH-slab": "#2ca02c",
    "CoNiOH2S-noH-slab-flip": "#ff7f0e", "CoMnH2CO5-slab": "#d62728",
    "hetero_intrinsic": "#1f77b4", "hetero_s_doped": "#2ca02c", "hetero_s_exposed": "#9467bd",
}

SYSTEM_GROUPS = {"Bulk": SYSTEMS_BULK, "Slab": SYSTEMS_SLAB, "Hetero": SYSTEMS_HETERO}

ELEM_PDOS_CONFIG = {
    "Co": ("d", "#2171b5", "#6baed6"),
    "Ni": ("d", "#238b45", "#74c476"),
    "O":  ("p", "#cb181d", "#fb6a4a"),
    "S":  ("p", "#d94801", "#fd8d3c"),
}


# ═══════════════════════════════════════════════
# Loading
# ═══════════════════════════════════════════════

def load_dos(name: str) -> CompleteDos | None:
    path = DATA_DIR / name / "dos_out.json"
    if not path.exists():
        return None
    try:
        with open(path) as f:
            raw = json.load(f)
        dos_dict = raw.get("vasp_objects", {}).get("dos")
        if dos_dict is None:
            return None
        return MontyDecoder().process_decoded(dos_dict)
    except Exception as e:
        log.warning("Failed to load DOS for %s: %s", name, e)
        return None


# ═══════════════════════════════════════════════
# Analysis Functions
# ═══════════════════════════════════════════════

def analyze_system(name: str) -> dict:
    """提取一个体系的 DOS 信息."""
    dos = load_dos(name)
    if dos is None:
        return {"name": name, "has_dos": False}

    result = {"name": name, "has_dos": True, "dos": dos}
    result["efermi"] = dos.efermi
    result["bandgap"] = dos.get_gap()

    # d-band center for Co/Ni
    for elem in ["Co", "Ni"]:
        try:
            result[f"d_center_{elem}"] = dos.get_band_center(
                band=OrbitalType.d, elements=[elem], erange=(-10, 5)
            )
        except Exception:
            result[f"d_center_{elem}"] = None

    # PDOS at E_F
    result["pdos_ef"] = _pdos_at_fermi(dos)

    # Spin polarization
    result["spin_pol"] = _spin_polarization(dos)

    # Site-projected S-3p
    result["has_S"] = any(s.species_string == "S" for s in dos.structure)

    return result


def _pdos_at_fermi(dos: CompleteDos, window: float = 0.3) -> dict:
    energies = dos.energies - dos.efermi
    mask = np.abs(energies) <= window
    pd = {"S_p": 0.0, "Co_d": 0.0, "Ni_d": 0.0, "total_up": 0.0, "total_down": 0.0}

    for site in dos.structure:
        sym = site.species_string
        if sym not in ("S", "Co", "Ni"):
            continue
        spd = dos.get_site_spd_dos(site)
        if sym == "S" and OrbitalType.p in spd:
            for spin in (Spin.up, Spin.down):
                if spin in spd[OrbitalType.p].densities:
                    pd["S_p"] += float(np.mean(spd[OrbitalType.p].densities[spin][mask]))
        elif sym in ("Co", "Ni"):
            key = f"{sym}_d"
            for spin in (Spin.up, Spin.down):
                if spin in spd[OrbitalType.d].densities:
                    pd[key] += float(np.mean(spd[OrbitalType.d].densities[spin][mask]))

    for spin in (Spin.up, Spin.down):
        if spin in dos.densities:
            pd[f"total_{'up' if spin == Spin.up else 'down'}"] = float(
                np.mean(dos.get_densities(spin)[mask])
            )
    return pd


def _spin_polarization(dos: CompleteDos, e_window=0.05, threshold=1e-3):
    try:
        energies = dos.energies - dos.efermi
        mask = np.abs(energies) <= e_window
        dos_up = dos.get_densities(Spin.up)[mask]
        if Spin.down not in dos.densities:
            return None
        dos_down = dos.get_densities(Spin.down)[mask]
        if np.mean(dos_up + dos_down) < threshold:
            return None
        n = dos.get_interpolated_value(dos.efermi)
        nu, nd = n.get(Spin.up, 0), n.get(Spin.down, 0)
        if nu + nd == 0:
            return None
        return abs((nu - nd) / (nu + nd))
    except Exception:
        return None


# ═══════════════════════════════════════════════
# Print Summary
# ═══════════════════════════════════════════════

def print_summary(results: dict):
    print("\n" + "=" * 120)
    print("DOS/PDOS SUMMARY")
    print("=" * 120)
    print(f"  {'System':<28s}  {'Gap':>8s}  {'E_Fermi':>10s}  {'d-ctr Co':>10s}  "
          f"{'d-ctr Ni':>10s}  {'Spin Pol':>10s}  {'S-3p@Ef':>10s}")
    print(f"  {'─'*28}  {'─'*8}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*10}")

    for name in SYSTEMS_ALL:
        r = results.get(name, {})
        if not r.get("has_dos"):
            continue
        gap_s = f"{r['bandgap']:.3f}" if r["bandgap"] else "—"
        ef_s = f"{r['efermi']:.3f}"
        dco = r.get("d_center_Co")
        dni = r.get("d_center_Ni")
        dco_s = f"{dco:.4f}" if dco else "—"
        dni_s = f"{dni:.4f}" if dni else "—"
        sp = r.get("spin_pol")
        sp_s = f"{sp:.3f}" if sp else "—"
        s3p = r.get("pdos_ef", {}).get("S_p", 0)
        s3p_s = f"{s3p:.4f}" if s3p else "—"
        print(f"  {SYSTEM_LABELS.get(name, name):<28s}  {gap_s:>8s}  {ef_s:>10s}  "
              f"{dco_s:>10s}  {dni_s:>10s}  {sp_s:>10s}  {s3p_s:>10s}")

    # Key contrast: bandgap collapse
    bg_p = results.get("CoNiOH2", {}).get("bandgap")
    bg_d = results.get("CoNiOH2S-noH", {}).get("bandgap")
    if bg_p and bg_d:
        print(f"\n  ★ S doping bandgap collapse (Bulk): {bg_p:.3f} → {bg_d:.3f} eV  "
              f"(Δ = {bg_d - bg_p:+.3f})")
    for name in SYSTEMS_HETERO:
        bg = results.get(name, {}).get("bandgap")
        if bg is not None:
            print(f"  ★ {SYSTEM_LABELS[name]} gap = {bg:.3f} eV")


# ═══════════════════════════════════════════════
# Plots
# ═══════════════════════════════════════════════

def plot_tdos_grouped(results: dict):
    """三组 TDOS 对比: Bulk / Slab / Hetero."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    e_range = (-6, 4)

    for ax_idx, (gname, systems) in enumerate(SYSTEM_GROUPS.items()):
        ax = axes[ax_idx]

        for name in systems:
            r = results.get(name)
            if not r or not r.get("dos"):
                continue
            dos = r["dos"]
            energies = dos.energies - dos.efermi
            mask = (energies >= e_range[0]) & (energies <= e_range[1])
            e_plot = energies[mask]
            up = dos.get_densities(Spin.up)[mask]
            down = dos.get_densities(Spin.down)[mask] if Spin.down in dos.densities else np.zeros_like(up)
            c = SYSTEM_COLORS[name]
            ax.plot(e_plot, up, color=c, lw=1.0, label=SYSTEM_LABELS[name])
            ax.plot(e_plot, -down, color=c, lw=0.7, alpha=0.5)
            ax.fill_between(e_plot, 0, up, color=c, alpha=0.12)
            gap = r.get("bandgap")
            if gap and gap > 0.05:
                ax.axvspan(0, gap, color="gray", alpha=0.08, lw=0)
                ax.axvline(0, color="gray", ls=":", lw=0.4, alpha=0.5)
                ax.axvline(gap, color="gray", ls=":", lw=0.4, alpha=0.5)

        ax.axvline(0, color="k", ls="--", lw=0.5)
        ax.set_xlim(e_range)
        ax.set_title(f"{gname} — TDOS", fontsize=12, fontweight="bold")
        ax.legend(fontsize=7, loc="upper right", framealpha=0.8)
        ax.grid(True, alpha=0.2)
        ax.set_xlabel("E − E$_F$ (eV)")

    axes[0].set_ylabel("DOS (states/eV)")
    fig.suptitle("Total DOS Comparison Across Layers", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUTPUT_DIR / "tdos_grouped.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "tdos_grouped.pdf", bbox_inches="tight")
    plt.close(fig)
    log.info("Plot → output/tdos_grouped.{png,pdf}")


def plot_pdos_elemental(results: dict):
    """元素分辨 PDOS (Co-3d, Ni-3d, O-2p, S-3p)."""
    focus = ["CoNiOH2", "CoNiOH2S-noH", "CoNiOH2-slab", "CoNiOH2S-noH-slab",
             "hetero_intrinsic", "hetero_s_doped"]

    fig, axes = plt.subplots(len(focus), 1, figsize=(12, 3 * len(focus) + 1), sharex=True)
    e_range = (-6, 4)

    for i, name in enumerate(focus):
        ax = axes[i]
        r = results.get(name)
        if not r or not r.get("dos"):
            ax.text(0.5, 0.5, "No DOS", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(SYSTEM_LABELS[name])
            continue

        dos = r["dos"]
        energies = dos.energies - dos.efermi
        mask = (energies >= e_range[0]) & (energies <= e_range[1])
        e_plot = energies[mask]

        for elem, (ot_str, c_up, c_down) in ELEM_PDOS_CONFIG.items():
            if elem == "S" and not r.get("has_S"):
                continue
            try:
                spd = dos.get_element_spd_dos(elem)
                ot = OrbitalType.d if ot_str == "d" else OrbitalType.p
                if ot not in spd:
                    continue
                scale = 5.0 if elem == "S" else 1.0
                for spin, is_up, color, ls in [
                    (Spin.up, True, c_up, "-"), (Spin.down, False, c_down, "--")
                ]:
                    if spin in spd[ot].densities:
                        d = spd[ot].densities[spin][mask] * scale
                        ax.plot(e_plot, d if is_up else -d, color=color, ls=ls, lw=0.9, alpha=0.85)
            except Exception:
                pass

        gap = r.get("bandgap")
        if gap and gap > 0.05:
            ax.axvspan(0, gap, color="gray", alpha=0.08, lw=0)
        ax.axvline(0, color="k", ls="--", lw=0.4)
        ax.set_xlim(e_range)
        gs = f" (gap={gap:.3f}eV)" if gap else ""
        ax.set_title(f"{SYSTEM_LABELS[name]}{gs}", fontsize=10, fontweight="bold")
        ax.set_ylabel("PDOS")
        ax.grid(True, alpha=0.2)

    axes[-1].set_xlabel("E − E$_F$ (eV)")

    legend = [
        Line2D([0], [0], color="#2171b5", lw=1.5, label="Co-3d (up)"),
        Line2D([0], [0], color="#2171b5", lw=1.5, ls="--", alpha=0.5, label="Co-3d (down)"),
        Line2D([0], [0], color="#238b45", lw=1.5, label="Ni-3d (up)"),
        Line2D([0], [0], color="#238b45", lw=1.5, ls="--", alpha=0.5, label="Ni-3d (down)"),
        Line2D([0], [0], color="#cb181d", lw=1.5, label="O-2p (up)"),
        Line2D([0], [0], color="#d94801", lw=1.5, label="S-3p (×5, up)"),
    ]
    fig.legend(handles=legend, fontsize=8, loc="upper right", framealpha=0.9, ncol=2)
    fig.suptitle("Element-Resolved PDOS", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUTPUT_DIR / "pdos_elemental.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "pdos_elemental.pdf", bbox_inches="tight")
    plt.close(fig)
    log.info("Plot → output/pdos_elemental.{png,pdf}")


def plot_dcenter_scatter(results: dict):
    """d-band center vs bandgap scatter."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax_idx, elem in enumerate(["Co", "Ni"]):
        ax = axes[ax_idx]
        xv, yv, lb, cl = [], [], [], []
        for name in SYSTEMS_ALL:
            r = results.get(name, {})
            gap = r.get("bandgap")
            dc = r.get(f"d_center_{elem}")
            if gap is not None and dc is not None:
                xv.append(gap)
                yv.append(dc)
                lb.append(SYSTEM_LABELS.get(name, name))
                cl.append(SYSTEM_COLORS.get(name, "#333"))
        ax.scatter(xv, yv, c=cl, s=50, edgecolors="k", linewidth=0.5, zorder=5)
        for xi, yi, nm in zip(xv, yv, lb):
            ax.annotate(nm, (xi, yi), textcoords="offset points", xytext=(4, 4),
                        fontsize=6, alpha=0.85,
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", alpha=0.7, ec="none"))
        ax.set_xlabel("Band Gap (eV)")
        ax.set_ylabel(f"{elem} d-band center (eV)")
        ax.set_title(f"{elem} d-center vs Band Gap", fontsize=11, fontweight="bold")
        ax.grid(True, alpha=0.2)
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "dcenter_vs_gap.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "dcenter_vs_gap.pdf", bbox_inches="tight")
    plt.close(fig)
    log.info("Plot → output/dcenter_vs_gap.{png,pdf}")


# ═══════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="DOS/PDOS analysis")
    parser.add_argument("--systems", type=str, nargs="+", help="体系名（默认全部）")
    parser.add_argument("--skip-plots", action="store_true", help="跳过绘图")
    args = parser.parse_args()

    systems = args.systems or SYSTEMS_ALL

    log.info("DOS/PDOS analysis for %d systems", len(systems))
    results = {}
    for name in systems:
        r = analyze_system(name)
        results[name] = r
        if r.get("has_dos"):
            log.info("  ✓ %s: gap=%.3f, d-Co=%.3f", name, r["bandgap"], r.get("d_center_Co", 0))
        else:
            log.warning("  ✗ %s: no DOS data", name)

    print_summary(results)

    if not args.skip_plots:
        log.info("Generating plots...")
        plot_tdos_grouped(results)
        plot_pdos_elemental(results)
        plot_dcenter_scatter(results)

    log.info("Done. Output → %s", OUTPUT_DIR)


if __name__ == "__main__":
    main()
