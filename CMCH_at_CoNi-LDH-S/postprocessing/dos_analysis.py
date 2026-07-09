"""
DOS/PDOS 分析 — 10 体系电子结构对比 (Level 2)

功能:
  1. TDOS 对比图 (Bulk / Slab / Hetero 三组)
  2. PDOS 元素分辨 (Co-3d, Ni-3d, O-2p, S-3p)
  3. d-band center 汇总
  4. PDOS at E_F (费米能级处轨道贡献)
  5. Key contrast: 本征 vs S 掺杂带隙崩塌

体系:
  Bulk (3):  LDH_bulk, LDH_S_bulk, CMCH_bulk
  Slab (8):  LDH_slab, LDH_S_slab, LDH_S_flip_slab, CMCH_slab,
             CMCH_strained, LDH_strained, LDH_S_strained, LDH_S_flip_strained
  Hetero (3): hetero_intrinsic, hetero_s_doped, hetero_s_exposed

用法:
  python dos_analysis.py                          # 全量分析 (Level 2)
  python dos_analysis.py --systems LDH_strained   # 指定体系
"""

import argparse
import json
import logging
import os
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

plt.rcParams.update({"font.family": "sans-serif", "font.size": 11, "axes.linewidth": 1.2})

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parent / "data"
OUTPUT_DIR = SCRIPT_DIR / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 系统命名 → 数据目录映射（解耦命名与目录结构）
DATA_DIR_MAP = {
    "LDH_bulk": "CoNiOH2",
    "LDH_S_bulk": "CoNiOH2S-noH",
    "CMCH_bulk": "CMCH_bulk",
    "LDH_slab": "CoNiOH2-slab",
    "LDH_S_slab": "CoNiOH2S-noH-slab",
    "LDH_S_flip_slab": "CoNiOH2S-noH-slab-flip",
    "CMCH_slab": "CoMnH2CO5-slab",
    "CMCH_strained": "CMCH_strained",
    "LDH_strained": "LDH_strained",
    "LDH_S_strained": "LDH_S_strained",
    "LDH_S_flip_strained": "LDH_S_flip_strained",
    "hetero_intrinsic": "hetero_intrinsic",
    "hetero_s_doped": "hetero_s_doped",
    "hetero_s_exposed": "hetero_s_exposed",
}

SYSTEMS_BULK = ["LDH_bulk", "LDH_S_bulk", "CMCH_bulk"]
SYSTEMS_SLAB = ["LDH_slab", "LDH_S_slab", "LDH_S_flip_slab", "CMCH_slab",
                 "CMCH_strained", "LDH_strained", "LDH_S_strained", "LDH_S_flip_strained"]
SYSTEMS_HETERO = ["hetero_intrinsic", "hetero_s_doped", "hetero_s_exposed"]
SYSTEMS_ALL = SYSTEMS_BULK + SYSTEMS_SLAB + SYSTEMS_HETERO

SYSTEM_LABELS = {
    "LDH_bulk": "LDH (Bulk)", "LDH_S_bulk": "LDH+S (Bulk)", "CMCH_bulk": "CMCH (Bulk)",
    "LDH_slab": "LDH (Slab)", "LDH_S_slab": "LDH+S (Slab)",
    "LDH_S_flip_slab": "LDH+S (Slab-flip)", "CMCH_slab": "CMCH (Slab)",
    "CMCH_strained": "CMCH (Strained)", "LDH_strained": "LDH (Strained)",
    "LDH_S_strained": "LDH+S (Strained)", "LDH_S_flip_strained": "LDH+S (Strained-flip)",
    "hetero_intrinsic": "Intrinsic Het.", "hetero_s_doped": "S-Doped Het.", "hetero_s_exposed": "S-Exposed Het.",
}

SYSTEM_COLORS = {
    "LDH_bulk": "#1f77b4", "LDH_S_bulk": "#2ca02c", "CMCH_bulk": "#1f77b4",
    "LDH_slab": "#1f77b4", "LDH_S_slab": "#2ca02c",
    "LDH_S_flip_slab": "#ff7f0e", "CMCH_slab": "#d62728",
    "CMCH_strained": "#1f77b4", "LDH_strained": "#2ca02c",
    "LDH_S_strained": "#d62728", "LDH_S_flip_strained": "#ff7f0e",
    "hetero_intrinsic": "#333333", "hetero_s_doped": "#e67e22", "hetero_s_exposed": "#9467bd",
}

SYSTEM_GROUPS = {"Bulk": SYSTEMS_BULK, "Slab": SYSTEMS_SLAB, "Hetero": SYSTEMS_HETERO}

ELEM_PDOS_CONFIG = {
    "Co": ("d", "#2171b5", "#6baed6"),
    "Ni": ("d", "#238b45", "#74c476"),
    "Mn": ("d", "#6a51a3", "#9e9ac8"),
    "O":  ("p", "#cb181d", "#fb6a4a"),
    "S":  ("p", "#d94801", "#fd8d3c"),
}


# ═══════════════════════════════════════════════
# Loading
# ═══════════════════════════════════════════════

def _dir(name: str) -> Path:
    """返回系统名对应的数据目录（通过 DATA_DIR_MAP 解耦命名与目录路径）。"""
    return DATA_DIR / DATA_DIR_MAP.get(name, name)


def load_dos(name: str) -> CompleteDos | None:
    path = _dir(name) / "dos_out.json"
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

    # Spin-resolved band gaps
    try:
        up_dos = dos.get_densities(Spin.up)
        down_dos = dos.get_densities(Spin.down) if Spin.down in dos.densities else None
    except Exception:
        up_dos = down_dos = None

    result["gap_up"] = None
    result["gap_down"] = None
    if up_dos is not None:
        try:
            result["gap_up"] = dos.get_gap(spin=Spin.up)
        except Exception:
            pass
    if down_dos is not None:
        try:
            result["gap_down"] = dos.get_gap(spin=Spin.down)
        except Exception:
            pass

    # d-band center for Co, Ni, Mn
    for elem in ["Co", "Ni", "Mn"]:
        try:
            result[f"d_center_{elem}"] = dos.get_band_center(
                band=OrbitalType.d, elements=[elem], erange=(-10, 5)
            )
        except Exception:
            result[f"d_center_{elem}"] = None

    # Total magnetization (per formula unit estimate)
    try:
        spin_up = np.trapz(dos.get_densities(Spin.up), dos.energies)
        if Spin.down in dos.densities:
            spin_down = np.trapz(dos.get_densities(Spin.down), dos.energies)
            result["mag_moment"] = float(spin_up - spin_down)
        else:
            result["mag_moment"] = None
    except Exception:
        result["mag_moment"] = None

    # PDOS at E_F
    result["pdos_ef"] = _pdos_at_fermi(dos)

    # Spin polarization
    result["spin_pol"] = _spin_polarization(dos)

    # Site-projected S-3p
    result["has_S"] = any(s.species_string == "S" for s in dos.structure)
    result["has_Mn"] = any(s.species_string == "Mn" for s in dos.structure)

    return result


def _pdos_at_fermi(dos: CompleteDos, window: float = 0.3) -> dict:
    energies = dos.energies - dos.efermi
    mask = np.abs(energies) <= window
    pd = {"S_p": 0.0, "Co_d": 0.0, "Ni_d": 0.0, "Mn_d": 0.0, "total_up": 0.0, "total_down": 0.0}

    for site in dos.structure:
        sym = site.species_string
        if sym not in ("S", "Co", "Ni", "Mn"):
            continue
        spd = dos.get_site_spd_dos(site)
        if sym == "S" and OrbitalType.p in spd:
            for spin in (Spin.up, Spin.down):
                if spin in spd[OrbitalType.p].densities:
                    pd["S_p"] += float(np.mean(spd[OrbitalType.p].densities[spin][mask]))
        elif sym in ("Co", "Ni", "Mn"):
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
    print("\n" + "=" * 140)
    print("DOS/PDOS SUMMARY")
    print("=" * 140)
    print(f"  {'System':<25s}  {'Gap':>7s}  {'Gap↑':>7s}  {'Gap↓':>7s}  {'E_F':>8s}  {'d-Co':>8s}  "
          f"{'d-Ni':>8s}  {'d-Mn':>8s}  {'Mag':>7s}  {'SpinPol':>8s}  {'S@Ef':>7s}")
    print(f"  {'─'*25}  {'─'*7}  {'─'*7}  {'─'*7}  {'─'*8}  {'─'*8}  {'─'*8}  {'─'*8}  {'─'*7}  {'─'*8}  {'─'*7}")

    for name in SYSTEMS_ALL:
        r = results.get(name, {})
        if not r.get("has_dos"):
            continue
        def _f(v, fmt=".3f", fallback="—"):
            return f"{v:{fmt}}" if v is not None else fallback
        gap_s     = _f(r.get("bandgap"))
        gap_up_s  = _f(r.get("gap_up"))
        gap_dn_s  = _f(r.get("gap_down"))
        ef_s      = _f(r.get("efermi"), ".3f")
        dco_s     = _f(r.get("d_center_Co"), ".4f")
        dni_s     = _f(r.get("d_center_Ni"), ".4f")
        dmn_s     = _f(r.get("d_center_Mn"), ".4f")
        mag_s     = _f(r.get("mag_moment"), ".2f")
        sp_s      = _f(r.get("spin_pol"), ".3f")
        s3p       = r.get("pdos_ef", {}).get("S_p", 0)
        s3p_s     = f"{s3p:.3f}" if s3p else "—"
        print(f"  {SYSTEM_LABELS.get(name, name):<25s}  {gap_s:>7s}  {gap_up_s:>7s}  {gap_dn_s:>7s}  "
              f"{ef_s:>8s}  {dco_s:>8s}  {dni_s:>8s}  {dmn_s:>8s}  {mag_s:>7s}  {sp_s:>8s}  {s3p_s:>7s}")

    # Key contrast: bandgap collapse
    bg_p = results.get("LDH_slab", {}).get("bandgap")
    for n in SYSTEMS_SLAB:
        if "S" in n:
            bg_d = results.get(n, {}).get("bandgap")
            break
    if bg_p and bg_d is not None:
        print(f"\n  ★ S doping bandgap collapse (Slab): {bg_p:.3f} → {bg_d:.3f} eV  "
              f"(Δ = {bg_d - bg_p:+.3f})")


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
    fig.savefig(OUTPUT_DIR / "tdos_grouped.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("Plot → output/tdos_grouped.png")


def plot_pdos_elemental(results: dict):
    """元素分辨 PDOS (Co-3d, Ni-3d, O-2p, S-3p)."""
    focus = ["LDH_strained", "LDH_S_strained",
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
        Line2D([0], [0], color="#6a51a3", lw=1.5, label="Mn-3d (up)"),
        Line2D([0], [0], color="#6a51a3", lw=1.5, ls="--", alpha=0.5, label="Mn-3d (down)"),
        Line2D([0], [0], color="#cb181d", lw=1.5, label="O-2p (up)"),
        Line2D([0], [0], color="#d94801", lw=1.5, label="S-3p (×5, up)"),
    ]
    fig.legend(handles=legend, fontsize=8, loc="upper right", framealpha=0.9, ncol=2)
    fig.suptitle("Element-Resolved PDOS", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUTPUT_DIR / "pdos_elemental.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("Plot → output/pdos_elemental.png")


def plot_dcenter_scatter(results: dict):
    """d-band center vs bandgap scatter for Co, Ni, Mn."""
    elems = [e for e in ["Co", "Ni", "Mn"] if any(
        results.get(n, {}).get(f"d_center_{e}") is not None for n in SYSTEMS_ALL
    )]
    fig, axes = plt.subplots(1, len(elems), figsize=(5 * len(elems), 4))
    if len(elems) == 1:
        axes = [axes]
    for ax_idx, elem in enumerate(elems):
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
    fig.savefig(OUTPUT_DIR / "dcenter_vs_gap.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("Plot → output/dcenter_vs_gap.png")


# ═══════════════════════════════════════════════
# Export / Load Utilities
# ═══════════════════════════════════════════════

def export_summary_csv(results: dict):
    """Export DOS summary to CSV."""
    rows = []
    for name in SYSTEMS_ALL:
        r = results.get(name, {})
        if not r.get("has_dos"):
            continue
        rows.append({
            "system": name,
            "label": SYSTEM_LABELS.get(name, name),
            "bandgap": r.get("bandgap"),
            "gap_up": r.get("gap_up"),
            "gap_down": r.get("gap_down"),
            "efermi": r.get("efermi"),
            "mag_moment": r.get("mag_moment"),
            "d_center_Co": r.get("d_center_Co"),
            "d_center_Ni": r.get("d_center_Ni"),
            "d_center_Mn": r.get("d_center_Mn"),
            "spin_pol": r.get("spin_pol"),
            "S_p_at_Ef": r.get("pdos_ef", {}).get("S_p"),
            "Co_d_at_Ef": r.get("pdos_ef", {}).get("Co_d"),
            "Ni_d_at_Ef": r.get("pdos_ef", {}).get("Ni_d"),
            "Mn_d_at_Ef": r.get("pdos_ef", {}).get("Mn_d"),
        })
    import csv
    path = OUTPUT_DIR / "dos_summary.csv"
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader()
        w.writerows(rows)
    log.info("CSV → output/dos_summary.csv")


def load_results() -> dict | None:
    path = OUTPUT_DIR / "dos.json"
    if not path.exists():
        return None
    try:
        with open(path) as f:
            return json.load(f)
    except Exception:
        return None


def load_dos_object(name: str) -> CompleteDos | None:
    path = _dir(name) / "dos_out.json"
    if not path.exists():
        return None
    try:
        with open(path) as f:
            raw = json.load(f)
        dos_dict = raw.get("vasp_objects", {}).get("dos")
        if dos_dict is None:
            return None
        return MontyDecoder().process_decoded(dos_dict)
    except Exception:
        return None


def export_dos_csv(results: dict):
    export_summary_csv(results)


def export_results_json(results: dict):
    clean = {}
    for name, r in results.items():
        clean[name] = {k: v for k, v in r.items() if k != "dos"}
    with open(OUTPUT_DIR / "dos.json", "w") as f:
        json.dump(clean, f, indent=2, ensure_ascii=False)
    log.info("Data → output/dos.json")


def export_dos_curves(results: dict):
    """导出 DOS 曲线逐点数据 → output/dos_curves/dos_{name}.csv

    每体系输出: energy, tdos_up, tdos_down, Co_d, Ni_d, Mn_d, O_p, S_p
    """
    out_dir = OUTPUT_DIR / "dos_curves"
    out_dir.mkdir(parents=True, exist_ok=True)
    count = 0
    for name in SYSTEMS_ALL:
        r = results.get(name, {})
        if not r.get("has_dos"):
            continue
        dos = r.get("dos")  # could be None if loaded from JSON
        if dos is None:
            dos = load_dos_object(name)
        if dos is None:
            continue

        energies = dos.energies - dos.efermi
        tdos_up = dos.get_densities(Spin.up)
        tdos_down = dos.get_densities(Spin.down) if Spin.down in dos.densities else np.zeros_like(tdos_up)

        # Element/orbital resolved PDOS
        pdos_cols = {}
        for elem, (ot_str, _, _) in ELEM_PDOS_CONFIG.items():
            ot = OrbitalType.d if ot_str == "d" else OrbitalType.p
            try:
                spd = dos.get_element_spd_dos(elem)
                if ot in spd:
                    d_up = spd[ot].densities.get(Spin.up, np.zeros_like(energies))
                    d_dn = spd[ot].densities.get(Spin.down, np.zeros_like(energies))
                    pdos_cols[f"{elem}_{ot_str}"] = (d_up + d_dn) / 2  # average spin
                else:
                    pdos_cols[f"{elem}_{ot_str}"] = np.zeros_like(energies)
            except Exception:
                pdos_cols[f"{elem}_{ot_str}"] = np.zeros_like(energies)

        header = ["energy_eV", "tdos_up", "tdos_down"]
        col_data = [energies, tdos_up, tdos_down]
        for elem_label in ["Co_d", "Ni_d", "Mn_d", "O_p", "S_p"]:
            if elem_label in pdos_cols:
                header.append(elem_label)
                col_data.append(pdos_cols[elem_label])

        data = np.column_stack(col_data)
        fpath = out_dir / f"dos_{name}.csv"
        np.savetxt(fpath, data, delimiter=",", header=",".join(header), comments="")
        count += 1

    log.info("DOS curves → output/dos_curves/ (%d systems)", count)


# ═══════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="DOS/PDOS analysis")
    parser.add_argument("--systems", type=str, nargs="+", help="体系名（默认全部）")
    parser.add_argument("--no-plot", action="store_true", help="仅导出数据，不生成图片")
    parser.add_argument("--force", action="store_true", help="强制重新计算（忽略已有 JSON）")
    args = parser.parse_args()

    systems = args.systems or SYSTEMS_ALL

    # 优先从 JSON 恢复（标量指标）
    results = None if args.force else load_results()

    if results is None:
        log.info("DOS/PDOS analysis for %d systems", len(systems))
        results = {}
        for name in systems:
            r = analyze_system(name)
            results[name] = r
            if r.get("has_dos"):
                log.info("  ✓ %s: gap=%.3f, d-Co=%.3f", name, r["bandgap"], r.get("d_center_Co", 0))
            else:
                log.warning("  ✗ %s: no DOS data", name)
    else:
        log.info("Loaded %d system results from JSON", len(results))
        # 绘图时需补充 CompleteDos 对象（JSON 中未存储）
        if not args.no_plot:
            for name, r in results.items():
                if r.get("has_dos") and "dos" not in r:
                    dos = load_dos_object(name)
                    if dos is not None:
                        r["dos"] = dos

    print_summary(results)

    # Data export
    log.info("Exporting data...")
    export_dos_csv(results)
    export_results_json(results)
    export_dos_curves(results)

    if not args.no_plot:
        log.info("Generating plots...")
        plot_tdos_grouped(results)
        plot_pdos_elemental(results)
        plot_dcenter_scatter(results)

    log.info("Done. Output → %s", OUTPUT_DIR)


if __name__ == "__main__":
    main()
