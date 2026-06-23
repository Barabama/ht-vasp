#!/usr/bin/env python3
"""
CMCH@CoNi-LDH-S 异质结多体系综合分析

10 体系跨层对比：
  Bulk(3): CoNiOH2, CoNiOH2S-noH, CoMnH2CO5
  Slab(4): CoNiOH2-slab, CoNiOH2S-noH-slab, CoNiOH2S-noH-slab-flip, CoMnH2CO5-slab
  Hetero(3): hetero_intrinsic, hetero_s_doped, hetero_s_exposed
"""

import json
import logging
import os
import gzip
import shutil
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
from monty.json import MontyDecoder
from pymatgen.core import Structure
from pymatgen.electronic_structure.dos import CompleteDos, Spin, OrbitalType
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.electronic_structure.plotter import BSPlotter
from pymatgen.io.vasp.outputs import Outcar, Locpot

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.linewidth": 1.2, "figure.dpi": 100})

# ═══════════════════════════════════════════════
# CONFIG
# ═══════════════════════════════════════════════

DATA_DIR = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S/data")
OUTPUT_DIR = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S/output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ZVAL for charge transfer reference
ZVAL = {"Co": 9.0, "Ni": 10.0, "O": 6.0, "H": 1.0, "C": 4.0, "Mn": 7.0, "S": 6.0}

SYSTEMS_BULK = ["CoNiOH2", "CoNiOH2S-noH", "CoMnH2CO5"]
SYSTEMS_SLAB = ["CoNiOH2-slab", "CoNiOH2S-noH-slab", "CoNiOH2S-noH-slab-flip", "CoMnH2CO5-slab"]
SYSTEMS_HETERO = ["hetero_intrinsic", "hetero_s_doped", "hetero_s_exposed"]
SYSTEMS_ALL = SYSTEMS_BULK + SYSTEMS_SLAB + SYSTEMS_HETERO

SYSTEM_LABELS = {
    "CoNiOH2": "CoNiOH2 (Bulk)",
    "CoNiOH2S-noH": "CoNiOH2S (Bulk)",
    "CoMnH2CO5": "CoMnH2CO5 (Bulk)",
    "CoNiOH2-slab": "CoNiOH2 (Slab)",
    "CoNiOH2S-noH-slab": "CoNiOH2S (Slab)",
    "CoNiOH2S-noH-slab-flip": "CoNiOH2S (Slab-flip)",
    "CoMnH2CO5-slab": "CoMnH2CO5 (Slab)",
    "hetero_intrinsic": "Intrinsic Het.",
    "hetero_s_doped": "S-Doped Het.",
    "hetero_s_exposed": "S-Exposed Het.",
}

SYSTEM_COLORS = {
    "CoNiOH2": "#1f77b4",
    "CoNiOH2S-noH": "#2ca02c",
    "CoMnH2CO5": "#d62728",
    "CoNiOH2-slab": "#1f77b4",
    "CoNiOH2S-noH-slab": "#2ca02c",
    "CoNiOH2S-noH-slab-flip": "#ff7f0e",
    "CoMnH2CO5-slab": "#d62728",
    "hetero_intrinsic": "#1f77b4",
    "hetero_s_doped": "#2ca02c",
    "hetero_s_exposed": "#9467bd",
}

SYSTEM_GROUPS = {
    "Bulk": SYSTEMS_BULK,
    "Slab": SYSTEMS_SLAB,
    "Hetero": SYSTEMS_HETERO,
}


# ═══════════════════════════════════════════════
# Data Loading
# ═══════════════════════════════════════════════

def load_json(path: Path):
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


def load_dos(path: Path) -> CompleteDos | None:
    raw = load_json(path)
    if raw is None:
        return None
    try:
        dos_dict = raw.get("vasp_objects", {}).get("dos")
        if dos_dict is None:
            return None
        return MontyDecoder().process_decoded(dos_dict)
    except Exception as e:
        log.warning("DOS parse error: %s", e)
        return None


def load_band(path: Path) -> BandStructureSymmLine | None:
    raw = load_json(path)
    if raw is None:
        return None
    try:
        band_dict = raw.get("vasp_objects", {}).get("bandstructure")
        if band_dict is None:
            return None
        return MontyDecoder().process_decoded(band_dict)
    except Exception as e:
        log.warning("Band parse error: %s", e)
        return None


def load_structure(contcar: Path) -> Structure | None:
    if not contcar.exists():
        return None
    return Structure.from_file(contcar)


# ═══════════════════════════════════════════════
# Info Extraction
# ═══════════════════════════════════════════════

def extract_static_info(store_path: Path) -> dict:
    data = load_json(store_path)
    if not data:
        return {}
    for job in data:
        if job.get("name") == "static":
            output = job.get("output", {})
            inner = output.get("output", {})
            return {
                "energy": inner.get("energy"),
                "energy_per_atom": inner.get("energy_per_atom"),
                "bandgap": inner.get("bandgap"),
                "nsites": output.get("nsites"),
                "formula": output.get("formula_pretty"),
                "total_mag": output.get("total_magnetization"),
            }
    return {}


def extract_bandgap_dos(dos: CompleteDos) -> float | None:
    try:
        return dos.get_gap()
    except Exception:
        return None


def d_band_center(dos: CompleteDos, element: str, erange=(-10, 5)) -> float | None:
    try:
        return dos.get_band_center(band=OrbitalType.d, elements=[element], erange=erange)
    except Exception:
        return None


def vbm_cbm_from_band(band) -> dict:
    try:
        vbm = band.get_vbm()["energy"]
        cbm = band.get_cbm()["energy"]
        return {"vbm": vbm, "cbm": cbm, "bandgap_band": cbm - vbm}
    except Exception:
        return {"vbm": None, "cbm": None, "bandgap_band": None}


def extract_mag(outcar_path: Path) -> dict:
    try:
        o = Outcar(outcar_path)
        atomic = [m["tot"] for m in o.magnetization] if o.magnetization else []
        return {"total_mag": o.total_mag, "atomic_mags": atomic}
    except Exception:
        return {"total_mag": None, "atomic_mags": []}


def extract_work_function(locpot_path: Path, structure: Structure, outcar_path: Path = None) -> dict:
    """Extract work function from LOCPOT using pymatgen's WorkFunctionAnalyzer.

    Φ = V_vacuum - E_Fermi
    Automatically detects slab/vacuum regions via bond-length analysis.
    """
    try:
        # Decompress .gz to temp if needed
        tmp_path = None
        if str(locpot_path).endswith(".gz"):
            tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".LOCPOT")
            with gzip.open(locpot_path, "rb") as f_in:
                shutil.copyfileobj(f_in, tmp)
            tmp_path = tmp.name
            tmp.close()
        else:
            tmp_path = str(locpot_path)

        locpot = Locpot.from_file(tmp_path)
        avg_potential = locpot.get_average_along_axis(2)

        if tmp_path and str(locpot_path).endswith(".gz"):
            os.unlink(tmp_path)

        # E_Fermi from OUTCAR
        efermi = None
        if outcar_path and outcar_path.exists():
            o = Outcar(outcar_path)
            efermi = o.efermi

        if efermi is None:
            return {"error": "No E_Fermi available"}

        # Use pymatgen's built-in WorkFunctionAnalyzer
        from pymatgen.analysis.surface_analysis import WorkFunctionAnalyzer
        wfa = WorkFunctionAnalyzer(
            structure=structure,
            locpot_along_c=avg_potential,
            efermi=efermi,
        )

        return {
            "locpot_avg": avg_potential.tolist(),
            "z_coords": np.linspace(0, structure.lattice.c, len(avg_potential)).tolist(),
            "vacuum_level": wfa.vacuum_locpot,
            "efermi": efermi,
            "work_function": wfa.work_function,
        }
    except Exception as e:
        return {"error": str(e)}


# ═══════════════════════════════════════════════
# System Analysis
# ═══════════════════════════════════════════════

def analyze_system(name: str) -> dict:
    base = DATA_DIR / name
    result = {"name": name, "label": SYSTEM_LABELS.get(name, name)}

    # Static info from store.json
    store = base / "store.json"
    result.update(extract_static_info(store))

    # DOS
    dos_path = base / "dos_out.json"
    dos = load_dos(dos_path)
    if dos:
        result["dos"] = dos
        result["efermi"] = dos.efermi
        result["bandgap_dos"] = extract_bandgap_dos(dos)
        result["d_center_Co"] = d_band_center(dos, "Co")
        result["d_center_Ni"] = d_band_center(dos, "Ni")

    # Band
    band_path = base / "band_out.json"
    band = load_band(band_path)
    if band:
        result["band"] = band
        result["is_line_mode"] = isinstance(band, BandStructureSymmLine)
        result.update(vbm_cbm_from_band(band))

    # Structure from static step
    static_dir = base / "3-static"
    contcar = static_dir / "CONTCAR.gz"
    if not contcar.exists():
        contcar = static_dir / "CONTCAR"
    struct = load_structure(contcar)
    if struct:
        result["structure"] = struct
        lat = struct.lattice
        result["lattice"] = {"a": lat.a, "b": lat.b, "c": lat.c, "volume": lat.volume}
        # S-M bonds if present
        s_idx = next((i for i, s in enumerate(struct) if s.species_string == "S"), None)
        if s_idx is not None:
            result["has_S"] = True
            s_site = struct[s_idx]
            nns = struct.get_neighbors(s_site, r=3.5)
            metal_nn = [(n.index, n.species_string, n.nn_distance) for n in nns if n.species_string in ("Co", "Ni", "Mn")]
            metal_nn.sort(key=lambda x: x[2])
            result["S_bonds"] = {
                "indices": [x[0] for x in metal_nn[:4]],
                "elements": [x[1] for x in metal_nn[:4]],
                "lengths": [x[2] for x in metal_nn[:4]],
                "avg_length": np.mean([x[2] for x in metal_nn[:4]]),
            }

    # Magnetic moments
    outcar = static_dir / "OUTCAR.gz"
    if outcar.exists():
        mag = extract_mag(outcar)
        result["total_mag"] = mag["total_mag"]
        result["atomic_mags"] = mag["atomic_mags"]

    # Work function from LOCPOT (prefer dipole-corrected version)
    locpot_dipole = static_dir / "LOCPOT_dipole.gz"
    locpot_orig = static_dir / "LOCPOT.gz"
    outcar_dipole = static_dir / "OUTCAR_dipole.gz"
    outcar_orig = static_dir / "OUTCAR.gz"

    # Try dipole-corrected first, fall back to original
    if locpot_dipole.exists() and outcar_dipole.exists():
        locpot = locpot_dipole
        outcar_wf = outcar_dipole
    else:
        locpot = locpot_orig if locpot_orig.exists() else static_dir / "LOCPOT"
        outcar_wf = outcar_orig if outcar_orig.exists() else None

    if locpot.exists() and struct:
        wf = extract_work_function(locpot, struct, outcar_path=outcar_wf if outcar_wf and outcar_wf.exists() else None)
        result["work_function_data"] = wf
        result["work_function"] = wf.get("work_function")

    return result


# ═══════════════════════════════════════════════
# 1. Overview Table
# ═══════════════════════════════════════════════

def make_overview(results: dict) -> pd.DataFrame:
    rows = []
    props = [
        ("Layer", None), ("Energy (eV)", "energy"), ("E/atom (eV)", "energy_per_atom"),
        ("Band Gap (eV)", "bandgap_dos"), ("VBM (eV)", "vbm"), ("CBM (eV)", "cbm"),
        ("E_Fermi (eV)", "efermi"), ("Total Mag (μB)", "total_mag"),
        ("d-ctr Co (eV)", "d_center_Co"), ("d-ctr Ni (eV)", "d_center_Ni"),
        ("a (Å)", None), ("c (Å)", None), ("Vol (Å³)", None),
        ("Work Fn (eV)", "work_function"),
    ]

    # Identify layer for each system
    def get_layer(name):
        if name in SYSTEMS_BULK: return "Bulk"
        if name in SYSTEMS_SLAB: return "Slab"
        return "Hetero"

    for prop_name, key in props:
        row = {"Property": prop_name}
        for name in SYSTEMS_ALL:
            r = results[name]
            if prop_name == "Layer":
                row[name] = get_layer(name)
            elif prop_name == "a (Å)":
                lat = r.get("lattice", {})
                row[name] = f"{lat.get('a', '—'):.4f}" if lat.get('a') else "—"
            elif prop_name == "c (Å)":
                lat = r.get("lattice", {})
                row[name] = f"{lat.get('c', '—'):.4f}" if lat.get('c') else "—"
            elif prop_name == "Vol (Å³)":
                lat = r.get("lattice", {})
                row[name] = f"{lat.get('volume', '—'):.2f}" if lat.get('volume') else "—"
            else:
                val = r.get(key)
                if val is None:
                    row[name] = "—"
                else:
                    row[name] = f"{val:.4f}" if isinstance(val, float) else str(val)
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_DIR / "overview_10_systems.csv", index=False)
    log.info("Overview → output/overview_10_systems.csv")
    return df


def print_overview(results: dict):
    df = make_overview(results)
    print("\n" + "=" * 120)
    print("OVERVIEW: 10 Systems")
    print("=" * 120)
    print(df.to_string(index=False))

    # Key contrasts
    print("\n── Key Electronic Contrasts ──")
    b_gap = results["CoNiOH2"].get("bandgap_dos")
    s_gap = results["CoNiOH2S-noH"].get("bandgap_dos")
    if b_gap and s_gap:
        print(f"  S doping bandgap collapse: {b_gap:.3f} → {s_gap:.3f} eV  (Δ = {s_gap - b_gap:+.3f})")

    for het in SYSTEMS_HETERO:
        gap = results[het].get("bandgap_dos")
        if gap is not None:
            print(f"  {SYSTEM_LABELS[het]}: gap = {gap:.3f} eV")


# ═══════════════════════════════════════════════
# 2. DOS Comparison Plots
# ═══════════════════════════════════════════════

def plot_dos_grouped(results: dict):
    """Three panels: Bulk (3), Slab (4), Hetero (3) — TDOS comparison."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=False)
    e_range = (-6, 4)

    for ax_idx, (group_name, systems) in enumerate(SYSTEM_GROUPS.items()):
        ax = axes[ax_idx]

        for name in systems:
            dos = results[name].get("dos")
            if not dos:
                continue
            energies = dos.energies - dos.efermi
            mask = (energies >= e_range[0]) & (energies <= e_range[1])
            e_plot = energies[mask]
            dos_up = dos.get_densities(Spin.up)[mask]
            dos_down = dos.get_densities(Spin.down)[mask] if Spin.down in dos.densities else np.zeros_like(dos_up)

            color = SYSTEM_COLORS[name]
            label = SYSTEM_LABELS[name]
            gap = results[name].get("bandgap_dos")

            ax.plot(e_plot, dos_up, color=color, lw=1.0, label=label)
            ax.plot(e_plot, -dos_down, color=color, lw=0.7, alpha=0.5)
            ax.fill_between(e_plot, 0, dos_up, color=color, alpha=0.12)

            # Annotate gap if significant
            if gap and gap > 0.05:
                ax.axvspan(0, gap, color="gray", alpha=0.08, lw=0)
                ax.axvline(0, color="gray", ls=":", lw=0.4, alpha=0.5)
                ax.axvline(gap, color="gray", ls=":", lw=0.4, alpha=0.5)

        ax.axvline(0, color="k", ls="--", lw=0.5)
        ax.set_xlim(e_range)
        ax.set_title(f"{group_name} Systems — TDOS", fontsize=12, fontweight="bold")
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
    """Element-resolved PDOS: Co-3d, Ni-3d, O-2p, S-3p for selected systems."""
    focus_systems = ["CoNiOH2", "CoNiOH2S-noH", "CoNiOH2-slab", "CoNiOH2S-noH-slab",
                     "hetero_intrinsic", "hetero_s_doped"]

    n = len(focus_systems)
    fig, axes = plt.subplots(n, 1, figsize=(12, 3*n+1), sharex=True)
    e_range = (-6, 4)

    elem_style = {
        "Co": ("d", "#2171b5", "#6baed6"),
        "Ni": ("d", "#238b45", "#74c476"),
        "O":  ("p", "#cb181d", "#fb6a4a"),
        "S":  ("p", "#d94801", "#fd8d3c"),
    }

    for i, name in enumerate(focus_systems):
        ax = axes[i]
        dos = results[name].get("dos")
        if not dos:
            ax.text(0.5, 0.5, "No DOS", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(SYSTEM_LABELS[name])
            continue

        energies = dos.energies - dos.efermi
        mask = (energies >= e_range[0]) & (energies <= e_range[1])
        e_plot = energies[mask]

        has_s = any(s.species_string == "S" for s in dos.structure)

        for elem, (orb_type, c_up, c_down) in elem_style.items():
            if elem == "S" and not has_s:
                continue
            try:
                spd = dos.get_element_spd_dos(elem)
                ot = OrbitalType.d if orb_type == "d" else OrbitalType.p
                if ot not in spd:
                    continue
                edos = spd[ot]
                scale = 5.0 if elem == "S" else 1.0
                for spin, spin_up_flag, color, ls in [
                    (Spin.up, True, c_up, "-"),
                    (Spin.down, False, c_down, "--"),
                ]:
                    if spin in edos.densities:
                        d = edos.densities[spin][mask] * scale
                        ax.plot(e_plot, d if spin_up_flag else -d, color=color, ls=ls, lw=0.9, alpha=0.85)
            except Exception:
                pass

        gap = results[name].get("bandgap_dos")
        if gap and gap > 0.05:
            ax.axvspan(0, gap, color="gray", alpha=0.08, lw=0)
        ax.axvline(0, color="k", ls="--", lw=0.4)
        ax.set_xlim(e_range)
        ax.set_ylabel("PDOS")
        gap_str = f" (gap={gap:.3f}eV)" if gap else ""
        ax.set_title(f"{SYSTEM_LABELS[name]}{gap_str}", fontsize=10, fontweight="bold")
        ax.grid(True, alpha=0.2)

    axes[-1].set_xlabel("E − E$_F$ (eV)")

    # Legend
    legend_elements = [
        Line2D([0], [0], color="#2171b5", lw=1.5, label="Co-3d (up)"),
        Line2D([0], [0], color="#2171b5", lw=1.5, ls="--", alpha=0.5, label="Co-3d (down)"),
        Line2D([0], [0], color="#238b45", lw=1.5, label="Ni-3d (up)"),
        Line2D([0], [0], color="#238b45", lw=1.5, ls="--", alpha=0.5, label="Ni-3d (down)"),
        Line2D([0], [0], color="#cb181d", lw=1.5, label="O-2p (up)"),
        Line2D([0], [0], color="#d94801", lw=1.5, label="S-3p (×5, up)"),
    ]
    fig.legend(handles=legend_elements, fontsize=8, loc="upper right",
               framealpha=0.9, ncol=2, columnspacing=0.5)

    fig.suptitle("Element-Resolved PDOS", fontsize=14, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUTPUT_DIR / "pdos_elemental.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "pdos_elemental.pdf", bbox_inches="tight")
    plt.close(fig)
    log.info("Plot → output/pdos_elemental.{png,pdf}")


# ═══════════════════════════════════════════════
# 3. Band Structure Comparison
# ═══════════════════════════════════════════════

def plot_band_comparison(results: dict):
    """Band structures for key systems in 2×3 grid."""
    key_systems = [
        ("Bulk", ["CoNiOH2", "CoNiOH2S-noH", "CoMnH2CO5"]),
        ("Slab", ["CoNiOH2-slab", "CoNiOH2S-noH-slab", "CoMnH2CO5-slab"]),
    ]

    for group_name, systems in key_systems:
        n = len(systems)
        fig, axes = plt.subplots(1, n, figsize=(5*n, 5))
        if n == 1:
            axes = [axes]

        for i, name in enumerate(systems):
            ax = axes[i]
            band = results[name].get("band")
            is_line = results[name].get("is_line_mode")

            if band and is_line:
                plotter = BSPlotter(band)
                pd_data = plotter.bs_plot_data(zero_to_efermi=True)

                for bidx, dists in enumerate(pd_data["distances"]):
                    for spin_key in ("1", "-1"):
                        if spin_key in pd_data["energy"] and bidx < len(pd_data["energy"][spin_key]):
                            for ib in range(len(pd_data["energy"][spin_key][bidx])):
                                e_vals = pd_data["energy"][spin_key][bidx][ib]
                                color = "#2171b5" if spin_key == "1" else "#cb181d"
                                alpha = 0.7 if spin_key == "1" else 0.4
                                ax.plot(dists, e_vals, color=color, lw=0.5, alpha=alpha)

                # Bandgap
                vbm = results[name].get("vbm")
                cbm = results[name].get("cbm")
                if vbm and cbm and (cbm - vbm) > 0.05:
                    ax.axhspan(0, cbm - vbm, color="gray", alpha=0.08, lw=0)
                    ax.axhline(0, color="gray", ls=":", lw=0.4)
                    ax.axhline(cbm - vbm, color="gray", ls=":", lw=0.4)

                ax.axhline(0, color="k", ls="--", lw=0.4)
                ticks = pd_data["ticks"]
                for dist, label in zip(ticks["distance"], ticks["label"]):
                    ax.axvline(dist, color="k", lw=0.3, alpha=0.2)
                ax.set_xticks(ticks["distance"])
                ax.set_xticklabels([l.replace("$\\Gamma$", "Γ").replace("$", "") for l in ticks["label"]], fontsize=7)
                ax.set_xlim(ticks["distance"][0], ticks["distance"][-1])

                # Get y-range from data
                all_e = []
                for s_key in ("1", "-1"):
                    if s_key in pd_data["energy"]:
                        for bd in pd_data["energy"][s_key]:
                            for ib in bd:
                                all_e.extend(ib)
                if all_e:
                    m = (max(all_e) - min(all_e)) * 0.1
                    ax.set_ylim(min(all_e) - m, max(all_e) + m)
            else:
                ax.text(0.5, 0.5, "No band data", ha="center", va="center", transform=ax.transAxes)

            gap = results[name].get("bandgap_dos")
            gap_str = f", gap={gap:.3f}eV" if gap else ""
            ax.set_title(f"{SYSTEM_LABELS[name]}{gap_str}", fontsize=10, fontweight="bold")
            ax.grid(True, alpha=0.2)

        fig.suptitle(f"Band Structure — {group_name}", fontsize=13, fontweight="bold")
        fig.text(0.5, 0.01, "k-path", ha="center", fontsize=10)
        fig.text(0.01, 0.5, "E − E$_F$ (eV)", va="center", rotation="vertical", fontsize=10)
        plt.tight_layout(rect=[0.03, 0.03, 1, 0.95])
        fig.savefig(OUTPUT_DIR / f"band_{group_name.lower()}.png", dpi=200, bbox_inches="tight")
        fig.savefig(OUTPUT_DIR / f"band_{group_name.lower()}.pdf", bbox_inches="tight")
        plt.close(fig)
        log.info(f"Plot → output/band_{group_name.lower()}.{{png,pdf}}")


# ═══════════════════════════════════════════════
# 4. Interface Binding Energy (Hetero)
# ═══════════════════════════════════════════════

def analyze_interfaces(results: dict):
    print("\n" + "=" * 80)
    print("INTERFACE BINDING ENERGY")
    print("=" * 80)

    # Slab references: which LDH slab and CMCH slab correspond to each hetero
    # hetero_intrinsic   → CoNiOH2-slab + CoMnH2CO5-slab
    # hetero_s_doped     → CoNiOH2S-noH-slab + CoMnH2CO5-slab
    # hetero_s_exposed   → CoNiOH2S-noH-slab-flip + CoMnH2CO5-slab (closer to interface)

    ref_map = {
        "hetero_intrinsic": ("CoNiOH2-slab", "CoMnH2CO5-slab"),
        "hetero_s_doped": ("CoNiOH2S-noH-slab", "CoMnH2CO5-slab"),
        "hetero_s_exposed": ("CoNiOH2S-noH-slab-flip", "CoMnH2CO5-slab"),
    }

    print(f"  {'Hetero':<25s}  {'E_hetero':>12s}  {'E_LDH':>12s}  {'E_CMCH':>12s}  {'E_bind(eV)':>12s}  {'E_bind/area':>14s}")
    print(f"  {'─'*25}  {'─'*12}  {'─'*12}  {'─'*12}  {'─'*12}  {'─'*14}")

    for het_name, (ldh_name, cmch_name) in ref_map.items():
        eh = results[het_name].get("energy")
        el = results[ldh_name].get("energy")
        ec = results[cmch_name].get("energy")

        if eh is None or el is None or ec is None:
            print(f"  {het_name:<25s}  ENERGY DATA MISSING")
            continue

        e_bind = eh - el - ec

        # Estimate interface area from hetero structure (a*b in Å²)
        struct = results[het_name].get("structure")
        if struct:
            a, b = struct.lattice.a, struct.lattice.b
            area = a * b
            e_bind_per_area = e_bind / area * 1000  # meV/Å²
            area_str = f"{e_bind_per_area:+.2f} meV/Å²"
        else:
            area_str = "N/A"

        print(f"  {het_name:<25s}  {eh:12.4f}  {el:12.4f}  {ec:12.4f}  {e_bind:+12.4f}  {area_str:>14s}")


# ═══════════════════════════════════════════════
# 5. S Doping Effect — Bulk Analysis
# ═══════════════════════════════════════════════

def analyze_doping_effect(results: dict):
    print("\n" + "=" * 80)
    print("S DOPING EFFECT: CoNiOH2 → CoNiOH2S (Bulk)")
    print("=" * 80)

    p = results["CoNiOH2"]
    d = results["CoNiOH2S-noH"]

    # Energy
    e_p, e_d = p.get("energy"), d.get("energy")
    if e_p and e_d:
        print(f"\n  Substitution energy: ΔE = {e_d - e_p:+.4f} eV")

    # Lattice
    lat_p, lat_d = p.get("lattice", {}), d.get("lattice", {})
    if lat_p and lat_d:
        print(f"\n  Lattice change:")
        for axis in ("a", "c", "volume"):
            vp, vd = lat_p.get(axis), lat_d.get(axis)
            if vp and vd:
                print(f"    {axis}: {vp:.4f} → {vd:.4f}  (Δ = {vd-vp:+.4f}, {((vd-vp)/vp*100):+.2f}%)")

    # Bandgap — KEY finding
    gp, gd = p.get("bandgap_dos"), d.get("bandgap_dos")
    if gp is not None and gd is not None:
        print(f"\n  ★ Bandgap: {gp:.4f} → {gd:.4f} eV  (collapse of {gp-gd:.4f} eV)")

    # Magnetism
    mag_p, mag_d = p.get("total_mag"), d.get("total_mag")
    if mag_p is not None and mag_d is not None:
        print(f"\n  Total magnetic moment: {mag_p:.2f} → {mag_d:.2f} μB")

    # d-band centers
    for elem in ["Co", "Ni"]:
        dc_p = p.get(f"d_center_{elem}")
        dc_d = d.get(f"d_center_{elem}")
        if dc_p is not None and dc_d is not None:
            print(f"  {elem} d-band center: {dc_p:.4f} → {dc_d:.4f} eV  (Δ = {dc_d-dc_p:+.4f})")

    # S-M bonds
    sb = d.get("S_bonds", {})
    if sb:
        print(f"\n  S-M bonds:")
        for elem, bl in zip(sb["elements"], sb["lengths"]):
            print(f"    {elem}: {bl:.4f} Å")


# ═══════════════════════════════════════════════
# 6. Magnetism Analysis
# ═══════════════════════════════════════════════

def analyze_magnetism(results: dict):
    print("\n" + "=" * 80)
    print("MAGNETISM ACROSS ALL SYSTEMS")
    print("=" * 80)

    print(f"  {'System':<28s}  {'Total Mag (μB)':>16s}  {'Δ from bulk/slab':>18s}")
    print(f"  {'─'*28}  {'─'*16}  {'─'*18}")

    for name in SYSTEMS_ALL:
        mag = results[name].get("total_mag")
        if mag is None:
            continue

        # Find corresponding baseline
        base = None
        if name in SYSTEMS_BULK:
            base = 0.0
        elif name in SYSTEMS_SLAB:
            base_name = name.replace("-slab", "")
            base = results.get(base_name, {}).get("total_mag")
            if base is None:
                base_name = name.replace("-slab", "").replace("-flip", "")
                base = results.get(base_name, {}).get("total_mag")
        elif name in SYSTEMS_HETERO:
            base = 0.0

        delta = mag - base if base is not None else None
        delta_str = f"{delta:+.2f}" if delta is not None else "—"
        print(f"  {SYSTEM_LABELS[name]:<28s}  {mag:>16.4f}  {delta_str:>18s}")


# ═══════════════════════════════════════════════
# 7. Bandgap Trends Summary
# ═══════════════════════════════════════════════

def bandgap_trend(results: dict):
    print("\n" + "=" * 80)
    print("BANDGAP TREND SUMMARY")
    print("=" * 80)

    # Group by layer
    for group, systems in SYSTEM_GROUPS.items():
        print(f"\n  {group}:")
        for name in systems:
            gap = results[name].get("bandgap_dos")
            gap_str = f"{gap:.4f} eV" if gap else "—"
            has_s = "(S)" if "S" in name else ""
            print(f"    {name:<30s}  {gap_str:>12s}  {has_s}")


# ═══════════════════════════════════════════════
# 8. Work Function Summary
# ═══════════════════════════════════════════════

def work_function_summary(results: dict):
    print("\n" + "=" * 80)
    print("WORK FUNCTION (Φ)")
    print("=" * 80)
    print(f"  {'System':<30s}  {'Φ (eV)':>10s}")
    print(f"  {'─'*30}  {'─'*10}")

    for name in SYSTEMS_ALL:
        wf = results[name].get("work_function")
        if wf is not None:
            print(f"  {SYSTEM_LABELS[name]:<30s}  {wf:>10.4f}")


# ═══════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════

def main():
    print("=" * 80)
    print("CMCH@CoNi-LDH-S: COMPREHENSIVE ANALYSIS (10 Systems)")
    print("=" * 80)
    print("\nAnalyzing systems ...")

    # Analyze all systems
    results = {}
    for name in SYSTEMS_ALL:
        try:
            results[name] = analyze_system(name)
            gap = results[name].get("bandgap_dos")
            mag = results[name].get("total_mag")
            gap_str = f", gap={gap:.3f}" if gap else ""
            mag_str = f", mag={mag:.2f}μB" if mag else ""
            log.info(f"  ✓ {name}{gap_str}{mag_str}")
        except Exception as e:
            log.error(f"  ✗ {name}: {e}")
            results[name] = {"name": name, "label": SYSTEM_LABELS.get(name, name)}

    # === PRINTOUTS ===
    print_overview(results)
    analyze_doping_effect(results)
    bandgap_trend(results)
    analyze_interfaces(results)
    analyze_magnetism(results)
    work_function_summary(results)

    # === PLOTS ===
    print("\n" + "=" * 80)
    print("GENERATING FIGURES")
    print("=" * 80)

    plot_dos_grouped(results)
    plot_pdos_elemental(results)
    plot_band_comparison(results)

    # Simple scatter: bandgap vs d-band center
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax_idx, elem in enumerate(["Co", "Ni"]):
        ax = axes[ax_idx]
        xv, yv, lb, cl = [], [], [], []
        for name in SYSTEMS_ALL:
            gap = results[name].get("bandgap_dos")
            dc = results[name].get(f"d_center_{elem}")
            if gap is not None and dc is not None:
                xv.append(gap); yv.append(dc)
                lb.append(SYSTEM_LABELS[name]); cl.append(SYSTEM_COLORS[name])
        ax.scatter(xv, yv, c=cl, s=50, edgecolors="k", linewidth=0.5, zorder=5)
        for xi, yi, nm in zip(xv, yv, lb):
            ax.annotate(nm, (xi, yi), textcoords="offset points", xytext=(4, 4),
                        fontsize=6, alpha=0.85,
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", alpha=0.7, ec="none"))
        ax.set_xlabel("Band Gap (eV)", fontsize=10)
        ax.set_ylabel(f"{elem} d-band center (eV)", fontsize=10)
        ax.set_title(f"{elem} d-center vs Band Gap", fontsize=11, fontweight="bold")
        ax.grid(True, alpha=0.2)
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "dcenter_vs_gap.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "dcenter_vs_gap.pdf", bbox_inches="tight")
    plt.close(fig)
    log.info("Plot → output/dcenter_vs_gap.{png,pdf}")

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 80)


if __name__ == "__main__":
    main()
