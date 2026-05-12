"""
Bader charge analysis for CoNiHOS-Co3-noH.

Runs Bader analysis on existing static calculation or optionally
re-runs static with correct MAGMOM + LAECHG.

Usage:
    python bader_charge_analysis.py              # Run on existing 3-static data
    python bader_charge_analysis.py --redo       # Re-run static with corrected settings then Bader
    python bader_charge_analysis.py --no-spin    # Skip spin density analysis
"""

import argparse
import json
import logging
import os
import shutil
from pathlib import Path

from pymatgen.analysis.local_env import CrystalNN
from pymatgen.command_line.bader_caller import bader_analysis_from_path
from pymatgen.core import Structure

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

# ── Paths ──────────────────────────────────────────────────────────
STATIC_DIR = Path("data/CoNiHOS-Co3-noH/3-static")
OUT_DIR = Path("data/CoNiHOS-Co3-noH/bader_analysis")
PROJECT_DIR = Path("/workspace/gaominliang/ht-vasp")


def find_s_neighbors(structure: Structure, dist_cutoff: float = 3.0) -> dict:
    """Find S atom and its neighboring Co atoms using CrystalNN and distance cutoff."""
    s_sites = [i for i, s in enumerate(structure) if s.species_string == "S"]
    if not s_sites:
        raise ValueError("No S atom found in structure")

    s_idx = s_sites[0]

    nn_engine = CrystalNN()
    neighbors = nn_engine.get_nn_info(structure, s_idx)

    co_neighbors = [n for n in neighbors if n["site"].species_string == "Co"]

    dist_data = {}
    for i, site in enumerate(structure):
        if site.species_string == "Co":
            d = structure[s_idx].distance(site)
            dist_data[i] = d

    dist_sorted = sorted(dist_data.items(), key=lambda x: x[1])

    return {
        "s_idx": s_idx,
        "nn_co": [(n["site_index"], n["site"].species_string)
                  for n in co_neighbors],
        "all_co_dists": dist_sorted,
        "s_species": str(structure[s_idx].species),
    }


def run_bader_on_static(static_dir: Path, output_dir: Path):
    """Run Bader charge + spin analysis on a static calculation directory."""
    output_dir.mkdir(parents=True, exist_ok=True)

    log.info(f"Running Bader analysis on {static_dir}")
    summary = bader_analysis_from_path(str(static_dir))

    # Supplement with per-atom charge transfer
    log.info("Bader analysis completed. Summary:")
    log.info(f"  Vacuum charge: {summary['vacuum_charge']:.4f}")
    log.info(f"  Vacuum volume: {summary['vacuum_volume']:.4f}")
    log.info(f"  Bader version: {summary['bader_version']}")
    log.info(f"  Reference used (AECCAR0+AECCAR2): {summary['reference_used']}")
    log.info(f"  Charge transfer available: {'charge_transfer' in summary}")
    if "magmom" in summary:
        log.info(f"  Spin density analysis: available")

    return summary


def get_atom_data(summary: dict, structure: Structure, potcar_path: Path | None = None):
    """Extract per-atom Bader data, merging with structure and magnetic info."""
    charges = summary["charge"]
    min_dists = summary["min_dist"]
    atomic_vols = summary["atomic_volume"]
    charge_transfer = summary.get("charge_transfer", [])
    magmoms = summary.get("magmom", [])

    atoms = []
    for i in range(len(charges)):
        atom = {
            "index": i,
            "element": str(structure[i].species),
            "bader_charge": charges[i],
            "min_dist": min_dists[i] if i < len(min_dists) else None,
            "atomic_vol": atomic_vols[i] if i < len(atomic_vols) else None,
            "charge_transfer": charge_transfer[i] if i < len(charge_transfer) else None,
            "magmom": magmoms[i] if i < len(magmoms) else None,
        }
        atoms.append(atom)
    return atoms


def print_s_environment(s_info: dict, atoms: list[dict]):
    """Pretty-print S and its neighboring Co atoms."""
    s_idx = s_info["s_idx"]
    s_atom = atoms[s_idx]

    print("\n" + "=" * 70)
    print(f"  S ATOM (site {s_idx}): Bader charge = {s_atom['bader_charge']:.4f} e")
    if s_atom["charge_transfer"] is not None:
        sign = "+" if s_atom["charge_transfer"] > 0 else ""
        print(f"    Charge transfer: {sign}{s_atom['charge_transfer']:.4f} e  "
              f"(positive = gained e⁻, anionic)")
    if s_atom["magmom"] is not None:
        print(f"    Magnetic moment: {s_atom['magmom']:.4f} μB")
    print(f"    Atomic volume:   {s_atom['atomic_vol']:.2f} Å³")

    print(f"\n  S NEIGHBORS (Co within coordination shell):")
    nn_sites = [idx for idx, _ in s_info["nn_co"]]
    for idx, elem in s_info["nn_co"]:
        a = atoms[idx]
        print(f"    Co site {idx}: Bader charge = {a['bader_charge']:.4f} e,  "
              f"transfer = {a['charge_transfer']:+.4f} e,  "
              f"magmom = {a['magmom']:.4f} μB")

    print(f"\n  ALL Co ATOMS (sorted by distance from S):")
    print(f"    {'Site':>5s}  {'Dist(Å)':>8s}  {'Bader(e)':>9s}  {'Transfer(e)':>11s}  {'Mag(μB)':>8s}  {'Vol(Å³)':>8s}")
    print(f"    {'─'*5}  {'─'*8}  {'─'*9}  {'─'*11}  {'─'*8}  {'─'*8}")
    for idx, dist in s_info["all_co_dists"]:
        a = atoms[idx]
        nn_marker = " ◀" if idx in nn_sites else ""
        ct_str = f"{a['charge_transfer']:+.4f}" if a["charge_transfer"] is not None else "N/A"
        mm_str = f"{a['magmom']:.4f}" if a['magmom'] is not None else "N/A"
        print(f"    {idx:5d}  {dist:8.3f}  {a['bader_charge']:9.4f}  {ct_str:>11s}  {mm_str:>8s}  {a['atomic_vol']:8.2f}{nn_marker}")

    print(f"\n  SUMMARY:")
    nn_charges = [atoms[idx]["charge_transfer"] for idx, _ in s_info["nn_co"]]
    nn_bader = [atoms[idx]["bader_charge"] for idx, _ in s_info["nn_co"]]
    nn_mag = [atoms[idx]["magmom"] for idx, _ in s_info["nn_co"]]
    print(f"    S charge transfer:            {s_atom['charge_transfer']:+.4f} e")
    print(f"    3 Co-NN avg transfer:          {sum(nn_charges)/3:+.4f} e")
    print(f"    3 Co-NN avg Bader charge:      {sum(nn_bader)/3:.4f} e")
    print(f"    3 Co-NN avg magnetic moment:   {sum(nn_mag)/3:.4f} μB")
    print("=" * 70)


def save_results(atoms: list[dict], s_info: dict, summary: dict, output_dir: Path):
    """Save full Bader results as JSON."""
    result = {
        "s_info": {
            "s_idx": s_info["s_idx"],
            "nn_co_sites": [idx for idx, _ in s_info["nn_co"]],
        },
        "vacuum_charge": summary["vacuum_charge"],
        "vacuum_volume": summary["vacuum_volume"],
        "reference_used": summary["reference_used"],
        "bader_version": summary["bader_version"],
        "atoms": atoms,
    }

    json_path = output_dir / "bader_results.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, default=str)
    log.info(f"Full results saved to {json_path}")


def run_redo_static_and_bader():
    """Re-run static calculation with corrected MAGMOM, then Bader."""
    import sys
    sys.path.insert(0, str(PROJECT_DIR))

    from atomate2.vasp.jobs.core import StaticMaker
    from atomate2.vasp.sets.core import StaticSetGenerator
    from atomate2.vasp.flows.core import DoubleRelaxMaker
    from jobflow import Flow
    from htvasp.workflows.base import Worker
    from htvasp.utils.local import run_locally_custom
    from jobflow import JobStore
    from maggma.stores import JSONStore, MemoryStore

    # Read relaxed structure from CONTCAR of last relax
    relax2_dir = Path("data/CoNiHOS-Co3-noH/2-relax_2")
    structure = Structure.from_file(relax2_dir / "CONTCAR.gz")

    # Set MAGMOM explicitly in INCAR: Co=5.0, Ni=2.0, H/O/S=0.6
    # VASP orders by element type from POSCAR/POTCAR
    global_incar = {
        "ENCUT": 520,
        "PREC": "Accurate",
        "ALGO": "Fast",
        "NELM": 200,
        "EDIFF": 1e-6,
        "ISPIN": 2,
        "MAGMOM": {"Co": 5.0, "Ni": 2.0, "H": 0.6, "O": 0.6, "S": 0.6},
        "AMIX": 0.1,
        "BMIX": 1e-4,
        "AMIX_MAG": 0.4,
        "BMIX_MAG": 1e-4,
        "LREAL": "Auto",
        "KPAR": 4,
        "NCORE": 4,
        "GGA": "PE",
        "IVDW": 12,
        "LDAU": True,
        "LDAUTYPE": 2,
        "LDAUPRINT": 1,
        "LASPH": True,
        "LMAXMIX": 4,
        "LDAUL": {"Co": 2, "Ni": 2, "H": -1, "O": -1, "S": -1},
        "LDAUU": {"Co": 3.32, "Ni": 6.20},
        "LDAUJ": {"Co": 0.0, "Ni": 0.0},
    }

    vasp_args = {
        "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_std'",
        "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_gam'",
    }

    static_maker = StaticMaker(
        run_vasp_kwargs={
            "handlers": [],
            "custodian_kwargs": {"max_errors_per_job": 3, "gzipped_output": False},
            **vasp_args,
        },
        stop_children_kwargs={"handle_unsuccessful": False},
        copy_vasp_kwargs={"additional_vasp_files": ("WAVECAR",)},
        input_set_generator=StaticSetGenerator(
            user_potcar_functional="PBE_64",
            user_incar_settings={
                **global_incar,
                "ISTART": 1,
                "ICHARG": 11,  # Read charge density from previous CHGCAR
                "IBRION": -1,
                "ISIF": 2,
                "NSW": 0,
                "LWAVE": True,
                "LCHARG": True,
                "LAECHG": True,
                "LORBIT": 11,
                "ISMEAR": -5,
                "ALGO": "Normal",
            },
        ),
    )

    static_job = static_maker.make(structure, prev_dir=str(STATIC_DIR))

    flow = Flow([static_job], output=static_job.output)

    work_dir = Path("data/CoNiHOS-Co3-noH/bader_static_redo")
    store = JobStore(
        JSONStore(str(work_dir / "store.json"), read_only=False),
        additional_stores={"data": MemoryStore()},
    )

    log.info(f"Running redo static in {work_dir}")
    run_locally_custom(flow, store=store, root_dir=work_dir,
                       ensure_success=True, raise_immediately=False, resume=False)

    # After flow completes, run Bader on the new static directory
    static_dirs = sorted(work_dir.glob("1-*"))
    if static_dirs:
        new_static_dir = static_dirs[0]
        log.info(f"Running Bader on new static: {new_static_dir}")
        return run_bader_on_static(new_static_dir, OUT_DIR)
    else:
        raise RuntimeError(f"No static dir found in {work_dir}")


def main():
    parser = argparse.ArgumentParser(description="Bader charge analysis for CoNiHOS-Co3-noH")
    parser.add_argument("--redo", action="store_true",
                        help="Re-run static with corrected MAGMOM+ICHARG=11 before Bader")
    parser.add_argument("--no-spin", action="store_true",
                        help="Skip spin density (magnetization) analysis")
    args = parser.parse_args()

    # 1. Run Bader analysis
    if args.redo:
        summary = run_redo_static_and_bader()
    else:
        summary = run_bader_on_static(STATIC_DIR, OUT_DIR)

    # 2. Read structure from CONTCAR
    structure = Structure.from_file(STATIC_DIR / "CONTCAR.gz")

    # 3. Find S and its Co neighbors
    s_info = find_s_neighbors(structure)

    # 4. Get per-atom data
    atoms = get_atom_data(summary, structure)

    # 5. Print results
    print_s_environment(s_info, atoms)

    # 6. Save full results
    save_results(atoms, s_info, summary, OUT_DIR)

    log.info("Done.")


if __name__ == "__main__":
    main()
