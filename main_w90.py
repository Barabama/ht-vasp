#!/usr/bin/env python3
"""HT-VASP TB2J Workflow - Wannier90-based magnetic exchange (J parameters).

Endmember portal script for the Co-Fe-Mn-Ni TB2J pipeline::

    endmember POSCAR -> R3 relax (DoubleRelax) -> SCF + Wannier90 -> TB2J
                                                          (exchange.out J params)

Mirrors main_qha.py: GLOBAL_INCAR + run_tick + argparse (--slurm/--batch/--dry-run).
TC / BMAGN post-processing is left to em-tdb; this script only produces the J
parameters (stored as ``<name>-tb2j.json``).

Usage:
    python main_w90.py --name SER-Co                 # single structure locally
    python main_w90.py --batch                       # all structures locally
    python main_w90.py --dry-run SER-Co              # build structure+worker+flow, no VASP
    python main_w90.py --slurm                       # submit all to Slurm
    python main_w90.py --check                       # print J summaries of completed runs

Parallelism (run-critical): the W90 SCF must run with NCORE=1 / KPAR=1
(Wannier90/PEAD interface), so NPAR = NTASKS and NBANDS is pre-rounded to a
multiple of NTASKS (Tb2jInputSetGenerator.get_num_bands). Tb2jWorker defaults
to ntasks=32, and submit_jobs() submits with ``--ntasks 32`` to match. A
different Slurm --ntasks would let VASP silently change NBANDS and desync the
INCAR from wannier90.win num_bands.
"""

import json
import shutil
import logging
import argparse
import shlex
import os
from pathlib import Path
from htvasp.slurm import (
    artifact_marker_matches,
    current_artifact_matches,
    mark_artifact_from_ledger,
    new_run_id,
    payload_fingerprint,
    stable_run_id,
)

from htvasp.model import Endmember
from htvasp.workflows import Tb2jWorker
from htvasp.slurm import SlurmJobManager
from htvasp.tb2j import Tb2jResult
from htvasp.tb2j.context import build_magnetic_context, context_summary

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

root_dir = Path(__file__).resolve().parent

# VASP configuration (module purge + vasp-cpu, verified for TB2J runs)
VASP_ARGS = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-cpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-cpu && srun vasp_gam'",
}

# Initial magnetic moments per element (muB) - from verified submit_tb2j.sh
# 13-element initial guess table (restored; shared with main_qha.py)
MAGMOM = {
    "Al": 1.0,
    "Co": 3.0,
    "Cr": 3.0,
    "Cu": 5.0,
    "Fe": 5.0,
    "Mn": 5.0,
    "Nb": 3.0,
    "Ni": 2.0,
    "Ta": 3.0,
    "Ti": 3.0,
    "V": 3.0,
    "W": 3.0,
    "Zr": 1.0,
}

# Global INCAR for the TB2J pipeline. Values follow main_qha.py structure but
# take the physics from the verified submit_tb2j.sh (ENCUT=520, PREC=Accurate,
# LREAL=.FALSE., KPAR=1, ISPIN=2 + element MAGMOM dict, ISMEAR=1/SIGMA=0.2).
# MAGMOM is an element->moment dict: atomate2 expands it per-site at input-set
# generation time. Ionic keys apply to the R3 relax; the worker forces the W90
# SCF to a static run (NSW=0, IBRION=-1, KPAR=1, EDIFF=1e-6, LORBIT=11).
#
# LASPH=True and ENAUG=1360 are intentionally NOT listed below: atomate2's base
# VASP set (_BASE_VASP_SET['INCAR']) injects them into every generated INCAR,
# and they are deliberately kept. LASPH (aspherical corrections) and the higher
# ENAUG (augmentation grid) are generally more accurate for magnetic transition
# metals, and TB2J post-processing is insensitive to them. The verified
# tb2j_test baseline (submit_tb2j.sh) ran WITHOUT these two keys, so the small
# quantitative shift vs that baseline (e.g. BCC J +4.8%) is a known, accepted
# decision - not a bug. Do not pin LASPH=False / lower ENAUG to match it.
GLOBAL_INCAR = {
    # Basis / precision (verified for TB2J)
    "ENCUT": 520,
    "ISTART": 0,
    "ICHARG": 2,

    "ISMEAR": 1,
    "SIGMA": 0.2,
    "ALGO": "Normal",
    "NELM": 120,
    "NELMIN": 6,
    "NELMDL": -6,
    # Ionic (R3 relax; W90 overridden by the worker)
    "IBRION": 2,
    "ISIF": 3,
    "NSW": 80,
    "POTIM": 0.2,
    "EDIFF": 1e-5,          # R3 convergence; W90 forced to 1e-6 by the worker
    "EDIFFG": -0.02,
    # Magnetism
    "ISPIN": 2,
    "MAGMOM": MAGMOM,
    "ISYM": 0,
    "LREAL": False,         # verified for TB2J (submit_tb2j.sh uses .FALSE.)
    "PREC": "Accurate",
    "SYMPREC": 1e-5,
    # Output
    "LWAVE": False,
    "LCHARG": False,
    "LORBIT": 11,
    "LOPTICS": False,
    "LVTOT": False,
    "GGA": "PE",
    # Parallelisation: Wannier90 interface requires KPAR=1 (keep R3 consistent)
    "KPAR": 1,
    "NCORE": 1,
    # Magnetic mixing deliberately left at VASP defaults (no AMIX/BMIX/AMIX_MAG/BMIX_MAG; 
    # NELMIN/NELMDL set to VASP defaults above). 
    # Custom mixing (AMIX=0.2, BMIX=1e-4, AMIX_MAG=0.8, BMIX_MAG=1e-4) drove the W90 SCF 
    # of antiferromagnetic FCC-Fe-Mn into a non-magnetic minimum (mag -> 0.001 muB); 
    # VASP defaults (AMIX=0.4, BMIX=1.0, AMIX_MAG=1.6, BMIX_MAG=1.0) converge the 
    # AFM state (Fe -1.48 / Mn +0.41 muB, matching the tb2j_test baseline).
    # Ferromagnetic systems (Co/Ni) are insensitive to this and remain stable.
    # "AMIX": 0.2,
    # "BMIX": 1e-4,
    # "AMIX_MAG": 0.8,
    # "BMIX_MAG": 1e-4,
}

# R3 relax overrides (on top of GLOBAL_INCAR) - matches submit_tb2j.sh INCAR_R3
RELAX_INCAR = {
    "ISIF": 3,               # full structural relaxation (also set by the worker)
    "NSW": 80,
    "POTIM": 0.2,
    "EDIFFG": -0.02,
    "NELM": 120,
    "LORBIT": 11,
}

# Co-Fe-Mn-Ni full endmember set (SER pure elements + BCC/FCC binary pairs).
# BCC-X-Y / FCC-X-Y: phase template with X on sublattice 1, Y on sublattice 2.
STRUCTURE_NAMES = [
    # Al-Nb-Ti-V
    "SER-Al",
    "SER-Nb",
    "SER-Ti",
    "SER-V",
    "BCC-Al-Al",
    "BCC-Al-Nb",
    "BCC-Al-Ti",
    "BCC-Al-V",
    "BCC-Nb-Nb",
    "BCC-Nb-Ti",
    "BCC-Nb-V",
    "BCC-Ti-Ti",
    "BCC-Ti-V",
    "BCC-V-V",

    # Co-Fe-Mn-Ni
    "SER-Co",
    "SER-Fe",
    "SER-Mn",
    "SER-Ni",
    "BCC-Co-Co",
    "BCC-Co-Fe",
    "BCC-Co-Mn",
    "BCC-Co-Ni",
    "BCC-Fe-Fe",
    "BCC-Fe-Mn",
    "BCC-Fe-Ni",
    "BCC-Mn-Mn",
    "BCC-Mn-Ni",
    "BCC-Ni-Ni",
    "FCC-Co-Co",
    "FCC-Co-Fe",
    "FCC-Co-Mn",
    "FCC-Co-Ni",
    "FCC-Fe-Co",
    "FCC-Fe-Fe",
    "FCC-Fe-Mn",
    "FCC-Fe-Ni",
    "FCC-Mn-Co",
    "FCC-Mn-Fe",
    "FCC-Mn-Mn",
    "FCC-Mn-Ni",
    "FCC-Ni-Co",
    "FCC-Ni-Fe",
    "FCC-Ni-Mn",
    "FCC-Ni-Ni",
]


def structure_elements(structure) -> list[str]:
    """Distinct element symbols in POSCAR site order (mirrors Tb2jWorker)."""
    seen: list[str] = []
    for site in structure:
        sym = site.species_string
        if sym not in seen:
            seen.append(sym)
    return seen


def _dry_run(name: str, structure, magnetic_context: dict | None = None) -> None:
    """Build the structure + Tb2jWorker + flow and print flow info (no VASP)."""
    magnetic_context = magnetic_context or build_magnetic_context(name, structure)
    log.info(f"  Structure {name}: formula={structure.formula}, atoms={len(structure)}")
    log.info(f"  Lattice abc: {structure.lattice.abc}  angles: {structure.lattice.angles}")
    log.info(f"  Species (site order): {structure_elements(structure)}")
    log.info(f"  Magnetic context: {context_summary(magnetic_context)}")

    worker = Tb2jWorker(
        vasp_args=VASP_ARGS,
        potcar_functional="PBE",
        global_incar=GLOBAL_INCAR,
        relax_incar=RELAX_INCAR,
        kmesh=(9, 9, 9),
        magnetic_context=magnetic_context,
    )
    flow = worker._make_flow(structure)
    # Elements resolved by the worker at flow-build time
    elements = worker.tb2j_maker.elements or worker.tb2j_maker.input_set_generator.elements
    num_wann = worker.tb2j_maker.input_set_generator.get_num_wann(structure)

    log.info(f"  Tb2jWorker: elements={elements}, num_wann={num_wann}, kmesh={worker.kmesh}")
    log.info(f"  Flow: {flow.name}")
    for job in flow.jobs:
        log.info(f"    job: {job.name}")
    log.info(f"  Dry-run OK (no VASP launched)")


def run_tick(name: str, force: bool = False, dry_run: bool = False, run_id: str | None = None):
    """Run a single endmember structure through the TB2J workflow."""
    endmember = Endmember()
    flow_dir = Path("/nfs_ssd/tmp")
    if force:
        run_id = run_id or new_run_id("w90")
    if run_id and Path(run_id).name != run_id:
        raise ValueError(f"Invalid run_id: {run_id!r}")
    store_base = root_dir / "data" / "endmembers" / name / "tb2jflow"
    store_dir = store_base / "runs" / run_id if run_id else store_base
    json_path = store_dir / f"{name}-tb2j.json" if run_id else store_dir.parent / f"{name}-tb2j.json"

    # Dry-run is always an inspection request: it must build and print the
    # current context even when a previous artifact already exists.
    if dry_run:
        struct = endmember.get_poscar(name, root_dir / "data" / "poscars")
        magnetic_context = build_magnetic_context(name, struct)
        _dry_run(name, struct, magnetic_context=magnetic_context)
        return

    # Skip if already done
    if not force and json_path.exists() and (
        (run_id is None and not os.environ.get("HTVASP_INTENT_KEY"))
        or current_artifact_matches(store_base / "slurm_submissions.json", json_path)
    ):
        log.info(f"Structure {name} already done (use --rerun to re-run)")
        return

    struct = endmember.get_poscar(name, root_dir / "data" / "poscars")
    magnetic_context = build_magnetic_context(name, struct)

    log.info(f"Structure {name} start")

    try:
        worker = Tb2jWorker(
            vasp_args=VASP_ARGS,
            potcar_functional="PBE",
            global_incar=GLOBAL_INCAR,
            relax_incar=RELAX_INCAR,
            # elements=None -> auto-derived from the structure site order
            kmesh=(9, 9, 9),
            magnetic_context=magnetic_context,
        )
        worker.run_flow(
            name=name,
            structure=struct,
            flow_dir=flow_dir,
            store_dir=store_dir,
            resume=not force,
        )
        output = worker.get_result()
        worker.write_result(data=output, json_path=json_path)
        mark_artifact_from_ledger(store_base / "slurm_submissions.json", json_path)
        if output:
            _report_result(name, output)
    except Exception as e:
        log.error(f"Structure {name} failed: {e}")

    log.info(f"Structure {name} done")


def _report_result(name: str, output: dict) -> None:
    """Print a short J summary from a fresh TB2J output dict."""
    try:
        summary = Tb2jResult.from_dict(output).to_summary()
        log.info(
            f"  {name}: {summary['num_j_pairs']} J pairs, {summary['num_shells']} shells"
        )
        if summary.get("shells"):
            nn = summary["shells"][0]
            log.info(
                f"  NN shell: d={nn['distance']:.3f} A, count={nn['count']}, "
                f"mean J={nn['mean_J_meV']} meV"
            )
        if summary.get("max_abs_J_meV") is not None:
            log.info(f"  max|J|: {summary['max_abs_J_meV']} meV")
    except Exception as e:
        log.warning(f"  {name}: failed to summarise result: {e}")


def check_jobs() -> None:
    """Scan STRUCTURE_NAMES and print a J-parameter summary per completed run."""
    for name in STRUCTURE_NAMES:
        json_path = root_dir / "data" / "endmembers" / name / f"{name}-tb2j.json"
        if not json_path.exists():
            log.info(f"{name}: not run")
            continue
        try:
            with open(json_path, encoding="utf-8") as f:
                data = json.load(f)
            summary = Tb2jResult.from_dict(data).to_summary()

            moments = data.get("atom_moments", [])
            m_str = (
                ", ".join(f"{m['atom']}={m['w_magmom']:.3f}" for m in moments)
                if moments
                else "n/a"
            )

            shells = summary.get("shells", [])
            nn_str = (
                f"d={shells[0]['distance']:.3f}A J={shells[0]['mean_J_meV']}meV"
                if shells
                else "n/a"
            )

            log.info(
                f"{name}: {summary['num_j_pairs']} J pairs | "
                f"{summary['num_shells']} shells | NN {nn_str} | "
                f"max|J| {summary.get('max_abs_J_meV', 'n/a')} meV | "
                f"moments {m_str}"
            )
        except Exception as e:
            log.warning(f"{name}: failed to parse result: {e}")


def run_batch(force: bool = False, dry_run: bool = False):
    """Run all structures locally."""
    for name in STRUCTURE_NAMES:
        run_tick(name, force=force, dry_run=dry_run)


def submit_jobs(force: bool = False) -> None:
    """Submit every structure to Slurm via SlurmJobManager."""
    manager = SlurmJobManager()
    for name in STRUCTURE_NAMES:
        # ntasks must match Tb2jWorker.ntasks (default 32): NBANDS is
        # pre-rounded to a multiple of NPAR (= NTASKS with KPAR=NCORE=1) in
        # Tb2jInputSetGenerator.get_num_bands(), so a Slurm --ntasks different
        # from the worker's ntasks would let VASP silently change NBANDS and
        # desync it from wannier90.win num_bands.
        config = manager.get_cpu_config(
            job_name=f"{name}-tb2j",
            output_log=f"logs/{name}-tb2j.log",
            ntasks=32,
            memory="64G",
        )
        run_payload = {
            "system": name,
            "stage": "tb2j",
            "recipe": "main_w90-v2",
            "vasp": VASP_ARGS,
            "incar": GLOBAL_INCAR,
            "relax_incar": RELAX_INCAR,
            "kmesh": (9, 9, 9),
        }
        run_id = new_run_id("w90") if force else stable_run_id("w90", run_payload)
        command = f"{shlex.quote(str(root_dir / '.conda' / 'bin' / 'python'))} {shlex.quote(str(Path(__file__).resolve()))} --name {shlex.quote(name)} --run-id {shlex.quote(run_id)}"
        if force:
            command += " --rerun"
        artifact_path = root_dir / "data" / "endmembers" / name / f"tb2jflow/runs/{run_id}/{name}-tb2j.json"
        result = manager.submit_intent(
            system_name=name,
            stage="tb2j",
            recipe_version="main_w90-v2",
            settings_fingerprint=payload_fingerprint({"vasp": VASP_ARGS, "incar": GLOBAL_INCAR, "relax": RELAX_INCAR, "kmesh": (9, 9, 9)}),
            manifest_fingerprint=payload_fingerprint({"system": name, "stage": "tb2j", "recipe": "main_w90-v2"}),
            command=command,
            config=config,
            ledger_path=root_dir / "data" / "endmembers" / name / "tb2jflow" / "slurm_submissions.json",
            conda_env="htvasp",
            work_dir=root_dir,
            run_id=run_id,
            retry=force,
            artifact_path=artifact_path,
            artifact_checker=lambda entry: artifact_marker_matches(artifact_path, entry),
            artifact_shortcut=True,
        )
        if not result.ok:
            log.error(f"Failed to submit job for {name}")
        else:
            log.info(f"Submitted job for {name}: {result.state} ({result.job_id})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="TB2J magnetic exchange workflow (Wannier90-based)"
    )
    parser.add_argument("--name", "-n", type=str, help="Run single structure by name")
    parser.add_argument("--batch", "-b", action="store_true", help="Run all structures locally")
    parser.add_argument("--slurm", "-s", action="store_true", help="Submit all structures to Slurm")
    parser.add_argument("--rerun", "-r", action="store_true", help="Force re-run (overwrite existing results)")
    parser.add_argument(
        "--dry-run",
        type=str,
        nargs="?",
        const="",
        metavar="NAME",
        help="Build structure + worker + flow and print (no VASP). "
        "With a NAME for one structure; alone or with --batch for all.",
    )
    parser.add_argument("--check", action="store_true", help="Print J summaries of completed runs")
    parser.add_argument("--run-id", help="Operate on an explicit immutable run")
    args = parser.parse_args()

    if args.check:
        check_jobs()
    elif args.dry_run is not None:
        names = (
            [args.dry_run]
            if args.dry_run
            else ([args.name] if args.name else STRUCTURE_NAMES)
        )
        for name in names:
            run_tick(name, force=args.rerun, dry_run=True, run_id=args.run_id)
    elif args.name:
        run_tick(args.name, force=args.rerun, run_id=args.run_id)
    elif args.batch:
        run_batch(force=args.rerun)
    elif args.slurm:
        submit_jobs(force=args.rerun)
    else:
        parser.print_help()
