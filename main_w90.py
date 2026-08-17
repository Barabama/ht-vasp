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
from pathlib import Path

from htvasp.model import Endmember
from htvasp.workflows import Tb2jWorker
from htvasp.slurm import SlurmJobManager
from htvasp.tb2j import Tb2jResult

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

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
    "PREC": "Accurate",
    "ISTART": 0,
    "ICHARG": 2,
    "EDIFF": 1e-5,           # R3 convergence; W90 forced to 1e-6 by the worker
    "EDIFFG": -0.02,
    # Magnetism
    "ISPIN": 2,
    "MAGMOM": MAGMOM,
    "ISMEAR": 1,
    "SIGMA": 0.2,
    # LMIXTAU=False / LMAXMIX=2 (VASP defaults) - explicitly override atomate2's
    # defaults (LMIXTAU=True, LMAXMIX=4) so the W90 INCAR matches the verified
    # tb2j_test baseline (no LMIXTAU/LMAXMIX keys) that converges the AFM state.
    "LMIXTAU": False,
    "LMAXMIX": 2,
    "ALGO": "Fast",
    "NELM": 120,
    # Matching VASP defaults (NELMIN=2, NELMDL=-5). These explicitly override
    # the Worker base-class defaults (NELMIN=6, NELMDL=-6); kept here so the
    # base defaults cannot re-inject the custom electronic-step settings that
    # (with the custom mixing) drove the AFM W90 SCF into a non-magnetic state.
    "NELMIN": 2,
    "NELMDL": -5,
    # Ionic (R3 relax; W90 overridden by the worker)
    "IBRION": 2,
    "ISIF": 3,
    "NSW": 80,
    "POTIM": 0.2,
    "LREAL": False,          # verified for TB2J (submit_tb2j.sh uses .FALSE.)
    "ISYM": 0,
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
    "NCORE": 2,
    # Magnetic mixing deliberately left at VASP defaults (no AMIX/BMIX/
    # AMIX_MAG/BMIX_MAG; NELMIN/NELMDL set to VASP defaults above). Custom
    # mixing (AMIX=0.2, BMIX=1e-4, AMIX_MAG=0.8, BMIX_MAG=1e-4, NELMIN=6,
    # NELMDL=-6) drove the W90 SCF of antiferromagnetic FCC-Fe-Mn into a
    # non-magnetic minimum (mag -> 0.001 muB); VASP defaults (AMIX=0.4,
    # BMIX=1.0, AMIX_MAG=1.6, BMIX_MAG=1.0, NELMIN=2, NELMDL=-5) converge the
    # AFM state (Fe -1.48 / Mn +0.41 muB, matching the tb2j_test baseline).
    # Ferromagnetic systems (Co/Ni) are insensitive to this and remain stable.
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
    # SER pure elements
    "SER-Co",
    "SER-Fe",
    "SER-Mn",
    "SER-Ni",
    # BCC endmembers (4 elements -> 10 binary pairs)
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
    # FCC endmembers (all 16 directional pairs)
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


def _dry_run(name: str, structure) -> None:
    """Build the structure + Tb2jWorker + flow and print flow info (no VASP)."""
    log.info(f"  Structure {name}: formula={structure.formula}, atoms={len(structure)}")
    log.info(f"  Lattice abc: {structure.lattice.abc}  angles: {structure.lattice.angles}")
    log.info(f"  Species (site order): {structure_elements(structure)}")

    worker = Tb2jWorker(
        vasp_args=VASP_ARGS,
        potcar_functional="PBE",
        global_incar=GLOBAL_INCAR,
        relax_incar=RELAX_INCAR,
        kmesh=(9, 9, 9),
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


def run_tick(name: str, force: bool = False, dry_run: bool = False):
    """Run a single endmember structure through the TB2J workflow."""
    endmember = Endmember()
    flow_dir = Path("/nfs_ssd/tmp")
    store_dir = Path("data/endmembers") / name / "tb2jflow"
    json_path = store_dir / f"{name}-tb2j.json"

    # Skip if already done
    if not force and json_path.exists():
        log.info(f"Structure {name} already done (use --rerun to re-run)")
        return

    struct = endmember.get_poscar(name, Path("data/poscars"))

    if dry_run:
        # Read-only: build the structure + worker + flow and print, never touch
        # the store. --dry-run --rerun must NOT rmtree the completed store_dir
        # (dry-run is a pure dry-run of the run path).
        _dry_run(name, struct)
        return

    log.info(f"Structure {name} start")

    if force and store_dir.exists():
        shutil.rmtree(store_dir)

    try:
        worker = Tb2jWorker(
            vasp_args=VASP_ARGS,
            potcar_functional="PBE",
            global_incar=GLOBAL_INCAR,
            relax_incar=RELAX_INCAR,
            # elements=None -> auto-derived from the structure site order
            kmesh=(9, 9, 9),
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
        json_path = Path("data/endmembers") / name / "tb2jflow" / f"{name}-tb2j.json"
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
        job_id = manager.submit_command(
            command=f"python {__file__} --name {name} {'--rerun' if force else ''}",
            config=config,
            conda_env="htvasp",
            work_dir=".",
        )
        if not job_id:
            log.error(f"Failed to submit job for {name}")
        else:
            log.info(f"Submitted job for {name} with ID {job_id}")


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
            run_tick(name, force=args.rerun, dry_run=True)
    elif args.name:
        run_tick(args.name, force=args.rerun)
    elif args.batch:
        run_batch(force=args.rerun)
    elif args.slurm:
        submit_jobs(force=args.rerun)
    else:
        parser.print_help()
