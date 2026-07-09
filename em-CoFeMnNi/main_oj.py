#!/usr/bin/env python3
"""HT-VASP OJ Workflow - Magnetic exchange calculation using OstravaJ.

Usage:
    python main_oj.py --oj SER-Fe                 # single structure oj
    python main_oj.py --r8 SER-Fe                 # single structure isif8
    python main_oj.py --dry-run SER-Fe            # dry-run: just verify config generation
    python main_oj.py --batch                     # all structures locally
    python main_oj.py --slurm                     # submit to Slurm
"""

import json
import logging
import shutil
import argparse
from pathlib import Path
from pymatgen.core import Structure
from htvasp.model import Endmember
from htvasp.workflows import OJWorker, R8Worker
from htvasp.slurm import SlurmJobManager

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

# VASP configuration
VASP_ARGS = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_gam'",
}

# Static INCAR for magnetic exchange calculations
# OstravaJ requires ISPIN=2 WITHOUT MAGMOM in INCAR (it handles it)
# MUST be static (IBRION=-1) — relaxation breaks the Heisenberg model
GLOBAL_INCAR = {
    "ENCUT": 400,
    "ISTART": 0,
    "ICHARG": 2,
    "ISMEAR": 1,
    "SIGMA": 0.2,
    "ALGO": "Normal",           # "Fast" unstable for magnetic systems
    "NELM": 120,
    "NELMIN": 6,
    "NELMDL": -6,
    "IBRION": -1,               # CRITICAL: static, not relaxation
    "NSW": 0,                   # CRITICAL: no ionic steps
    "EDIFF": 1e-6,
    "ISPIN": 2,
    "ISYM": 0,
    "LREAL": "Auto",
    "PREC": "Normal",
    "SYMPREC": 1e-5,
    "LWAVE": False,
    "LCHARG": False,
    "LORBIT": 10,               # required by OstravaJ for reading vasprun.xml
    "GGA": "PE",
    # Mixing parameters tuned for magnetic systems
    "AMIX": 0.1,
    "BMIX": 1e-4,
    "AMIX_MAG": 0.4,
    "BMIX_MAG": 1e-4,
    "KPAR": 2,
    "NCORE": 4,
}
RELAX_INCAR = {
    "ISIF": 8,
    "NELM": 100,
    "NSW": 100,
    "EDIFF": 1e-6,
    "EDIFFG": -0.02,
    "KPAR": 1,
    "NCORE": 1,
}
# Structure names for CoFeMnNi magnetic exchange
# B2/endmember naming: BCC-X-Y means BCC structure with X on sublattice 1, Y on sublattice 2
STRUCTURE_NAMES = [
    # Pure elements (SER)
    "SER-Co",
    "SER-Fe",
    # "SER-Mn",
    # "SER-Ni",
    # BCC endmembers (4 elements = 10 binary pairs)
    # "BCC-Co-Co",
    # "BCC-Co-Fe",
    "BCC-Co-Mn",
    # "BCC-Co-Ni",
    # "BCC-Fe-Fe",
    # "BCC-Fe-Mn",
    "BCC-Fe-Ni",
    # "BCC-Mn-Mn",
    # "BCC-Mn-Ni",
    # "BCC-Ni-Ni",
    # FCC endmembers
    # "FCC-Co-Co",
    "FCC-Co-Fe",
    # "FCC-Co-Mn",
    # "FCC-Co-Ni",
    # "FCC-Fe-Fe",
    "FCC-Fe-Ni",
    # "FCC-Fe-Mn",
    # "FCC-Mn-Mn",
    # "FCC-Mn-Ni",
    # "FCC-Ni-Ni",
    # HCP endmembers
    # "HCP-Co-Co",
    # "HCP-Co-Fe",
    "HCP-Co-Mn",
    "HCP-Co-Ni",
]

# Magnetic species: set base_spin per element (μB)
# From OstravaJ paper / experimental values
MAGNETIC_SPECIES = ["Co", "Fe", "Mn", "Ni"]
BASE_SPIN_MAP = {
    "Co": 2.0,
    "Fe": 2.5,
    "Mn": 2.5,
    "Ni": 1.0,
}

root_dir = Path(__file__).parent.resolve()
flow_dir = Path("/nfs_ssd/tmp")
poscar_dir = root_dir / "data" / "poscars"
conda_env = root_dir.parent / ".conda"

def detect_magnetic_ions(stem: str) -> list[str]:
    """Detect which elements are magnetic from structure name.
    Preserves order from STRUCT_NAMES entry, deduplicates.
    """
    seen = set()
    ions = []
    parts = stem.replace("-", " ").split()
    for p in parts:
        if p in MAGNETIC_SPECIES and p not in seen:
            seen.add(p)
            ions.append(p)
    return ions


def detect_base_spin(magnetic_ions: list[str], default: float = 2.0) -> list[float]:
    """Return base_spin per magnetic ion, in the same order as magnetic_ions list.

    OJ's OJ.conf requires base_spin to be a space-separated list whose
    order matches magnetic_ion_types exactly, e.g.:
        magnetic_ion_types Fe Co
        base_spin 2.5 2.0
    """
    spins = []
    for el in magnetic_ions:
        spin = BASE_SPIN_MAP.get(el, default)
        spins.append(spin)
    return spins

def run_relax(name: str, rerun: bool = False):
    store_dir = root_dir / "data" / name / "ojflow"
    json_path = store_dir / f"{name}-r8.json"

    # Skip if already done
    if not rerun and json_path.exists():
        log.info(f"R8Worker for Structure {name} already done (use --rerun to re-run)")
        return 0
    log.info(f"R8Worker for Structure {name} started")

    endmember = Endmember()
    structure = endmember.get_poscar(name, poscar_dir)
    try:
        r8_worker = R8Worker(
            vasp_args=VASP_ARGS,
            global_incar=GLOBAL_INCAR,
            r8_incar=RELAX_INCAR,
        )
        r8_worker.run_flow(
            name=name,
            structure=structure,
            flow_dir=flow_dir,
            store_dir=store_dir,
            resume=not rerun,
        )
        output = r8_worker.get_result()
        r8_worker.write_result(data=output, json_path=json_path)
        log.info(str(output.get))
    except Exception as e:
        log.error(f"R8Worker for Structure {name} failed")
        log.error(str(e))
        return 1
    log.info(f"R8Worker for Structure {name} done")
    return 0

def run_oj(name: str, rerun: bool = False, dry_run: bool = False):
    """Run a single structure locally."""
    store_dir = root_dir / "data" / name / "ojflow"
    json_path = store_dir / f"{name}-oj.json"

    # Skip if already done
    if not rerun and json_path.exists():
        log.info(f"OJWorker for Structure {name} already done (use --rerun to re-run)")
        return 0

    log.info(f"OJWorker for Structure {name} started")

    endmember = Endmember()
    structure = endmember.get_poscar(name, poscar_dir)

    # Detect magnetic ions and base_spin from structure name
    magnetic_types = detect_magnetic_ions(name)
    base_spin = detect_base_spin(magnetic_types)
    # If no magnetic species detected, OJ defaults to all atoms with base_spin=1.0
    if not magnetic_types:
        base_spin = [1.0]
    log.info(f"  Magnetic ions: {magnetic_types}, base_spin: {base_spin}")
    log.info(f"  Magnetic ions: {magnetic_types}, base_spin: {base_spin}")
    if dry_run:
        _dry_run_oj(name, structure, magnetic_types, base_spin)
        return 0

    try:
        worker = OJWorker(
            vasp_args=VASP_ARGS,
            global_incar=GLOBAL_INCAR,
            dist_cutoff=6.0,                # ~6 nearest neighbors in BCC
            extend_poscar=(4, 4, 2),         # 64 atoms for BCC 2-atom cell
            magnetic_ion_types=magnetic_types,
            base_spin=base_spin,
        )
        worker.run_flow(
            name=name,
            structure=structure,
            flow_dir=flow_dir,
            store_dir=store_dir,
            resume=not rerun,
        )
        output = worker.get_result()
        worker.write_result(data=output, json_path=json_path)
        _report_result(output)
    except Exception as e:
        log.error(f"OJWorker for Structure {name} failed")
        log.error(e)
        return 1
    log.info(f"OJWorker for Structure {name} done")
    return 0


def _dry_run_oj(name: str, structure, magnetic_types: list[str], base_spin: float | list[float]):
    """Run OstravaJ generate with --dry-run to test configuration."""
    from htvasp.oj.input_set import OJInputSetGenerator, write_oj_input_set
    import subprocess
    from htvasp.oj.jobs import _get_oj_script_path

    tmp_dir = Path(f"/tmp/oj_dryrun_{name}")
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True)

    input_gen = OJInputSetGenerator(
        dist_cutoff=6.0,
        extend_poscar=(4, 4, 2),
        magnetic_ion_types=magnetic_types,
        base_spin=base_spin,
    )

    write_oj_input_set(structure, tmp_dir, input_gen, include_oj_conf=True)

    # Print OJ.conf
    log.info(f"  OJ.conf:\n{input_gen.get_oj_conf_string()}")

    # Run dry-run
    cmd = [str(_get_oj_script_path()), "generate", "-d", "-i", str(tmp_dir), "-r", str(tmp_dir)]
    log.info(f"  Run: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

    if result.returncode == 0:
        log.info(f"  ✅ Dry-run SUCCESS")
        # OJ prints system info to stdout; config count may appear there
        for line in result.stdout.split("\n"):
            if any(kw in line for kw in ["Found", "Success", "matrix", "config", "J=", "J("]):
                log.info(f"  {line.strip()}")
    else:
        log.warning(f"  ❌ Dry-run FAILED")
        for line in result.stderr.split("\n"):
            log.warning(f"  {line.strip()}")

    # Cleanup
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)


def _report_result(output: dict | None):
    """Print a short summary of the OJ result."""
    if not output:
        log.warning("  No output returned")
        return
    if "error" in output:
        log.warning(f"  Error: {output['error']}")
        return

    log.info(f"  Js found: {output.get('J_reprs', 'N/A')}")
    log.info(f"  Tc_MFA: {output.get('Tc_MFA', 'N/A'):.1f} K" if output.get('Tc_MFA') else "  Tc_MFA: N/A")
    log.info(f"  Tc_RPA: {output.get('Tc_RPA', 'N/A'):.1f} K" if output.get('Tc_RPA') else "  Tc_RPA: N/A")
    log.info(f"  Num configs: {output.get('num_configs', 'N/A')}")
    log.info(f"  Condition number: {output.get('condition_number', 'N/A')}")
    log.info(f"  Rank: {output.get('rank', 'N/A')}")
    log.info(f"  Magnetic moment (avg): {output.get('avg_magnetic_moment', 'N/A')} μB")


def run_batch(rerun: bool = False, dry_run: bool = False):
    """Run all structures locally."""
    for name in STRUCTURE_NAMES:
        run_oj(name, rerun=rerun, dry_run=dry_run)


def submit_jobs(rerun: bool = False) -> None:
    manager = SlurmJobManager()
    for name in STRUCTURE_NAMES:
        jid = manager.submit_command(
            command=f"python {__file__} --relax {name} {'--rerun' if rerun else ''}",
            config=manager.get_cpu_config(
                job_name=f"{name}-oj",
                output_log=f"logs/{name}-oj.log",
                ntasks=8,
                memory="8G",
            ),
            conda_env=str(conda_env),
            work_dir=root_dir,
        )
        if not jid:
            log.error(f"Failed to submit oj job for Structure {name}")
        else:
            log.info(f"Submitted oj job for Structure {name} with ID {jid}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="OJ magnetic exchange workflow")
    parser.add_argument("--oj", type=str, help="Run single structure")
    parser.add_argument("--relax", type=str, help="Run single structure")
    parser.add_argument("--dry-run", type=str, help="Dry-run: test config generation only")
    parser.add_argument("--batch", action="store_true", help="Run all locally")
    parser.add_argument("--dry-batch", action="store_true", help="Dry-run all structures")
    parser.add_argument("--slurm", action="store_true", help="Submit to Slurm")
    parser.add_argument("--rerun", action="store_true", help="rerun re-run")
    args = parser.parse_args()

    if args.oj:
        run_oj(args.oj, rerun=args.rerun)
    elif args.relax:
        run_relax(args.relax, rerun=args.rerun)
    elif args.dry_run:
        run_oj(args.dry_run, dry_run=True)
    elif args.batch:
        run_batch(rerun=args.rerun)
    elif args.dry_batch:
        run_batch(dry_run=True)
    elif args.slurm:
        submit_jobs(rerun=args.rerun)
    else:
        parser.print_help()
