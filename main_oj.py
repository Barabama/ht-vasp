#!/usr/bin/env python3
"""HT-VASP OJ Workflow - Magnetic exchange calculation using OstravaJ.

Usage:
    python main_oj.py --tick SER-Fe              # single structure
    python main_oj.py --dry-run SER-Fe            # dry-run: just verify config generation
    python main_oj.py --batch                     # all structures locally
    python main_oj.py --slurm                     # submit to Slurm
"""

import json
import logging
import shutil
import argparse
from pathlib import Path

from htvasp.workflows import OJWorker
from htvasp.slurm import SlurmJobManager

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)


def get_structure(name: str) -> "Structure":
    """Read POSCAR directly from data/poscars/, bypassing MPRester dependency."""
    from pymatgen.core import Structure
    poscar_path = Path("data/poscars") / f"{name}.vasp"
    if poscar_path.exists():
        return Structure.from_file(poscar_path)
    raise FileNotFoundError(f"POSCAR not found: {poscar_path}")

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
    "ENCUT": 500,
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
    "EDIFF": 1e-7,
    "ISYM": 0,
    "LREAL": "Auto",
    "PREC": "Normal",
    "SYMPREC": 1e-5,
    "LWAVE": False,
    "LCHARG": False,
    "GGA": "PE",
    "KPAR": 2,
    "NCORE": 4,
    "LORBIT": 10,               # required by OstravaJ for reading vasprun.xml
    # Mixing parameters tuned for magnetic systems
    "AMIX": 0.1,
    "BMIX": 0.0001,
    "AMIX_MAG": 0.4,
    "BMIX_MAG": 0.0001,
}

# Structure names for CoFeMnNi magnetic exchange
# B2/endmember naming: BCC-X-Y means BCC structure with X on sublattice 1, Y on sublattice 2
STRUCT_NAMES = [
    # Pure elements (SER)
    "SER-Co",
    "SER-Fe",
    # "SER-Mn",
    # "SER-Ni",
    # BCC endmembers (4 elements = 10 binary pairs)
    # "BCC-Co-Co",
    # "BCC-Co-Fe",
    # "BCC-Co-Mn",
    # "BCC-Co-Ni",
    # "BCC-Fe-Fe",
    # "BCC-Fe-Mn",
    # "BCC-Fe-Ni",
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
]

# Magnetic species: set base_spin per element (μB)
# From OstravaJ paper / experimental values
MAGNETIC_SPECIES = ["Co", "Fe", "Mn", "Ni"]
BASE_SPIN_MAP = {
    "Co": 2.0,    # FCC Co ~1.7 μB, BCC ~1.5 μB
    "Fe": 2.5,    # BCC Fe ~2.2 μB, FCC ~2.5 μB
    "Mn": 2.5,    # ~2.4 μB
    "Ni": 1.0,    # FCC Ni ~0.6 μB
}


def detect_magnetic_ions(stem: str) -> list[str]:
    """Detect which elements are magnetic from structure name."""
    ions = set()
    parts = stem.replace("-", " ").split()
    for p in parts:
        if p in MAGNETIC_SPECIES:
            ions.add(p)
    return sorted(ions)


def detect_base_spin(stem: str, default: float = 2.0) -> float | list[float]:
    """Use highest base_spin among detected magnetic elements."""
    max_spin = default
    for el, spin in BASE_SPIN_MAP.items():
        if el in stem:
            max_spin = max(max_spin, spin)
    return max_spin


def run_tick(name: str, force: bool = False, dry_run: bool = False):
    """Run a single structure locally."""
    flow_dir = Path("/nfs_ssd/tmp")
    store_dir = Path("data/endmembers") / name / "ojflow"
    json_path = store_dir / f"{name}-oj.json"

    # Skip if already done
    if not force and json_path.exists():
        log.info(f"Structure {name} already done (use --force to re-run)")
        return

    log.info(f"Structure {name} start")

    if force and store_dir.exists():
        shutil.rmtree(store_dir)

    struct = get_structure(name)

    # Detect magnetic ions and base_spin from structure name
    magnetic_types = detect_magnetic_ions(name)
    base_spin = detect_base_spin(name)
    log.info(f"  Magnetic ions: {magnetic_types}, base_spin: {base_spin}")

    if dry_run:
        _dry_run_oj(name, struct, magnetic_types, base_spin)
        return

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
            structure=struct,
            flow_dir=flow_dir,
            store_dir=store_dir,
            resume=not force,
        )
        output = worker.get_result()
        worker.write_result(data=output, json_path=json_path)
        _report_result(output)
    except Exception as e:
        log.error(f"Structure {name} failed: {e}")


def _dry_run_oj(name: str, structure, magnetic_types: list[str], base_spin: float):
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


def run_batch(force: bool = False, dry_run: bool = False):
    """Run all structures locally."""
    for name in STRUCT_NAMES:
        run_tick(name, force=force, dry_run=dry_run)


def submit_jobs(force: bool = False) -> None:
    manager = SlurmJobManager()
    for name in STRUCT_NAMES:
        config = manager.get_cpu_config(
            job_name=f"{name}-oj",
            output_log=f"logs/{name}-oj.log",
            ntasks=32,
            memory="20G",
        )
        job_id = manager.submit_command(
            command=f"python {__file__} --tick {name} {'--force' if force else ''}",
            config=config,
            conda_env="htvasp",
            work_dir=".",
        )
        if not job_id:
            log.error(f"Failed to submit job for {name}")
        else:
            log.info(f"Submitted job for {name} with ID {job_id}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="OJ magnetic exchange workflow")
    parser.add_argument("--tick", type=str, help="Run single structure")
    parser.add_argument("--dry-run", type=str, help="Dry-run: test config generation only")
    parser.add_argument("--batch", action="store_true", help="Run all locally")
    parser.add_argument("--dry-batch", action="store_true", help="Dry-run all structures")
    parser.add_argument("--slurm", action="store_true", help="Submit to Slurm")
    parser.add_argument("--force", action="store_true", help="Force re-run")
    args = parser.parse_args()

    if args.tick:
        run_tick(args.tick, force=args.force)
    elif args.dry_run:
        run_tick(args.dry_run, dry_run=True)
    elif args.batch:
        run_batch(force=args.force)
    elif args.dry_batch:
        run_batch(dry_run=True)
    elif args.slurm:
        submit_jobs(force=args.force)
    else:
        parser.print_help()
