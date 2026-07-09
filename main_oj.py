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
from htvasp.oj.supercell_suggest import suggest_supercell
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


def suggest_supercell(
    struct,
    dist_cutoff: float = 6.0,
    atoms_per_cell: int | None = None,
) -> tuple[tuple[int, int, int], str]:
    """Suggest supercell size for OstravaJ based on structure + dist_cutoff.

    Physics:
      OJ needs the supercell to contain each J-pair's distance vector fully within
      the box.  The Nyquist-like criterion is:
          n_i * a_i  >  2 * d_k    for each direction i and each J-pair k.
      The tightest single constraint uses the shortest lattice vector:
          n_i  >  int(2 * d_max / a_min) + 1

      However, the *total number of atoms* in the supercell also matters.
      OJ searches over 2^{N_atoms} spin configurations and picks only those
      that satisfy the "identical magnetic surrounding" criterion (an
      exponentially hard constraint).  Empirically (from our dry-runs):

          BCC 2-atom:  sum(n_i)  ≥ 10  → enough independent configs
                        sum(n_i)   = 9  →  failure (3×3×3 = 54 atoms)

          HCP 2-atom:  (4,4,2) works,  sum = 10
          FCC 4-atom:  (2,2,2) → 32 atoms  should be tried first
                       (4,4,2) → 128 atoms → OJ search space too large

    Returns:
      (nx, ny, nz), reason  —  recommended supercell and a human-readable
                                justification or warning.
    """
    import math
    from pymatgen.core import Structure

    if atoms_per_cell is None:
        atoms_per_cell = len(struct)

    # Lattice vector lengths
    a, b, c = struct.lattice.abc
    a_min = min(a, b, c)

    # Minimum cells per direction to contain 2 × dist_cutoff
    n_req = int(math.ceil(2.0 * dist_cutoff / a_min)) + 1
    nx, ny, nz = max(n_req, 2), max(n_req, 2), max(n_req, 2)

    # Total atoms in supercell
    n_total = atoms_per_cell * nx * ny * nz
    search_space = 2 ** n_total  # 2^N possible spin configurations

    reasons = []

    # --- Adjust: search space must not explode ---
    # OJ's parallel generator can handle ~2^{64} in reasonable time.
    # >2^{100} is almost certainly too large.
    if search_space > 2 ** 100:
        # Try reducing the non-optimal directions first
        if nx == ny and nx > nz:
            # Cubic or near-cubic  →  prioritise the shortest direction for reduction
            pass

        # Strategy: keep the primary directions full, shrink less-important ones
        while search_space > 2 ** 80 and nx > 2 and ny > 2 and nz > 1:
            if nz > 1:
                nz -= 1
            elif ny > nx and ny > 2:
                ny -= 1
            elif nx > ny and nx > 2:
                nx -= 1
            else:
                ny -= 1
            n_total = atoms_per_cell * nx * ny * nz
            search_space = 2 ** n_total

        reasons.append(
            f"search space {math.log10(search_space):.0f} decibans "
            f"(trying small supercell to keep OJ tractable)"
        )

    # --- Adjust: sum(n_i) heuristic for independent configs ---
    # Empirically, BCC 2-atom needs sum(n_i) ≥ 10
    n_sum = nx + ny + nz
    while n_sum < 10 and nx < 6 and ny < 6 and nz < 6:
        # Grow the smallest dimension
        if nz < nx and nz < ny:
            nz += 1
        elif ny < nx:
            ny += 1
        else:
            nx += 1
        n_sum = nx + ny + nz
        n_total = atoms_per_cell * nx * ny * nz

    if n_sum >= 10:
        reasons.append(
            f"supercell sum({nx},{ny},{nz})={n_sum} ≥ 10  "
            f"(empirical threshold for enough independent configs)"
        )
    else:
        reasons.append(
            f"⚠  sum({nx},{ny},{nz})={n_sum} < 10 — may fail to find "
            f"enough independent magnetic configurations"
        )

    # --- Verify supercell contains full J-pair distances ---
    # OJ needs: for each J-pair vector d_k, its image in the supercell is unique.
    # A sufficient condition: n_i * a_i > max_k |d_k| for each direction i.
    # The maximum J-pair distance is dist_cutoff.
    ok = True
    for n_i, a_i in [(nx, a), (ny, b), (nz, c)]:
        if n_i * a_i <= dist_cutoff:
            ok = False
            break
    if ok:
        reasons.append(
            f"extent ≥ {dist_cutoff}Å in all directions  (covers all J-pairs)"
        )
    else:
        # Try growing the cell if it does NOT cover the distance
        if nx * a <= dist_cutoff:
            nx = int(ceil(dist_cutoff / a)) + 1
        if ny * b <= dist_cutoff:
            ny = int(ceil(dist_cutoff / b)) + 1
        if nz * c <= dist_cutoff:
            nz = int(ceil(dist_cutoff / c)) + 1
        n_total = atoms_per_cell * nx * ny * nz
        reasons.append(
            f"⚠ grew cell to cover J-pair distances → ({nx},{ny},{nz})"
        )

    # Final check: OJ parallelism uses multiprocessing with cpu_count()
    # Large supercell config generation is memory-intensive
    total_atoms = atoms_per_cell * nx * ny * nz
    reasons.append(f"{total_atoms} atoms, 2^{total_atoms} ≈ "
                   f"10^{total_atoms * math.log10(2):.0f} config space")

    return (nx, ny, nz), "; ".join(reasons)
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
    base_spin = detect_base_spin(magnetic_types)
    # If no magnetic species detected, OJ defaults to all atoms with base_spin=1.0
    if not magnetic_types:
        base_spin = [1.0]
    log.info(f"  Magnetic ions: {magnetic_types}, base_spin: {base_spin}")
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
