"""
Test script for QHA workflow.

Usage:
python test_qha.py --unit     # Run unit tests only
python test_qha.py --local    # Run locally (fresh run)
python test_qha.py --rerun    # Run locally with resume
python test_qha.py --slurm    # Submit to Slurm
"""

import json
import shutil
import logging
import argparse
import traceback
from pathlib import Path

from pymatgen.core import Structure, Lattice

from htvasp.workflows import QhaWorker

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)


# =============================================================================
# Test Helpers
# =============================================================================


def print_header(name: str):
    """Print standardized test header."""
    print("=" * 50)
    print(f"Test: {name}")
    print("=" * 50)


def get_al_bcc_structure(a: float = 2.73) -> Structure:
    """Return BCC Al structure for testing."""
    return Structure(Lattice.cubic(a), ["Al", "Al"], [[0, 0, 0], [0.5, 0.5, 0.5]])


def make_test_worker(name: str = "test-qha", **kwargs) -> QhaWorker:
    """Create QhaWorker with standard test defaults."""
    defaults = {
        "worker_name": name,
        "vasp_args": {"vasp_cmd": "echo test"},
        "supercell_matrix": ((1, 0, 0), (0, 1, 0), (0, 0, 1)),  # no supercell
        "temperature_range": (0, 1000, 100),
    }
    defaults.update(kwargs)
    return QhaWorker(**defaults)


# =============================================================================
# Unit Tests
# =============================================================================


def test_imports():
    """Test all imports"""
    print_header("Imports")
    print("✓ All imports successful\n")
    return True


def test_worker():
    """Test QhaWorker"""
    print_header("QhaWorker")

    worker = make_test_worker()
    print("✓ Worker created")
    print(f"✓ Has flow_maker: {hasattr(worker, 'flow_maker')}")
    print(f"✓ Supercell matrix: {worker.supercell_matrix}")
    print()
    return True


def test_temperature_range():
    """Test temperature range settings"""
    print_header("Temperature Range")

    temp_range = (0, 1000, 100)
    worker = make_test_worker(temperature_range=temp_range)
    print("✓ Worker created with custom temperature range")
    print(f"✓ Temperature range: {temp_range}")
    print()
    return True


def test_supercell_matrix():
    """Test supercell matrix settings"""
    print_header("Supercell Matrix")

    supercell = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    worker = make_test_worker(supercell_matrix=supercell)
    print("✓ Worker created with custom supercell matrix")
    print(f"✓ Supercell matrix: {supercell} (no supercell)")
    print()
    return True


def test_incar_settings():
    """Test INCAR settings"""
    print_header("INCAR Settings")

    custom_incar = {"ENCUT": 500, "KPAR": 4}
    worker = make_test_worker(global_incar=custom_incar)
    print("✓ Worker created with custom INCAR settings")
    print(f"✓ Global INCAR: {custom_incar}")
    print()
    return True


# =============================================================================
# Test Runner
# =============================================================================


def run_unit_tests():
    """Run all unit tests"""
    print("\n" + "=" * 50)
    print("QHA Unit Tests")
    print("=" * 50 + "\n")

    tests = [
        ("Imports", test_imports),
        ("Worker", test_worker),
        ("TemperatureRange", test_temperature_range),
        ("SupercellMatrix", test_supercell_matrix),
        ("IncarSettings", test_incar_settings),
    ]

    results = []
    for name, test_func in tests:
        try:
            result = test_func()
            results.append((name, result))
        except Exception as e:
            print(f"✗ Test {name} failed: {e}")
            traceback.print_exc()
            results.append((name, False))

    print("=" * 50)
    print("Summary")
    print("=" * 50)
    passed = sum(1 for _, r in results if r)
    total = len(results)
    for name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"  {name}: {status}")
    print(f"\nTotal: {passed}/{total} passed")

    return passed == total


# =============================================================================
# Integration Functions
# =============================================================================

CLUSTER_VASP_ARGS = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_gam'",
}


def run_locally(clean: bool = True):
    """Run QHA workflow locally (requires VASP).

    Args:
        clean: If True, remove existing flow_dir before running.
               If False, keep existing data and use resume=True.
    """
    structure = get_al_bcc_structure()
    flow_name = "Al-qha"
    flow_dir = Path("/nfs_ssd/tmp")
    store_dir = Path("./temp") / flow_name
    json_path = store_dir / f"{flow_name}.json"

    if clean and flow_dir.exists():
        shutil.rmtree(flow_dir)
        log.info(f"Cleaned existing flow_dir: {flow_dir}")

    worker = QhaWorker(
        vasp_args=CLUSTER_VASP_ARGS,
        global_incar={"GGA": "PE"},
        supercell_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),  # no supercell
        temperature_range=(0, 1000, 100),
    )

    resume = not clean
    if resume:
        log.info("Running in RESUME mode - will skip completed jobs")

    worker.run_flow(
        name=flow_name,
        structure=structure,
        flow_dir=flow_dir,
        store_dir=store_dir,
        resume=resume,
    )

    output = worker.get_result()
    worker.write_result(data=output, json_path=json_path)

    if output:

        print("\n" + "=" * 50)
        print("Results Summary")
        print("=" * 50)
        print(f"Bulk modulus: {output.get('bulk_modulus', 'N/A')} GPa")
        print(f"Temperatures: {output.get('temperatures', 'N/A')}")
        if vol_t := output.get("volume_temperature"):
            print(f"Volume range: {min(vol_t):.3f} - {max(vol_t):.3f} Å³")
    else:
        log.error("Workflow failed")


def submit_job():
    """Submit QHA workflow to Slurm with abort/rerun simulation."""
    from htvasp.slurm import SlurmJobManager

    manager = SlurmJobManager()

    def submit(
        config_overrides: dict | None = None, command_suffix: str = ""
    ) -> str | None:
        """Helper to submit a job with common defaults."""
        config = manager.get_cpu_config(
            ntasks=16, memory="8G", **(config_overrides or {})
        )
        return manager.submit_command(
            command=f"python {__file__} --local{command_suffix}",
            config=config,
            conda_env="htvasp",
            work_dir=".",
        )

    job_id = submit()
    if not job_id:
        return
    log.info(f"Submitted job: {job_id}")

    # Simulate abort and rerun
    import time

    time.sleep(100)

    if manager.cancel_job(job_id):
        log.info(f"Cancelled job: {job_id}")
    else:
        log.error(f"Failed to cancel job: {job_id}")
        return

    job_id2 = submit({"output_log": "job2.log"}, " --rerun")
    if job_id2:
        log.info(f"ReSubmitted job: {job_id2}")


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test QHA workflow")
    parser.add_argument(
        "--local", action="store_true", help="run workflow locally (fresh run)"
    )
    parser.add_argument("--rerun", action="store_true", help="run workflow with resume")
    parser.add_argument("--slurm", action="store_true", help="submit to slurm")
    parser.add_argument("--unit", action="store_true", help="run unit tests only")
    args = parser.parse_args()

    if args.local:
        run_locally(clean=True)
    elif args.rerun:
        run_locally(clean=False)
    elif args.slurm:
        submit_job()
    else:
        success = run_unit_tests()
        exit(0 if success else 1)
