"""
Test script for Relax workflow.

Usage:
    python test_relax.py --local    # Run locally
    python test_relax.py --slurm    # Submit to Slurm
    python test_relax.py --unit     # Run unit tests only
"""

import json
import shutil
import logging
import argparse
import traceback
from pathlib import Path

from pymatgen.core import Structure, Lattice

from htvasp.workflows import RelaxWorker

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


def make_test_worker(name: str = "test-relax", **kwargs) -> RelaxWorker:
    """Create RelaxWorker with standard test defaults."""
    defaults = {
        "worker_name": name,
        "vasp_args": {"vasp_cmd": "echo test"},
    }
    defaults.update(kwargs)
    return RelaxWorker(**defaults)


# =============================================================================
# Unit Tests
# =============================================================================


def test_imports():
    """Test all imports"""
    print_header("Imports")
    print("✓ All imports successful\n")
    return True


def test_worker():
    """Test RelaxWorker"""
    print_header("RelaxWorker")

    worker = make_test_worker()
    print("✓ Worker created")
    print(f"✓ Has flow_maker: {hasattr(worker, 'flow_maker')}")
    print()
    return True


def test_incar_settings():
    """Test INCAR settings"""
    print_header("INCAR Settings")

    custom_incar = {"ENCUT": 520, "KPAR": 4}
    worker = make_test_worker(global_incar=custom_incar)

    print("✓ Worker created with custom INCAR settings")
    print(f"✓ Custom settings: {custom_incar}")
    print()
    return True


def test_potcar_functional():
    """Test POTCAR functional selection"""
    print_header("POTCAR Functional")

    worker = make_test_worker(potcar_functional="PBE_64")
    print("✓ Worker created with PBE_64 functional")
    print()
    return True


# =============================================================================
# Test Runner
# =============================================================================


def run_unit_tests():
    """Run all unit tests"""
    print("\n" + "=" * 50)
    print("Relax Unit Tests")
    print("=" * 50 + "\n")

    tests = [
        ("Imports", test_imports),
        ("Worker", test_worker),
        ("IncarSettings", test_incar_settings),
        ("PotcarFunctional", test_potcar_functional),
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


def run_locally():
    """Run Relax workflow locally (requires VASP)"""
    structure = get_al_bcc_structure()
    flow_name = "Al-relax"
    flow_dir = Path("/nfs_ssd/tmp")
    store_dir = Path("./temp") / flow_name
    json_path = store_dir / f"{flow_name}.json"

    if flow_dir.exists():
        shutil.rmtree(flow_dir)

    worker = RelaxWorker(
        vasp_args=CLUSTER_VASP_ARGS,
        global_incar={"GGA": "PE"},
    )

    worker.run_flow(
        name=flow_name,
        structure=structure,
        flow_dir=flow_dir,
        store_dir=store_dir,
    )

    output = worker.get_result()
    worker.write_result(data=output, json_path=json_path)


def submit_job():
    """Submit Relax workflow to Slurm"""
    from htvasp.slurm import SlurmJobManager

    manager = SlurmJobManager()
    config = manager.get_cpu_config(ntasks=8, memory="4G")
    job_id = manager.submit_command(
        command=f"python {__file__} --local",
        config=config,
        conda_env="htvasp",
        work_dir=".",
    )
    if job_id:
        log.info(f"Submitted job: {job_id}")


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test Relax workflow")
    parser.add_argument("--local", action="store_true", help="run workflow locally")
    parser.add_argument("--slurm", action="store_true", help="submit to slurm")
    parser.add_argument("--unit", action="store_true", help="run unit tests only")
    args = parser.parse_args()

    if args.local:
        run_locally()
    elif args.slurm:
        submit_job()
    else:
        success = run_unit_tests()
        exit(0 if success else 1)
