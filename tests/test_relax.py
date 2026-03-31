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
from pathlib import Path
from datetime import datetime

from pymatgen.core import Structure, Lattice

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)


class DateTimeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, datetime):
            return obj.isoformat()
        return super().default(obj)


def test_imports():
    """Test all imports"""
    print("=" * 50)
    print("Test: Imports")
    print("=" * 50)

    from htvasp.workflows import RelaxWorker

    print("✓ All imports successful\n")
    return True


def test_worker():
    """Test RelaxWorker"""
    print("=" * 50)
    print("Test: RelaxWorker")
    print("=" * 50)

    from htvasp.workflows import RelaxWorker

    vasp_args = {
        "vasp_cmd": "echo test",
        "vasp_gamma_cmd": "echo test",
    }

    worker = RelaxWorker(
        worker_name="test-relax",
        vasp_args=vasp_args,
    )

    print(f"✓ Worker created")
    print(f"✓ Has relax_flow: {hasattr(worker, 'relax_flow')}")

    lattice = Lattice.cubic(2.85)
    structure = Structure(lattice, ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    print()
    return True


def test_incar_settings():
    """Test INCAR settings"""
    print("=" * 50)
    print("Test: INCAR Settings")
    print("=" * 50)

    from htvasp.workflows import RelaxWorker

    vasp_args = {
        "vasp_cmd": "echo test",
    }

    custom_incar = {
        "ENCUT": 500,
        "KPAR": 4,
    }

    worker = RelaxWorker(
        worker_name="test-relax",
        vasp_args=vasp_args,
        global_incar=custom_incar,
    )

    print(f"✓ Worker created with custom INCAR settings")
    print(f"✓ Custom settings: {custom_incar}")

    print()
    return True


def test_potcar_functional():
    """Test POTCAR functional selection"""
    print("=" * 50)
    print("Test: POTCAR Functional")
    print("=" * 50)

    from htvasp.workflows import RelaxWorker

    vasp_args = {
        "vasp_cmd": "echo test",
    }

    worker = RelaxWorker(
        worker_name="test-relax",
        vasp_args=vasp_args,
        potcar_functional="PBE_64",
    )

    print(f"✓ Worker created with PBE_64 functional")

    print()
    return True


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
            import traceback
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


def run_locally():
    """Run Relax workflow locally (requires VASP)"""
    from htvasp.workflows import RelaxWorker

    struct = Structure(
        lattice=[[2.73, 0, 0], [0, 2.73, 0], [0, 0, 2.73]],
        species=["Al", "Al"],
        coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
    )

    flow_dir = Path("temp/relax-Al")
    json_path = flow_dir.joinpath("relax_Al.json")

    # if flow_dir.exists():
    #     shutil.rmtree(flow_dir)

    vasp_args = {
        "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_std'",
        "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_gam'",
    }

    worker = RelaxWorker(
        worker_name="relax-Al",
        vasp_args=vasp_args,
        global_incar={
            "GGA": "PE",
        },
    )

    output = worker.run_flow(
        name="Al",
        structure=struct,
        flow_dir=flow_dir,
    )

    if output:
        with open(json_path, "w", encoding="utf-8") as jf:
            json.dump(output, jf, indent=2, cls=DateTimeEncoder)
        log.info(f"Output saved to {json_path}")
    else:
        log.error("Workflow failed")


def submit_job():
    """Submit Relax workflow to Slurm"""
    from htvasp.slurm import SlurmJobManager

    manager = SlurmJobManager()
    config = manager.get_cpu_config(
        ntasks=8,
        memory="4G",
    )
    job_id = manager.submit_command(
        command=f"python {__file__} --local",
        config=config,
        conda_env="htvasp",
        workdir=".",
    )
    log.info(f"Submitted job: {job_id}")


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
    elif args.unit:
        success = run_unit_tests()
        exit(0 if success else 1)
    else:
        success = run_unit_tests()
        exit(0 if success else 1)
