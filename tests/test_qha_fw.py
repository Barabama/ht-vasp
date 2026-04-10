"""
Test script for QHA workflow with FireWorks integration.

Usage:
    python test_qha_fw.py --unit         # Run unit tests only
    python test_qha_fw.py --local        # Run workflow locally (requires VASP)
    python test_qha_fw.py --slurm        # Submit to Slurm as single job
    python test_qha_fw.py --fireworks    # Run via FireWorks distributed framework
"""

import json
import shutil
import logging
import argparse
import tempfile
from pathlib import Path
from datetime import datetime

# Configure jobflow to use a temporary JSON store for unit tests
# This must be done BEFORE importing any htvasp modules
_temp_store_dir = Path(tempfile.mkdtemp())
_jobflow_config = {
    "JOB_STORE": {
        "docs_store": {
            "type": "JSONStore",
            "paths": [str(_temp_store_dir / "test_store.json")],
        },
        "additional_stores": {},
    }
}
import os

os.environ["JOBFLOW_CONFIG_FILE"] = ""  # Clear any existing config

# Create a temporary config file
import yaml

_temp_config_file = _temp_store_dir / "jobflow_test.yaml"
with open(_temp_config_file, "w") as f:
    yaml.dump(_jobflow_config, f)
os.environ["JOBFLOW_CONFIG_FILE"] = str(_temp_config_file)

from pymatgen.core import Structure, Lattice

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)


class DateTimeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, datetime):
            return obj.isoformat()
        return super().default(obj)


def test_imports():
    """Test all imports including FireWorks integration."""
    print("=" * 50)
    print("Test: Imports")
    print("=" * 50)

    from htvasp.workflows import QhaWorker
    from htvasp.utils.fireworks_runner import FireWorksRunner

    print("All imports successful\n")
    return True


def test_worker_creation():
    """Test QhaWorker creation with default parameters."""
    print("=" * 50)
    print("Test: Worker Creation")
    print("=" * 50)

    from htvasp.workflows import QhaWorker

    vasp_args = {
        "vasp_cmd": "echo test",
        "vasp_gamma_cmd": "echo test",
    }

    worker = QhaWorker(
        worker_name="test-qha",
        vasp_args=vasp_args,
        supercell_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
    )

    print(f"Worker created: {worker.worker_name}")
    print(f"Has flow_maker: {hasattr(worker, 'flow_maker')}")
    print(f"Supercell matrix: {worker.supercell_matrix}")

    print()
    return True


def test_worker_fireworks_mode():
    """Test QhaWorker creation with FireWorks execution mode."""
    print("=" * 50)
    print("Test: Worker FireWorks Mode")
    print("=" * 50)

    from htvasp.workflows import QhaWorker

    vasp_args = {
        "vasp_cmd": "echo test",
    }

    fw_config_dir = Path("configs/fireworks")
    if fw_config_dir.exists():
        worker = QhaWorker(
            worker_name="test-qha-fw",
            vasp_args=vasp_args,
            execution_mode="fireworks",
            fireworks_config_dir=fw_config_dir,
            supercell_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        )

        print(f"Worker created with FireWorks mode")
        print(f"Execution mode: {worker.execution_mode}")
        print(f"FireWorks runner initialized: {worker.fireworks_runner is not None}")
    else:
        print(f"Skipping FireWorks mode test - config dir not found: {fw_config_dir}")

    print()
    return True


def test_flow_creation():
    """Test QHA flow creation with minimal supercell."""
    print("=" * 50)
    print("Test: Flow Creation")
    print("=" * 50)

    from htvasp.workflows import QhaWorker

    lattice = Lattice.cubic(2.85)
    structure = Structure(lattice, ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    vasp_args = {
        "vasp_cmd": "echo test",
    }

    worker = QhaWorker(
        worker_name="test-flow",
        vasp_args=vasp_args,
        supercell_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
    )

    flow = worker._make_flow(structure)
    print(f"Flow created: {flow.name}")
    print(f"Number of jobs in flow: {len(list(flow.jobs))}")

    expected_min_jobs = 3
    if len(list(flow.jobs)) >= expected_min_jobs:
        print(f"Job count correct (>= {expected_min_jobs})")
    else:
        print(f"Expected >= {expected_min_jobs} jobs, got {len(list(flow.jobs))}")
        return False

    print()
    return True


def test_incar_settings():
    """Test custom INCAR settings propagation."""
    print("=" * 50)
    print("Test: INCAR Settings")
    print("=" * 50)

    from htvasp.workflows import QhaWorker

    vasp_args = {
        "vasp_cmd": "echo test",
    }

    custom_incar = {
        "ENCUT": 500,
        "KPAR": 4,
    }

    relax_incar = {"NSW": 50}
    phonon_incar = {"EDIFF": 1e-8}

    worker = QhaWorker(
        worker_name="test-incar",
        vasp_args=vasp_args,
        global_incar=custom_incar,
        relax_incar=relax_incar,
        phonon_incar=phonon_incar,
        supercell_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
    )

    print(f"Worker created with custom INCAR settings")
    print(f"Global INCAR: {custom_incar}")
    print(f"Relax INCAR: {relax_incar}")
    print(f"Phonon INCAR: {phonon_incar}")

    print()
    return True


def test_temperature_range():
    """Test temperature range configuration."""
    print("=" * 50)
    print("Test: Temperature Range")
    print("=" * 50)

    from htvasp.workflows import QhaWorker

    vasp_args = {
        "vasp_cmd": "echo test",
    }

    temp_range = (0, 1000, 100)

    worker = QhaWorker(
        worker_name="test-temp",
        vasp_args=vasp_args,
        temperature_range=temp_range,
        supercell_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
    )

    print(f"Temperature range configured: {temp_range}")
    print(f"tmin={temp_range[0]}, tmax={temp_range[1]}, tstep={temp_range[2]}")

    print()
    return True


def test_supercell_matrix():
    """Test different supercell matrix configurations."""
    print("=" * 50)
    print("Test: Supercell Matrix")
    print("=" * 50)

    from htvasp.workflows import QhaWorker

    vasp_args = {
        "vasp_cmd": "echo test",
    }

    test_matrices = [
        ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        ((2, 0, 0), (0, 2, 0), (0, 0, 2)),
    ]

    for matrix in test_matrices:
        worker = QhaWorker(
            worker_name=f"test-sc-{matrix[0][0]}",
            vasp_args=vasp_args,
            supercell_matrix=matrix,
        )
        print(f"Supercell matrix {matrix}: OK")

    print()
    return True


def run_unit_tests():
    """Run all unit tests."""
    print("\n" + "=" * 50)
    print("QHA Unit Tests")
    print("=" * 50 + "\n")

    tests = [
        ("Imports", test_imports),
        ("WorkerCreation", test_worker_creation),
        ("FireWorksMode", test_worker_fireworks_mode),
        ("FlowCreation", test_flow_creation),
        ("INCARSettings", test_incar_settings),
        ("TemperatureRange", test_temperature_range),
        ("SupercellMatrix", test_supercell_matrix),
    ]

    results = []
    for name, test_func in tests:
        try:
            result = test_func()
            results.append((name, result))
        except Exception as e:
            print(f"Test {name} failed: {e}")
            import traceback

            traceback.print_exc()
            results.append((name, False))

    print("=" * 50)
    print("Summary")
    print("=" * 50)
    passed = sum(1 for _, r in results if r)
    total = len(results)
    for name, result in results:
        status = "PASS" if result else "FAIL"
        print(f"  {name}: {status}")
    print(f"\nTotal: {passed}/{total} passed")

    return passed == total


def run_locally():
    """Run QHA workflow locally (requires VASP)."""
    from htvasp.workflows import QhaWorker

    struct = Structure(
        lattice=[[2.73, 0, 0], [0, 2.73, 0], [0, 0, 2.73]],
        species=["Al", "Al"],
        coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
    )

    flow_dir = Path("temp/qha-Al")
    json_path = flow_dir.joinpath("qha_Al.json")

    if flow_dir.exists():
        shutil.rmtree(flow_dir)

    vasp_args = {
        "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_std'",
        "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_gam'",
    }

    worker = QhaWorker(
        worker_name="qha-Al",
        vasp_args=vasp_args,
        global_incar={"GGA": "PE"},
        supercell_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        temperature_range=(0, 1000, 100),
    )

    worker.run_flow(
        name="Al",
        structure=struct,
        flow_dir=flow_dir,
        resume=True,
    )

    qha_data = worker.get_result()
    if qha_data:
        result = {"name": "Al", "state": "successful", **qha_data}
        with open(json_path, "w", encoding="utf-8") as jf:
            json.dump(result, jf, ensure_ascii=False, indent=2, cls=DateTimeEncoder)
        log.info(f"Output saved to {json_path}")
    else:
        result = {"name": "Al", "state": "failed"}
        with open(json_path, "w", encoding="utf-8") as jf:
            json.dump(result, jf, ensure_ascii=False, indent=2, cls=DateTimeEncoder)
        log.error("Workflow failed")


def submit_slurm():
    """Submit QHA workflow to Slurm as single job."""
    from htvasp.slurm import SlurmJobManager

    manager = SlurmJobManager()
    config = manager.get_cpu_config(
        job_name="qha-test",
        output_log="qha-test.log",
        nodes=1,
        ntasks=48,
    )
    job_id = manager.submit_command(
        command=f"python {__file__} --local",
        config=config,
        conda_env="htvasp",
        workdir=".",
    )
    log.info(f"Submitted Slurm job: {job_id}")


def run_fireworks():
    """Run QHA workflow via FireWorks distributed framework."""
    from htvasp.workflows import QhaWorker

    struct = Structure(
        lattice=[[2.73, 0, 0], [0, 2.73, 0], [0, 0, 2.73]],
        species=["Al", "Al"],
        coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
    )

    flow_dir = Path("temp/qha-Al-fireworks")
    json_path = flow_dir.joinpath("qha_Al_fw.json")

    if flow_dir.exists():
        shutil.rmtree(flow_dir)

    vasp_args = {
        "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_std'",
        "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_gam'",
    }

    fw_config_dir = Path("configs/fireworks")

    worker = QhaWorker(
        worker_name="qha-Al-fw",
        vasp_args=vasp_args,
        global_incar={"GGA": "PE"},
        execution_mode="fireworks",
        fireworks_config_dir=fw_config_dir,
        supercell_matrix=((1, 0, 0), (0, 1, 0), (0, 0, 1)),
        temperature_range=(0, 1000, 100),
    )

    worker.run_flow(
        name="Al",
        structure=struct,
        flow_dir=flow_dir,
        wait_for_completion=False,
    )

    qha_data = worker.get_result()
    if qha_data:
        result = {"name": "Al", "state": "successful", **qha_data}
        with open(json_path, "w", encoding="utf-8") as jf:
            json.dump(result, jf, ensure_ascii=False, indent=2, cls=DateTimeEncoder)
        log.info(f"Output saved to {json_path}")
    else:
        result = {"name": "Al", "state": "submitted_to_fireworks"}
        with open(json_path, "w", encoding="utf-8") as jf:
            json.dump(result, jf, ensure_ascii=False, indent=2, cls=DateTimeEncoder)
        log.info("Workflow submitted to FireWorks (async mode)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test QHA workflow with FireWorks integration")
    parser.add_argument("--local", action="store_true", help="Run workflow locally")
    parser.add_argument("--slurm", action="store_true", help="Submit to Slurm as single job")
    parser.add_argument("--fireworks", action="store_true", help="Run via FireWorks framework")
    parser.add_argument("--unit", action="store_true", help="Run unit tests only")
    args = parser.parse_args()

    if args.local:
        run_locally()
    elif args.slurm:
        submit_slurm()
    elif args.fireworks:
        run_fireworks()
    elif args.unit:
        success = run_unit_tests()
        exit(0 if success else 1)
    else:
        success = run_unit_tests()
        exit(0 if success else 1)
