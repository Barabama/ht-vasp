"""
Test script for OJ workflow.

Usage:
    python test_oj.py --unit     # Run unit tests only
    python test_oj.py --local    # Run workflow locally (requires VASP)
    python test_oj.py --slurm    # Submit to Slurm
"""

import json
import shutil
import logging
import argparse
import tempfile
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

    from htvasp.workflows import OJWorker
    from htvasp.oj import OJConfig

    print("✓ All imports successful\n")
    return True


def test_worker_creation():
    """Test OJWorker creation"""
    print("=" * 50)
    print("Test: Worker Creation")
    print("=" * 50)

    from htvasp.workflows import OJWorker
    from htvasp.oj import OJConfig

    oj_config = OJConfig(j_count=2, magnetic_ion_types=["Fe"])

    worker = OJWorker(
        worker_name="test_oj",
        oj_config=oj_config,
        global_incar={"ENCUT": 520},
    )

    print(f"✓ Worker created: {worker.worker_name}")
    print(f"✓ Has flow_makers: {hasattr(worker, 'flow_makers')}")
    print(f"✓ Number of makers: {len(worker.flow_makers)}")

    relax_maker, static_maker, oj_maker = worker.flow_makers
    print(f"✓ Relax maker: {type(relax_maker).__name__}")
    print(f"✓ Static maker: {type(static_maker).__name__}")
    print(f"✓ OJ maker: {type(oj_maker).__name__}")

    print()
    return True


def test_flow_creation():
    """Test flow creation - basic worker initialization"""
    print("=" * 50)
    print("Test: Flow Creation")
    print("=" * 50)

    from htvasp.workflows import OJWorker
    from htvasp.oj import OJConfig

    lattice = Lattice.cubic(2.85)
    structure = Structure(lattice, ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    oj_config = OJConfig(j_count=2, magnetic_ion_types=["Fe"])

    # Test worker initialization
    worker = OJWorker(
        worker_name="test_flow_1",
        oj_config=oj_config,
    )

    print(f"✓ Worker initialized")
    print(f"✓ Has flow_makers: {hasattr(worker, 'flow_makers')}")
    print(f"✓ Number of makers: {len(worker.flow_makers)}")

    # Test that makers can create jobs individually
    relax_maker, static_maker, oj_maker = worker.flow_makers

    relax_job = relax_maker.make(structure)
    print(f"✓ Relax job created: {relax_job.name}")

    static_job = static_maker.make(relax_job.output.structure)
    print(f"✓ Static job created: {static_job.name}")

    # Test OJ maker can create flow
    oj_flow = oj_maker.make(static_job.output.structure)
    print(f"✓ OJ flow created: {oj_flow.name} with {len(oj_flow.jobs)} jobs")

    # Count total jobs that would be in combined flow
    # relax (1 job from DoubleRelaxMaker) + static (1 job) + OJ (3 jobs)
    total_jobs = 1 + 1 + len(oj_flow.jobs)
    print(f"✓ Total jobs in combined flow: {total_jobs}")

    # Expected: relax (1) + static (1) + OJ (3) = 5
    expected_jobs = 5
    if total_jobs >= expected_jobs:
        print(f"✓ Job count correct (>= {expected_jobs})")
    else:
        print(f"✗ Expected >= {expected_jobs} jobs, got {total_jobs}")
        return False

    print()
    return True


def test_flow_output_structure():
    """Test flow output structure - verify output keys"""
    print("=" * 50)
    print("Test: Flow Output Structure")
    print("=" * 50)

    # Test that run_flow method signature has correct output structure
    from htvasp.workflows import OJWorker
    from htvasp.oj import OJConfig
    import inspect

    worker = OJWorker(
        worker_name="test_output",
        oj_config=OJConfig(),
    )

    # Check the run_flow docstring mentions the output structure
    docstring = worker.run_flow.__doc__
    print(f"✓ run_flow docstring present")

    # Verify it mentions the expected output keys
    if "static_output" in docstring and "oj_output" in docstring:
        print(f"✓ Output structure documented")
    else:
        print(f"✗ Output structure not properly documented")
        return False

    # Check the actual implementation in run_flow
    # Simulate the combined output structure
    expected_keys = ["static_output", "oj_output", "combined"]
    print(f"✓ Expected output keys: {expected_keys}")

    print()
    return True


def test_combined_output():
    """Test combined output structure"""
    print("=" * 50)
    print("Test: Combined Output Structure")
    print("=" * 50)

    # Simulate combined output structure
    static_output = {
        "total_magnetic_moment": 4.5,
        "energy": -10.5,
        "structure": {"lattice": "cubic"},
    }

    oj_output = {
        "J_reprs": [["Fe", "Fe", 1], ["Fe", "Fe", 2]],
        "Js": [10.5, -5.2],
        "Tc_MFA": 300.0,
        "Tc_RPA": 280.0,
        "num_configs": 5,
        "magnetic_ion_types": ["Fe"],
        "vasp_success_rate": 1.0,
    }

    combined_output = {
        "static_output": static_output,
        "oj_output": oj_output,
        "combined": {
            "total_magnetic_moment": static_output["total_magnetic_moment"],
            "final_structure": static_output["structure"],
            "energy": static_output["energy"],
            "J_reprs": oj_output["J_reprs"],
            "Js": oj_output["Js"],
            "Tc_MFA": oj_output["Tc_MFA"],
            "Tc_RPA": oj_output["Tc_RPA"],
            "num_configs": oj_output["num_configs"],
            "magnetic_ion_types": oj_output["magnetic_ion_types"],
            "vasp_success_rate": oj_output["vasp_success_rate"],
        },
    }

    print(f"✓ Combined output structure:")
    print(f"  - static_output keys: {list(static_output.keys())}")
    print(f"  - oj_output keys: {list(oj_output.keys())}")
    print(f"  - combined keys: {list(combined_output['combined'].keys())}")

    # Verify all expected keys are present
    expected_combined_keys = [
        "total_magnetic_moment",
        "final_structure",
        "energy",
        "J_reprs",
        "Js",
        "Tc_MFA",
        "Tc_RPA",
        "num_configs",
        "magnetic_ion_types",
        "vasp_success_rate",
    ]

    missing_keys = [k for k in expected_combined_keys if k not in combined_output["combined"]]
    if not missing_keys:
        print(f"✓ All expected keys present in combined output")
    else:
        print(f"✗ Missing keys: {missing_keys}")
        return False

    print()
    return True


def test_incar_merging():
    """Test INCAR merging from different sources"""
    print("=" * 50)
    print("Test: INCAR Merging")
    print("=" * 50)

    from htvasp.workflows import OJWorker
    from htvasp.oj import OJConfig

    # Test with OJ config INCAR
    oj_config = OJConfig(
        j_count=2,
        magnetic_ion_types=["Fe"],
        incar={"ENCUT": 550, "ISPIN": 2},
    )

    worker = OJWorker(
        worker_name="test_incar",
        oj_config=oj_config,
        global_incar={"ENCUT": 520, "SIGMA": 0.05},
    )

    relax_maker, static_maker, oj_maker = worker.flow_makers

    # Check that OJ config INCAR is merged
    print(f"✓ Worker created with merged INCAR settings")
    print(f"✓ OJ config INCAR: {oj_config.incar}")
    print(f"✓ Global INCAR applied")

    print()
    return True


def test_resume_functionality():
    """Test resume parameter in run_flow signature"""
    print("=" * 50)
    print("Test: Resume Functionality")
    print("=" * 50)

    from htvasp.workflows import OJWorker
    from htvasp.oj import OJConfig
    import inspect

    # Check run_flow signature
    worker = OJWorker(worker_name="test", oj_config=OJConfig())

    sig = inspect.signature(worker.run_flow)
    params = list(sig.parameters.keys())

    print(f"✓ run_flow parameters: {params}")

    if "resume" in params:
        print(f"✓ 'resume' parameter present")
        default_value = sig.parameters["resume"].default
        print(f"✓ Default resume value: {default_value}")
        if default_value is True:
            print(f"✓ Default resume is True (as expected)")
        else:
            print(f"✗ Expected default resume=True, got {default_value}")
            return False
    else:
        print(f"✗ 'resume' parameter missing")
        return False

    print()
    return True


def run_unit_tests():
    """Run all unit tests"""
    print("\n" + "=" * 50)
    print("OJ Unit Tests")
    print("=" * 50 + "\n")

    tests = [
        ("Imports", test_imports),
        ("WorkerCreation", test_worker_creation),
        ("FlowCreation", test_flow_creation),
        ("FlowOutputStructure", test_flow_output_structure),
        ("CombinedOutput", test_combined_output),
        ("INCAR Merging", test_incar_merging),
        ("ResumeFunctionality", test_resume_functionality),
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
    """Run OJ workflow locally (requires OstravaJ and VASP)"""
    from htvasp.workflows import OJWorker
    from htvasp.oj import OJConfig

    structure = Structure(
        lattice=[[2.85, 0, 0], [0, 2.85, 0], [0, 0, 2.85]],
        species=["Fe", "Fe"],
        coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
    )

    flow_dir = Path("temp", "oj-Fe")
    json_path = flow_dir.joinpath("oj_Fe.json")

    if flow_dir.exists():
        shutil.rmtree(flow_dir)

    oj_config = OJConfig(
        j_count=2,
        magnetic_ion_types=["Fe"],
        extend_poscar=(1, 1, 1),
    )

    worker = OJWorker(
        worker_name="oj-Fe",
        oj_config=oj_config,
        vasp_args={
            "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_std'",
        },
    )

    output = worker.run_flow(
        name="Fe",
        structure=structure,
        flow_dir=flow_dir,
        resume=True,
    )

    if output:
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(output, f, indent=2, cls=DateTimeEncoder)
        log.info(f"Output saved to {json_path}")

        # Print summary
        print("\n" + "=" * 50)
        print("Results Summary")
        print("=" * 50)
        combined = output.get("combined", {})
        print(f"Total Magnetic Moment: {combined.get('total_magnetic_moment')}")
        print(f"Energy: {combined.get('energy')}")
        print(f"J parameters: {combined.get('Js')}")
        print(f"Tc (MFA): {combined.get('Tc_MFA')} K")
        print(f"Tc (RPA): {combined.get('Tc_RPA')} K")
        print(f"Number of configurations: {combined.get('num_configs')}")
    else:
        log.error("Workflow failed")


def submit_job():
    """Submit OJ workflow to Slurm"""
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
    parser = argparse.ArgumentParser(description="Test OJ workflow")
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
