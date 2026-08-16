"""
Test script for OJ workflow.

Usage:
python test_oj.py --unit     # Run unit tests only
python test_oj.py --local    # Run locally (fresh run)
python test_oj.py --rerun    # Run locally with resume
python test_oj.py --slurm    # Submit to Slurm
"""

import json
import shutil
import logging
import argparse
import inspect
import traceback
from pathlib import Path

from pymatgen.core import Structure, Lattice

from htvasp.workflows import OJWorker
from htvasp.oj import OJInputSetGenerator
from htvasp.oj.task_doc import OJResult

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


def get_fe_bcc_structure(a: float = 2.85) -> Structure:
    """Return BCC Fe structure for testing."""
    return Structure(Lattice.cubic(a), ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])


def make_test_worker(name: str = "test_oj", **kwargs) -> OJWorker:
    """Create OJWorker with standard test defaults."""
    defaults = {
        "worker_name": name,
        "vasp_args": {"vasp_cmd": "vasp_std"},
        "j_count": 2,
        "magnetic_ion_types": ["Fe"],
    }
    defaults.update(kwargs)
    return OJWorker(**defaults)


# =============================================================================
# Unit Tests
# =============================================================================


def test_imports():
    """Test all imports"""
    print_header("Imports")
    print("✓ All imports successful\n")
    return True


def test_worker_creation():
    """Test OJWorker creation"""
    print_header("Worker Creation")

    worker = make_test_worker(user_incar_settings={"ENCUT": 520})
    print(f"✓ Worker created: {worker.worker_name}")
    print(f"✓ Has flow_maker: {hasattr(worker, 'flow_maker')}")
    print(f"✓ OJ maker: {type(worker.flow_maker).__name__}")
    print()
    return True


def test_flow_creation():
    """Test flow creation - OJ workflow"""
    print_header("Flow Creation")

    worker = make_test_worker(name="test_flow_1")
    print("✓ Worker initialized")
    print(f"✓ Has flow_maker: {hasattr(worker, 'flow_maker')}")

    structure = get_fe_bcc_structure()
    oj_flow = worker.flow_maker.make(structure)
    print(f"✓ OJ flow created: {oj_flow.name} with {len(oj_flow.jobs)} jobs")

    expected_jobs = 3
    if len(oj_flow.jobs) >= expected_jobs:
        print(f"✓ Job count correct (>= {expected_jobs})")
    else:
        print(f"✗ Expected >= {expected_jobs} jobs, got {len(oj_flow.jobs)}")
        return False

    print()
    return True


def test_oj_result_model():
    """Test OJResult pydantic model with new fields"""
    print_header("OJResult Model (new fields)")

    # Minimal test data - only fields that are actually verified
    solution = {
        "J_reprs": [["Fe", "Fe", 1], ["Fe", "Fe", 2]],
        "Js": [10.5, -5.2],
        "Tc_MFA": 300.0,
        "Tc_RPA": -1,
        "num_configs": 5,
        "magnetic_ion_types": ["Fe"],
        "vasp_success_rate": 1.0,
        "E_DLM": -25.6,
        "D_stiff": 0.00012,
        "is_complete": True,
        "rank": 4,
        "resid": 0.0,
        "condition_number": 12.0,
        "singular_values": [192.0, 96.0, 48.0, 16.0],
        "avg_magnetic_moment": 1.75,
        "total_magnetic_moment": 3.50,
        "ground_state_moments": [1.765, -1.740],
        "atom_types": ["Fe", "Fe"],
        "configs_data": [
            {
                "flip_name": "flip00000",
                "total_energy": -409.15,
                "free_energy": -409.16,
                "energy_wo_entrp": -409.14,
                "fermi_energy": 5.23,
                "initial_moments": [1.0, -1.0],
                "final_moments": [1.765, -1.740],
                "positions": [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
                "basis": [[2.87, 0.0, 0.0], [0.0, 2.87, 0.0], [0.0, 0.0, 2.87]],
            }
        ],
        "energies": [-409.15],
        "J_mat": [[64, -48]],
    }

    result = OJResult.from_solution(solution)

    print(f"✓ J_reprs: {result.J_reprs}")
    print(f"✓ Js: {result.Js}")
    print(f"✓ Tc_MFA: {result.Tc_MFA}")
    print(f"✓ Tc_RPA: {result.Tc_RPA}")
    print(f"✓ E_DLM: {result.E_DLM}")
    print(f"✓ D_stiff: {result.D_stiff}")
    print(f"✓ is_complete: {result.is_complete}")
    print(f"✓ rank: {result.rank}")
    print(f"✓ condition_number: {result.condition_number}")
    print(f"✓ singular_values: {result.singular_values}")
    print(f"✓ avg_magnetic_moment: {result.avg_magnetic_moment}")
    print(f"✓ total_magnetic_moment: {result.total_magnetic_moment}")
    print(f"✓ ground_state_moments: {result.ground_state_moments}")
    print(f"✓ configs_data count: {len(result.configs_data)}")
    print(f"✓ energies: {result.energies}")

    # Verify configs_data structure
    if result.configs_data:
        cfg = result.configs_data[0]
        expected_keys = [
            "flip_name",
            "total_energy",
            "free_energy",
            "energy_wo_entrp",
            "fermi_energy",
            "initial_moments",
            "final_moments",
            "positions",
            "basis",
        ]
        missing = [k for k in expected_keys if k not in cfg]
        if missing:
            print(f"✗ Missing keys in configs_data[0]: {missing}")
            return False
        print("✓ configs_data[0] has all expected keys")

    print()
    return True


def test_summary_with_warnings():
    """Test to_summary() with warnings for new fields"""
    print_header("to_summary() with warnings")

    # Minimal OJResult - only fields needed to trigger warnings
    result = OJResult(
        J_reprs=[["Fe", "Fe", 1]],
        Js=[10.5],
        Tc_MFA=300.0,
        Tc_RPA=-1.0,
        is_complete=True,
        condition_number=1e12,
        E_DLM=-25.6,
        avg_magnetic_moment=1.75,
        total_magnetic_moment=3.5,
        num_configs=1,
        vasp_success_rate=1.0,
        D_stiff=0.00012,
        rank=1,
        resid=0.0,
    )

    summary = result.to_summary()
    print(f"✓ Summary keys: {list(summary.keys())}")
    print(f"✓ E_DLM in summary: {summary.get('E_DLM')}")
    print(f"✓ condition_number in summary: {summary.get('condition_number')}")
    print(f"✓ avg_magnetic_moment in summary: {summary.get('avg_magnetic_moment')}")

    if "warnings" not in summary:
        print("✗ Expected warnings but none found")
        return False

    print(f"✓ Warnings: {summary['warnings']}")
    assert any("condition" in w.lower() for w in summary["warnings"])
    assert any("Tc_RPA" in w or "negative" in w.lower() for w in summary["warnings"])
    print("✓ Warnings correctly generated")
    print()
    return True


def test_incar_merging():
    """Test INCAR settings via user_incar_settings"""
    print_header("INCAR Settings")

    user_incar_settings = {"ENCUT": 520, "SIGMA": 0.05}
    worker = make_test_worker(name="test_incar", user_incar_settings=user_incar_settings)

    print(f"✓ Worker created with user_incar_settings: {user_incar_settings}")

    gen = worker.flow_maker.input_set_generator
    print(f"✓ Input set generator: {type(gen).__name__}")
    print(f"✓ user_incar_settings in generator: {gen.user_incar_settings}")
    print()
    return True


def test_input_set_generator():
    """Test OJInputSetGenerator directly"""
    print_header("OJInputSetGenerator")

    gen = OJInputSetGenerator(
        j_count=3,
        magnetic_ion_types=["Fe", "Co"],
        user_incar_settings={"ENCUT": 520},
    )

    print(f"✓ j_count: {gen.j_count}")
    print(f"✓ magnetic_ion_types: {gen.magnetic_ion_types}")
    print(f"✓ user_incar_settings: {gen.user_incar_settings}")

    oj_conf = gen.get_oj_conf_string()
    print("✓ OJ.conf string generated:")
    for line in oj_conf.strip().split("\n"):
        print(f"  {line}")
    print()
    return True


def test_resume_functionality():
    """Test resume parameter in run_flow signature"""
    print_header("Resume Functionality")

    worker = make_test_worker()
    sig = inspect.signature(worker.run_flow)
    params = list(sig.parameters.keys())

    print(f"✓ run_flow parameters: {params}")

    if "resume" not in params:
        print("✗ 'resume' parameter missing")
        return False

    print("✓ 'resume' parameter present")
    default_value = sig.parameters["resume"].default
    print(f"✓ Default resume value: {default_value}")

    if default_value is not True:
        print(f"✗ Expected default resume=True, got {default_value}")
        return False

    print("✓ Default resume is True (as expected)")
    print()
    return True


# =============================================================================
# Test Runner
# =============================================================================


def run_unit_tests():
    """Run all unit tests"""
    print("\n" + "=" * 50)
    print("OJ Unit Tests")
    print("=" * 50 + "\n")

    tests = [
        ("Imports", test_imports),
        ("WorkerCreation", test_worker_creation),
        ("FlowCreation", test_flow_creation),
        ("OJResultModel", test_oj_result_model),
        ("SummaryWithWarnings", test_summary_with_warnings),
        ("INCAR Settings", test_incar_merging),
        ("InputSetGenerator", test_input_set_generator),
        ("ResumeFunctionality", test_resume_functionality),
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
        print(f" {name}: {status}")
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
    """Run OJ workflow locally (requires OstravaJ and VASP).

    Args:
        clean: If True, remove existing flow_dir before running.
               If False, keep existing data and use resume=True.
    """
    structure = get_fe_bcc_structure()
    flow_name = "Fe-oj"
    flow_dir = Path("/nfs_ssd/tmp")
    store_dir = Path("./temp") / flow_name
    json_path = store_dir / f"{flow_name}.json"

    if clean and flow_dir.exists():
        shutil.rmtree(flow_dir)
        log.info(f"Cleaned existing flow_dir: {flow_dir}")

    worker = OJWorker(
        vasp_args=CLUSTER_VASP_ARGS,
        global_incar={"GGA": "PE"},
        j_count=2,
        extend_poscar=(1, 1, 1),
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
        oj_result = OJResult.from_solution(output)
        summary = oj_result.to_summary()

        print("\n" + "=" * 50)
        print("Results Summary")
        print("=" * 50)
        print(f"J parameters: {summary.get('J_pairs')}")
        print(f"Tc (MFA): {summary.get('Tc_MFA')} K")
        print(f"Tc (RPA): {summary.get('Tc_RPA')} K")
        print(f"E_DLM: {summary.get('E_DLM')} eV/atom")
        print(f"avg magnetic moment: {summary.get('avg_magnetic_moment')} mu_B")
        if cond := summary.get("condition_number"):
            print(f"condition number: {cond:.2e}")
        if "warnings" in summary:
            print(f"Warnings: {summary['warnings']}")
    else:
        log.error("Workflow failed")


def submit_job():
    """Submit OJ workflow to Slurm with abort/rerun simulation."""
    from htvasp.slurm import SlurmJobManager

    manager = SlurmJobManager()

    def submit(config_overrides: dict | None = None, command_suffix: str = "") -> str | None:
        """Helper to submit a job with common defaults."""
        config = manager.get_cpu_config(ntasks=16, memory="8G", **(config_overrides or {}))
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

    time.sleep(30)

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
    parser = argparse.ArgumentParser(description="Test OJ workflow")
    parser.add_argument("--local", action="store_true", help="run workflow locally (fresh run)")
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
