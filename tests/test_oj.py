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
    from htvasp.oj import OJInputSetGenerator
    from htvasp.oj.task_doc import OJResult

    print("✓ All imports successful\n")
    return True


def test_worker_creation():
    """Test OJWorker creation"""
    print("=" * 50)
    print("Test: Worker Creation")
    print("=" * 50)

    from htvasp.workflows import OJWorker

    worker = OJWorker(
        worker_name="test_oj",
        vasp_args={"vasp_cmd": "vasp_std"},
        j_count=2,
        magnetic_ion_types=["Fe"],
        user_incar_settings={"ENCUT": 520},
    )

    print(f"✓ Worker created: {worker.worker_name}")
    print(f"✓ Has oj_maker: {hasattr(worker, 'oj_maker')}")

    oj_maker = worker.oj_maker
    print(f"✓ OJ maker: {type(oj_maker).__name__}")

    print()
    return True


def test_flow_creation():
    """Test flow creation - OJ workflow"""
    print("=" * 50)
    print("Test: Flow Creation")
    print("=" * 50)

    from htvasp.workflows import OJWorker

    lattice = Lattice.cubic(2.85)
    structure = Structure(lattice, ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    worker = OJWorker(
        worker_name="test_flow_1",
        vasp_args={"vasp_cmd": "vasp_std"},
        j_count=2,
        magnetic_ion_types=["Fe"],
    )

    print(f"✓ Worker initialized")
    print(f"✓ Has oj_maker: {hasattr(worker, 'oj_maker')}")

    oj_maker = worker.oj_maker
    oj_flow = oj_maker.make(structure)
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
    print("=" * 50)
    print("Test: OJResult Model (new fields)")
    print("=" * 50)

    from htvasp.oj.task_doc import OJResult

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
        if not missing:
            print(f"✓ configs_data[0] has all expected keys")
        else:
            print(f"✗ Missing keys in configs_data[0]: {missing}")
            return False

    print()
    return True


def test_summary_with_warnings():
    """Test to_summary() with warnings for new fields"""
    print("=" * 50)
    print("Test: to_summary() with warnings")
    print("=" * 50)

    from htvasp.oj.task_doc import OJResult

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

    if "warnings" in summary:
        print(f"✓ Warnings: {summary['warnings']}")
        assert any(
            "condition" in w.lower() for w in summary["warnings"]
        ), "Expected high condition number warning"
        assert any(
            "Tc_RPA" in w or "negative" in w.lower() for w in summary["warnings"]
        ), "Expected Tc_RPA warning"
        print(f"✓ Warnings correctly generated")
    else:
        print("✗ Expected warnings but none found")
        return False

    print()
    return True


def test_incar_merging():
    """Test INCAR settings via user_incar_settings"""
    print("=" * 50)
    print("Test: INCAR Settings")
    print("=" * 50)

    from htvasp.workflows import OJWorker
    from htvasp.oj import OJInputSetGenerator

    user_incar_settings = {"ENCUT": 520, "SIGMA": 0.05}

    worker = OJWorker(
        worker_name="test_incar",
        vasp_args={"vasp_cmd": "vasp_std"},
        j_count=2,
        magnetic_ion_types=["Fe"],
        user_incar_settings=user_incar_settings,
    )

    print(f"✓ Worker created with user_incar_settings")
    print(f"✓ user_incar_settings: {user_incar_settings}")
    print(f"✓ OJ maker created")

    gen = worker.oj_maker.input_set_generator
    print(f"✓ Input set generator: {type(gen).__name__}")
    print(f"✓ user_incar_settings in generator: {gen.user_incar_settings}")

    print()
    return True


def test_input_set_generator():
    """Test OJInputSetGenerator directly"""
    print("=" * 50)
    print("Test: OJInputSetGenerator")
    print("=" * 50)

    from htvasp.oj import OJInputSetGenerator

    gen = OJInputSetGenerator(
        j_count=3,
        magnetic_ion_types=["Fe", "Co"],
        user_incar_settings={"ENCUT": 520},
    )

    print(f"✓ j_count: {gen.j_count}")
    print(f"✓ magnetic_ion_types: {gen.magnetic_ion_types}")
    print(f"✓ user_incar_settings: {gen.user_incar_settings}")

    oj_conf = gen.get_oj_conf_string()
    print(f"✓ OJ.conf string generated:")
    for line in oj_conf.strip().split("\n"):
        print(f"  {line}")

    print()
    return True


def test_resume_functionality():
    """Test resume parameter in run_flow signature"""
    print("=" * 50)
    print("Test: Resume Functionality")
    print("=" * 50)

    from htvasp.workflows import OJWorker
    import inspect

    worker = OJWorker(worker_name="test", vasp_args={"vasp_cmd": "vasp_std"})

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
        print(f" {name}: {status}")
    print(f"\nTotal: {passed}/{total} passed")

    return passed == total


def run_locally():
    """Run OJ workflow locally (requires OstravaJ and VASP)"""
    from htvasp.workflows import OJWorker

    structure = Structure(
        lattice=[[2.85, 0, 0], [0, 2.85, 0], [0, 0, 2.85]],
        species=["Fe", "Fe"],
        coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
    )

    flow_dir = Path("temp", "oj-Fe")
    json_path = flow_dir.joinpath("oj_Fe.json")

    # if flow_dir.exists():
    #     shutil.rmtree(flow_dir)

    vasp_args = {
        "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_std'",
        "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_gam'",
    }
    worker = OJWorker(
        worker_name="oj-Fe",
        vasp_args=vasp_args,
        global_incar={"GGA": "PE"},
        j_count=2,
        extend_poscar=(1, 1, 1),
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

        from htvasp.oj.task_doc import OJResult

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
        print(
            f"condition number: {summary.get('condition_number'):.2e}"
            if summary.get("condition_number")
            else "condition number: N/A"
        )
        if "warnings" in summary:
            print(f"Warnings: {summary['warnings']}")
    else:
        log.error("Workflow failed")


def submit_job():
    """Submit OJ workflow to Slurm"""
    from htvasp.slurm import SlurmJobManager

    manager = SlurmJobManager()
    config = manager.get_cpu_config(
        ntasks=16,
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
