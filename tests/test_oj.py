"""
Test script for OJ workflow.

Usage:
    python test_oj.py --local    # Run locally
    python test_oj.py --slurm    # Submit to Slurm
    python test_oj.py --unit     # Run unit tests only
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
    from htvasp.oj import (
        OJConfig,
        OJInputSetGenerator,
        write_oj_input_set,
        OJMaker,
        OJSimpleMaker,
        OJRelaxMaker,
        OJResult,
        oj_generate,
        oj_vasp,
        oj_solve,
        oj_workflow,
    )

    print("✓ All imports successful\n")
    return True


def test_config():
    """Test OJConfig"""
    print("=" * 50)
    print("Test: OJConfig")
    print("=" * 50)

    from htvasp.oj import OJConfig

    config = OJConfig(j_count=3, magnetic_ion_types=["Fe", "Co"])
    print(f"✓ Created: {config.model_dump()}")

    conf_str = config.to_oj_conf_string()
    print(f"✓ OJ.conf:\n{conf_str}")

    incar_str = config.to_incar_string()
    print(f"✓ INCAR (first 200 chars):\n{incar_str[:200]}...")

    default_config = OJConfig()
    print(f"✓ Default config created: {default_config.model_dump()}")

    config_updated = config.update_incar(ENCUT=520, NEW_PARAM=123)
    print(f"✓ Updated INCAR ENCUT: {config_updated.incar['ENCUT']}")
    print(f"✓ Added NEW_PARAM: {config_updated.incar.get('NEW_PARAM')}")

    config_removed = config.update_incar(ISPIN=None)
    print(f"✓ Removed ISPIN: {'ISPIN' not in config_removed.incar}")

    config_merged = config.merge_incar({"ENCUT": 600, "KPAR": 4})
    print(f"✓ Merged INCAR ENCUT: {config_merged.incar['ENCUT']}")
    print(f"✓ Merged INCAR KPAR: {config_merged.incar['KPAR']}")

    print()
    return True


def test_input_set_generator():
    """Test OJInputSetGenerator"""
    print("=" * 50)
    print("Test: OJInputSetGenerator")
    print("=" * 50)

    from htvasp.oj import OJInputSetGenerator

    lattice = Lattice.cubic(2.85)
    structure = Structure(lattice, ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    generator = OJInputSetGenerator(
        j_count=2,
        magnetic_ion_types=["Fe"],
        reciprocal_density=100,
        reciprocal_density_metal=400,
    )

    print(f"✓ Generator created")
    print(f"  j_count: {generator.j_count}")
    print(f"  reciprocal_density: {generator.reciprocal_density}")
    print(f"  auto_metal_kpoints: {generator.auto_metal_kpoints}")

    input_set = generator.get_input_set(structure)
    print(f"✓ Input set generated: {type(input_set).__name__}")
    print(f"  INCAR keys: {list(input_set.incar.keys())[:5]}...")
    print(f"  KPOINTS: {input_set.kpoints}")

    oj_conf = generator.get_oj_conf_string()
    print(f"✓ OJ.conf generated:\n{oj_conf}")

    print()
    return True


def test_write_input_set():
    """Test write_oj_input_set"""
    print("=" * 50)
    print("Test: write_oj_input_set")
    print("=" * 50)

    from htvasp.oj import OJInputSetGenerator, write_oj_input_set

    lattice = Lattice.cubic(2.85)
    structure = Structure(lattice, ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    generator = OJInputSetGenerator(
        j_count=2,
        magnetic_ion_types=["Fe"],
        reciprocal_density=100,
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        write_oj_input_set(structure, tmpdir, generator)

        for fname in ["POSCAR", "INCAR", "KPOINTS", "POTCAR", "OJ.conf"]:
            fpath = Path(tmpdir).joinpath(fname)
            if fpath.exists():
                print(f"✓ {fname} exists")
            else:
                print(f"✗ {fname} missing")
                return False

    print()
    return True


def test_jobs():
    """Test job functions"""
    print("=" * 50)
    print("Test: Job Functions")
    print("=" * 50)

    from htvasp.oj import oj_generate, oj_vasp, oj_solve

    lattice = Lattice.cubic(2.85)
    structure = Structure(lattice, ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    gen_job = oj_generate(
        structure,
        j_count=2,
        magnetic_ion_types=["Fe"],
        reciprocal_density=100,
    )
    print(f"✓ oj_generate job created: {gen_job.name}")

    vasp_job = oj_vasp(["/tmp/flip1", "/tmp/flip2"], "vasp_std")
    print(f"✓ oj_vasp job created: {vasp_job.name}")

    solve_job = oj_solve("/tmp/run_dir")
    print(f"✓ oj_solve job created: {solve_job.name}")

    print()
    return True


def test_maker():
    """Test OJMaker"""
    print("=" * 50)
    print("Test: OJMaker")
    print("=" * 50)

    from htvasp.oj import OJMaker

    lattice = Lattice.cubic(2.85)
    structure = Structure(lattice, ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    maker = OJMaker(
        name="test_exchange",
        j_count=2,
        magnetic_ion_types=["Fe"],
        vasp_cmd="echo test",
        reciprocal_density=100,
    )

    flow = maker.make(structure)
    print(f"✓ Flow created: {flow.name}")
    print(f"✓ Flow jobs: {len(flow.jobs)}")

    job_names = [job.name for job in flow.jobs]
    print(f"✓ Job names: {job_names}")

    expected = ["test_exchange_generate", "test_exchange_vasp", "test_exchange_solve"]
    if job_names == expected:
        print("✓ Job chain correct")
    else:
        print(f"✗ Expected {expected}, got {job_names}")
        return False

    print()
    return True


def test_maker_with_config():
    """Test OJMaker with OJConfig (backward compatibility)"""
    print("=" * 50)
    print("Test: OJMaker with OJConfig")
    print("=" * 50)

    from htvasp.oj import OJConfig, OJMaker

    lattice = Lattice.cubic(2.85)
    structure = Structure(lattice, ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    config = OJConfig(j_count=2, magnetic_ion_types=["Fe"])
    maker = OJMaker(config=config, vasp_cmd="echo test")

    flow = maker.make(structure)
    print(f"✓ Flow created with config: {flow.name}")
    print(f"✓ j_count from config: {maker.j_count}")
    print(f"✓ magnetic_ion_types from config: {maker.magnetic_ion_types}")

    print()
    return True


def test_maker_with_user_incar():
    """Test OJMaker with user_incar_settings (atomate2 style)"""
    print("=" * 50)
    print("Test: OJMaker with user_incar_settings")
    print("=" * 50)

    from htvasp.oj import OJMaker

    lattice = Lattice.cubic(2.85)
    structure = Structure(lattice, ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    maker = OJMaker(
        j_count=2,
        magnetic_ion_types=["Fe"],
        user_incar_settings={"ENCUT": 520, "ISPIN": 2},
        reciprocal_density=100,
    )

    print(f"✓ Maker created")
    print(f"✓ user_incar_settings: {maker.user_incar_settings}")
    print(f"✓ ENCUT in settings: {maker.user_incar_settings.get('ENCUT')}")

    flow = maker.make(structure)
    print(f"✓ Flow created: {flow.name}")

    print()
    return True


def test_result():
    """Test OJResult"""
    print("=" * 50)
    print("Test: OJResult")
    print("=" * 50)

    from htvasp.oj import OJResult
    from htvasp.oj.task_doc import format_j_repr

    solution = {
        "J_reprs": [["Fe", "Fe", 1], ["Fe", "Fe", 2]],
        "Js": [10.5, -5.2],
        "Tc_MFA": 300.0,
        "Tc_RPA": 280.0,
        "num_configs": 5,
        "magnetic_ion_types": ["Fe"],
    }

    result = OJResult.from_solution(solution)
    print(f"✓ OJResult created: {result}")

    formatted = format_j_repr(["Fe", "Fe", 1])
    print(f"✓ format_j_repr: {formatted}")
    assert formatted == "Fe-Fe-1", f"Expected 'Fe-Fe-1', got {formatted}"

    j_pairs = result.get_j_pairs_dict()
    print(f"✓ get_j_pairs_dict: {j_pairs}")
    assert "Fe-Fe-1" in j_pairs, "Expected 'Fe-Fe-1' in J_pairs"
    assert j_pairs["Fe-Fe-1"] == 10.5, f"Expected 10.5, got {j_pairs['Fe-Fe-1']}"

    summary = result.to_summary()
    print(f"✓ to_summary: {summary}")
    assert "J_pairs" in summary, "Expected 'J_pairs' in summary"
    assert summary["Tc_MFA"] == 300.0, f"Expected 300.0, got {summary['Tc_MFA']}"

    print("✓ OJResult all methods tested successfully")

    print()
    return True


def test_worker():
    """Test OJWorker"""
    print("=" * 50)
    print("Test: OJWorker")
    print("=" * 50)

    from htvasp.workflows import OJWorker
    from htvasp.oj import OJConfig

    worker = OJWorker(worker_name="test", config=OJConfig())

    print(f"✓ Worker created")
    print(f"✓ Has _maker: {hasattr(worker, '_maker')}")

    lattice = Lattice.cubic(2.85)
    structure = Structure(lattice, ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    flow, maker = worker.create_flow(structure, "Fe")
    print(f"✓ Flow created: {flow.name}, jobs: {len(flow.jobs)}")

    print()
    return True


def test_atomate2_integration():
    """Test integration with atomate2 StaticMaker"""
    print("=" * 50)
    print("Test: atomate2 Integration")
    print("=" * 50)

    from jobflow import Flow
    from atomate2.vasp.jobs.core import StaticMaker
    from htvasp.oj import OJMaker

    lattice = Lattice.cubic(2.85)
    structure = Structure(lattice, ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    static_maker = StaticMaker()
    static_job = static_maker.make(structure)
    print(f"✓ StaticMaker job created: {static_job.name}")

    oj_maker = OJMaker(
        j_count=2,
        magnetic_ion_types=["Fe"],
        reciprocal_density=100,
        user_incar_settings={"ENCUT": 520},
    )
    oj_flow = oj_maker.make(structure)
    print(f"✓ OJMaker flow created: {oj_flow.name}")
    
    # Create combined flow using job references (not the flow directly)
    # This demonstrates that OJ Maker jobs can be combined with atomate2 jobs
    combined_jobs = [static_job] + list(oj_flow.jobs)
    print(f"✓ Combined jobs count: {len(combined_jobs)}")
    print(f"✓ Job types: {[type(j).__name__ for j in combined_jobs]}")

    print()
    return True


def run_unit_tests():
    """Run all unit tests"""
    print("\n" + "=" * 50)
    print("OJ Unit Tests")
    print("=" * 50 + "\n")

    tests = [
        ("Imports", test_imports),
        ("Config", test_config),
        ("InputSetGenerator", test_input_set_generator),
        ("WriteInputSet", test_write_input_set),
        ("Jobs", test_jobs),
        ("Maker", test_maker),
        ("MakerWithConfig", test_maker_with_config),
        ("MakerWithUserIncar", test_maker_with_user_incar),
        ("Result", test_result),
        ("Worker", test_worker),
        ("Atomate2Integration", test_atomate2_integration),
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

    worker = OJWorker(
        worker_name="oj-Fe",
        config=OJConfig(extend_poscar=[1, 1, 1]),
        vasp_args={
            "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_std'",
        },
    )

    output = worker.run_flow(
        name="Fe",
        structure=structure,
        flow_dir=flow_dir,
    )

    if output:
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(output, f, indent=2, cls=DateTimeEncoder)
        log.info(f"Output saved to {json_path}")
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
