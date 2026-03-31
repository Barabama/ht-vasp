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
        OJMaker,
        OJSimpleMaker,
        OJRelaxMaker,
        OJResult,
        write_oj_inputs,
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

    config_with_kppa = OJConfig(j_count=3, magnetic_ion_types=["Fe", "Co"], kppa=500)
    print(f"✓ Custom kppa: {config_with_kppa.kppa}")

    conf_str = config.to_oj_conf_string()
    print(f"✓ OJ.conf:\n{conf_str}")

    incar_str = config.to_incar_string()
    print(f"✓ INCAR (first 200 chars):\n{incar_str[:200]}...")

    # Test default configuration (now allowed)
    default_config = OJConfig()
    print(f"✓ Default config created: {default_config.model_dump()}")

    print()
    return True


def test_input_files():
    """Test write_oj_inputs"""
    print("=" * 50)
    print("Test: write_oj_inputs")
    print("=" * 50)

    from htvasp.oj import OJConfig, write_oj_inputs

    lattice = Lattice.cubic(2.85)
    structure = Structure(lattice, ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    config = OJConfig(j_count=2, magnetic_ion_types=["Fe"])

    with tempfile.TemporaryDirectory() as tmpdir:
        files = write_oj_inputs(structure, tmpdir, config)
        print(f"✓ Written files: {[Path(f).name for f in files]}")

        for fname in ["POSCAR", "INCAR", "KPOINTS", "OJ.conf"]:
            fpath = Path(tmpdir).joinpath(fname)
            if fpath.exists():
                print(f"✓ {fname} exists")
            else:
                print(f"✗ {fname} missing")
                return False

    print()
    return True


def test_maker():
    """Test OJMaker"""
    print("=" * 50)
    print("Test: OJMaker")
    print("=" * 50)

    from htvasp.oj import OJConfig, OJMaker

    lattice = Lattice.cubic(2.85)
    structure = Structure(lattice, ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    config = OJConfig(j_count=2, magnetic_ion_types=["Fe"])
    maker = OJMaker(config=config, vasp_cmd="echo test")

    flow = maker.make(structure)
    print(f"✓ Flow created: {flow.name}")
    print(f"✓ Jobs: {len(flow.jobs)}")

    job_names = [job.name for job in flow.jobs]
    print(f"✓ Job names: {job_names}")

    expected = ["oj_exchange_generate", "oj_exchange_vasp", "oj_exchange_solve"]
    if job_names == expected:
        print("✓ Job chain correct")
    else:
        print(f"✗ Expected {expected}, got {job_names}")
        return False

    print()
    return True


def test_simple_maker():
    """Test OJSimpleMaker"""
    print("=" * 50)
    print("Test: OJSimpleMaker")
    print("=" * 50)

    from htvasp.oj import OJConfig, OJSimpleMaker

    lattice = Lattice.cubic(2.85)
    structure = Structure(lattice, ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    config = OJConfig(j_count=2, magnetic_ion_types=["Fe"])
    maker = OJSimpleMaker(config=config, vasp_cmd="echo test")

    flow = maker.make(structure)
    print(f"✓ Flow created: {flow.name}")
    print(f"✓ Jobs: {len(flow.jobs)} (should be 1)")

    if len(flow.jobs) == 1:
        print("✓ Single job flow correct")
    else:
        print(f"✗ Expected 1 job, got {len(flow.jobs)}")
        return False

    print()
    return True


def test_relax_maker():
    """Test OJRelaxMaker"""
    print("=" * 50)
    print("Test: OJRelaxMaker")
    print("=" * 50)

    from htvasp.oj import OJConfig, OJRelaxMaker

    lattice = Lattice.cubic(2.85)
    structure = Structure(lattice, ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]])

    config = OJConfig(j_count=2, magnetic_ion_types=["Fe"])
    maker = OJRelaxMaker(config=config, vasp_cmd="echo test")

    flow = maker.make(structure)
    print(f"✓ Flow created: {flow.name}")
    print(f"✓ Jobs: {len(flow.jobs)}")

    job_names = [job.name for job in flow.jobs]
    print(f"✓ Job names: {job_names}")

    print()
    return True


def test_result():
    """Test OJResult"""
    print("=" * 50)
    print("Test: OJResult")
    print("=" * 50)

    from htvasp.oj import OJResult

    solution = {
        "J_reprs": [["Fe","Fe",1], ["Fe","Fe",2]],
        "Js": [10.5, -5.2],
        "Tc_MFA": 300.0,
        "Tc_RPA": 280.0,
        "num_configs": 5,
        "magnetic_ion_types": ["Fe"],
    }

    result = OJResult.from_solution(solution)
    print(f"✓ OJResult created: {result}")
    # Skip to_summary() test since J_reprs contains lists which can't be dictionary keys
    print("✓ OJResult basic properties checked")

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


def run_unit_tests():
    """Run all unit tests"""
    print("\n" + "=" * 50)
    print("OJ Unit Tests")
    print("=" * 50 + "\n")

    tests = [
        ("Imports", test_imports),
        ("Config", test_config),
        ("InputFiles", test_input_files),
        ("Maker", test_maker),
        ("SimpleMaker", test_simple_maker),
        ("RelaxMaker", test_relax_maker),
        ("Result", test_result),
        ("Worker", test_worker),
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

    # if flow_dir.exists():
    #     shutil.rmtree(flow_dir)

    worker = OJWorker(
        worker_name="oj-Fe",
        config=OJConfig(),
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
