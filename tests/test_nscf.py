"""
Test script for NscfWorker workflow.

Usage:
    python test_nscf.py --local    # Run locally
    python test_nscf.py --slurm    # Submit to Slurm
    python test_nscf.py --unit     # Run unit tests only
"""

import json
import shutil
import logging
import argparse
import traceback
from pathlib import Path

from pymatgen.core import Structure, Lattice

from htvasp.workflows import NscfWorker

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


def make_test_worker(name: str = "test-nscf", **kwargs) -> NscfWorker:
    """Create NscfWorker with standard test defaults."""
    defaults = {
        "worker_name": name,
        "vasp_args": {"vasp_cmd": "echo test"},
    }
    defaults.update(kwargs)
    return NscfWorker(**defaults)


# =============================================================================
# Unit Tests
# =============================================================================


def test_imports():
    """Test all imports"""
    print_header("Imports")
    print("✓ All imports successful\n")
    return True


def test_worker():
    """Test NscfWorker creation"""
    print_header("NscfWorker")

    worker = make_test_worker()
    print("✓ Worker created")
    print(f"✓ Has flow_makers: {hasattr(worker, 'flow_makers')}")
    print(f"✓ compute_band: {worker.compute_band}")
    print()
    return True


def test_worker_without_band():
    """Test NscfWorker without band calculation"""
    print_header("NscfWorker (no band)")

    worker = make_test_worker(compute_band=False)
    print("✓ Worker created with compute_band=False")
    print(f"✓ compute_band: {worker.compute_band}")
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


def test_dos_parameters():
    """Test DOS-specific parameters"""
    print_header("DOS Parameters")

    worker = make_test_worker(
        reciprocal_density=300,
        dedos=0.01,
    )
    print()
    return True


def test_band_parameters():
    """Test Band-specific parameters"""
    print_header("Band Parameters")

    worker = make_test_worker(line_density=30)
    print()
    return True


def test_flow_makers():
    """Test flow makers"""
    print_header("Flow Makers")

    worker = make_test_worker()
    relax_maker, static_maker, nscf_dos_maker, nscf_band_maker = worker.flow_makers
    print(f"✓ Relax maker: {relax_maker.name}")
    print(f"✓ Static maker: {static_maker.name}")
    print(f"✓ NSCF-DOS maker: {nscf_dos_maker.name}")
    print(f"✓ NSCF-Band maker: {nscf_band_maker.name}")
    print()
    return True


def test_flow_creation():
    """Test Flow creation"""
    print_header("Flow Creation")

    worker = make_test_worker()
    structure = get_al_bcc_structure()
    flow = worker._make_flow(structure)

    print(f"✓ Flow created with {len(flow.jobs)} jobs")
    job_names = [job.name for job in flow.jobs]
    print(f"✓ Job names: {job_names}")
    print()
    return True


def test_flow_creation_no_band():
    """Test Flow creation without band"""
    print_header("Flow Creation (no band)")

    worker = make_test_worker(compute_band=False)
    structure = get_al_bcc_structure()
    flow = worker._make_flow(structure)

    print(f"✓ Flow created with {len(flow.jobs)} jobs")
    job_names = [job.name for job in flow.jobs]
    print(f"✓ Job names: {job_names}")
    print()
    return True


# =============================================================================
# Test Runner
# =============================================================================


def run_unit_tests():
    """Run all unit tests"""
    print("\n" + "=" * 50)
    print("NscfWorker Unit Tests")
    print("=" * 50 + "\n")

    tests = [
        ("Imports", test_imports),
        ("Worker", test_worker),
        ("WorkerNoBand", test_worker_without_band),
        ("IncarSettings", test_incar_settings),
        ("DosParameters", test_dos_parameters),
        ("BandParameters", test_band_parameters),
        ("FlowMakers", test_flow_makers),
        ("FlowCreation", test_flow_creation),
        ("FlowCreationNoBand", test_flow_creation_no_band),
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
    """Run NSCF workflow locally (requires VASP)"""
    structure = get_al_bcc_structure()
    flow_name = "Al-nscf"
    flow_dir = Path("/nfs_ssd/tmp")
    store_dir = Path("./temp") / flow_name
    json_path = store_dir / f"{flow_name}.json"

    if flow_dir.exists():
        shutil.rmtree(flow_dir)

    worker = NscfWorker(
        vasp_args=CLUSTER_VASP_ARGS,
        global_incar={"GGA": "PE"},
    )

    worker.run_flow(
        name=flow_name,
        structure=structure,
        flow_dir=flow_dir,
        store_dir=store_dir,
    )

    dos_output = worker.get_result("nscf uniform")
    band_output = worker.get_result("nscf line") if worker.compute_band else None

    results = {"dos": dos_output}
    if band_output:
        results["band"] = band_output

    worker.write_result(data=results, json_path=json_path)


def submit_job():
    """Submit NSCF workflow to Slurm"""
    from htvasp.slurm import SlurmJobManager

    manager = SlurmJobManager()
    config = manager.get_cpu_config(ntasks=8, memory="4G")
    job_id = manager.submit_command(
        command=f"python {__file__} --local",
        config=config,
        conda_env="htvasp",
        workdir=".",
    )
    if job_id:
        log.info(f"Submitted job: {job_id}")


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test NscfWorker workflow")
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
