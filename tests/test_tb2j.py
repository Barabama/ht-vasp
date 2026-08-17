"""
Test script for the TB2J workflow (Wannier90-based magnetic exchange).

Usage:
python test_tb2j.py --unit     # Run unit tests only (no VASP)
python test_tb2j.py --local    # Run locally (fresh run; requires VASP + TB2J)
python test_tb2j.py --rerun    # Run locally with resume
python test_tb2j.py --slurm    # Submit to Slurm
"""

import os
import sys
import json
import shutil
import logging
import argparse
import tempfile
import traceback
import warnings
from pathlib import Path

from pymatgen.core import Structure, Lattice

from htvasp.tb2j import Tb2jInputSetGenerator, ORBITALS_PER_ATOM, Tb2jResult
from htvasp.tb2j.jobs import parse_exchange_out
from htvasp.workflows import Tb2jWorker

# Make the repo root importable so `import main_w90` works regardless of cwd.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

# Root of the verified TB2J test runs. exchange.out files here are read-only
# reference data used by the unit tests that validate the real parser output.
TB2J_TEST_ROOT = Path("/nfs_hdd/2025/gaominliang/tb2j_test")

# Cosmetic: silence POTCAR-functional / low-spin warnings emitted by atomate2
# and pymatgen while building input sets with simple test structures.
warnings.filterwarnings(
    "ignore",
    message=".*Overriding the POTCAR functional.*",
    category=UserWarning,
)
warnings.filterwarnings(
    "ignore",
    message=".*oxidation state.*low spin.*",
    category=UserWarning,
)


# =============================================================================
# Test Helpers
# =============================================================================


def print_header(name: str):
    """Print standardized test header."""
    print("=" * 50)
    print(f"Test: {name}")
    print("=" * 50)


def get_test_structure(elements=("Fe", "Co"), a: float = 2.85) -> Structure:
    """Return a simple 2-atom BCC test cell with the given species."""
    return Structure(Lattice.cubic(a), list(elements), [[0, 0, 0], [0.5, 0.5, 0.5]])


def make_test_worker(name: str = "test_tb2j", **kwargs) -> Tb2jWorker:
    """Create a Tb2jWorker with standard test defaults."""
    defaults = {
        "worker_name": name,
        "vasp_args": {"vasp_cmd": "vasp_std"},
        "global_incar": {
            "ENCUT": 520,
            "ISPIN": 2,
            "MAGMOM": {"Fe": 5.0, "Co": 3.0},
        },
    }
    defaults.update(kwargs)
    return Tb2jWorker(**defaults)


def get_exchange_out_path(name: str) -> Path:
    """Return the path to a verified TB2J exchange.out file (read-only)."""
    path = TB2J_TEST_ROOT / name / "TB2J" / "exchange.out"
    if not path.exists():
        raise FileNotFoundError(f"exchange.out not found (read-only reference): {path}")
    return path


def parse_reference(name: str) -> dict:
    """Parse a reference exchange.out once and cache it for reuse."""
    return parse_exchange_out(get_exchange_out_path(name))


# =============================================================================
# Unit Tests
# =============================================================================


def test_imports():
    """Test all imports."""
    print_header("Imports")
    from htvasp.tb2j.input_set import W90_INCAR_DEFAULTS  # noqa: F401
    from htvasp.tb2j.jobs import tb2j_solve  # noqa: F401
    from htvasp.tb2j.maker import Tb2jMaker  # noqa: F401
    from htvasp.workflows import Tb2jWorker  # noqa: F401
    print("✓ All imports successful\n")
    return True


def test_worker_creation():
    """Test Tb2jWorker creation and maker wiring."""
    print_header("Worker Creation")

    worker = make_test_worker(elements=["Fe", "Co"])
    print(f"✓ Worker created: {worker.worker_name}")
    print(f"✓ Relax maker: {type(worker.relax_maker).__name__}")
    print(f"✓ TB2J maker: {type(worker.tb2j_maker).__name__}")

    if not hasattr(worker, "relax_maker") or not hasattr(worker, "tb2j_maker"):
        print("✗ Missing relax_maker/tb2j_maker attributes")
        return False
    if type(worker.relax_maker).__name__ != "RelaxMaker":
        print(f"✗ Expected RelaxMaker, got {type(worker.relax_maker).__name__}")
        return False
    if type(worker.tb2j_maker).__name__ != "Tb2jMaker":
        print(f"✗ Expected Tb2jMaker, got {type(worker.tb2j_maker).__name__}")
        return False
    if worker.kmesh != (9, 9, 9):
        print(f"✗ Default kmesh should be (9,9,9), got {worker.kmesh}")
        return False
    print("✓ Relax = RelaxMaker, exchange = Tb2jMaker, kmesh=(9,9,9)")
    print()
    return True


def test_input_set_generator():
    """Test Tb2jInputSetGenerator attributes and num_wann derivation."""
    print_header("Input Set Generator")

    gen = Tb2jInputSetGenerator(
        elements=["Co", "Ni"],
        num_wann=18,
        user_incar_settings={"ENCUT": 520},
    )
    print(f"✓ elements: {gen.elements}")
    print(f"✓ num_wann explicit: {gen.get_num_wann()}")

    if gen.elements != ["Co", "Ni"]:
        print(f"✗ elements not stored, got {gen.elements}")
        return False
    if gen.get_num_wann() != 18:
        print(f"✗ Explicit num_wann should win, got {gen.get_num_wann()}")
        return False

    # Auto-derivation from a structure: n_atoms * 9 orbitals.
    struct = get_test_structure(("Co", "Ni"))
    gen2 = Tb2jInputSetGenerator(elements=["Co", "Ni"])
    if gen2.get_num_wann(struct) != 18:
        print(f"✗ Auto num_wann = n_atoms*9 failed, got {gen2.get_num_wann(struct)}")
        return False
    print(f"✓ Auto num_wann from 2-atom cell: {gen2.get_num_wann(struct)} (n_atoms*{ORBITALS_PER_ATOM})")

    # No structure and no explicit num_wann must raise.
    gen3 = Tb2jInputSetGenerator(elements=["Co"])
    try:
        gen3.get_num_wann()
        print("✗ Expected ValueError for get_num_wann() without structure")
        return False
    except ValueError:
        print("✓ get_num_wann() raises ValueError without structure/num_wann")

    print()
    return True


def test_win_string_content():
    """Test wannier90.win string content: num_wann, projections, write_hr."""
    print_header("wannier90.win String")

    struct = get_test_structure(("Co", "Ni"))
    gen = Tb2jInputSetGenerator(elements=["Co", "Ni"])
    win = gen.get_wannier90_win_string(struct, num_wann=18)
    print("✓ win string generated:")
    for line in win.strip().split("\n"):
        print(f"  {line}")

    checks = {
        "num_wann line": "num_wann = 18" in win,
        # num_bands = ceil((18+20)/32)*32 = 64 (NPAR=ntasks=32, KPAR=NCORE=1)
        "num_bands rounded to NPAR": "num_bands = 64" in win,
        "begin projections": "begin projections" in win,
        "end projections": "end projections" in win,
        "Co projection": "Co: s,p,d" in win,
        "Ni projection": "Ni: s,p,d" in win,
        "write_hr": "write_hr = .true." in win,
        "write_xyz": "write_xyz = .true." in win,
    }
    for label, ok in checks.items():
        print(f"  ✓ {label}: {'OK' if ok else 'MISSING'}")
        if not ok:
            print(f"✗ {label} missing from win string")
            return False

    # Projection order follows structure site order.
    if win.index("Co:") > win.index("Ni:"):
        print("✗ Projections not in structure site order (Co before Ni)")
        return False
    print("✓ Projection order follows structure site order")

    print()
    return True


def test_get_input_set_side_effect():
    """Test get_input_set(): writes wannier90.win and merges the INCAR."""
    print_header("get_input_set() side effect")

    struct = get_test_structure()
    gen = Tb2jInputSetGenerator(
        elements=["Fe", "Co"],
        user_incar_settings={"ENCUT": 520},
    )

    cwd_backup = os.getcwd()
    try:
        with tempfile.TemporaryDirectory() as td:
            os.chdir(td)
            vis = gen.get_input_set(struct)

            # Side effect: wannier90.win written to the job directory (cwd).
            win_path = Path("wannier90.win")
            if not win_path.exists():
                print("✗ wannier90.win was not written to cwd")
                return False
            win = win_path.read_text(encoding="utf-8")
            print(f"✓ wannier90.win written to {win_path.resolve()}")
            if "num_wann = 18" not in win or "begin projections" not in win:
                print("✗ wannier90.win content incorrect")
                return False
            print("✓ wannier90.win content has num_wann + projections")

            incar = vis.incar
            w90_checks = {
                "NUM_WANN": incar.get("NUM_WANN") == 18,
                # NBANDS = ceil((18+20)/32)*32 = 64 (NPAR=ntasks=32)
                "NBANDS": incar.get("NBANDS") == 64,
                "LWANNIER90": incar.get("LWANNIER90") is True,
                "LWANNIER90_RUN": incar.get("LWANNIER90_RUN") is True,
                "LWRITE_MMN_AMN": incar.get("LWRITE_MMN_AMN") is True,
                "KPAR": incar.get("KPAR") == 1,
                "LREAL": incar.get("LREAL") is False,
                "NSW": incar.get("NSW") == 0,
                "IBRION": incar.get("IBRION") == -1,
                "NELM": incar.get("NELM") == 200,
                "LORBIT": incar.get("LORBIT") == 11,
                "user ENCUT": incar.get("ENCUT") == 520,
            }
            all_ok = True
            for label, ok in w90_checks.items():
                val = incar.get(label.split()[-1])
                print(f"  ✓ {label}: {val}" if ok else f"  ✗ {label}: {val}")
                all_ok = all_ok and ok
            if not all_ok:
                print("✗ Wannier90-interface INCAR merge incorrect")
                return False
            print("✓ W90 INCAR keys merged (W90 defaults + user ENCUT=520)")

            if vis.get("POTCAR") is None:
                print("✗ POTCAR missing from input set")
                return False
            if vis.get("KPOINTS") is None:
                print("✗ KPOINTS missing from input set")
                return False
            print("✓ POTCAR and KPOINTS present")
    finally:
        os.chdir(cwd_backup)

    print()
    return True


def test_incar_merging():
    """Test user_incar_settings override W90 defaults in the merged INCAR."""
    print_header("INCAR Merging (user overrides)")

    # A user may override any W90 default (e.g., KPAR) on top of the defaults.
    gen = Tb2jInputSetGenerator(
        elements=["Fe", "Co"],
        user_incar_settings={"KPAR": 2, "NELM": 300, "ENCUT": 600},
    )

    cwd_backup = os.getcwd()
    try:
        with tempfile.TemporaryDirectory() as td:
            os.chdir(td)
            vis = gen.get_input_set(get_test_structure())
            incar = vis.incar
            print(f"✓ KPAR (user): {incar.get('KPAR')}")
            print(f"✓ NELM (user): {incar.get('NELM')}")
            print(f"✓ ENCUT (user): {incar.get('ENCUT')}")
            print(f"✓ LWANNIER90 (W90 default preserved): {incar.get('LWANNIER90')}")

            if incar.get("KPAR") != 2:
                print("✗ User KPAR=2 should override W90 default KPAR=1")
                return False
            if incar.get("NELM") != 300:
                print(f"✗ User NELM=300 should win, got {incar.get('NELM')}")
                return False
            if incar.get("ENCUT") != 600:
                print(f"✗ User ENCUT=600 should win, got {incar.get('ENCUT')}")
                return False
            if incar.get("LWANNIER90") is not True:
                print("✗ W90 default LWANNIER90 should be preserved")
                return False
    finally:
        os.chdir(cwd_backup)

    print("✓ user_incar_settings win over W90 defaults; non-overridden defaults kept")
    print()
    return True


def test_get_num_bands():
    """Test get_num_bands(): NBANDS pre-rounded to a multiple of NPAR (= ntasks).

    VASP silently bumps NBANDS up to the next multiple of NPAR
    (``number of bands changed``); pre-rounding keeps the INCAR NBANDS equal to
    wannier90.win ``num_bands``. NPAR = NTASKS since KPAR=1 and NCORE=1.
    """
    print_header("get_num_bands() adaptive NBANDS")

    # num_wann=18 -> nbands_min = 38.  NBANDS = ceil(38 / ntasks) * ntasks.
    gen = Tb2jInputSetGenerator(num_wann=18)
    for ntasks, expected in [
        (16, 48),  # ceil(38/16)=3 -> 48
        (32, 64),  # ceil(38/32)=2 -> 64  (default ntasks=32)
        (48, 48),  # ceil(38/48)=1 -> 48  (never over-allocates above need)
    ]:
        gen.ntasks = ntasks
        got = gen.get_num_bands()
        print(f"✓ ntasks={ntasks}: NBANDS={got}")
        if got != expected:
            print(f"✗ ntasks={ntasks}: expected NBANDS={expected}, got {got}")
            return False
        if got % ntasks != 0:
            print(f"✗ ntasks={ntasks}: NBANDS={got} not divisible by NPAR")
            return False

    # Larger system (FCC-Fe-Mn, num_wann=36 -> nbands_min=56) at ntasks=48:
    # ceil(56/48)=2 -> 96. This is the former 48/96 parallelisation pair.
    gen2 = Tb2jInputSetGenerator(num_wann=36, ntasks=48)
    if gen2.get_num_bands() != 96:
        print(f"✗ num_wann=36, ntasks=48: expected 96, got {gen2.get_num_bands()}")
        return False
    print("✓ num_wann=36, ntasks=48: NBANDS=96 (ceil(56/48)=2*48)")

    # Explicit num_wann argument beats the instance num_wann.
    gen3 = Tb2jInputSetGenerator(num_wann=18, ntasks=32)
    if gen3.get_num_bands(num_wann=36) != 64:  # ceil(56/32)=2 -> 64
        print(f"✗ explicit num_wann arg should win, got {gen3.get_num_bands(36)}")
        return False
    print("✓ explicit num_wann arg wins over instance num_wann")

    # NBANDS must be >= num_wann + num_bands_extra.
    if gen.get_num_bands() < 18 + gen.num_bands_extra:
        print("✗ NBANDS < num_wann + num_bands_extra")
        return False
    print("✓ NBANDS >= num_wann + num_bands_extra (20)")

    print()
    return True


def test_parse_exchange_out():
    """Test parse_exchange_out() on real exchange.out files (read-only)."""
    print_header("parse_exchange_out (real files)")

    # (name, expected J-pair count, expected atom-moment count)
    references = [
        ("SER-Co", 2914, 2),
        ("BCC-Co-Ni", 2914, 2),
        ("FCC-Fe-Mn", 11660, 4),
    ]
    for name, exp_pairs, exp_moments in references:
        parsed = parse_reference(name)
        print(f"  {name}: {parsed['num_j_pairs']} J pairs, "
              f"{len(parsed['atom_moments'])} moments")

        if parsed["num_j_pairs"] != exp_pairs:
            print(f"✗ {name}: expected {exp_pairs} J pairs, got {parsed['num_j_pairs']}")
            return False
        if len(parsed["atom_moments"]) != exp_moments:
            print(f"✗ {name}: expected {exp_moments} moments, got {len(parsed['atom_moments'])}")
            return False
        if len(parsed["j_tensors"]) != exp_pairs:
            print("✗ j_tensors must be parallel to j_pairs")
            return False
        if parsed["cell"] is None or len(parsed["cell"]) != 3:
            print(f"✗ {name}: cell not parsed")
            return False

        pair = parsed["j_pairs"][0]
        for key in ("i", "j", "R", "J_iso_meV", "vector", "distance"):
            if key not in pair:
                print(f"✗ {name}: pair missing key '{key}'")
                return False

        tensor = parsed["j_tensors"][0]
        if len(tensor) != 3 or any(len(row) != 3 for row in tensor):
            print(f"✗ {name}: first tensor not 3x3: {tensor}")
            return False

        moment = parsed["atom_moments"][0]
        for key in ("atom", "xyz", "w_charge", "w_magmom"):
            if key not in moment:
                print(f"✗ {name}: moment missing key '{key}'")
                return False

    print("✓ j_pairs/j_tensors/atom_moments/cell parsed from all 3 reference files")
    print()
    return True


def test_tb2j_result_from_dict():
    """Test Tb2jResult.from_dict() on real data + validate BCC Co-Ni NN J."""
    print_header("Tb2jResult.from_dict (real data)")

    parsed = parse_reference("BCC-Co-Ni")
    result = Tb2jResult.from_dict(
        {
            "elements": ["Co", "Ni"],
            "num_atoms": 2,
            "efermi": 5.02625994,
            "num_wann": 18,
            "kmesh": [9, 9, 9],
            "j_pairs": parsed["j_pairs"],
            "j_tensors": parsed["j_tensors"],
            "atom_moments": parsed["atom_moments"],
            "exchange_out_path": str(get_exchange_out_path("BCC-Co-Ni")),
            "warnings": [],
        }
    )
    print(f"✓ elements: {result.elements}")
    print(f"✓ num_atoms: {result.num_atoms}, num_wann: {result.num_wann}")
    print(f"✓ efermi: {result.efermi} eV, kmesh: {result.kmesh}")
    print(f"✓ j_pairs: {len(result.j_pairs)}, j_tensors: {len(result.j_tensors)}")

    if result.elements != ["Co", "Ni"]:
        print("✗ elements not round-tripped")
        return False
    if result.num_atoms != 2:
        print("✗ num_atoms should be 2")
        return False
    if result.num_wann != 18:
        print(f"✗ num_wann should be 18, got {result.num_wann}")
        return False
    if len(result.j_pairs) != 2914:
        print(f"✗ expected 2914 j_pairs, got {len(result.j_pairs)}")
        return False
    if len(result.j_tensors) != 2914:
        print(f"✗ expected 2914 j_tensors, got {len(result.j_tensors)}")
        return False

    # Nearest-neighbour Co-Ni shell of BCC (a=2.81): distance 2.435 A, J ~ +4.82 meV.
    nn_pairs = [p for p in result.j_pairs if p["distance"] == 2.435]
    if len(nn_pairs) != 16:
        print(f"✗ BCC Co-Ni NN shell should have 16 pairs, got {len(nn_pairs)}")
        return False
    mean_nn_j = sum(p["J_iso_meV"] for p in nn_pairs) / len(nn_pairs)
    print(f"✓ NN shell: {len(nn_pairs)} pairs at d=2.435 A, mean J = {mean_nn_j:.4f} meV")
    if not 4.82 <= mean_nn_j <= 4.83:
        print(f"✗ BCC Co-Ni NN J expected ≈ +4.82 meV, got {mean_nn_j:.4f}")
        return False

    # First tensor is 3x3 and (for BCC Co-Ni) diagonal, trace ≈ 3 * J_iso.
    nn0 = min(result.j_pairs, key=lambda p: p["distance"])
    idx = result.j_pairs.index(nn0)
    tensor = result.j_tensors[idx]
    if len(tensor) != 3 or any(len(row) != 3 for row in tensor):
        print(f"✗ NN tensor not 3x3: {tensor}")
        return False
    trace = sum(tensor[i][i] for i in range(3))
    print(f"✓ NN tensor 3x3, trace={trace:.3f} ≈ 3*J_iso={3*nn0['J_iso_meV']:.3f}")
    if abs(trace - 3 * nn0["J_iso_meV"]) > 0.05:
        print(f"✗ NN tensor trace should track 3*J_iso")
        return False

    print()
    return True


def test_to_summary():
    """Test Tb2jResult.to_summary() output structure and shell stats."""
    print_header("to_summary()")

    # SER-Co (HCP Co): 2 atoms, 2914 J pairs.
    parsed = parse_reference("SER-Co")
    result = Tb2jResult.from_dict(
        {
            "elements": ["Co"],
            "num_atoms": 2,
            "j_pairs": parsed["j_pairs"],
            "j_tensors": parsed["j_tensors"],
            "warnings": ["test warning"],
        }
    )
    summary = result.to_summary()
    print(f"✓ Summary keys: {list(summary.keys())}")

    if summary["num_j_pairs"] != 2914:
        print(f"✗ num_j_pairs in summary wrong: {summary['num_j_pairs']}")
        return False
    if summary["num_j_tensors"] != 2914:
        print(f"✗ num_j_tensors in summary wrong: {summary['num_j_tensors']}")
        return False
    if summary["num_shells"] <= 0:
        print("✗ num_shells should be > 0")
        return False
    if len(summary["shells"]) > 10:
        print("✗ summary should cap shells at 10")
        return False

    shell0 = summary["shells"][0]
    for key in ("distance", "count", "mean_J_meV", "max_abs_J_meV"):
        if key not in shell0:
            print(f"✗ shell entry missing '{key}'")
            return False
    print(f"✓ NN shell: d={shell0['distance']} A, count={shell0['count']}, "
          f"mean J={shell0['mean_J_meV']} meV")

    if summary.get("max_abs_J_meV") is None:
        print("✗ max_abs_J_meV missing from summary")
        return False
    if "strongest_pair" not in summary:
        print("✗ strongest_pair missing from summary")
        return False
    print(f"✓ max|J| = {summary['max_abs_J_meV']} meV")

    if summary.get("warnings") != ["test warning"]:
        print("✗ warnings not propagated to summary")
        return False
    print("✓ warnings propagated")

    print()
    return True


def test_flow_creation():
    """Test Tb2jWorker._make_flow(): R3 relax + W90/TB2J chain, output = tb2j solve."""
    print_header("Flow Creation (R3 relax -> W90 -> TB2J)")

    worker = make_test_worker()  # elements=None -> auto-derived from structure
    structure = get_test_structure(("Fe", "Co"))
    flow = worker._make_flow(structure)

    print(f"✓ Flow: {flow.name} with {len(flow.jobs)} top-level jobs")

    if len(flow.jobs) != 2:
        print(f"✗ Expected 2 top-level jobs, got {len(flow.jobs)}")
        return False

    relax_flow = flow.jobs[0]
    if relax_flow.name != "r3 relax":
        print(f"✗ First job should be 'r3 relax', got {relax_flow.name}")
        return False
    relax_jobs = list(relax_flow.jobs) if hasattr(relax_flow, "jobs") else [relax_flow]
    print(f"✓ relax flow: {relax_flow.name} -> {[j.name for j in relax_jobs]}")
    if len(relax_jobs) != 1:
        print("✗ Relax should contain a single R3 relax job")
        return False

    tb2j_flow = flow.jobs[1]
    tb2j_jobs = list(tb2j_flow.jobs)
    job_names = [j.name for j in tb2j_jobs]
    print(f"✓ tb2j flow: {tb2j_flow.name} -> {job_names}")
    if job_names != ["w90 static", "tb2j solve"]:
        print(f"✗ Expected ['w90 static', 'tb2j solve'], got {job_names}")
        return False

    # Flow output must point at the tb2j solve job (output reference).
    solve_job = tb2j_jobs[1]
    if flow.output != solve_job.output:
        print("✗ Flow output does not point to the tb2j solve output")
        return False
    if flow.output.uuid != solve_job.uuid:
        print("✗ Flow output uuid mismatch")
        return False
    print("✓ Flow output points to tb2j solve")

    # num_wann auto-derived from the input structure (invariant under relaxation).
    if worker.tb2j_maker.input_set_generator.num_wann != 18:
        print(f"✗ num_wann should be 18 (2 atoms * 9), got "
              f"{worker.tb2j_maker.input_set_generator.num_wann}")
        return False
    if worker.tb2j_maker.elements != ["Fe", "Co"]:
        print(f"✗ elements auto-derived from structure expected ['Fe','Co'], "
              f"got {worker.tb2j_maker.elements}")
        return False
    print("✓ elements + num_wann auto-derived from structure")

    # ntasks threads Worker -> Maker -> InputSetGenerator (default 32).
    gen = worker.tb2j_maker.input_set_generator
    if worker.ntasks != 32:
        print(f"✗ worker.ntasks default should be 32, got {worker.ntasks}")
        return False
    if gen.ntasks != 32:
        print(f"✗ input_set_generator.ntasks should be 32, got {gen.ntasks}")
        return False
    if worker.ntasks != gen.ntasks:
        print("✗ ntasks not threaded Worker -> InputSetGenerator")
        return False
    print("✓ ntasks=32 threaded Worker -> Maker -> InputSetGenerator")

    print()
    return True


def test_structure_names():
    """Test main_w90.STRUCTURE_NAMES: 30 endmembers in SER/BCC/FCC format."""
    print_header("main_w90 STRUCTURE_NAMES")

    import main_w90

    names = main_w90.STRUCTURE_NAMES
    print(f"✓ {len(names)} structure names")

    if len(names) != 30:
        print(f"✗ Expected 30 structure names, got {len(names)}")
        return False
    if len(set(names)) != 30:
        print("✗ STRUCTURE_NAMES contains duplicates")
        return False

    if any(not n.startswith(("SER-", "BCC-", "FCC-")) for n in names):
        print("✗ All names must start with SER-/BCC-/FCC-")
        return False

    ser = [n for n in names if n.startswith("SER-")]
    bcc = [n for n in names if n.startswith("BCC-")]
    fcc = [n for n in names if n.startswith("FCC-")]
    print(f"✓ SER={len(ser)}  BCC={len(bcc)}  FCC={len(fcc)}")
    if len(ser) != 4 or len(bcc) != 10 or len(fcc) != 16:
        print("✗ Expected SER=4, BCC=10, FCC=16")
        return False

    expected_ser = {"SER-Co", "SER-Fe", "SER-Mn", "SER-Ni"}
    if set(ser) != expected_ser:
        print(f"✗ SER names wrong: {ser}")
        return False

    # Format check: prefix + capitalized element symbol tokens.
    for n in names:
        prefix, *elems = n.split("-")
        if prefix not in ("SER", "BCC", "FCC"):
            print(f"✗ Bad prefix in {n}")
            return False
        if not elems or any(not (e.isalpha() and e[0].isupper()) for e in elems):
            print(f"✗ Bad element token in {n}")
            return False
    print("✓ Name format valid (PREFIX-Elem[-Elem])")

    print()
    return True


# =============================================================================
# Test Runner
# =============================================================================


def run_unit_tests():
    """Run all unit tests."""
    print("\n" + "=" * 50)
    print("TB2J Unit Tests")
    print("=" * 50 + "\n")

    tests = [
        ("Imports", test_imports),
        ("WorkerCreation", test_worker_creation),
        ("InputSetGenerator", test_input_set_generator),
        ("WinStringContent", test_win_string_content),
        ("GetInputSetSideEffect", test_get_input_set_side_effect),
        ("IncarMerging", test_incar_merging),
        ("GetNumBands", test_get_num_bands),
        ("ParseExchangeOut", test_parse_exchange_out),
        ("Tb2jResultFromDict", test_tb2j_result_from_dict),
        ("ToSummary", test_to_summary),
        ("FlowCreation", test_flow_creation),
        ("StructureNames", test_structure_names),
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

TB2J_GLOBAL_INCAR = {
    "ENCUT": 520,
    "PREC": "Accurate",
    "EDIFF": 1e-5,
    "EDIFFG": -0.02,
    "ISPIN": 2,
    "MAGMOM": {"Co": 3.0, "Fe": 5.0, "Mn": 5.0, "Ni": 2.0},
    "ISMEAR": 1,
    "SIGMA": 0.2,
    "ALGO": "Fast",
    "LREAL": False,
    "ISYM": 0,
    "GGA": "PE",
}


def run_locally(clean: bool = True):
    """Run the TB2J workflow locally (requires VASP + TB2J + Wannier90).

    Args:
        clean: If True, remove existing flow_dir before running.
               If False, keep existing data and use resume=True.
    """
    structure = get_test_structure(("Fe", "Co"))
    flow_name = "FeCo-tb2j"
    flow_dir = Path("/nfs_ssd/tmp")
    store_dir = Path("./temp") / flow_name
    json_path = store_dir / f"{flow_name}.json"

    if clean and flow_dir.exists():
        shutil.rmtree(flow_dir)
        log.info(f"Cleaned existing flow_dir: {flow_dir}")

    worker = Tb2jWorker(
        vasp_args=CLUSTER_VASP_ARGS,
        global_incar=TB2J_GLOBAL_INCAR,
        elements=["Fe", "Co"],
        kmesh=(9, 9, 9),
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
        summary = Tb2jResult.from_dict(output).to_summary()

        print("\n" + "=" * 50)
        print("TB2J Results Summary")
        print("=" * 50)
        print(f"J pairs: {summary['num_j_pairs']}")
        print(f"Shells: {summary['num_shells']}")
        if summary.get("shells"):
            nn = summary["shells"][0]
            print(f"NN shell: d={nn['distance']:.3f} A, count={nn['count']}, "
                  f"mean J={nn['mean_J_meV']} meV")
        if summary.get("max_abs_J_meV") is not None:
            print(f"max|J|: {summary['max_abs_J_meV']} meV")
        if "warnings" in summary:
            print(f"Warnings: {summary['warnings']}")
    else:
        log.error("Workflow failed")


def submit_job():
    """Submit the TB2J workflow to Slurm with abort/rerun simulation."""
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
    parser = argparse.ArgumentParser(description="Test TB2J workflow")
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
