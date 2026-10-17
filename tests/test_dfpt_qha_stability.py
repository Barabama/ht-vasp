"""Tier-A and tier-B diagnostics on real DFPT force constants."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar

from htvasp.dfpt_qha.phonons import phonopy_for
from htvasp.dfpt_qha.stability import (
    commensurate_qpoints,
    diagnose_stability,
    force_constants_checks,
    force_constants_sha256,
    reciprocal_check,
)
from htvasp.dfpt_qha.jobs import force_constants_from_dynmat, parse_dynmat

BASE = Path("/nfs_ssd/tmp/BCC-Al-Ta-2-5e3275de/1-DFPT_volume_1/5")
POOL = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0")
MATRIX = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]


def _fixture():
    if not (BASE / "DYNMAT").is_file():
        pytest.skip("real DFPT fixture unavailable")
    import gzip

    with gzip.open(BASE / "POSCAR.gz", "rt") as handle:
        supercell = Poscar.from_str(handle.read()).structure
    cell = Structure(
        np.array(supercell.lattice.matrix) / 2,
        ["Al", "Ta"],
        [[0, 0, 0], [0.5, 0.5, 0.5]],
    )
    fc = force_constants_from_dynmat(
        cell, supercell, parse_dynmat(BASE / "DYNMAT")[0], supercell_matrix=MATRIX
    )
    recipe = json.loads((POOL / "pool_recipes.json").read_text())["endpoints"][
        "BCC-Al-Ta-2"
    ]["recipe"]
    ph = phonopy_for(cell, MATRIX)
    ph.force_constants = fc
    return cell, fc, ph, recipe


def test_commensurate_points_for_two_by_two_by_two():
    points = commensurate_qpoints(MATRIX)
    assert len(points) == 8
    assert (0.5, 0.5, 0.0) in points
    assert (0.5, 0.5, 0.5) in points


def test_force_constant_checks_real_fixture():
    _, fc, _, _ = _fixture()
    checks = force_constants_checks(fc)
    assert checks["shape_ok"] is True
    assert checks["asr_ok"] is True
    assert checks["asr_rel"] < 1e-3
    assert checks["permutation_ok"] is True
    assert len(checks["force_constants_sha256"]) == 64


def test_reciprocal_identity_real_fixture():
    _, _, ph, _ = _fixture()
    assert reciprocal_check(ph)["ok"] is True


def test_real_alta_is_tier_b_finite_q_soft_mode():
    cell, fc, ph, recipe = _fixture()
    report = diagnose_stability(
        ph, fc, structure=cell, recipe=recipe, supercell_matrix=MATRIX, mesh=(12, 12, 12)
    )
    assert report["status"] == "finite_q_soft_mode_diagnostic"
    assert report["qha_allowed"] is False
    soft = [q for q in report["commensurate_qpoints"] if q["n_imaginary"]]
    assert any(q["q_fractional"] == [0.0, 0.5, 0.5] for q in soft)
    assert report["checks"]["gamma"]["noise_ok"] is True
    assert report["mesh_diagnostic"]["note"].startswith("mesh minimum")
    json.dumps(report)


def test_tier_b_hash_is_stable():
    _, fc, _, _ = _fixture()
    assert force_constants_sha256(fc) == force_constants_sha256(np.array(fc, dtype="<f8"))


def test_shape_or_asr_failure_is_not_a_physical_soft_mode():
    _, fc, _, recipe = _fixture()
    bad = np.array(fc)
    bad[0, 0, 0, 0] += 100.0
    checks = force_constants_checks(bad)
    assert checks["asr_ok"] is False
    # diagnosis receives the same Phonopy object only to exercise the API; the
    # classification must be numerical_unresolved, not a soft-mode conclusion.
    cell, _, ph, _ = _fixture()
    report = diagnose_stability(
        ph, bad, structure=cell, recipe=recipe, supercell_matrix=MATRIX, mesh=(4, 4, 4)
    )
    assert report["status"] == "numerical_unresolved"
    assert report["qha_allowed"] is False
