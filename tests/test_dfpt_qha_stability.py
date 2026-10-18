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
    DfptStabilityReport,
    commensurate_qpoints,
    diagnose_stability,
    force_constants_checks,
    force_constants_sha256,
    gamma_acoustic_report,
    reciprocal_check,
    unresolved_report,
)
from htvasp.dfpt_qha.flow_jobs import classify_dfpt_stability
from htvasp.dfpt_qha.jobs import (
    collect_dfpt_task_doc,
    force_constants_from_dynmat,
    parse_dynmat,
)

BASE = Path(__file__).parent / "fixtures" / "BCC-Al-Ta-2-dfpt-volume-1"
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
    cell, _, ph, _ = _fixture()
    points = commensurate_qpoints(ph)
    assert len(points) == 8
    assert (0.5, 0.5, 0.0) in points
    assert (0.5, 0.5, 0.5) in points


def test_commensurate_points_follow_the_primitive_cell_not_the_recipe_cell():
    """fcc A=A: a 4-atom conventional cell with a 1-atom primitive.

    The 2x2x2 recipe supercell holds 32 primitive cells, so 32 q points are
    exactly representable.  Computing them from the recipe matrix instead gives
    8, and the other 24 -- including (1/4,1/4,1/4) -- are never tested.
    """

    fcc = Structure(
        [[4.0, 0, 0], [0, 4.0, 0], [0, 0, 4.0]], ["Al"] * 4,
        [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
    )
    phonon = phonopy_for(fcc, MATRIX)
    points = commensurate_qpoints(phonon)
    assert len(points) == 32
    # a quarter of the way along the primitive reciprocal axes: representable
    # by 32 primitive cells, invisible to a {0, 1/2}^3 grid
    assert any(abs(x[0]) < 1e-9 and abs(x[1] - 0.25) < 1e-9 for x in points)
    cell, _, b2, _ = _fixture()
    assert len(commensurate_qpoints(b2)) == 8


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
    assert report["status"] == "finite_q_soft_mode"
    assert report["qha_allowed"] is False
    soft = [q for q in report["commensurate_qpoints"] if q["n_imaginary"]]
    assert any(q["q_fractional"] == [0.0, 0.5, 0.5] for q in soft)
    assert report["checks"]["gamma"]["noise_ok"] is True
    assert report["mesh_diagnostic"]["note"].startswith("interpolated diagnostic")
    # -3.59 THz is far past the borderline threshold: no larger supercell will
    # turn this into a truncation artifact
    assert report["confirmation_required"] is False
    assert report["supercell_convergence"] == "not_checked"
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


def test_report_model_rejects_a_contradictory_verdict():
    """``qha_allowed`` is derived from ``status``, not carried alongside it."""

    report = {
        "schema_version": "htvasp-dfpt-stability/v1",
        "status": "finite_q_soft_mode",
        "qha_allowed": True,
        "reason": "contradiction",
        "checks": {},
        "commensurate_qpoints": [],
        "mesh_diagnostic": {},
        "force_constants_sha256": "a" * 64,
    }
    with pytest.raises(ValueError, match="contradicts status"):
        DfptStabilityReport.model_validate(report)


def test_unresolved_report_is_a_usable_artifact():
    report = unresolved_report("phonopy would not build", force_constants=np.zeros((2,)))
    parsed = DfptStabilityReport.model_validate(report)
    assert parsed.status == "numerical_unresolved"
    assert parsed.qha_allowed is False
    assert len(parsed.force_constants_sha256) == 64
    json.dumps(parsed.model_dump(mode="json"))


def test_a_borderline_soft_mode_asks_for_a_larger_supercell():
    """The agreed policy: strong soft modes die, borderline ones get confirmed.

    Scaling the force constants by s scales every frequency by sqrt(s).  The
    real BCC-Al-Ta-2 N-point mode at -3.59 THz becomes -0.5 THz, which is
    inside the 1 THz confirmation window -- same physics, different verdict.
    """

    cell, fc, ph, recipe = _fixture()
    scale = (0.5 / 3.59) ** 2
    softened = np.array(fc) * scale
    phonon = phonopy_for(cell, MATRIX)
    phonon.force_constants = softened
    report = diagnose_stability(
        phonon, softened, structure=cell, recipe=recipe,
        supercell_matrix=MATRIX, mesh=(4, 4, 4),
    )
    assert report["status"] == "finite_q_soft_mode"
    assert report["qha_allowed"] is False
    assert report["confirmation_required"] is True
    assert "larger supercell" in report["confirmation_reason"]


def test_classify_endpoint_blocks_a_soft_volume_and_keeps_the_artifact():
    """One soft volume fails the whole endpoint, and says which one."""

    cell, fc, _, recipe = _fixture()
    collected = {
        "volumes": [
            {
                "index": 0,
                "cell": cell.as_dict(),
                "supercell": phonopy_supercell(cell).as_dict(),
                "volume_angstrom3": float(cell.volume),
                "static_energy_eV": -10.0,
                "force_constants": np.array(fc).tolist(),
                "n_atoms_supercell": len(fc),
                "ordering_fingerprint": "x" * 64,
            }
        ],
        "rejected_volume_indices": [],
        "n_requested": 1,
        "n_validated": 1,
    }
    result = classify_dfpt_stability(
        collected, recipe=recipe, supercell_matrix=MATRIX, mesh=(4, 4, 4)
    ).run(_Store())
    output = result.output
    assert result.stop_children is True
    assert output["ok"] is False
    assert output["blocked_volume_indices"] == [0]
    assert output["reports"][0]["status"] == "finite_q_soft_mode"
    assert output["reports"][0]["commensurate_qpoints"]
    json.dumps(output)


@pytest.mark.skipif(
    not Path(
        "/nfs_hdd/2025/gaominliang/ht-vasp/data/endpoints/BCC-Nb-Nb-2"
        "/dfpt-qhaflow/runs/dfpt-qha-f23783fe35518cb42436/store.json"
    ).is_file(),
    reason="archived BCC-Nb-Nb-2 run unavailable",
)
def test_archived_nb_nb_passes_a_gate_that_used_to_reject_it():
    """The regression the three-tier rule exists for.

    BCC-Nb-Nb-2 volume 1 has a 12x12x12 mesh minimum of -0.0238 THz at
    q = (-1/24, 1/24, 1/24) while every exact commensurate q point is stable.
    The old mesh gate rejected the endpoint; the interpolation artifact must now
    be recorded as a warning and the endpoint must pass.
    """

    import json as _json

    run = Path(
        "/nfs_hdd/2025/gaominliang/ht-vasp/data/endpoints/BCC-Nb-Nb-2"
        "/dfpt-qhaflow/runs/dfpt-qha-f23783fe35518cb42436"
    )
    documents = _json.loads((run / "store.json").read_text())
    collected = next(
        d for d in documents if d.get("name") == "collect validated volumes"
    )["output"]
    recipe = _json.loads((POOL / "pool_recipes.json").read_text())["endpoints"][
        "BCC-Nb-Nb-2"
    ]["recipe"]

    result = classify_dfpt_stability(
        collected, recipe=recipe, supercell_matrix=MATRIX, mesh=(12, 12, 12)
    ).run(_Store())
    output = result.output
    assert result.stop_children is False
    assert output["ok"] is True
    assert output["blocked_volume_indices"] == []
    shallow = [
        r for r in output["reports"]
        if (r["mesh_diagnostic"] or {}).get("mesh_min_frequency_THz", 0.0) < 0.0
    ]
    assert shallow, "the archived run no longer reproduces the artifact"
    for report in shallow:
        assert report["status"] == "harmonic_QHA_applicable"
        assert any("interpolation artifact" in w for w in report["warnings"])


def phonopy_supercell(cell):
    from htvasp.dfpt_qha.phonons import supercell_structure

    return supercell_structure(cell, MATRIX)


def _Store():
    from jobflow import JobStore
    from maggma.stores import MemoryStore

    docs = MemoryStore()
    docs.connect()
    return JobStore(docs)


def _dimer_model(cell, phonon, spring):
    """A force-constant matrix with a flat optical branch of controlled sign.

    Every atom couples to every atom of the other sublattice and to nothing
    else, so the acoustic branches are exactly flat at zero (the coupling is
    q-independent) and the optical branches have no dispersion.  With
    ``spring > 0`` the optical branch is imaginary at Gamma and at every other
    q -- which is what makes it a test of the *tier order*: Gamma must win.
    """

    mapping = np.array(phonon.primitive.s2p_map)
    parts = {value: index for index, value in enumerate(np.unique(mapping))}
    sublattice = np.array([parts[value] for value in mapping])
    n = len(sublattice)
    n_other = {
        0: int((sublattice == 1).sum()), 1: int((sublattice == 0).sum())
    }
    fc = np.zeros((n, n, 3, 3))
    for i in range(n):
        part = sublattice[i]
        for j in range(n):
            if sublattice[j] != part:
                fc[i, j] = spring * np.eye(3)
        fc[i, i] -= n_other[part] * spring * np.eye(3)
    return fc


def test_gamma_optical_soft_mode_is_not_swallowed_by_the_acoustic_slice():
    """phonopy sorts frequencies ascending, so an imaginary optical branch
    sorts BELOW the three ~0 THz acoustic modes.

    Slicing ``frequencies[:3]`` as "acoustic" therefore hands the instability
    the acoustic noise allowance and reports the volume as usable.  The
    branches are identified by eigenvector character instead.
    """

    cell, _, phonon, recipe = _fixture()
    model = _dimer_model(cell, phonon, spring=1.0)
    phonon.force_constants = model
    gamma = gamma_acoustic_report(phonon)
    frequencies = np.array(gamma["frequency_THz"])
    # the imaginary optical branch sorts first ...
    assert frequencies[0] < -1.0
    # ... but it is NOT one of the acoustic branches
    assert 0 not in gamma["acoustic_bands"]
    assert max(abs(x) for x in gamma["acoustic_frequencies_THz"]) < 0.1
    assert gamma["optical_soft_THz"]

    report = diagnose_stability(
        phonon, model, structure=cell, recipe=recipe,
        supercell_matrix=MATRIX, mesh=(4, 4, 4),
    )
    assert report["status"] == "gamma_soft_mode"
    assert report["qha_allowed"] is False


def test_a_deep_mesh_minimum_away_from_every_exact_q_is_suspected():
    """The tier that catches an instability the supercell cannot represent.

    A 2x2x2 supercell represents a finite set of q points.  A branch far below
    zero on a dense mesh at a q no commensurate point covers is neither a
    truncation artifact nor an exact-q instability: the cell is too small to
    decide, and QHA must not proceed.  The matrix here is stable at every
    exact q (the dimer model's flat optical branch is positive, its acoustic
    branches are exactly zero), so only the patched mesh minimum can fire.
    """

    from unittest.mock import patch

    cell, _, phonon, recipe = _fixture()
    model = _dimer_model(cell, phonon, spring=-1.0)
    phonon.force_constants = model
    fake = {
        "mesh": [12, 12, 12], "mesh_min_qpoint": [0.4167, 0.0, 0.0],
        "mesh_min_frequency_THz": -1.5, "n_mesh_imaginary_modes": 3,
        "n_modes": 1728 * 3, "note": "interpolated diagnostic only",
    }
    with patch("htvasp.dfpt_qha.stability.mesh_diagnostic", return_value=fake):
        report = diagnose_stability(
            phonon, model, structure=cell, recipe=recipe,
            supercell_matrix=MATRIX, mesh=(12, 12, 12),
        )
    assert report["status"] == "incommensurate_soft_mode_suspected"
    assert report["qha_allowed"] is False
    assert "no commensurate q point" in report["reason"]


def test_a_shallow_mesh_minimum_is_recorded_and_does_not_block():
    """The BCC-Nb-Nb-2 shape: -0.024 THz at q=(-1/24,1/24,1/24), stable at
    every exact q point.  The endpoint must pass with a warning."""

    from unittest.mock import patch

    cell, _, phonon, recipe = _fixture()
    model = _dimer_model(cell, phonon, spring=-1.0)
    phonon.force_constants = model
    fake = {
        "mesh": [12, 12, 12], "mesh_min_qpoint": [-0.0417, 0.0417, 0.0417],
        "mesh_min_frequency_THz": -0.0238, "n_mesh_imaginary_modes": 1,
        "n_modes": 1728 * 3, "note": "interpolated diagnostic only",
    }
    with patch("htvasp.dfpt_qha.stability.mesh_diagnostic", return_value=fake):
        report = diagnose_stability(
            phonon, model, structure=cell, recipe=recipe,
            supercell_matrix=MATRIX, mesh=(12, 12, 12),
        )
    assert report["status"] == "harmonic_QHA_applicable"
    assert report["qha_allowed"] is True
    assert any("interpolation artifact" in w for w in report["warnings"])
