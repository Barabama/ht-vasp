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
    reciprocal_check,
    unresolved_report,
)
from htvasp.dfpt_qha.flow_jobs import classify_dfpt_stability
from htvasp.dfpt_qha.jobs import (
    collect_dfpt_task_doc,
    force_constants_from_dynmat,
    parse_dynmat,
)

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
