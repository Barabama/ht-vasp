"""Offline tests for the canonical TB2J magnetic contract."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from htvasp.tb2j.jobs import canonicalize_tb2j_result, parse_exchange_out
from htvasp.tb2j.magnetic_contract import (
    K_B_MEV_PER_K,
    MAGNETIC_CONTRACT_SCHEMA_VERSION,
)
from htvasp.tb2j.task_doc import Tb2jResult


FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "data/endmembers/BCC-Al-Nb/tb2jflow/3-tb2j_solve/exchange.out"
)


def _context() -> dict:
    return {
        "phase": "B2",
        "type_character": "AB",
        "sublattice_ratios": [1, 1],
        "site_to_sublattice": {"A1": 0, "B1": 1},
        "endmember_components": [["A"], ["B"]],
        "order_num": 0,
        "structure_role": "ordered_endmember",
        "tc_model": {
            "method": "mfa_sublattice",
            "transition_type": "auto",
            "cutoff": "all_pairs",
            "k_B": K_B_MEV_PER_K,
        },
        "bmagn": {
            "source": "atom_moments.w_magmom",
            "unit": "mu_B",
            "aggregation": "mean_abs",
            "denominator": "structure_sites",
        },
    }


def _payload(j: float, moments=(2.0, -3.0)) -> dict:
    return {
        "elements": ["A", "B"],
        "num_atoms": 2,
        "structure": {
            "lattice": {"matrix": [[3, 0, 0], [0, 3, 0], [0, 0, 3]]},
            "sites": [
                {"label": "A1", "species": [{"element": "A", "occu": 1}]},
                {"label": "B1", "species": [{"element": "B", "occu": 1}]},
            ],
        },
        "atom_moments": [
            {"atom": "A1", "w_magmom": moments[0]},
            {"atom": "B1", "w_magmom": moments[1]},
        ],
        "j_pairs": [
            {"i": "A1", "j": "B1", "R": [0, 0, 0], "J_iso_meV": j, "vector": [1, 1, 1], "distance": 1.0},
            {"i": "B1", "j": "A1", "R": [0, 0, 0], "J_iso_meV": j, "vector": [-1, -1, -1], "distance": 1.0},
        ],
        "j_tensors": [
            [[j, 0, 0], [0, j, 0], [0, 0, j]],
            [[j, 0, 0], [0, j, 0], [0, 0, j]],
        ],
        "warnings": [],
    }


def test_missing_context_is_needs_input_and_not_production_ready():
    result = canonicalize_tb2j_result(_payload(2.0))
    assert result["magnetic_contract"]["schema_version"] == MAGNETIC_CONTRACT_SCHEMA_VERSION
    assert result["magnetic_observables"]["status"] == "needs_input"
    assert result["audit"]["production_ready"] is False
    assert result["audit"]["contract_status"] == "raw_reader_only"
    assert result["audit"]["review_required"] is True


def test_b2_two_sublattice_fm_uses_mfa_tc():
    result = canonicalize_tb2j_result(_payload(2.0), magnetic_context=_context())
    observables = result["magnetic_observables"]
    assert observables["status"] == "computed"
    assert observables["transition_type"] == "Curie"
    assert observables["TC"] == pytest.approx((2.0 / 3.0) * 2.0 / K_B_MEV_PER_K)
    assert observables["TN"] is None
    assert observables["exchange_matrix"] == [[0.0, 2.0], [2.0, 0.0]]
    assert max(observables["eigenvalues"]) == pytest.approx(2.0)
    assert observables["conversion_trace"]["pair_count_used"] == 2


def test_b2_two_sublattice_afm_is_neel_and_keeps_negative_j():
    result = canonicalize_tb2j_result(_payload(-2.0), magnetic_context=_context())
    observables = result["magnetic_observables"]
    assert observables["status"] == "computed"
    assert observables["transition_type"] == "Neel"
    assert observables["TC"] is None
    assert observables["TN"] == pytest.approx((2.0 / 3.0) * 2.0 / K_B_MEV_PER_K)
    assert observables["exchange_matrix"] == [[0.0, -2.0], [-2.0, 0.0]]


def test_zero_local_moments_are_nonmagnetic():
    result = canonicalize_tb2j_result(
        _payload(2.0, moments=(0.0, 0.0)), magnetic_context=_context()
    )
    assert result["magnetic_observables"]["status"] == "nonmagnetic"
    assert result["magnetic_observables"]["TC"] is None
    assert result["audit"]["production_ready"] is True


def test_bmagn_uses_explicit_structure_site_denominator():
    result = canonicalize_tb2j_result(_payload(2.0), magnetic_context=_context())
    assert result["magnetic_observables"]["BMAGN"] == pytest.approx(2.5)
    assert result["magnetic_contract"]["bmagn"]["denominator"] == "structure_sites"
    assert "BMAGN=mean(abs(w_magmom))" in " ".join(result["audit"]["verified_checks"])


def test_json_round_trip_and_contract_provenance():
    result = canonicalize_tb2j_result(_payload(2.0), magnetic_context=_context())
    restored = Tb2jResult.model_validate(json.loads(json.dumps(result)))
    assert restored.magnetic_contract["schema_version"] == MAGNETIC_CONTRACT_SCHEMA_VERSION
    assert restored.magnetic_contract["provenance"]["context_provided"] is True
    assert restored.magnetic_observables["TC"] == pytest.approx(result["magnetic_observables"]["TC"])
    assert restored.audit["production_ready"] is True


def test_existing_parser_fixture_does_not_regress():
    if not FIXTURE.is_file():
        pytest.skip(f"fixture not available: {FIXTURE}")
    parsed = parse_exchange_out(FIXTURE)
    assert parsed["num_j_pairs"] > 0
    assert len(parsed["j_pairs"]) == len(parsed["j_tensors"])
    assert parsed["units"]["J_iso"] == "meV"
    assert parsed["units"]["w_magmom"] == "mu_B"
    assert parsed["audit"]["production_ready"] is False


def test_structure_and_exchange_atom_order_can_differ():
    payload = _payload(2.0)
    payload["structure"]["sites"] = [
        {"label": "B1", "species": [{"element": "B", "occu": 1}]},
        {"label": "A1", "species": [{"element": "A", "occu": 1}]},
    ]
    context = _context()
    result = canonicalize_tb2j_result(payload, magnetic_context=context)
    assert result["magnetic_observables"]["status"] == "computed"
    assert [site["label"] for site in result["structure"]["sites"]] == ["A1", "B1"]


def test_real_bcc_fixture_context_is_complete_without_vasp(tmp_path):
    from pymatgen.core import Structure
    from htvasp.tb2j.context import build_magnetic_context
    from htvasp.tb2j.jobs import parse_exchange_out

    root = Path(__file__).resolve().parents[1]
    structure = Structure.from_file(root / "data/poscars/BCC-Al-Nb.vasp")
    parsed = parse_exchange_out(
        root / "data/endmembers/BCC-Al-Nb/tb2jflow/3-tb2j_solve/exchange.out"
    )
    result = canonicalize_tb2j_result(
        {**parsed, "structure": structure.as_dict(), "elements": ["Al", "Nb"]},
        output_path=tmp_path / "BCC-Al-Nb-tb2j.json",
        magnetic_context=build_magnetic_context("BCC-Al-Nb", structure),
    )
    assert result["magnetic_contract"]["calphad"]["sublattice_ratios"] == [0.5, 0.5]
    assert result["magnetic_contract"]["calphad"]["structure_role"] == "ordered_bcc_endmember"
    assert result["magnetic_observables"]["status"] in {"computed", "nonmagnetic", "unsupported"}
