"""Offline regression tests for QHA/TB2J production JSON metadata."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from htvasp.tb2j.jobs import canonicalize_tb2j_result, parse_exchange_out
from htvasp.tb2j.task_doc import Tb2jResult
from htvasp.workflows.qha import QhaWorker, canonicalize_qha_result


ROOT = Path(__file__).resolve().parents[1]
EXCHANGE_FIXTURE = ROOT / "data" / "endmembers" / "BCC-Al-Nb" / "tb2jflow" / "3-tb2j_solve" / "exchange.out"


def _qha_payload() -> dict:
    return {
        "structure": {
            "@module": "pymatgen.core.structure",
            "@class": "Structure",
            "lattice": {"matrix": [[2.8, 0, 0], [0, 2.8, 0], [0, 0, 2.8]]},
            "sites": [
                {"species": [{"element": "Al", "occu": 1}], "abc": [0, 0, 0]},
                {"species": [{"element": "Nb", "occu": 1}], "abc": [0.5, 0.5, 0.5]},
            ],
        },
        "nsites": 2,
        "composition_reduced": {"Al": 1.0, "Nb": 1.0},
        "temperatures": [0.0, 50.0, 100.0],
        "gibbs_temperature": [-1.0, -0.9],
        "bulk_modulus": 1.0,
        "volumes": [20.0],
    }


def test_qha_worker_write_result_emits_canonical_contract(tmp_path):
    output = tmp_path / "BCC-Al-Nb-qha.json"
    worker = QhaWorker(
        vasp_args={"vasp_cmd": "echo test"},
        reference_name="BCC-Al-Nb",
        temperature_range=(0, 100, 50),
    )

    worker.write_result(_qha_payload(), output)
    document = json.loads(output.read_text(encoding="utf-8"))

    assert document["schema_version"] == "qha-v1"
    assert document["units"] == {
        "temperature": "K",
        "gibbs_temperature": "eV/reduced_formula_unit",
    }
    assert document["reference_context"] == {
        "role": "DFT_QHA_ORDERED",
        "role_status": "explicit",
        "role_source": "htvasp_workflow_structure_name",
    }
    assert len(document["temperatures"]) == len(document["gibbs_temperature"]) == 2
    assert document["provenance"]["input_temperature_count"] == 3
    assert document["provenance"]["temperature_alignment_operation"] == "trailing_temperature_dropped"
    assert document["provenance"]["source_sha256"]


def test_qha_worker_supports_fcc_roles_without_calculation(tmp_path):
    output = tmp_path / "FCC-Co-Co-qha.json"
    worker = QhaWorker(
        vasp_args={"vasp_cmd": "echo test"},
        reference_name="FCC-Co-Co",
        temperature_range=(0, 50, 50),
    )

    worker.write_result(
        {
            **_qha_payload(),
            "structure": {
                **_qha_payload()["structure"],
                "sites": [
                    {"species": [{"element": "Co", "occu": 1}], "abc": [0, 0, 0]},
                    {"species": [{"element": "Co", "occu": 1}], "abc": [0.5, 0.5, 0.5]},
                ],
            },
            "composition_reduced": {"Co": 1.0},
            "temperatures": [0.0, 50.0],
            "gibbs_temperature": [-1.0, -0.9],
        },
        output,
    )
    document = json.loads(output.read_text(encoding="utf-8"))
    assert document["reference_context"]["role"] == "DFT_QHA_FCC_PURE"


def test_tb2j_parser_emits_contract_and_worker_result_is_canonical(tmp_path):
    parsed = parse_exchange_out(EXCHANGE_FIXTURE)
    assert parsed["schema_version"] == "emtdb-magnetic-dataset/v1"
    assert parsed["units"]["J_iso"] == "meV"
    assert parsed["units"]["w_magmom"] == "mu_B"
    assert parsed["exchange_convention"]["pair_counting"] == "outgoing_per_site"
    assert parsed["provenance"]["reverse_pair_validation"]["reverse_consistency_checked"] is True

    persistent = tmp_path / "BCC-Al-Nb" / "tb2jflow" / "3-tb2j_solve" / "exchange.out"
    persistent.parent.mkdir(parents=True)
    shutil.copy2(EXCHANGE_FIXTURE, persistent)
    output = tmp_path / "BCC-Al-Nb" / "BCC-Al-Nb-tb2j.json"
    payload = {
        **parsed,
        "structure": {
            "lattice": {"matrix": parsed["cell"]},
            "sites": [{"label": moment["atom"]} for moment in parsed["atom_moments"]],
        },
        "elements": ["Al", "Nb"],
        "efermi": 1.0,
        "kmesh": [9, 9, 9],
        "num_wann": 18,
    }
    canonical = canonicalize_tb2j_result(payload, output_path=output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(canonical), encoding="utf-8")
    result = Tb2jResult.model_validate(json.loads(output.read_text(encoding="utf-8")))

    assert result.schema_version == "emtdb-magnetic-dataset/v1"
    assert result.source["sha256"]
    assert result.units["J_iso"] == "meV"
    assert result.exchange_convention["pair_counting"] == "outgoing_per_site"
    # Parser-level metadata is not a physical CALPHAD approval.  A producer
    # context is required before TC/BMAGN can be marked production-ready.
    assert result.audit["production_ready"] is False
    assert result.magnetic_observables["status"] == "needs_input"


def test_tb2j_worker_rejects_site_moment_element_mismatch(tmp_path):
    parsed = parse_exchange_out(EXCHANGE_FIXTURE)
    output = tmp_path / "BCC-Al-Nb" / "BCC-Al-Nb-tb2j.json"
    with pytest.raises(ValueError, match="does not match"):
        canonicalize_tb2j_result(
            {
                **parsed,
                "structure": {
                    "lattice": {"matrix": parsed["cell"]},
                    "sites": [
                        {"species": [{"element": "Al", "occu": 1}]},
                        {"species": [{"element": "Al", "occu": 1}]},
                    ],
                },
                "elements": ["Al", "Nb"],
                "efermi": 1.0,
                "kmesh": [9, 9, 9],
                "num_wann": 18,
            },
            output_path=output,
        )
