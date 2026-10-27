"""Regression tests for B2 reverse-endmember QHA targets and plans."""

from __future__ import annotations

import json
from pathlib import Path

import main_qha
from htvasp.slurm import SlurmConfig
from htvasp.workflows.qha import canonicalize_qha_result


REVERSE = [
    "BCC-Nb-Al",
    "BCC-Ti-Al",
    "BCC-V-Al",
    "BCC-Ti-Nb",
    "BCC-V-Nb",
    "BCC-V-Ti",
]


class NoSubmitManager:
    def __init__(self):
        self.submit_calls = 0

    def get_cpu_config(self, **kwargs):
        return SlurmConfig(**kwargs)

    def submit_intent(self, **kwargs):  # pragma: no cover - must not be called
        self.submit_calls += 1
        raise AssertionError("dry-run must not call submit_intent")


def _payload():
    return {
        "structure": {
            "lattice": {"matrix": [[2.8, 0, 0], [0, 2.8, 0], [0, 0, 2.8]]},
            "sites": [
                {"species": [{"element": "Nb", "occu": 1}], "abc": [0, 0, 0]},
                {"species": [{"element": "Al", "occu": 1}], "abc": [0.5, 0.5, 0.5]},
            ],
        },
        "nsites": 2,
        "composition_reduced": {"Al": 1.0, "Nb": 1.0},
        "temperatures": [0.0, 50.0],
        "gibbs_temperature": [-1.0, -0.9],
    }


def test_reverse_manifest_has_six_unique_targets_and_hashes():
    path = Path("data/poscars/b2_reverse_manifest.json")
    manifest = json.loads(path.read_text(encoding="utf-8"))
    records = manifest["targets"]
    assert [record["reverse_name"] for record in records] == REVERSE
    assert len({record["reverse_name"] for record in records}) == 6
    for record in records:
        assert record["canonical_phase"] == "B2_BCC"
        assert record["ordered"] is True
        assert len(record["source_sha256"]) == 64
        assert len(record["derived_sha256"]) == 64
        assert record["site_mapping"]["sublattice_1"]["fractional_coordinate"] == [0.0, 0.0, 0.0]
        assert record["site_mapping"]["sublattice_2"]["fractional_coordinate"] == [0.5, 0.5, 0.5]


def test_reverse_dry_run_is_single_target_and_does_not_submit():
    manager = NoSubmitManager()
    plan = main_qha.submit_job("BCC-Nb-Al", "qha", manager, dry_run=True)
    assert plan["name"] == "BCC-Nb-Al"
    assert plan["phase_context"]["canonical_phase"] == "B2_BCC"
    assert plan["phase_context"]["ordered"] is True
    assert plan["phase_context"]["reverse_of"] == "BCC-Al-Nb"
    assert plan["config_dict"]["partition"] == "partCPU"
    assert plan["config_dict"]["ntasks"] == 32
    assert plan["config_dict"]["memory"] == "64G"
    assert plan["config_dict"]["module_name"] == "vasp-cpu"
    assert manager.submit_calls == 0
    assert len(main_qha.resolve_targets(name="BCC-Nb-Al")) == 1
    assert main_qha.resolve_targets(targets=REVERSE) == REVERSE


def test_qha_result_carries_reverse_phase_context(tmp_path):
    context = main_qha._b2_phase_context("BCC-Nb-Al")
    output = canonicalize_qha_result(
        _payload(),
        structure_name="BCC-Nb-Al",
        output_path=tmp_path / "BCC-Nb-Al-qha.json",
        phase_context=context,
    )
    assert output["phase_context"] == context
    assert output["structure_context"]["phase_context"] == context
    assert output["provenance"]["phase_context"] == context
    assert output["reference_context"]["role"] == "DFT_QHA_ORDERED"
