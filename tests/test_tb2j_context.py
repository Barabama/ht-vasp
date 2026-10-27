"""Offline tests for producer-side TB2J magnetic context handling."""

from __future__ import annotations

import json

from jobflow import job
from pymatgen.core import Lattice, Structure

from htvasp.tb2j.context import (
    build_magnetic_context,
    deserialize_magnetic_context,
    serialize_magnetic_context,
)
from htvasp.tb2j.maker import Tb2jMaker
from htvasp.workflows.tb2j import Tb2jWorker


def _bcc_reversed() -> Structure:
    return Structure(
        Lattice.cubic(3.0),
        ["Nb", "Al"],
        [[0.5, 0.5, 0.5], [0.0, 0.0, 0.0]],
    )


def test_b2_context_maps_fractional_sites_and_actual_occupation_direction():
    context = build_magnetic_context("BCC-Al-Nb", _bcc_reversed())

    assert context["source"] == "htvasp_workflow_context"
    assert context["context_version"]
    assert context["producer_recipe_version"]
    assert len(context["structure_checksum"]) == 64
    assert context["calphad"]["phase"] == "B2_BCC"
    assert context["calphad"]["type_character"] == "&"
    assert context["calphad"]["ratios"] == [0.5, 0.5]
    assert context["calphad"]["site_to_sublattice"] == {"Al1": 0, "Nb1": 1}
    assert context["calphad"]["endmember_components"] == ["AL", "NB"]
    assert context["calphad"]["order_num"] == 0
    assert context["mapping_rationale"]["corner_site_id"] == "Al1"
    assert context["mapping_rationale"]["body_center_site_id"] == "Nb1"
    assert [site["id"] for site in context["sites"]] == ["Nb1", "Al1"]
    assert "initial VASP MAGMOM is not an atom moment" in context["magnetic_moment_source"]


def test_ser_context_is_not_applicable_and_does_not_fabricate_b2():
    structure = Structure(Lattice.cubic(3.0), ["Al"], [[0, 0, 0]])
    context = build_magnetic_context("SER-Al", structure)

    assert context["structure_role"] == "unary_reference"
    assert context["calphad"]["status"] == "not_applicable"
    assert context["calphad"]["phase"] is None
    assert context["calphad"]["site_to_sublattice"] == {}


def test_ser_context_marks_magnetic_observable_not_applicable():
    from htvasp.tb2j.jobs import canonicalize_tb2j_result

    context = build_magnetic_context(
        "SER-Al", Structure(Lattice.cubic(3.0), ["Al"], [[0, 0, 0]])
    )
    payload = {
        "structure": {
            "lattice": {"matrix": [[3, 0, 0], [0, 3, 0], [0, 0, 3]]},
            "sites": [{"label": "Al1", "species": [{"element": "Al", "occu": 1}]}],
        },
        "num_atoms": 1,
        "atom_moments": [{"atom": "Al1", "w_magmom": 0.0}],
        "j_pairs": [{"i": "Al1", "j": "Al1", "R": [1, 0, 0], "J_iso_meV": 0.0, "vector": [1, 0, 0], "distance": 3.0}],
        "j_tensors": [[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]],
    }
    result = canonicalize_tb2j_result(payload, magnetic_context=context)
    assert result["magnetic_observables"]["status"] == "not_applicable"
    assert result["audit"]["production_ready"] is True


def test_fcc_context_is_reviewed_until_1a_3c_is_explicit():
    structure = Structure(
        Lattice.cubic(3.0), ["Co", "Ni"], [[0, 0, 0], [0.5, 0.5, 0.5]]
    )
    context = build_magnetic_context("FCC-Co-Ni", structure)

    assert context["structure_role"] == "ordered_fcc_endmember"
    assert context["calphad"]["phase"] == "L12_FCC"
    assert context["calphad"]["status"] == "needs_review"
    assert context["calphad"]["production_ready"] is False
    assert context["calphad"]["site_to_sublattice"] == {}


def test_unknown_role_requires_input():
    structure = Structure(Lattice.cubic(3.0), ["Al"], [[0, 0, 0]])
    context = build_magnetic_context("UNKNOWN-Al", structure)

    assert context["structure_role"] == "unknown"
    assert context["calphad"]["status"] == "needs_input"
    assert context["calphad"]["phase"] is None


def test_context_round_trip_is_json_safe():
    original = build_magnetic_context("BCC-Al-Nb", _bcc_reversed())
    restored = deserialize_magnetic_context(serialize_magnetic_context(original))

    assert restored == original


def test_worker_and_maker_keep_context_without_running_vasp():
    context = build_magnetic_context("BCC-Al-Nb", _bcc_reversed())
    worker = Tb2jWorker(
        "context-worker",
        vasp_args={"vasp_cmd": "unused"},
        magnetic_context=context,
    )
    assert worker.magnetic_context == context
    assert worker.tb2j_maker.magnetic_context == context

    maker = Tb2jMaker(magnetic_context=context)
    assert maker.magnetic_context == context


def test_maker_passes_context_to_new_tb2j_job_without_execution(monkeypatch):
    context = build_magnetic_context("BCC-Al-Nb", _bcc_reversed())

    @job
    def fake_tb2j_solve(
        prev_dir, elements, kmesh=(9, 9, 9), num_wann=None, magnetic_context=None
    ):
        return {"magnetic_context": magnetic_context}

    monkeypatch.setattr("htvasp.tb2j.maker.tb2j_solve", fake_tb2j_solve)
    flow = Tb2jMaker(magnetic_context=context).make(_bcc_reversed())

    assert flow.jobs[-1].function_kwargs["magnetic_context"] == context
