"""Spin screening: reading TaskDocs and deciding a branch.

The screening's only job is to compare energies and moments from three
relaxations.  Both halves have a way to fail silently -- a TaskDoc that comes
back as a pydantic model has no ``.get``, and emmet keeps the site moments
outside ``TaskDoc.output`` -- so both are pinned here.
"""

from __future__ import annotations

import pytest
from jobflow import JobStore
from maggma.stores import MemoryStore
from monty.json import MSONable

import main_spin_screen as screen


class _FakeTaskDoc(MSONable):
    """The shape the screening reads, as monty would rebuild it."""

    def __init__(self, energy, moments, abc=None):
        self.energy = energy
        self.moments = moments
        self.abc = abc or [4.0, 4.0, 4.0]

    def as_dict(self):
        return {
            "@module": __name__, "@class": "_FakeTaskDoc",
            "output": {
                "output": {
                    "energy": self.energy,
                    "structure": {
                        "sites": [{} for _ in self.moments],
                        "lattice": {"abc": self.abc},
                    },
                },
            },
            "calcs_reversed": [{
                "output": {"outcar": {
                    "magnetization": [{"tot": value} for value in self.moments]
                }},
            }],
        }


def _store() -> JobStore:
    docs = MemoryStore()
    docs.connect()
    return JobStore(docs)


def test_fields_are_read_from_a_plain_dict_and_from_a_model():
    document = _FakeTaskDoc(-10.0, [3.0, -3.0]).as_dict()
    assert screen._taskdoc_fields(document) == (-10.0, 2, [3.0, -3.0], [4.0, 4.0, 4.0])
    # across a job boundary monty rebuilds the object; ``.get`` would raise
    assert screen._taskdoc_fields(_FakeTaskDoc(-10.0, [3.0, -3.0])) == (
        -10.0, 2, [3.0, -3.0], [4.0, 4.0, 4.0]
    )


def _summary(states):
    documents = {
        label: _FakeTaskDoc(energy, [moments, 0.0]).as_dict()
        for label, (energy, moments) in states.items()
    }
    return screen.summarize_states(
        "X", documents, threshold_mev_per_atom=5.0
    ).run(_store()).output


def test_a_magnetic_branch_that_wins_by_more_than_the_threshold_is_chosen():
    # two atoms, so a 283.7 meV/atom gap is 0.5674 eV in total
    output = _summary({
        "nm": (-20.0, 0.0),
        "fm": (-20.5674, 2.66),
    })
    assert output["ok"] is True
    assert output["best_state"] == "fm"
    assert output["decision"] == "confirmed_branch"
    assert output["magnetic_gain_meV_per_atom"] == pytest.approx(283.7, abs=0.1)


def test_a_margin_below_the_threshold_is_nonmagnetic():
    output = _summary({
        "nm": (-20.0, 0.0),
        "fm": (-20.002, 0.4),
    })
    assert output["decision"] == "nonmagnetic"
    assert output["best_state"] == "fm"


def test_the_lowest_state_wins_between_fm_and_afm():
    output = _summary({
        "nm": (-10.0, 0.0),
        "fm": (-10.10, 2.0),
        "afm": (-10.20, 0.1),
    })
    assert output["best_state"] == "afm"
    assert output["decision"] == "confirmed_branch"


def test_a_structurally_collapsed_state_is_excluded():
    """A staggered start can collapse the cell, not magnetise it.

    In L1_2 the AFM start drives cubic -> tetragonal by ~20%, the moments
    vanish, and the energy drops by hundreds of meV/atom.  That is a different
    phase; letting it win would write a magnetic policy for a structure the
    endpoint does not have.
    """

    documents = {
        "nm": _FakeTaskDoc(-40.0, [0.0, 0.0, 0.0, 0.0]).as_dict(),
        "fm": _FakeTaskDoc(-40.1, [1.0, 1.0, 1.0, 1.0]).as_dict(),
        "afm": _FakeTaskDoc(
            -42.0, [0.0, 0.0, 0.0, 0.0], abc=[4.8, 3.4, 4.8]
        ).as_dict(),
    }
    output = screen.summarize_states(
        "X", documents, threshold_mev_per_atom=5.0
    ).run(_store()).output
    assert output["excluded_states"] == ["afm"]
    assert output["best_state"] == "fm"
    assert output["decision"] == "confirmed_branch"
    assert set(output["delta_vs_nm_meV_per_atom"]) == {"fm"}
