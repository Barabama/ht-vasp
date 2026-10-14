"""Tests for the 84-endpoint pool inventory.

The pool is where a wrong lattice constant would be invisible: every endpoint
would still have a cell, still pass the gate, and still produce a phonon.  So
the checks here are about the two decisions the cells rest on -- Vegard rather
than the `endmember.py` template, and conventional rather than primitive -- and
about the four pairs that are one crystal with two CALPHAD roles.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from htvasp.dfpt_qha.pool import (
    GROUND_STATE,
    POOL_ELEMENTS,
    POLYMORPH_STATE,
    build_pool,
    cubic_a,
    enumerate_endpoints,
    physical_state_id,
    volume_per_atom,
)

CACHE = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0c/endmember_cache")

# The four A=A endpoints whose Vegard volume IS the element's own volume, so
# they coincide with an SER entry (Al fcc; Cr, Nb, Ta bcc).
DUPLICATE_PAIRS = {
    "fcc-Al": ("FCC-Al-Al-4", "SER-Al-1"),
    "bcc-Cr": ("BCC-Cr-Cr-2", "SER-Cr-1"),
    "bcc-Nb": ("BCC-Nb-Nb-2", "SER-Nb-1"),
    "bcc-Ta": ("BCC-Ta-Ta-2", "SER-Ta-1"),
}


def test_pool_is_28_49_7():
    specs = enumerate_endpoints()
    counts: dict[str, int] = {}
    for spec in specs:
        counts[spec["class"]] = counts.get(spec["class"], 0) + 1
    assert counts == {"BCC_B2": 28, "FCC_L12": 49, "SER": 7}, counts
    assert len(specs) == 84


def test_every_ground_state_element_is_in_the_pool():
    """A pool element with no ground state would KeyError deep inside a build."""

    assert set(POOL_ELEMENTS) == set(GROUND_STATE), set(POOL_ELEMENTS) ^ set(GROUND_STATE)


def test_ground_state_volumes_are_physical():
    for element, state in GROUND_STATE.items():
        assert 8.0 < state.volume_per_atom < 30.0, (element, state)
        assert state.phase in ("BCC", "FCC", "HCP"), (element, state)


@pytest.mark.parametrize(
    ("element", "phase"),
    [("Al", "FCC"), ("Ta", "BCC"), ("Zr", "HCP"), ("Al", "BCC"), ("Hf", "BCC")],
)
def test_volume_per_atom_provenance(element, phase):
    vpa, source, is_proxy = volume_per_atom(element, phase)
    assert vpa > 0
    known = POLYMORPH_STATE.get((element, phase))
    if known is not None:
        assert not is_proxy and str(known.mp_id) in source, (element, phase, source)
    elif GROUND_STATE[element].phase == phase:
        assert not is_proxy and GROUND_STATE[element].mp_id in source
    else:
        # a stand-in must SAY it is a stand-in
        assert is_proxy and source.startswith("PROXY:"), (element, phase, source)


def test_cubic_a_reproduces_the_volume():
    for vpa in (10.0, 16.4717, 23.4997):
        for n_atoms in (2, 4):
            a = cubic_a(vpa, n_atoms)
            assert math.isclose(a**3 / n_atoms, vpa, rel_tol=1e-12)


def test_physical_state_id_marks_only_the_four_real_duplicates():
    specs = {spec["name"]: spec for spec in enumerate_endpoints()}
    duplicates: dict[str, list[str]] = {}
    for name, spec in specs.items():
        state = physical_state_id(spec)
        if state is not None:
            duplicates.setdefault(state, []).append(name)
    collapsed = {state: names for state, names in duplicates.items() if len(names) > 1}
    assert set(collapsed) == set(DUPLICATE_PAIRS), collapsed
    for state, (binary, ser) in DUPLICATE_PAIRS.items():
        assert set(collapsed[state]) == {binary, ser}, state


def test_a_mismatched_a_a_endpoint_is_not_a_duplicate():
    """Hf is hcp, so an hcp-phase A=A endpoint would duplicate SER-Hf -- but the
    pool has no hcp binary endpoint, and a bcc one must NOT be called a
    duplicate."""

    assert physical_state_id(
        {"class": "BCC_B2", "elements": ["Hf", "Hf"], "phase": "BCC"}
    ) is None
    assert physical_state_id(
        {"class": "BCC_B2", "elements": ["Ta", "Ta"], "phase": "BCC"}
    ) == "bcc-Ta"
    # different elements are never the same crystal
    assert physical_state_id(
        {"class": "BCC_B2", "elements": ["Al", "Ta"], "phase": "BCC"}
    ) is None


@pytest.mark.skipif(not CACHE.is_dir(), reason="SER structure cache not present")
def test_built_pool_matches_the_manifest_summary():
    manifest = build_pool(CACHE)
    assert manifest.n_endpoints == 84
    assert manifest.by_class == {"BCC_B2": 28, "FCC_L12": 49, "SER": 7}
    assert manifest.decisions.template_lattice_used is False
    assert set(manifest.duplicate_states) == set(DUPLICATE_PAIRS)

    from pymatgen.core import Structure

    for name, endpoint in manifest.endpoints.items():
        structure = Structure.from_dict(endpoint.structure)
        assert len(structure) == endpoint.n_atoms, name
        assert math.isclose(
            structure.volume / len(structure), endpoint.volume_per_atom, rel_tol=1e-9
        ), name
        # a Vegard cell must be cubic; an SER cell is whatever its element is
        if endpoint.structure_class != "SER":
            assert math.isclose(
                structure.lattice.a, endpoint.a, rel_tol=1e-12
            ), name


@pytest.mark.skipif(not CACHE.is_dir(), reason="SER structure cache not present")
def test_duplicate_pairs_are_the_same_cell():
    """The whole reason they collapse to one calculation."""

    manifest = build_pool(CACHE)
    for state, names in manifest.duplicate_states.items():
        first, second = (manifest.endpoints[name] for name in names)
        assert math.isclose(first.a, second.a, rel_tol=2e-5), (
            state, first.a, second.a
        )
        assert first.n_atoms == second.n_atoms, state
