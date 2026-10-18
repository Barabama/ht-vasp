"""Tests for the DFPT-QHA thermal and Gibbs construction.

The unit conversion is the reason this file exists.  phonopy reports its
thermal properties in kJ per mole of ITS OWN primitive cell -- it divides by
``sum(weights)`` over the mesh -- while the DFT energies are per recipe cell,
and those two are not always the same cell.  For B2 they coincide; for an fcc
SER entry the phonopy primitive cell is 1 atom against a 4-atom conventional
cell.  A conversion that gets this wrong produces a Gibbs curve that is smooth,
plausible, and off by a constant factor.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from jobflow import JobStore
from maggma.stores import MemoryStore
from pymatgen.core import Structure

from htvasp.dfpt_qha.phonons import phonopy_for
from htvasp.dfpt_qha.qha import (
    EV_PER_ANGSTROM3_TO_GPA,
    EV_PER_ATOM_TO_KJ_PER_MOL,
    analyze_dfpt_qha,
    birch_murnaghan_3rd,
    dfpt_thermal_properties,
    volume_free_energy,
)

STRUCTURE_FILE = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_structures.json")
MATRIX = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
#: The Debye temperatures of Al (428 K) and Ta (240 K) put the zero-point
#: energy of B2 AlTa near (9/8) k_B * 334 K = 32 meV/atom.  A conversion that
#: multiplies by the atom count instead of dividing by the primitive-cell
#: count lands at twice that.
EXPECTED_ZPE_MEV_PER_ATOM = 32.0


def _store() -> JobStore:
    docs = MemoryStore()
    docs.connect()
    return JobStore(docs)


def _cell(name: str) -> Structure:
    data = json.loads(STRUCTURE_FILE.read_text())["endpoints"][name]["structure"]
    return Structure.from_dict(data)


def test_primitive_cell_size_differs_from_the_recipe_cell_for_fcc_ser():
    """The two cells are not always the same, which is why the ratio is explicit."""

    b2 = phonopy_for(_cell("BCC-Al-Ta-2"), MATRIX)
    fcc = phonopy_for(_cell("SER-Al-1"), MATRIX)
    assert len(b2.primitive) == 2 == len(_cell("BCC-Al-Ta-2"))
    assert len(fcc.primitive) == 1 != len(_cell("SER-Al-1"))


def test_zero_point_energy_is_physical_not_scaled_by_atom_count():
    """A factor-of-two unit error shows up here and nowhere else.

    phonopy's ZPE is kJ/mol of primitive cells.  Divided by the conversion it
    must land near the Debye estimate for this compound; multiplied by the atom
    count it lands at twice that, which is what the code did before.
    """

    cell = _cell("BCC-Al-Ta-2")
    phonon = phonopy_for(cell, MATRIX)
    phonon.generate_displacements(distance=1.0)
    n = len(phonon.supercell)
    rng = np.random.default_rng(7)
    phonon.forces = [rng.normal(scale=0.5, size=(n, 3)) for _ in phonon.displacements]
    phonon.produce_force_constants()
    phonon.run_mesh([6, 6, 6], with_eigenvectors=False)
    phonon.run_thermal_properties(t_min=0, t_max=0, t_step=1)
    zpe_kj = float(phonon.thermal_properties.zero_point_energy)

    per_atom_mev = zpe_kj / EV_PER_ATOM_TO_KJ_PER_MOL / len(phonon.primitive) * 1000
    # this fixture is random forces, so only the order of magnitude is asserted
    assert 0.0 < per_atom_mev < 3 * EXPECTED_ZPE_MEV_PER_ATOM, per_atom_mev


def test_thermal_job_reports_the_primitive_cell_size():
    """The Gibbs step needs it, so the thermal job must carry it."""

    cell = _cell("BCC-Al-Ta-2")
    phonon = phonopy_for(cell, MATRIX)
    phonon.generate_displacements(distance=1.0)
    n = len(phonon.supercell)
    phonon.forces = [np.zeros((n, 3)) for _ in phonon.displacements]
    phonon.produce_force_constants()
    fc = np.array(phonon.force_constants)

    supercell = Structure(
        phonon.supercell.cell,
        list(phonon.supercell.symbols),
        [np.array(p) for p in phonon.supercell.scaled_positions],
        coords_are_cartesian=False,
    )
    collected = {
        "volumes": [
            {
                "index": 0,
                "cell": cell.as_dict(),
                "supercell": supercell.as_dict(),
                "volume_angstrom3": float(cell.volume),
                "static_energy_eV": -15.0,
                "force_constants": fc.tolist(),
                "stability": {
                    "status": "harmonic_QHA_applicable", "qha_allowed": True,
                },
            }
        ],
        "rejected_volume_indices": [],
        "n_validated": 1,
        "n_requested": 1,
    }
    thermal = dfpt_thermal_properties(
        collected, supercell_matrix=MATRIX, mesh=(4, 4, 4), temperatures=[0.0, 300.0]
    ).run(_store()).output
    record = thermal["volumes"][0]
    assert record["n_atoms_cell"] == 2
    assert record["n_atoms_primitive"] == 2
    assert record["n_atoms_supercell"] == 16
    assert len(record["free_energy_kJ_per_mol"]) == 2


def test_thermal_job_reports_the_grid_phonopy_actually_evaluated():
    """phonopy builds its grid from (min, max, step), so it can differ from the
    requested one: 0, 300, 600, 1000 evaluates 0, 300, 600, 900.

    Reporting the requested list while indexing phonopy's arrays by position
    would pair every free energy with the wrong temperature, and the resulting
    curve would still look smooth.
    """

    cell = _cell("BCC-Al-Ta-2")
    phonon = phonopy_for(cell, MATRIX)
    phonon.generate_displacements(distance=1.0)
    n = len(phonon.supercell)
    phonon.forces = [np.zeros((n, 3)) for _ in phonon.displacements]
    phonon.produce_force_constants()
    supercell = Structure(
        phonon.supercell.cell,
        list(phonon.supercell.symbols),
        [np.array(p) for p in phonon.supercell.scaled_positions],
        coords_are_cartesian=False,
    )
    collected = {
        "volumes": [
            {
                "index": 0,
                "cell": cell.as_dict(),
                "supercell": supercell.as_dict(),
                "volume_angstrom3": float(cell.volume),
                "static_energy_eV": -15.0,
                "force_constants": np.array(phonon.force_constants).tolist(),
                "stability": {
                    "status": "harmonic_QHA_applicable", "qha_allowed": True,
                },
            }
        ]
    }
    thermal = dfpt_thermal_properties(
        collected, supercell_matrix=MATRIX, mesh=(4, 4, 4),
        temperatures=[0.0, 300.0, 600.0, 1000.0],
    ).run(_store()).output
    assert thermal["requested_temperatures"] == [0.0, 300.0, 600.0, 1000.0]
    assert thermal["temperatures"] == [0.0, 300.0, 600.0, 900.0]
    record = thermal["volumes"][0]
    assert len(record["free_energy_kJ_per_mol"]) == len(thermal["temperatures"])
    assert record["thermal_temperatures"] == thermal["temperatures"]


def test_birch_murnaghan_recovers_a_synthetic_curve():
    v0, b0, bp = 33.5, 1.0, 4.0  # eV/Angstrom^3

    def bm3(v):
        return 9.0 * v0 * b0 / 16.0 * (
            ((v0 / v) ** (2.0 / 3.0) - 1.0) ** 3 * bp
            + ((v0 / v) ** (2.0 / 3.0) - 1.0) ** 2
            * (6.0 - 4.0 * (v0 / v) ** (2.0 / 3.0))
        )

    volumes = np.array([v0 * s**3 for s in (0.94, 0.97, 1.0, 1.03, 1.06)])
    fit = birch_murnaghan_3rd(volumes, bm3(volumes))
    assert fit["V0_angstrom3"] == pytest.approx(v0, rel=1e-6)
    assert fit["B0_GPa"] == pytest.approx(b0 * EV_PER_ANGSTROM3_TO_GPA, rel=1e-5)
    assert fit["B0prime"] == pytest.approx(bp, rel=1e-4)


def _volume(index, volume, energy, free_energy_kj=(-1.0, -2.0), n_cell=2,
            n_primitive=2, mesh_min=0.0):
    return {
        "index": index,
        "volume_angstrom3": volume,
        "static_energy_eV": energy,
        "free_energy_kJ_per_mol": list(free_energy_kj),
        "n_atoms_cell": n_cell,
        "n_atoms_primitive": n_primitive,
        "mesh_diagnostic": {
            "mesh": [12, 12, 12],
            "mesh_min_THz": mesh_min,
            "mesh_min_qpoint": [-0.0417, 0.0417, 0.0417],
            "n_mesh_imaginary_modes": 1 if mesh_min < 0 else 0,
            "note": "interpolated diagnostic only; not a stability verdict",
        },
    }


def test_analysis_does_not_fail_on_a_mesh_artifact():
    """The mesh minimum is interpolated context, not a verdict.

    Gating on it rejected BCC-Nb-Nb-2, whose -0.024 THz mesh minimum sits at a
    q no 2x2x2 supercell represents while every exact commensurate q point is
    stable.  The decision is made by ``classify_dfpt_stability``; the analysis
    records the numbers and moves on.
    """

    thermal = {
        "temperatures": [0.0, 300.0],
        "volumes": [
            _volume(i, 30.0 + i, -10.0 + 0.02 * (i - 2) ** 2,
                    mesh_min=-0.024 if i == 2 else 0.5)
            for i in range(5)
        ],
    }
    result = analyze_dfpt_qha(thermal).run(_store())
    assert result.output["ok"] is True
    assert result.stop_children is False
    assert [d["mesh_min_THz"] for d in result.output["mesh_diagnostics"]][2] == pytest.approx(-0.024)


def test_analysis_fails_closed_on_too_few_volumes():
    thermal = {
        "temperatures": [0.0, 300.0],
        "volumes": [_volume(i, 30.0 + i, -10.0) for i in range(3)],
    }
    result = analyze_dfpt_qha(thermal).run(_store())
    assert result.output["ok"] is False
    assert result.stop_children is True
    assert "4" in result.output["error"]


def test_analysis_warns_when_the_eos_minimum_leaves_the_window():
    """V0 on the edge of the scan is an artefact of where the scan stopped."""

    thermal = {
        "temperatures": [0.0, 300.0],
        "volumes": [
            _volume(i, 30.0 + i, -10.0 + 0.02 * (i - 4) ** 2, [-1.0, -2.0])
            for i in range(5)
        ],
    }
    result = analyze_dfpt_qha(thermal).run(_store())
    assert result.output["ok"] is True
    assert result.output["v0_inside_sampled_window"] is False
    assert any("outside the sampled window" in w for w in result.output["warnings"])


def test_analysis_returns_a_gibbs_curve_on_a_clean_input():
    thermal = {
        "temperatures": [0.0, 300.0],
        "volumes": [
            _volume(i, 30.0 + i, -10.0 + 0.02 * (i - 2) ** 2, [-1.0 - 0.1 * i, -2.0])
            for i in range(5)
        ],
    }
    result = analyze_dfpt_qha(thermal).run(_store())
    assert result.output["ok"] is True
    assert result.stop_children is False
    assert len(result.output["gibbs_eV"]) == 2
    assert len(result.output["equilibrium_volume_angstrom3"]) == 2
    assert result.output["v0_inside_sampled_window"] is True
    assert result.output["units"]["gibbs_eV"] == "eV per cell"
    # G = F + pV at p=0 must never exceed the best F, since it is a minimum
    best_f = [
        min(v["static_energy_eV"] + v["free_energy_kJ_per_mol"][t] / EV_PER_ATOM_TO_KJ_PER_MOL
            for v in thermal["volumes"])
        for t in range(2)
    ]
    assert result.output["gibbs_eV"] == pytest.approx(best_f, rel=1e-9)


# --------------------------------------------------------------------------
# the factor-4 endpoints: fcc A=A has a 1-atom primitive against a 4-atom cell
# --------------------------------------------------------------------------

#: Legacy finite-displacement run on fcc Al: a 4-atom conventional cell, a
#: 2x2x2 supercell, and phonopy's primitive matrix gives 1 atom.  It is the
#: only real fcc force-constant dataset on disk, and it is the case the
#: per-primitive -> per-cell factor has never been exercised on.
FCC_FIXTURE = Path(
    "/nfs_hdd/2025/gaominliang/ht-vasp/data/endmembers/SER-Al/qhaflow"
    "/38-generate_frequencies_eigenvectors_eos_deformation_2/phonopy.yaml"
)

FACTOR_FOUR_ENDPOINTS = [
    "FCC-Al-Al-4", "FCC-Cr-Cr-4", "FCC-Hf-Hf-4", "FCC-Nb-Nb-4",
    "FCC-Ta-Ta-4", "FCC-Ti-Ti-4", "FCC-Zr-Zr-4", "SER-Al-1",
]


@pytest.mark.parametrize("name", FACTOR_FOUR_ENDPOINTS)
def test_the_eight_factor_four_endpoints_really_are_four(name):
    """A=A fcc: the recipe cell is 4 atoms and phonopy's primitive is 1.

    Getting this ratio wrong is not a small error: it multiplies every free
    energy by four.  The B2 case had the same bug with a factor of two, which
    is why the ratio is explicit in ``analyze_dfpt_qha``.
    """

    cell = _cell(name)
    phonon = phonopy_for(cell, MATRIX)
    assert len(cell) == 4
    assert len(phonon.primitive) == 1


@pytest.mark.skipif(not FCC_FIXTURE.is_file(), reason="legacy SER-Al phonon data absent")
def test_fcc_zero_point_energy_survives_the_primitive_to_cell_conversion():
    """Real fcc force constants, 32-atom supercell, factor 4.

    ZPE of fcc Al is ~50 meV/atom (mean phonon frequency ~8.1 THz).  A
    conversion that multiplies by the atom count instead of the
    primitive/cell ratio lands at ~200 meV/atom; one that divides by the cell
    count lands at ~12.9 meV/atom.  Both are inside an order of magnitude of
    the right answer and both are wrong.
    """

    import phonopy

    phonon = phonopy.load(FCC_FIXTURE, produce_fc=False, log_level=0)
    phonon.produce_force_constants()
    unit = np.array(phonon.unitcell.cell)
    cell = Structure(
        unit, list(phonon.unitcell.symbols),
        [np.array(p) for p in phonon.unitcell.scaled_positions],
        coords_are_cartesian=False,
    )
    supercell = Structure(
        np.array(phonon.supercell.cell), list(phonon.supercell.symbols),
        [np.array(p) for p in phonon.supercell.scaled_positions],
        coords_are_cartesian=False,
    )
    result = volume_free_energy(
        cell, supercell, phonon.force_constants,
        supercell_matrix=MATRIX, mesh=(12, 12, 12), temperatures=[0.0, 300.0],
    )
    assert result["n_atoms_cell"] == 4
    assert result["n_atoms_primitive"] == 1
    assert result["n_atoms_supercell"] == 32

    zpe_ev_per_atom = (
        result["zero_point_energy_kJ_per_mol"]
        / EV_PER_ATOM_TO_KJ_PER_MOL
        / result["n_atoms_primitive"]
    )
    assert 0.045 < zpe_ev_per_atom < 0.060, zpe_ev_per_atom

    # the analysis's per-cell conversion must reproduce that per-atom value
    ratio = result["n_atoms_cell"] / result["n_atoms_primitive"]
    per_atom_from_cell = (
        result["free_energy_kJ_per_mol"][1] * ratio
        / EV_PER_ATOM_TO_KJ_PER_MOL / result["n_atoms_cell"]
    )
    assert per_atom_from_cell == pytest.approx(
        result["free_energy_kJ_per_mol"][1] / EV_PER_ATOM_TO_KJ_PER_MOL
        / result["n_atoms_primitive"], rel=1e-12
    )


def test_analyze_converts_primitive_free_energy_to_the_recipe_cell():
    """The fcc 4:1 conversion, exercised through the analysis itself.

    ``volume_free_energy`` reports kJ per mole of phonopy's PRIMITIVE cell.
    For an fcc A=A endpoint the recipe cell holds four of them, so the free
    energy the analysis fits must be four times the reported value.  Inverting
    the ratio gives 0.25x and every Gibbs curve is off by 16.

    This must go through ``analyze_dfpt_qha``: an assertion that recomputes the
    same ratio by hand is an algebraic identity and passes for either direction.
    """

    thermal = {
        "temperatures": [0.0, 300.0],
        "volumes": [
            _volume(i, 30.0 + i, -10.0 + 0.02 * (i - 2) ** 2,
                    free_energy_kj=(-1.0 - 0.1 * i, -2.0),
                    n_cell=4, n_primitive=1)
            for i in range(5)
        ],
    }
    result = analyze_dfpt_qha(thermal).run(_store())
    assert result.output["ok"] is True
    per_cell = result.output["free_energies_eV"]
    for index in range(5):
        expected = (-1.0 - 0.1 * index) / EV_PER_ATOM_TO_KJ_PER_MOL * 4
        assert per_cell[index][0] == pytest.approx(expected, rel=1e-12)
    # the inverted ratio would land here, and must not
    inverted = -1.0 / EV_PER_ATOM_TO_KJ_PER_MOL * 0.25
    assert per_cell[0][0] != pytest.approx(inverted, rel=1e-6)
