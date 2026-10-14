"""Focused construction tests for the atomate2-style DFPT-QHA Maker."""

from __future__ import annotations

from jobflow import Flow
from pymatgen.core import Lattice, Structure

from htvasp.dfpt_qha.input_set import DfptSetGenerator
from htvasp.dfpt_qha.maker import DfptPhononMaker, DfptQhaMaker, DfptVaspMaker


STRUCTURE = Structure(
    Lattice.cubic(3.2), ["Al", "Ta"], [[0, 0, 0], [0.5, 0.5, 0.5]]
)


def test_dfpt_generator_uses_ibrion_8():
    generator = DfptSetGenerator()
    assert generator.incar_updates["IBRION"] == 8
    assert generator.incar_updates["NSW"] == 1
    assert generator.incar_updates["KPAR"] == 2
    assert "NCORE" not in generator.incar_updates


def test_dfpt_vasp_maker_is_base_vasp_maker():
    maker = DfptVaspMaker()
    assert isinstance(maker, DfptVaspMaker)
    assert maker.input_set_generator.incar_updates["IBRION"] == 8


def test_dfpt_phonon_maker_has_static_and_dfpt_children():
    maker = DfptPhononMaker()
    assert isinstance(maker.static_energy_maker, DfptVaspMaker)
    assert isinstance(maker.phonon_displacement_maker, DfptVaspMaker)
    assert maker.phonon_displacement_maker.input_set_generator.incar_updates["IBRION"] == 8


def test_dfpt_qha_flow_uses_official_three_job_lifecycle():
    flow = DfptQhaMaker().make(
        STRUCTURE,
        supercell_matrix=[[2, 0, 0], [0, 2, 0], [0, 0, 2]],
    )
    assert isinstance(flow, Flow)
    assert [job.name for job in flow.jobs] == [
        "EOS Maker",
        "get_phonon_jobs",
        "analyze_free_energy",
    ]
