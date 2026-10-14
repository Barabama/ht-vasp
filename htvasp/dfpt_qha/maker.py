"""Maker for the DFPT phonon stage and the official QHA lifecycle."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from atomate2.common.flows.phonons import BasePhononMaker
from atomate2.common.flows.qha import CommonQhaMaker
from atomate2.vasp.jobs.base import BaseVaspMaker
from jobflow import Flow

from htvasp.dfpt_qha.input_set import (
    DfptRelaxSetGenerator,
    DfptSetGenerator,
    DfptStaticSetGenerator,
)
from htvasp.dfpt_qha.jobs import collect_dfpt_phonons, validate_dfpt_output
from htvasp.model.dfpt_qha_recipe import Recipe


@dataclass
class DfptVaspMaker(BaseVaspMaker):
    """One atomate2 VASP job using the recipe-backed input generator."""

    recipe: Recipe | None = None
    stage: str = "dfpt"
    input_set_generator: Any = None

    def __post_init__(self) -> None:
        if self.input_set_generator is None:
            if self.stage == "dfpt":
                self.input_set_generator = DfptSetGenerator(recipe=self.recipe)
            elif self.stage == "relax":
                self.input_set_generator = DfptRelaxSetGenerator(recipe=self.recipe)
            else:
                self.input_set_generator = DfptStaticSetGenerator(recipe=self.recipe)
        # BaseVaspMaker is a dataclass with no post-init hook.


@dataclass
class DfptPhononMaker(BasePhononMaker):
    """PhononMaker-shaped flow whose force jobs use DFPT (IBRION=8)."""

    recipe: Recipe | None = None
    name: str = "DFPT phonon"
    sym_reduce: bool = False
    use_symmetrized_structure: str | None = "conventional"
    bulk_relax_maker: BaseVaspMaker | None = None
    static_energy_maker: BaseVaspMaker | None = None
    born_maker: BaseVaspMaker | None = None
    phonon_displacement_maker: BaseVaspMaker | None = None
    store_force_constants: bool = True

    def __post_init__(self) -> None:
        # CommonPhononMaker does not define a post-init hook in this atomate2
        # release; its dataclass fields are validated by the explicit checks
        # below and by the child BaseVaspMaker instances.
        if self.phonon_displacement_maker is None:
            self.phonon_displacement_maker = DfptVaspMaker(
                name="dfpt force",
                recipe=self.recipe,
                stage="dfpt",
                input_set_generator=DfptSetGenerator(recipe=self.recipe),
                stop_children_kwargs={"handle_unsuccessful": False},
            )
        if self.static_energy_maker is None:
            self.static_energy_maker = DfptVaspMaker(
                name="dfpt static energy",
                recipe=self.recipe,
                stage="static",
                input_set_generator=DfptStaticSetGenerator(recipe=self.recipe),
                stop_children_kwargs={"handle_unsuccessful": False},
            )

    @property
    def prev_calc_dir_argname(self) -> str:
        return "prev_dir"


@dataclass
class DfptQhaMaker(CommonQhaMaker):
    """Production DFPT-QHA Maker composed from atomate2's CommonQhaMaker.

    Only the phonon implementation is specialized.  EOS deformation, per-
    volume flow replacement, and QHA analysis remain the official atomate2
    lifecycle inherited through ``CommonQhaMaker``/``QhaMaker``.
    """

    name: str = "DFPT-QHA Maker"
    recipe: Recipe | None = None
    initial_relax_maker: BaseVaspMaker | None = None
    eos_relax_maker: BaseVaspMaker | None = None
    phonon_maker: DfptPhononMaker | None = None
    linear_strain: tuple[float, float] = (-0.05, 0.05)
    number_of_frames: int = 6
    ignore_imaginary_modes: bool = False
    skip_analysis: bool = False
    min_length: float | None = None

    def __post_init__(self) -> None:
        if self.initial_relax_maker is None:
            self.initial_relax_maker = DfptVaspMaker(
                name="DFPT-QHA initial relax",
                recipe=self.recipe,
                stage="relax",
                input_set_generator=DfptRelaxSetGenerator(recipe=self.recipe),
                stop_children_kwargs={"handle_unsuccessful": False},
            )
        if self.eos_relax_maker is None:
            self.eos_relax_maker = DfptVaspMaker(
                name="DFPT-QHA EOS relax",
                recipe=self.recipe,
                stage="static",
                input_set_generator=DfptStaticSetGenerator(recipe=self.recipe),
                stop_children_kwargs={"handle_unsuccessful": False},
            )
        if self.phonon_maker is None:
            self.phonon_maker = DfptPhononMaker(recipe=self.recipe)
        # CommonQhaMaker has no post-init implementation, so validate our
        # fixed-volume constraint locally before the official lifecycle runs.
        if self.phonon_maker.bulk_relax_maker is not None:
            raise ValueError("DFPT-QHA phonon maker must not re-relax the fixed-volume cell")


    @property
    def prev_calc_dir_argname(self) -> str:
        return "prev_dir"

    def make(self, structure: Any, supercell_matrix: Any = None, prev_dir: str | Path | None = None) -> Flow:
        """Use CommonQhaMaker's official EOS → phonon → analysis lifecycle."""

        return super().make(
            structure,
            supercell_matrix=supercell_matrix,
            prev_dir=prev_dir,
        )


__all__ = ["DfptPhononMaker", "DfptQhaMaker", "DfptVaspMaker"]
