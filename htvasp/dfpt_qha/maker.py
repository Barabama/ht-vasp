"""DFPT-QHA Maker using the official EOS contract and a DFPT-only volume flow."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from atomate2.common.flows.eos import CommonEosMaker
from atomate2.common.flows.qha import CommonQhaMaker
from atomate2.vasp.jobs.base import BaseVaspMaker
from jobflow import Flow

from htvasp.dfpt_qha.flow_jobs import spawn_dfpt_volume_jobs
from htvasp.dfpt_qha.input_set import DfptInputSetGenerator
from htvasp.dfpt_qha.jobs import DfptVaspMaker
from htvasp.model.dfpt_qha_recipe import Recipe


@dataclass
class DfptQhaMaker(CommonQhaMaker):
    """Build the fixed-volume EOS → DFPT → validation graph.

    This intentionally does not subclass or invoke ``BasePhononMaker.make``.
    That official flow is a finite-displacement workflow; this Maker's phonon
    data come from VASP's DFPT ``<dynmat>`` in each volume's vasprun.xml.

    The QHA thermodynamic analysis job is deliberately not appended yet: it
    will consume validated ``DfptTaskDoc`` documents once its conversion from
    VASP dynamical matrices to the common phonon/QHA schema is implemented.
    """

    name: str = "DFPT-QHA Maker"
    recipe: Recipe | None = None
    initial_relax_maker: BaseVaspMaker | None = None
    eos_relax_maker: BaseVaspMaker | None = None
    linear_strain: tuple[float, float] = (-0.05, 0.05)
    number_of_frames: int = 6
    min_length: float | None = 20.0
    max_length: float | None = None
    ignore_imaginary_modes: bool = False
    skip_analysis: bool = True

    def __post_init__(self) -> None:
        if self.recipe is None:
            raise ValueError("DfptQhaMaker requires an explicit Recipe")
        if self.initial_relax_maker is None:
            self.initial_relax_maker = DfptVaspMaker(
                name="DFPT-QHA initial relax", recipe=self.recipe,
                input_set_generator=DfptInputSetGenerator(recipe=self.recipe)
            )
        if self.eos_relax_maker is None:
            # This is a fixed-volume ionic relaxation at each EOS volume, not a
            # static job.  DfptVaspMaker stage is temporarily only the input
            # contract; the EOS child must be supplied by the production
            # EosRelaxMaker wiring before real execution.
            self.eos_relax_maker = DfptVaspMaker(
                name="DFPT-QHA EOS relax", recipe=self.recipe,
                input_set_generator=DfptInputSetGenerator(recipe=self.recipe)
            )

    @property
    def prev_calc_dir_argname(self) -> str:
        return "prev_dir"

    def make(
        self,
        structure: Any,
        *,
        supercell_matrix: Any | None = None,
        prev_dir: str | Path | None = None,
    ) -> Flow:
        """Return EOS followed by one DFPT/collector/validation chain per volume."""

        eos = CommonEosMaker(
            initial_relax_maker=self.initial_relax_maker,
            eos_relax_maker=self.eos_relax_maker,
            static_maker=None,
            postprocessor=None,
            linear_strain=self.linear_strain,
            number_of_frames=self.number_of_frames,
        )
        eos_job = eos.make(structure)
        dfpt_job = spawn_dfpt_volume_jobs(
            eos_job.output,
            recipe=self.recipe.model_dump(mode="json"),
            run_vasp_kwargs=getattr(self.initial_relax_maker, "run_vasp_kwargs", {}),
        )
        return Flow([eos_job, dfpt_job], output=dfpt_job.output, name=self.name)


__all__ = ["DfptQhaMaker"]
