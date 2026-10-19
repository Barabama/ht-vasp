"""DFPT-QHA Maker using the official EOS contract and a DFPT-only volume flow."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from atomate2.common.flows.eos import CommonEosMaker
from atomate2.common.flows.qha import CommonQhaMaker
from atomate2.vasp.jobs.base import BaseVaspMaker
from jobflow import Flow

from htvasp.dfpt_qha.flow_jobs import (
    DEFAULT_MESH,
    DEFAULT_TEMPERATURES,
    spawn_dfpt_volume_jobs,
)
from htvasp.dfpt_qha.jobs import DfptRelaxMaker
from htvasp.model.dfpt_qha_recipe import Recipe


@dataclass
class DfptQhaMaker(CommonQhaMaker):
    """Build the EOS -> DFPT -> validation -> thermal -> QHA graph.

    This intentionally does not subclass or invoke ``BasePhononMaker.make``.
    That official flow is a finite-displacement workflow; here the phonon data
    come from VASP's DFPT ``<dynmat>`` in each volume's ``vasprun.xml``, and the
    thermal properties are integrated from those dynamical matrices.

    ``run_vasp_kwargs`` is threaded to every child VASP job.  It is not
    optional in practice: without it atomate2 calls a bare ``vasp_std``, which
    on this cluster is not on PATH and is not launched under ``srun``.
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
    skip_analysis: bool = False
    run_vasp_kwargs: dict[str, Any] = field(default_factory=dict)
    #: element -> initial moment (muB); None uses the package table
    magmom: dict[str, float] | list[float] | None = None
    mesh: tuple[int, int, int] = DEFAULT_MESH
    temperatures: list[float] = field(default_factory=lambda: list(DEFAULT_TEMPERATURES))
    pressure_gpa: float = 0.0

    def __post_init__(self) -> None:
        if self.recipe is None:
            raise ValueError("DfptQhaMaker requires an explicit Recipe")
        if self.initial_relax_maker is None:
            self.initial_relax_maker = DfptRelaxMaker(
                name="DFPT-QHA initial relax",
                recipe=self.recipe,
                isif=3,
                run_vasp_kwargs=dict(self.run_vasp_kwargs),
                magmom=self.magmom,
            )
        if self.eos_relax_maker is None:
            # Fixed-volume ionic relaxation at each EOS volume.  The volume is
            # held by the strain applied by CommonEosMaker; ISIF=2 frees the
            # internal coordinates only.
            self.eos_relax_maker = DfptRelaxMaker(
                name="DFPT-QHA EOS relax",
                recipe=self.recipe,
                isif=2,
                run_vasp_kwargs=dict(self.run_vasp_kwargs),
                magmom=self.magmom,
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
        """Return EOS followed by the per-volume DFPT and QHA chain."""

        recipe = self.recipe
        if recipe is None:  # pragma: no cover - __post_init__ rejects this
            raise ValueError("DfptQhaMaker requires an explicit Recipe")

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
            recipe=recipe.model_dump(mode="json"),
            run_vasp_kwargs=dict(self.run_vasp_kwargs),
            magmom=self.magmom,
            supercell_matrix=supercell_matrix,
            mesh=self.mesh,
            temperatures=list(self.temperatures),
            pressure_gpa=self.pressure_gpa,
        )
        return Flow([eos_job, dfpt_job], output=dfpt_job.output, name=self.name)


__all__ = ["DfptQhaMaker"]
