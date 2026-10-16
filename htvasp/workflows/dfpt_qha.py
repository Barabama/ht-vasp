"""DFPT-QHA Worker: the production wrapper around :class:`DfptQhaMaker`.

The Worker owns the run lifecycle only -- it builds the Maker, runs the flow
through ``run_locally_custom`` (with resume), and returns the terminal analysis
job's output.  Everything scientific lives in the Maker and its children.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from jobflow import Flow
from pymatgen.core import Structure

from htvasp.dfpt_qha.maker import DfptQhaMaker
from htvasp.model.dfpt_qha_recipe import Recipe
from htvasp.workflows.base import Worker


class DfptQhaWorker(Worker):
    """Execute one endpoint's DFPT-QHA flow."""

    def __init__(
        self,
        recipe: Recipe,
        worker_name: str = "dfpt-qha-worker",
        vasp_args: dict[str, Any] | None = None,
        potcar_functional: str = "PBE",
        global_incar: dict[str, Any] | None = None,
        linear_strain: tuple[float, float] = (-0.05, 0.05),
        number_of_frames: int = 6,
        supercell_matrix: Any | None = None,
        mesh: tuple[int, int, int] = (12, 12, 12),
        temperatures: list[float] | None = None,
        pressure_gpa: float = 0.0,
        **kwargs: Any,
    ) -> None:
        super().__init__(
            worker_name=worker_name,
            vasp_args=vasp_args,
            potcar_functional=potcar_functional,
            global_incar=global_incar,
        )
        self.recipe = recipe
        self.supercell_matrix = supercell_matrix
        self.flow_maker = DfptQhaMaker(
            name=worker_name,
            recipe=recipe,
            linear_strain=linear_strain,
            number_of_frames=number_of_frames,
            run_vasp_kwargs=dict(self.run_vasp_kwargs),
            mesh=mesh,
            temperatures=temperatures,
            pressure_gpa=pressure_gpa,
        )

    def _make_flow(
        self, structure: Structure, prev_dir: Path | str | None = None
    ) -> Flow:
        return self.flow_maker.make(
            structure,
            supercell_matrix=self.supercell_matrix,
            prev_dir=prev_dir,
        )

    def get_result(self, output_job_name: str = "analyze DFPT-QHA") -> dict[str, Any] | None:
        return super().get_result(output_job_name)


__all__ = ["DfptQhaWorker"]
