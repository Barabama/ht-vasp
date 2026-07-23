"""HT-VASP OJ Workflow - Magnetic exchange calculation using total energy difference method."""

import logging
from pathlib import Path
from typing import Any

from jobflow import Flow
from pymatgen.core import Structure

from atomate2.vasp.jobs.core import RelaxMaker
from atomate2.vasp.sets.core import RelaxSetGenerator
from htvasp.oj.input_set import OJInputSetGenerator
from htvasp.oj.maker import OJMaker
from htvasp.workflows.base import Worker

log = logging.getLogger(__name__)


class OJWorker(Worker):
    """Worker for OJ magnetic exchange calculations.

    Performs OJ calculation using the total energy difference method:
    1. Generates magnetic configurations from input structure
    2. Runs VASP for each magnetic configuration
    3. Solves for exchange parameters J and Curie temperature

    Args:
        worker_name: Worker name
        vasp_args: VASP execution arguments
        potcar_functional: POTCAR functional
        global_incar: Global INCAR settings
        j_count: Number of J pairs to consider
        dist_cutoff: Distance cutoff for J pairs
        magnetic_ion_types: List of magnetic ion types
        noncollinear: Whether to use noncollinear magnetism
        base_spin: Base spin value
        extend_poscar: Supercell extension factors
    """

    def __init__(
        self,
        worker_name: str = "oj-worker",
        vasp_args: dict[str, Any] | None = None,
        potcar_functional: str = "PBE_64",
        global_incar: dict[str, Any] | None = None,
        j_count: int = 4,
        dist_cutoff: float | None = None,
        magnetic_ion_types: list[str] | None = None,
        noncollinear: bool = False,
        base_spin: float | list[float] = 1.0,
        extend_poscar: tuple[int, int, int] = (2, 2, 2),
        **kwargs,
    ):
        super().__init__(
            worker_name=worker_name,
            vasp_args=vasp_args,
            potcar_functional=potcar_functional,
            global_incar=global_incar,
        )

        self.flow_maker = OJMaker(
            run_vasp_kwargs=self.run_vasp_kwargs,
            # stop_children_kwargs={"handle_unsuccessful": False},
            input_set_generator=OJInputSetGenerator(
                user_potcar_functional=self.potcar_functional,
                user_incar_settings={**self.global_incar, "LORBIT": 10},
                j_count=j_count,
                dist_cutoff=dist_cutoff,
                magnetic_ion_types=magnetic_ion_types or [],
                noncollinear=noncollinear,
                base_spin=base_spin,
                extend_poscar=extend_poscar,
            ),
        )

    def _make_flow(self, structure: Structure, prev_dir: str | None = None) -> Flow:
        return self.flow_maker.make(structure)

    def get_result(self, output_job_name: str = "solve") -> dict[str, Any] | None:
        return super().get_result(output_job_name)


