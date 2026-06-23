"""
HT-VASP - Static Workflows

Structural relax and static calculation.
"""

import logging
from pathlib import Path
from typing import Any

from jobflow import Flow
from pymatgen.core import Structure
from atomate2.vasp.flows.core import DoubleRelaxMaker
from atomate2.vasp.jobs.core import RelaxMaker, StaticMaker
from atomate2.vasp.sets.core import RelaxSetGenerator, StaticSetGenerator

from htvasp.workflows.base import Worker

log = logging.getLogger(__name__)


class StaticWorker(Worker):
    """
    Worker for static calculation.
    """

    def __init__(
        self,
        worker_name: str = "static-worker",
        vasp_args: dict[str, Any] | None = None,
        potcar_functional: str = "PBE_64",
        global_incar: dict[str, Any] | None = None,
        relax_incar: dict[str, Any] | None = None,
        static_incar: dict[str, Any] | None = None,
        relax_reciprocal_density: int | None = None,
        static_reciprocal_density: int | None = None,
        **kwargs,
    ):
        """
        Initialize StaticWorker.

        Args:
            worker_name: Name of the worker
            vasp_args: VASP command and handler settings
            potcar_functional: POTCAR functional type
            global_incar: Global INCAR settings
            relax_incar: Relax-specific INCAR settings
            static_incar: Static-specific INCAR settings
            reciprocal_density: K-point density for DOS (uniform mode)
            relax_reciprocal_density: K-point density for relax step.
                Default None uses atomate2 default (64).
            static_reciprocal_density: K-point density for static step.
                Default None uses atomate2 default (64).
            **kwargs: Additional keyword arguments
        """
        # Initialize base Worker
        super().__init__(
            worker_name=worker_name,
            vasp_args=vasp_args,
            potcar_functional=potcar_functional,
            global_incar=global_incar,
        )

        relax_incar = relax_incar or {}
        static_incar = static_incar or {}

        # Structural relaxation (ISIF=3)
        relax_maker = DoubleRelaxMaker.from_relax_maker(
            RelaxMaker(
                run_vasp_kwargs=self.run_vasp_kwargs,
                stop_children_kwargs={"handle_unsuccessful": False},
                copy_vasp_kwargs={"additional_vasp_files": ("WAVECAR",)},
                input_set_generator=RelaxSetGenerator(
                    user_potcar_functional=self.potcar_functional,
                    user_kpoints_settings=(
                        {"reciprocal_density": relax_reciprocal_density}
                        if relax_reciprocal_density is not None
                        else {}
                    ),
                    user_incar_settings={
                        **self.global_incar,
                        "ISTART": 1,
                        "ISIF": 3,
                        "LWAVE": True,
                        **relax_incar,
                    },
                ),
            ),
        )
        # Static calculation
        static_maker = StaticMaker(
            run_vasp_kwargs=self.run_vasp_kwargs,
            stop_children_kwargs={"handle_unsuccessful": False},
            copy_vasp_kwargs={"additional_vasp_files": ("WAVECAR",)},
            input_set_generator=StaticSetGenerator(
                user_potcar_functional=self.potcar_functional,
                user_kpoints_settings=(
                    {"reciprocal_density": static_reciprocal_density}
                    if static_reciprocal_density is not None
                    else {}
                ),
                user_incar_settings={
                    **self.global_incar,
                    "ISTART": 1,
                    "IBRION": -1,
                    "ISIF": 2,
                    "NSW": 0,
                    "LWAVE": True,
                    "LCHARG": True,
                    **static_incar,
                },
            ),
        )

        # Static flow
        self.flow_makers: tuple[DoubleRelaxMaker, StaticMaker] = (
            relax_maker,
            static_maker,
        )

    def _make_flow(self, structure: Structure, prev_dir: Path | str | None = None) -> Flow:
        relax_maker, static_maker = self.flow_makers
        relax_job = relax_maker.make(structure, prev_dir)
        static_job = static_maker.make(
            relax_job.output.structure,
            prev_dir=relax_job.output.dir_name,
        )
        return Flow([relax_job, static_job], output=static_job.output)

    def get_result(self, output_job_name: str = "static") -> dict[str, Any] | None:
        return super().get_result(output_job_name)
