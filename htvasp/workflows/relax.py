"""
HT-VASP - Relax Workflows

Structural relax workflows using atomate2 DoubleRelaxMaker.
"""

import logging
from typing import Any

from jobflow import Flow
from pymatgen.core import Structure
from atomate2.vasp.flows.core import DoubleRelaxMaker
from atomate2.vasp.jobs.core import RelaxMaker
from atomate2.vasp.sets.core import RelaxSetGenerator

from htvasp.workflows.base import Worker

log = logging.getLogger(__name__)


class RelaxWorker(Worker):
    """
    Worker for structural relaxation using atomate2 DoubleRelaxMaker。
    Performs R7 volume relaxation (ISIF=7) and R3 full relaxation (ISIF=3).
    """

    def __init__(
        self,
        worker_name: str = "relax-worker",
        vasp_args: dict[str, Any] | None = None,
        potcar_functional="PBE_64",
        global_incar: dict[str, Any] | None = None,
        r7_incar: dict[str, Any] | None = None,
        r3_incar: dict[str, Any] | None = None,
        **kwargs,
    ):
        # Initialize base Worker
        super().__init__(
            worker_name=worker_name,
            vasp_args=vasp_args,
            potcar_functional=potcar_functional,
            global_incar=global_incar,
        )

        r7_incar = r7_incar or {}
        r3_incar = r3_incar or {}

        # R7 structural relaxation
        relax_r7_maker = DoubleRelaxMaker.from_relax_maker(
            RelaxMaker(
                name="r7 relax",
                run_vasp_kwargs=self.run_vasp_kwargs,
                stop_children_kwargs={"handle_unsuccessful": False},
                input_set_generator=RelaxSetGenerator(
                    user_potcar_functional=self.potcar_functional,
                    user_incar_settings={
                        **self.global_incar,
                        "ISIF": 7,
                        **r7_incar,
                    },
                ),
            ),
        )
        # R3 structural relaxation
        relax_r3_maker = DoubleRelaxMaker.from_relax_maker(
            RelaxMaker(
                name="r3 relax",
                run_vasp_kwargs=self.run_vasp_kwargs,
                stop_children_kwargs={"handle_unsuccessful": False},
                input_set_generator=RelaxSetGenerator(
                    user_potcar_functional=self.potcar_functional,
                    user_incar_settings={
                        **self.global_incar,
                        "ISIF": 3,
                        **r3_incar,
                    },
                ),
            ),
        )

        # Relax flow
        self.flow_maker = DoubleRelaxMaker(
            relax_maker1=relax_r7_maker,
            relax_maker2=relax_r3_maker,
        )

    def _make_flow(self, structure: Structure) -> Flow:
        return self.flow_maker.make(structure)

    def get_result(self, output_job_name: str = "r3 relax") -> dict[str, Any] | None:
        return super().get_result(output_job_name)
