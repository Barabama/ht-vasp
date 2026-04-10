"""
HT-VASP - OJ Workflows

OJ workflow: performs magnetic exchange calculation using the total energy difference method.
The OJ worker generates magnetic configurations, runs VASP, and solves for J parameters.
"""

import json
import logging
from datetime import datetime
from typing import Any, Literal
from pymatgen.core import Structure
from custodian.vasp.handlers import VaspErrorHandler
from jobflow import Flow

from htvasp.workflows.base import Worker
from htvasp.oj.input_set import OJInputSetGenerator
from htvasp.oj.maker import OJMaker

log = logging.getLogger(__name__)


class DateTimeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, datetime):
            return obj.isoformat()
        return super().default(obj)


class OJWorker(Worker):
    """Worker for OJ magnetic exchange calculations

    Performs OJ calculation using the total energy difference method:
    1. Generates magnetic configurations from input structure
    2. Runs VASP for each magnetic configuration
    3. Solves for exchange parameters J and Curie temperature

    The input structure is used as the base crystal structure.
    OJ internally generates supercell and various magnetic configurations (MAGMOM).
    No prior structure relaxation step is required by this worker
    (the ostravaj README does not mandate pre-relaxation).
    You may optionally pass a pre-relaxed structure for better accuracy.

    Args:
        worker_name: Worker name
        vasp_args: VASP execution arguments
        potcar_functional: POTCAR functional
        oj_incar: OJ-specific INCAR settings
        j_count: Number of J pairs to consider
        dist_cutoff: Distance cutoff for J pairs
        magnetic_ion_types: List of magnetic ion types
        noncollinear: Whether to use noncollinear magnetism
        base_spin: Base spin value
        extend_poscar: Supercell extension factors


    Example:
        >>> from pymatgen.core import Structure
        >>> from htvasp.workflows import OJWorker
        >>>
        >>> structure = Structure.from_file("POSCAR")
        >>> worker = OJWorker(
        ...     worker_name="fe_co_oj",
        ...     j_count=3,
        ...     user_incar_settings={"ENCUT": 520, "KPAR": 4},
        ... )
    """

    def __init__(
        self,
        worker_name: str,
        vasp_args: dict[str, Any],
        potcar_functional: Literal["PBE", "PBE_54", "PBE_64"] = "PBE_64",
        use_fireworks: bool = True,
        global_incar: dict[str, Any] | None = None,
        j_count: int = 4,
        dist_cutoff: float | None = None,
        magnetic_ion_types: list[str] = [],
        noncollinear: bool = False,
        base_spin: float | list[float] = 1.0,
        extend_poscar: tuple[int, int, int] = (2, 2, 2),
        **kwargs,
    ):
        super().__init__(
            worker_name=worker_name,
            vasp_args=vasp_args,
            potcar_functional=potcar_functional,
            use_fireworks=use_fireworks,
            global_incar=global_incar,
        )

        oj_maker = OJMaker(
            run_vasp_kwargs={"handlers": [VaspErrorHandler()], **vasp_args},
            stop_children_kwargs={"handle_unsuccessful": False},
            input_set_generator=OJInputSetGenerator(
                user_potcar_functional=potcar_functional,
                user_incar_settings={
                    **self.global_incar,
                },
                j_count=j_count,
                dist_cutoff=dist_cutoff,
                magnetic_ion_types=magnetic_ion_types,
                noncollinear=noncollinear,
                base_spin=base_spin,
                extend_poscar=extend_poscar,
            ),
        )

        self.flow_maker = oj_maker

    def _make_flow(self, structure: Structure) -> Flow:
        return self.flow_maker.make(structure)
    
    def get_result(self, output_job_name: str = "solve") -> dict[str, Any] | None:
        return super().get_result(output_job_name)
