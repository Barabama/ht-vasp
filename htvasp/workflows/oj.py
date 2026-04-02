"""
HT-VASP - OJ Workflows

OJ workflow: performs magnetic exchange calculation using the total energy difference method.
The OJ worker generates magnetic configurations, runs VASP, and solves for J parameters.
"""

import json
import logging
import traceback
from pathlib import Path
from datetime import datetime
from typing import Any, Literal

from pymatgen.core import Structure
from maggma.stores import JSONStore, MemoryStore
from custodian.vasp.handlers import VaspErrorHandler
from jobflow import Flow, JobStore

from htvasp.workflows.base import Worker
from htvasp.oj.input_set import OJInputSetGenerator
from htvasp.oj.maker import OJMaker
from htvasp.utils.run_locally import run_locally_custom

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
        oj_incar: dict[str, Any] | None = None,
        j_count: int = 4,
        dist_cutoff: float | None = None,
        magnetic_ion_types: list[str] = [],
        noncollinear: bool = False,
        base_spin: float | list[float] = 1.0,
        extend_poscar: tuple[int, int, int] = (2, 2, 2),
        **kwargs,
    ):
        self.store = None
        self.worker_name = worker_name

        default_incar = {
            "ENCUT": 500,
            "ISTART": 0,
            "ICHARG": 2,
            # Electronic
            "ISMEAR": 1,
            "SIGMA": 0.1,
            "ALGO": "Normal",
            "NELM": 100,
            "NELMIN": 6,
            "NELMDL": -6,
            # Ionic
            "IBRION": 2,
            "ISIF": 3,
            "NSW": 50,
            "POTIM": 0.2,
            "EDIFF": 1e-6,
            "EDIFFG": -0.01,
            # Magnetic
            "ISPIN": 2,
            # Precision
            "ISYM": 0,
            "LREAL": "Auto",
            "PREC": "Accurate",
            "SYMPREC": 1e-7,
            # Output
            "LWAVE": False,
            "LCHARG": False,
            "LORBIT": 11,
        }

        default_incar.update(oj_incar or {})

        oj_maker = OJMaker(
            run_vasp_kwargs={"handlers": [VaspErrorHandler()], **vasp_args},
            stop_children_kwargs={"handle_unsuccessful": False},
            input_set_generator=OJInputSetGenerator(
                user_potcar_functional=potcar_functional,
                user_incar_settings={
                    **default_incar,
                },
                j_count=j_count,
                dist_cutoff=dist_cutoff,
                magnetic_ion_types=magnetic_ion_types,
                noncollinear=noncollinear,
                base_spin=base_spin,
                extend_poscar=extend_poscar,
            ),
        )

        self.oj_maker = oj_maker

    def run_flow(
        self,
        name: str,
        structure: Structure,
        flow_dir: Path | str,
        dir_format: str = "{name}",
        store_path: Path | str = "",
        resume: bool = True,
    ) -> dict[str, Any] | None:
        """Run the OJ workflow

        Args:
            name: Structure name
            structure: Input structure
            flow_dir: Flow directory
            dir_format: Directory format for job folders
            store_path: Path to store results
            resume: Whether to resume from previously completed jobs (default: True)

        Returns:
            Dictionary containing OJ calculation results (J parameters, Tc, etc.)
        """
        flow_dir = Path(flow_dir)
        flow_dir.mkdir(parents=True, exist_ok=True)

        if not store_path:
            store_path = Path(flow_dir, "store.json").resolve()
        else:
            store_path = Path(store_path).resolve()

        self.store = JobStore(
            JSONStore(str(store_path), read_only=False),
            additional_stores={"data": MemoryStore()},
        )

        oj_flow = self.oj_maker.make(structure)
        flow = Flow([oj_flow], output=oj_flow.output, name=name)

        log.info(f"Running OJ flow for {name} in {flow_dir}")

        try:
            run_locally_custom(
                flow,
                store=self.store,
                root_dir=flow_dir,
                dir_format=dir_format,
                resume=resume,
            )

            self.store.connect()

            oj_job_doc = self.store.query_one(
                criteria={"name": {"$regex": "solve"}},
                properties=["uuid", "index", "name"],
                sort={"index": -1},
            )
            if not oj_job_doc:
                raise ValueError(f"No 'solve' job found in store {store_path}")

            oj_output = self.store.get_output(uuid=oj_job_doc["uuid"], which="last", load=True)

            log.info(f"OJ flow for {name} completed successfully")
            return oj_output

        except Exception as e:
            log.error(f"OJ flow for {name} failed: {e}")
            log.error(traceback.format_exc())
            return None

        finally:
            self.close()
