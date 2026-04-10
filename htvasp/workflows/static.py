"""
HT-VASP - Static Workflows

Structural relax and static calculation.
"""

import json
import logging
import traceback
from pathlib import Path
from datetime import datetime
from typing import Any, Literal

from pymatgen.core import Structure
from atomate2.vasp.flows.core import DoubleRelaxMaker
from atomate2.vasp.jobs.core import StaticMaker, TightRelaxMaker
from atomate2.vasp.sets.core import StaticSetGenerator, TightRelaxSetGenerator
from custodian.vasp.handlers import VaspErrorHandler
from maggma.stores import JSONStore, MemoryStore
from jobflow import Flow, JobStore

from htvasp.workflows.base import Worker
from htvasp.utils.run_locally import run_locally_custom

log = logging.getLogger(__name__)


class DateTimeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, datetime):
            return obj.isoformat()
        return super().default(obj)


class StaticWorker(Worker):
    """
    Worker for static calculation.
    """

    def __init__(
        self,
        worker_name: str,
        vasp_args: dict[str, Any],
        potcar_functional: Literal["PBE", "PBE54", "PBE_64"] = "PBE_64",
        global_incar: dict[str, Any] | None = None,
        relax_incar: dict[str, Any] | None = None,
        static_incar: dict[str, Any] | None = None,
        **kwargs,
    ):
        # Store will be initialized in run_flow to allow custom paths
        self.store = None
        self.worker_name = worker_name

        # Default INCAR settings
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
            "NSW": 100,
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
        }
        # Override with custom settings
        if global_incar:
            default_incar.update(global_incar)
        global_incar = default_incar
        relax_incar = relax_incar or {}
        static_incar = static_incar or {}

        # Structural relaxation
        relax_maker = DoubleRelaxMaker.from_relax_maker(
            TightRelaxMaker(
                name="r3_relax",
                run_vasp_kwargs={"handlers": [VaspErrorHandler()], **vasp_args},
                stop_children_kwargs={"handle_unsuccessful": False},
                input_set_generator=TightRelaxSetGenerator(
                    user_potcar_functional=potcar_functional,
                    user_incar_settings={
                        **global_incar,
                        "ISIF": 3,
                        "LWAVE": True,
                        "LCHARG": True,
                        **relax_incar,
                    },
                ),
            ),
        )
        # Static calculation
        static_maker = StaticMaker(
            run_vasp_kwargs={"handlers": [VaspErrorHandler()], **vasp_args},
            stop_children_kwargs={"handle_unsuccessful": False},
            input_set_generator=StaticSetGenerator(
                user_potcar_functional=potcar_functional,
                user_incar_settings={
                    **global_incar,
                    "ISTART": 1,
                    "ICHARG": 1,
                    "NELM": 200,
                    "IBRION": -1,
                    "ISIF": 2,
                    "NSW": 0,
                    "EDIFF": 1e-7,
                    "EDIFFG": 1e-6,
                    "LORBIT": 11,
                    **static_incar,
                },
            ),
        )

        # Static flow
        self.flow_makers: tuple[DoubleRelaxMaker, StaticMaker] = (relax_maker, static_maker)

    def run_flow(
        self,
        name: str,
        structure: Structure,
        flow_dir: Path | str,
        dir_format: str = "{name}",
        store_path: Path | str = "",
        fireworks: bool = True,
        resume: bool = True,
    ) -> dict[str, Any] | None:
        """
        Run the static workflow.

        Args:
            name: Structure name for identification
            structure: Input structure
            flow_dir: Directory to run the flow in
            dir_format: Directory format for job subdirectories
            store_path: Path to the JSON store for job tracking
            fireworks: Whether to use FireWorks for job management
            resume: Whether to resume from previous completed jobs

        Returns:
            Dictionary with static calculation result or None if failed
        """
        flow_dir = Path(flow_dir)
        flow_dir.mkdir(parents=True, exist_ok=True)

        # Initialize store
        if not store_path:
            store_path = Path(flow_dir, "store.json").resolve()
        else:
            store_path = Path(store_path).resolve()
        self.store = JobStore(
            JSONStore(str(store_path), read_only=False),
            additional_stores={"data": MemoryStore()},
        )

        # Build the workflow
        relax_maker, static_maker = self.flow_makers
        relax_job = relax_maker.make(structure)
        static_job = static_maker.make(
            relax_job.output.structure,
            prev_dir=relax_job.output.dir_name,
        )
        flow = Flow(
            [relax_job, static_job],
            output=static_job.output,
            name=name,
        )
        log.info(f"Running Static flow for struct {name} in {flow_dir}")

        try:
            run_locally_custom(
                flow,
                store=self.store,
                root_dir=flow_dir,
                dir_format=dir_format,
            )

            self.store.connect()
            job = self.store.query_one(
                criteria={"name": {"$regex": "static"}},
                properties=["uuid", "index", "name"],
                sort={"index": -1},
            )
            if not job:
                raise ValueError(f"No 'static' job found in store {store_path}")

            output = self.store.get_output(uuid=job["uuid"], which="last", load=True)
            log.info(f"Static flow for sturct {name} completed successfully")
            return output

        except Exception as e:
            log.error(f"Static flow for struct {name} failed: {e}")
            log.error(traceback.format_exc())
            return None

        finally:
            self.close()
