"""
HT-VASP - Static-OJ Workflows

Combined workflow: relax + static + OJ (magnetic exchange)
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
from jobflow.core.store import JobStore
from jobflow.core.flow import Flow

from htvasp.workflows.base import Worker
from htvasp.oj.config import OJConfig
from htvasp.oj.maker import OJMaker
from htvasp.utils.run_locally import run_locally_custom

log = logging.getLogger(__name__)


class DateTimeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, datetime):
            return obj.isoformat()
        return super().default(obj)


class StaticOJWorker(Worker):
    """Worker for combined relax + static + OJ calculations

    Performs:
    1. Structural relaxation (R3 full relaxation)
    2. Static calculation (to get total magnetic moment)
    3. OJ calculation (magnetic exchange interactions)

    Args:
        worker_name: Name for the worker
        oj_config: OJ configuration
        potcar_functional: POTCAR functional type
        global_incar: Custom INCAR settings to override defaults
        relax_incar: Custom INCAR settings for relaxation
        static_incar: Custom INCAR settings for static calculation
        vasp_args: VASP execution arguments

    Example:
        >>> from pymatgen.core import Structure
        >>> from htvasp.workflows import StaticOJWorker
        >>> from htvasp.oj import OJConfig
        >>>
        >>> structure = Structure.from_file("POSCAR")
        >>> oj_config = OJConfig(j_count=3, magnetic_ion_types=["Fe", "Co"])
        >>> worker = StaticOJWorker(
        ...     worker_name="fe_co_static_oj",
        ...     oj_config=oj_config,
        ...     global_incar={"ENCUT": 520},
        ... )
        >>> output = worker.run_flow("FeCo", structure, "./static_oj_work")
    """

    def __init__(
        self,
        worker_name: str,
        oj_config: OJConfig | None = None,
        potcar_functional: Literal["PBE", "PBE_54", "PBE_64"] = "PBE_64",
        global_incar: dict[str, Any] | None = None,
        relax_incar: dict[str, Any] | None = None,
        static_incar: dict[str, Any] | None = None,
        vasp_args: dict[str, Any] | None = None,
        **kwargs,
    ):
        self.store = None
        self.worker_name = worker_name
        self.oj_config = oj_config or OJConfig()
        self.potcar_functional = potcar_functional
        self.vasp_args = vasp_args or {}

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

        # Apply OJ config INCAR settings
        if self.oj_config.incar:
            global_incar.update(self.oj_config.incar)

        vasp_cmd = self.vasp_args.get("vasp_cmd", "vasp_std")

        # Structural relaxation (R3 full relaxation)
        relax_maker = DoubleRelaxMaker.from_relax_maker(
            TightRelaxMaker(
                name="r3_relax",
                run_vasp_kwargs={"handlers": [VaspErrorHandler()], "vasp_cmd": vasp_cmd},
                stop_children_kwargs={"handle_unsuccessful": False},
                input_set_generator=TightRelaxSetGenerator(
                    user_potcar_functional=potcar_functional,
                    user_incar_settings={
                        **global_incar,
                        "ISIF": 3,
                        "LWAVE": True,
                        "LCHARG": True,
                        **(relax_incar or {}),
                    },
                ),
            ),
        )

        # Static calculation
        static_maker = StaticMaker(
            run_vasp_kwargs={"handlers": [VaspErrorHandler()], "vasp_cmd": vasp_cmd},
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
                    **(static_incar or {}),
                },
            ),
        )

        # OJ calculation
        oj_maker = OJMaker(
            name=f"{worker_name}_oj",
            config=self.oj_config,
            vasp_cmd=vasp_cmd,
            potcar_functional=potcar_functional,
        )

        self.flow_makers: tuple[DoubleRelaxMaker, StaticMaker, OJMaker] = (
            relax_maker,
            static_maker,
            oj_maker,
        )

    def run_flow(
        self,
        name: str,
        structure: Structure,
        flow_dir: Path | str,
        dir_format: str = "{name}",
        store_path: Path | str = "",
        resume: bool = True,
    ) -> dict[str, Any] | None:
        """Run the combined relax + static + OJ workflow

        Args:
            name: Structure name
            structure: Input structure
            flow_dir: Flow directory
            dir_format: Directory format for job folders
            store_path: Path to store results
            resume: Whether to resume from previously completed jobs (default: True)

        Returns:
            Dictionary containing:
            - static_output: Static calculation results (total magnetic moment, etc.)
            - oj_output: OJ calculation results (J parameters, Tc, etc.)
            - combined: Combined results including both static and OJ data
        """
        flow_dir = Path(flow_dir)
        flow_dir.mkdir(parents=True, exist_ok=True)

        if not store_path:
            store_path = Path(flow_dir, "store.json").resolve()
        else:
            store_path = Path(store_path).resolve()

        self.store = JobStore(
            JSONStore(store_path, read_only=False),
            additional_stores={"data": MemoryStore()},
        )

        relax_maker, static_maker, oj_maker = self.flow_makers

        # Build the flow: relax -> static -> oj
        relax_job = relax_maker.make(structure)
        static_job = static_maker.make(
            relax_job.output.structure,
            prev_dir=relax_job.output.dir_name,
        )
        oj_flow = oj_maker.make(static_job.output.structure)

        # Combine all jobs into a single flow
        all_jobs = [relax_job, static_job, *oj_flow.jobs]
        flow = Flow(
            all_jobs,
            output={
                "static_output": static_job.output,
                "oj_output": oj_flow.output,
            },
            name=name,
        )

        log.info(f"Running Static-OJ flow for {name} in {flow_dir}")

        try:
            run_locally_custom(
                flow,
                store=self.store,
                root_dir=flow_dir,
                dir_format=dir_format,
                resume=resume,
            )

            self.store.connect()

            # Get static job output
            static_job_doc = self.store.query_one(
                criteria={"name": {"$regex": "static"}},
                properties=["uuid", "index", "name"],
                sort={"index": -1},
            )
            if not static_job_doc:
                raise ValueError(f"No 'static' job found in store {store_path}")

            static_output = self.store.get_output(
                uuid=static_job_doc["uuid"], which="last", load=True
            )

            # Get OJ solve job output
            oj_job_doc = self.store.query_one(
                criteria={"name": {"$regex": "solve"}},
                properties=["uuid", "index", "name"],
                sort={"index": -1},
            )
            if not oj_job_doc:
                raise ValueError(f"No 'solve' job found in store {store_path}")

            oj_output = self.store.get_output(uuid=oj_job_doc["uuid"], which="last", load=True)

            # Combine results
            combined_output = {
                **static_output,
                **oj_output,
            }

            log.info(f"Static-OJ flow for {name} completed successfully")
            return combined_output

        except Exception as e:
            log.error(f"Static-OJ flow for {name} failed: {e}")
            log.error(traceback.format_exc())
            return None

        finally:
            self.close()

    def create_flow(self, structure: Structure, name: str = "static_oj") -> tuple[Flow, tuple]:
        """Create a Flow for external execution

        Note: This method creates jobs that can only be used once.
        For multiple flow creations, use run_flow() instead.

        Returns:
            Tuple of (flow, makers) for use with run_locally or FireWorks
        """
        relax_maker, static_maker, oj_maker = self.flow_makers

        relax_job = relax_maker.make(structure)
        static_job = static_maker.make(
            relax_job.output.structure,
            prev_dir=relax_job.output.dir_name,
        )
        oj_flow = oj_maker.make(static_job.output.structure)

        # Extract jobs from oj_flow and add to combined flow
        # Note: This will fail if jobs were already added to another flow
        all_jobs = [relax_job, static_job] + list(oj_flow.jobs)
        flow = Flow(
            all_jobs,
            output={
                **static_job.output,
                **oj_flow.output,
            },
            name=name,
        )

        return flow, (relax_maker, static_maker, oj_maker)
