"""
HT-VASP - Relax Workflows

Structural relax workflows using atomate2 DoubleRelaxMaker.
"""

import logging
import traceback
from pathlib import Path
from typing import Any, Literal

from pymatgen.core import Structure
from atomate2.vasp.flows.core import DoubleRelaxMaker
from atomate2.vasp.jobs.core import TightRelaxMaker
from atomate2.vasp.sets.core import TightRelaxSetGenerator
from custodian.vasp.handlers import VaspErrorHandler
from maggma.stores import JSONStore, MemoryStore
from jobflow import JobStore

from htvasp.workflows.base import Worker
from htvasp.utils import run_locally_custom

log = logging.getLogger(__name__)


class RelaxWorker(Worker):
    """
    Worker for structural relaxation using atomate2 DoubleRelaxMaker。
    Performs R7 volume relaxation (ISIF=7) and R3 full relaxation (ISIF=3).
    """

    def __init__(
        self,
        worker_name: str,
        vasp_args: dict[str, Any],
        potcar_functional: Literal["PBE", "PBE_54", "PBE_64"] = "PBE_64",
        global_incar: dict[str, Any] | None = None,
        r7_incar: dict[str, Any] | None = None,
        r3_incar: dict[str, Any] | None = None,
        **kwargs,
    ):
        # Store will be initialized in run_flow to allow custom paths
        self.store = None
        self.worker_name = worker_name

        # Default INCAR settings
        default_incar = {
            "ENCUT": 400,
            "ISTART": 0,
            "ICHARG": 2,
            # Electronic
            "ISMEAR": 1,
            "SIGMA": 0.2,
            "ALGO": "Fast",
            "NELM": 100,
            "NELMIN": 4,
            "NELMDL": -12,
            # Ionic
            "IBRION": 2,
            "ISIF": 3,
            "NSW": 10,
            "POTIM": 0.2,
            "EDIFF": 1e-5,
            "EDIFFG": 1e-4,
            # Magnetic
            "ISPIN": 2,
            # Precision
            "ISYM": 0,
            "LREAL": "Auto",
            "PREC": "Normal",
            "SYMPREC": 1e-7,
            # Output
            "LWAVE": False,
            "LCHARG": False,
        }
        # Override with custom settings
        if global_incar:
            default_incar.update(global_incar)
        global_incar = default_incar
        r7_incar = r7_incar or {}
        r3_incar = r3_incar or {}

        # R7 structural relaxation
        relax_r7_maker = DoubleRelaxMaker.from_relax_maker(
            TightRelaxMaker(
                name="r7_relax",
                run_vasp_kwargs={"handlers": [VaspErrorHandler()], **vasp_args},
                stop_children_kwargs={"handle_unsuccessful": False},
                input_set_generator=TightRelaxSetGenerator(
                    user_potcar_functional=potcar_functional,
                    user_incar_settings={
                        **global_incar,
                        "NELM": 200,
                        "ISIF": 7,
                        "NSW": 20,
                        "EDIFF": 1e-5,
                        "EDIFFG": 1e-4,
                        **r7_incar,
                    },
                ),
            ),
        )
        # R3 structural relaxation
        relax_r3_maker = DoubleRelaxMaker.from_relax_maker(
            TightRelaxMaker(
                name="r3_relax",
                run_vasp_kwargs={"handlers": [VaspErrorHandler()], **vasp_args},
                stop_children_kwargs={"handle_unsuccessful": False},
                input_set_generator=TightRelaxSetGenerator(
                    user_potcar_functional=potcar_functional,
                    user_incar_settings={
                        **global_incar,
                        "NELM": 300,
                        "ISIF": 3,
                        "NSW": 50,
                        "EDIFF": 1e-5,
                        "EDIFFG": -0.05,
                        **r3_incar,
                    },
                ),
            ),
        )

        # Relax flow
        self.relax_flow = DoubleRelaxMaker(
            relax_maker1=relax_r7_maker,
            relax_maker2=relax_r3_maker,
        )

    def run_flow(
        self,
        name: str,
        structure: Structure,
        flow_dir: Path | str,
        dir_format: str = "{name}",
        store_path: Path | str = "",
    ) -> dict[str, Any] | None:
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

        flow = self.relax_flow.make(structure)
        log.info(f"Running relax flow for struct {name} in {flow_dir}")

        try:
            run_locally_custom(
                flow,
                store=self.store,
                root_dir=flow_dir,
                dir_format=dir_format,
            )

            self.store.connect()
            job = self.store.query_one(
                criteria={"name": {"$regex": "r3_relax"}},
                properties=["uuid", "index", "name"],
                sort={"index": -1},
            )
            if not job:
                raise ValueError(f"No 'r3_relax' job found in the store {store_path}")

            output = self.store.get_output(uuid=job["uuid"], which="last", load=True)
            log.info(f"Relax flow for struct {name} completed successfully")
            return output

        except Exception as e:
            log.error(f"Relax flow for struct {name} failed: {e}")
            log.error(traceback.format_exc())
            return None

        finally:
            self.close()
