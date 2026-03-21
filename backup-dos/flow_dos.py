# src/workflow/flow_dos.py

import logging
import traceback
from pathlib import Path
from typing import Any, TypedDict

from pymatgen.core import Structure
from atomate2.vasp.flows.core import DoubleRelaxMaker
from atomate2.vasp.jobs.core import RelaxMaker, StaticMaker, NonSCFMaker
from atomate2.vasp.sets.core import (
    RelaxSetGenerator,
    StaticSetGenerator,
    NonSCFSetGenerator,
)
from maggma.stores import MongoStore, MemoryStore
from jobflow.core.flow import Flow
from jobflow.core.store import JobStore
from jobflow.managers.local import run_locally

log = logging.getLogger(__name__)


class DosWorker:
    def __init__(
        self,
        vasp_args: dict[str, Any],
        potcar_functional: Any,
        incar_settings: dict[str, Any] = {},
    ):
        incar_settings = {
            "ISTART": 0,
            "ICHARG": 1,
            # Ionic relaxation
            "IBRION": 1,
            "ISIF": 3,
            "NSW": 50,
            # "EDIFFG": -0.02,
            # Electronic
            "ENCUT": 400,
            "EDIFF": 1e-5,
            "NELM": 200,
            "NELMIN": 4,
            "NELMDL": -12,
            "ALGO": "Fast",
            "PREC": "Normal",
            "LREAL": "Auto",
            # Magnetic
            "ISPIN": 2,
            "SYMPREC": 1e-5,
            # Parallel
            # "NCORE": 2,
            "KPAR": 4,
            # Output
            "LWAVE": False,
            "LCHARG": False,
            **incar_settings,
        }

        # Relaxation maker
        relax_maker = DoubleRelaxMaker.from_relax_maker(
            RelaxMaker(
                run_vasp_kwargs=vasp_args,
                stop_children_kwargs={"handle_unsuccessful": False},
                input_set_generator=RelaxSetGenerator(
                    user_potcar_functional=potcar_functional,
                    user_incar_settings={
                        **incar_settings,
                        "ISIF": 3,
                        "EDIFFG": -0.05,
                    },
                ),
            )
        )
        # Static maker
        static_maker = StaticMaker(
            run_vasp_kwargs=vasp_args,
            stop_children_kwargs={"handle_unsuccessful": False},
            input_set_generator=StaticSetGenerator(
                user_potcar_functional=potcar_functional,
                user_incar_settings={
                    **incar_settings,
                    "IBRION": -1,
                    "ISIF": 2,
                    "NSW": 0,
                    "ALGO": "Normal",
                    "LWAVE": True,
                    "LCHARG": True,
                },
            ),
        )
        # Non-SCF maker
        non_scf_maker = NonSCFMaker(
            run_vasp_kwargs=vasp_args,
            stop_children_kwargs={"handle_unsuccessful": False},
            input_set_generator=NonSCFSetGenerator(
                user_potcar_functional=potcar_functional,
                user_incar_settings={
                    **incar_settings,
                    "ICHARG": 11,
                    "IBRION": -1,
                    "NSW": 0,
                    # DOS
                    "LORBIT": 11,
                    # "NEDOS": 1000,
                    # "EMIN": -5,
                    # "EMAX": 5,
                    # "ISMEAR": -5,
                    # "NCORE": 1,
                    # "KPAR": 1,
                },
            ),
        )
        # DOS makers
        self.dos_makers: tuple[DoubleRelaxMaker, StaticMaker, NonSCFMaker] = (
            relax_maker,
            static_maker,
            non_scf_maker,
        )
        self.store = JobStore(MemoryStore(), additional_stores={"data": MemoryStore()})

    def run_dos(
        self,
        name: str,
        struct: Structure,
        flowdir: Path,
    ) -> dict[str, Any] | None:

        relax_maker, static_maker, non_scf_maker = self.dos_makers
        relax_job = relax_maker.make(struct)
        static_job = static_maker.make(
            relax_job.output.structure,
            prev_dir=relax_job.output.dir_name,
        )
        non_scf_job = non_scf_maker.make(
            static_job.output.structure,
            prev_dir=static_job.output.dir_name,
        )
        flow = Flow(
            [relax_job, static_job, non_scf_job],
            output=non_scf_job.output,
            name=name,
        )
        log.info(f"Running DOS flow for system {name}")
        try:
            run_locally(
                flow,
                store=self.store,
                create_folders=True,
                root_dir=flowdir,
                # ensure_success=True,
                # raise_immediately=True,
            )

            self.store.connect()

            # Debug: list all jobs
            log.info("\nAll jobs in store:")
            for doc in self.store.query(properties=["name", "uuid", "index"]):
                name = doc.get("name", "N/A")
                uuid = doc.get("uuid", "N/A")
                index = doc.get("index", "N/A")
                log.info(f"  Job Name: {name} | UUID: {uuid} | Index: {index}")

            job = self.store.query_one(
                criteria={"name": {"$regex": "non-scf"}},
                properties=["uuid", "index", "name"],
                sort={"index": -1},  # latest index first
            )
            if job is None:
                raise ValueError(f"No 'non-scf' job found for system {name}")
            output = self.store.get_output(uuid=job["uuid"], which="last", load=True)
            log.info(f"DOS flow of system {name} running successfully")
            self.store.close()
            return output
        except Exception as e:
            log.error(f"Error running DOS flow for system {name}: {e}")
            log.error(traceback.format_exc())
            self.store.close()
            return None
