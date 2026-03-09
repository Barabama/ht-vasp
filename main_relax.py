# %%
# src/workflow/flow_relax.py

import logging
import traceback
from pathlib import Path
from typing import Any, TypedDict

from pymatgen.core import Structure
from atomate2.vasp.flows.core import DoubleRelaxMaker
from atomate2.vasp.jobs.core import RelaxMaker
from atomate2.vasp.sets.core import RelaxSetGenerator
from maggma.stores import MongoStore, MemoryStore
from jobflow.core.store import JobStore
from jobflow.managers.local import run_locally

log = logging.getLogger(__name__)


class RelaxWorker:
    def __init__(
        self,
        vasp_args: dict[str, Any],
        potcar_functional: Any,
        incar_settings: dict[str, Any] = {},
    ):
        incar_settings = {
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
            "AMIX": 0.04,
            "BMIX": 1e-4,
            "AMIX_MAG": 0.8,
            "BMIX_MAG": 1e-4,
            # Precision
            "KPAR": 2,
            "NCORE": 1,
            "ISYM": 2,
            "LREAL": "Auto",
            "PREC": "Normal",
            "SYMPREC": 1e-7,
            # Output
            "LWAVE": False,
            "LCHARG": False,
            "LORBIT": 11,
            # Overrides
            **incar_settings,
        }
        # R7 structural relaxation
        relax_r7_maker = DoubleRelaxMaker.from_relax_maker(
            RelaxMaker(
                run_vasp_kwargs=vasp_args,
                stop_children_kwargs={"handle_unsuccessful": False},
                input_set_generator=RelaxSetGenerator(
                    user_potcar_functional=potcar_functional,
                    user_incar_settings={
                        **incar_settings,
                        "NELM": 200,
                        "ISIF": 7,
                        "NSW": 20,
                        "EDIFF": 1e-5,
                        "EDIFFG": 1e-4,
                    },
                ),
            ),
        )
        # R3 structural relaxation
        relax_r3_maker = DoubleRelaxMaker.from_relax_maker(
            RelaxMaker(
                run_vasp_kwargs=vasp_args,
                stop_children_kwargs={"handle_unsuccessful": False},
                input_set_generator=RelaxSetGenerator(
                    user_potcar_functional=potcar_functional,
                    user_incar_settings={
                        **incar_settings,
                        "NELM": 300,
                        "ISIF": 3,
                        "NSW": 50,
                        "EDIFF": 1e-5,
                        "EDIFFG": -0.05,
                    },
                ),
            ),
        )
        
        # Relax flow
        self.relax_flow = DoubleRelaxMaker(
            relax_maker1=relax_r7_maker,
            relax_maker2=relax_r3_maker,
        )


        self.store = JobStore(MemoryStore(), additional_stores={"data": MemoryStore()})

    def run_relax(
        self,
        name: str,
        struct: Structure,
        flowdir: Path,
    ) -> dict[str, Any] | None:

        flow = self.relax_flow.make(struct)
        log.info(f"Running relax flow for system {name}")
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
                log.info(f"  Name: {name} | Index: {index} | UUID: {uuid}")

            job = self.store.query_one(
                criteria={"name": {"$regex": "relax"}},
                properties=["uuid", "index", "name"],
                sort={"index": -1},  # latest index first
            )
            if job is None:
                raise ValueError("No 'relax' job found")
            output = self.store.get_output(uuid=job["uuid"], which="last", load=True)
            log.info(f"Relax flow of system {name} running successfully")
            self.store.close()
            return output
        except Exception as e:
            log.error(f"Error running relax flow for system {name}: {e}")
            log.error(traceback.format_exc())
            self.store.close()
            return None


# %%
# main_relax.py

import sys
import json
import shutil
import logging
from pathlib import Path
from datetime import datetime

from custodian.vasp.handlers import VaspErrorHandler
from pymatgen.core import Structure


# from src.workflow.flow_relax import RelaxWorker

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)


class DateTimeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, datetime):
            return obj.isoformat()
        return super().default(obj)


def main():
    vasp_args = {
        "handlers": [VaspErrorHandler()],
        "vasp_cmd": "/bin/bash -c 'module load vasp-cpu && mpirun vasp_std'",
        "vasp_gamma_cmd": "/bin/bash -c 'module load vasp-cpu && mpirun vasp_gam'",
    }
    # vasp_args = {
    #     "handlers": [VaspErrorHandler()],
    #     "vasp_cmd": "/bin/bash -c 'module load vasp-gpu && mpirun vasp_std'",
    #     "vasp_gamma_cmd": "/bin/bash -c 'module load vasp-gpu && mpirun vasp_gam'",
    # }
    incar_settings = {
        "KPAR": 4,
        "NCORE": 2,
        "GGA": "PE",
    }
    potcar_functional = "PBE_64"

    # Input settings
    name = sys.argv[1]
    rlxdir = Path(sys.argv[0]).parent.joinpath("data", name, "relax")
    poscar = Path(sys.argv[0]).parent.joinpath("data", "poscars", f"{name}.vasp")
    json_path = Path(sys.argv[0]).parent.joinpath("data", name, f"{name}_relax.json")

    # Skip if already done
    if json_path.exists():
        with open(json_path, "r", encoding="utf-8") as jf:
            result = json.load(jf)
        if result["state"] == "successful":
            log.info(f"System {name} already done")
            return
    
    # Clean storage
    if rlxdir.exists():
        shutil.rmtree(rlxdir)
    rlxdir.mkdir(parents=True, exist_ok=True)

    # Run relax
    log.info(f"Running relax for system {name}")
    struct = Structure.from_file(poscar)
    worker = RelaxWorker(vasp_args, potcar_functional, incar_settings)
    rlx_data = worker.run_relax(name, struct, rlxdir)
    result = (
        {"name": name, "state": "successful", **rlx_data}
        if rlx_data
        else {"name": name, "state": "failed"}
    )
    with open(json_path, "w", encoding="utf-8") as jf:
        json.dump(result, jf, indent=2, cls=DateTimeEncoder)
    log.info(f"System {name} done")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise ValueError("Please provide a system name")
    main()

# ln -fs /home/mcmf507/documents/potpaw_PBE54 /home/mcmf507/documents/POT_GGA_PAW_PBE_54
# ln -fs /home/mcmf507/documents/potpaw_PBE54 /home/mcmf507/documents/POT_PAW_PBE_54
# echo "PMG_VASP_PSP_DIR: /home/mcmf507/documents/" > ~/.pmgrc.yaml
