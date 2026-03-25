import json
import shutil
import logging
import argparse
from datetime import datetime
from pathlib import Path

from htvasp.model import Endmember
from htvasp.workflows import QhaWorker
from htvasp.slurm import SlurmConfig, SlurmJobManager

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)


class DateTimeEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, datetime):
            return obj.isoformat()
        return super().default(obj)


vasp_args = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_gam'",
}

incar_settings = {
    "KPAR": 2,
    "NCORE": 2,
    "GGA": "PE",
}

struct_names = [
    "SER-Co",
    "SER-Fe",
    "SER-Mn",
    "SER-Ni",
    # "BCC-Co-Co",
    # "BCC-Co-Fe",
    # "BCC-Co-Mn",
    # "BCC-Co-Ni",
    # "BCC-Fe-Fe",
    # "BCC-Fe-Mn",
    # "BCC-Fe-Ni",
    # "BCC-Mn-Mn",
    # "BCC-Mn-Ni",
    # "BCC-Ni-Ni",
    # "FCC-Co-Co",
    # "FCC-Co-Fe",
    # "FCC-Co-Mn",
    # "FCC-Co-Ni",
    # "FCC-Fe-Co",
    # "FCC-Fe-Fe",
    # "FCC-Fe-Mn",
    # "FCC-Fe-Ni",
    # "FCC-Mn-Co",
    # "FCC-Mn-Fe",
    # "FCC-Mn-Mn",
    # "FCC-Mn-Ni",
    # "FCC-Ni-Co",
    # "FCC-Ni-Fe",
    # "FCC-Ni-Mn",
    # "FCC-Ni-Ni",
]


def run_locally(force=False):
    """Run all structures locally."""
    endmember = Endmember()
    for name in struct_names:
        workdir = Path("data/endmembers").joinpath(name)
        posdir = Path("data/poscars")
        flowdir = workdir.joinpath("qhaflow")
        json_path = workdir.joinpath(f"{name}-qha.json")
        workdir.mkdir(parents=True, exist_ok=True)
        posdir.mkdir(parents=True, exist_ok=True)

        # Skip if already done (unless force is True)
        if not force and json_path.exists():
            with open(json_path, "r", encoding="utf-8") as jf:
                result = json.load(jf)
            if result.get("state") == "successful":
                log.info(f"System {name} already done")
                continue

        log.info(f"Structure {name} start")

        # Run QhaWorker
        struct = endmember.get_poscar(name, posdir)
        if flowdir.exists():
            shutil.rmtree(flowdir)
        flowdir.mkdir(parents=True, exist_ok=True)
        try:
            worker = QhaWorker(
                worker_name=f"{name}-qha",
                vasp_args=vasp_args,
                incar_settings=incar_settings,
            )
            qha_data = worker.run_flow(name, struct, flowdir)
            if not qha_data:
                result = {"name": name, "state": "failed", "struct": struct.as_dict()}
            result = {"name": name, "state": "successful", **qha_data}

        except Exception as e:
            log.error(f"Structure {name} failed: {e}")
            result = {"name": name, "state": "failed", "struct": struct.as_dict()}

        # Save result
        with open(json_path, "w", encoding="utf-8") as jf:
            json.dump(result, jf, ensure_ascii=False, indent=2, cls=DateTimeEncoder)

        log.info(f"Structure {name} done")


def submit_job(force=False):
    """Submit a single job to run all structures."""
    manager = SlurmJobManager()
    config = SlurmConfig(
        job_name="gml-qha",
        output_log="em-qha.log",
        nodes=1,
        ntasks=56,
    )
    job_id = manager.submit_command(
        command=f"python {__file__} --local {'--force' if force else ''}",
        config=config,
        conda_env="htvasp",
        workdir=".",
    )
    log.info(f"Submitted job with ID: {job_id}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="run VASP task")
    parser.add_argument("--local", "-l", action="store_true", help="run locally")
    parser.add_argument("--slurm", "-s", action="store_true", help="submit to slurm")
    parser.add_argument(
        "--force", "-f", action="store_true", help="force run even if already done"
    )
    args = parser.parse_args()

    if args.local:
        run_locally(args.force)
    elif args.slurm:
        submit_job(args.force)
    else:
        print("specify run mode: --local or --slurm")
