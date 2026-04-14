import json
import logging
import argparse
from pathlib import Path
from datetime import datetime

from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar

from sqsgenerator import StructureFormat, parse_config, optimize
from sqsgenerator.core import SqsConfiguration

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)


def run_locally():
    struct = Structure(
        lattice=[[3, 0, 0], [0, 3, 0], [0, 0, 3]],
        species=["Au", "Cu"],
        coords=[[0, 0, 0], [0.5, 0.5, 0.5]],
    )

    json_info = {
        "iterations": 10000,
        "composition": {"Co": 8, "Fe": 8},
        "structure": {
            "lattice": struct.lattice.matrix.tolist(),
            "coords": struct.frac_coords.tolist(),
            "species": [str(s) for s in struct.species],
            "supercell": [2, 2, 2],
        },
    }

    print(json_info)
    sqscfg = parse_config(dict(json_info))
    print(sqscfg)
    if not isinstance(sqscfg, SqsConfiguration):
        raise ValueError("sqscfg is not a SqsConfiguration")
    print("sqs generating")
    pack = optimize(sqscfg)
    if len(pack) < 1:
        raise ValueError("pack result empty")

    best = pack.best()
    sqs_struct = best.structure()
    poscar_str = sqs_struct.dump(StructureFormat.poscar)
    print(poscar_str)
    new_struct = Poscar.from_str(poscar_str).structure


def submit_job():
    from htvasp.slurm import SlurmJobManager

    manager = SlurmJobManager()
    config = manager.get_cpu_config(
        ntasks=8,
        memory="4G",
    )
    job_id = manager.submit_command(
        command=f"python {__file__} --local",
        config=config,
        conda_env="htvasp",
        workdir=".",
    )
    log.info(f"Submitted job: {job_id}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test OJ workflow")
    parser.add_argument("--local", action="store_true", help="run workflow locally")
    parser.add_argument("--slurm", action="store_true", help="submit to slurm")
    parser.add_argument("--unit", action="store_true", help="run unit tests only")
    args = parser.parse_args()

    if args.local:
        run_locally()
    elif args.slurm:
        submit_job()
    elif args.unit:
        success = run_unit_tests()
        exit(0 if success else 1)
    else:
        success = run_unit_tests()
        exit(0 if success else 1)
