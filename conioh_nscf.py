# %%

import json
import shutil
import logging
import argparse
from pathlib import Path


from htvasp.workflows import NscfWorker
from htvasp.slurm import SlurmJobManager


logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

VASP_ARGS = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-cpu && srun vasp_gam'",
}


GLOBAL_INCAR = {
    "ENCUT": 520,
    "PREC": "Accurate",
    "ALGO": "Fast",
    "NELM": 100,
    "EDIFF": 1e-6,
    "ISPIN": 2,
    "MAGMOM": {"Co": 5.0, "Ni": 2.0, "H": 0.6, "O": 0.6, "S": 0.6},
    "AMIX": 0.1,
    "BMIX": 1e-4,
    "AMIX_MAG": 0.4,
    "BMIX_MAG": 1e-4,
    "LREAL": "Auto",
    "KPAR": 4,
    "NCORE": 4,
    "GGA": "PE",
    "IVDW": 12,  # ← 新增
    "LDAU": True,
    "LDAUTYPE": 2,
    "LDAUPRINT": 1,
    "LASPH": True,
    "LMAXMIX": 4,
    "LDAUL": {"Co": 2, "Ni": 2, "H": -1, "O": -1, "S": -1},
    "LDAUU": {"Co": 3.32, "Ni": 6.20},
    "LDAUJ": {"Co": 0.0, "Ni": 0.0},
}

RELAX_INCAR = {
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "IBRION": 2,
    "ISIF": 3,
    "NSW": 50,
    "EDIFFG": -0.02,
}

STATIC_INCAR = {
    "ISMEAR": -5,
    "ISTART": 1,
    "IBRION": -1,
    "NSW": 0,
    "NELM": 200,
    "LWAVE": True,
    "LCHARG": True,
}

NSCF_INCAR = {
    "ISMEAR": -5,
    # "ISMEAR": 0,
    "IBRION": -1,
    "NSW": 0,
    "ICHARG": 11,
    "LORBIT": 11,
}

poscars_path = Path("data/poscars")
structs = {
    "CoNiHO": "data/poscars/CoNiHO.vasp",
    "CoNiHOS-Co": "data/poscars/CoNiHOS-Co.vasp",
    "CoNiHOS-Co-H": "data/poscars/CoNiHOS-Co-H.vasp",
    "CoNiHOS-Ni": "data/poscars/CoNiHOS-Ni.vasp",
    "CoNiHOS-Ni-H": "data/poscars/CoNiHOS-Ni-H.vasp",
}


# %%

from pymatgen.core import Structure


def run_tick(name: str, force: bool = False):
    flow_dir = Path("/nfs_ssd/tmp")
    store_dir = Path("data") / name
    json_path = store_dir / f"{name}.json"

    # Skip if already done
    if not force and json_path.exists():
        log.info(f"Structure {name} already done")
        return

    log.info(f"Structure {name} start")

    if force and store_dir.exists():
        shutil.rmtree(store_dir)

    structure = Structure.from_file(structs[name])

    try:
        worker = NscfWorker(
            vasp_args=VASP_ARGS,
            global_incar=GLOBAL_INCAR,
            relax_incar=RELAX_INCAR,
            static_incar=STATIC_INCAR,
            nscf_band_incar=NSCF_INCAR,
            nscf_dos_incar=NSCF_INCAR,
        )
        worker.run_flow(name, structure, flow_dir, store_dir, resume=not force)
        dos_output = worker.get_dos_result() or {}
        log.info(f"dos_output has keys: {list(dos_output.keys())}")
        worker.write_result(dos_output, store_dir / "dos_out.json")

        band_output = worker.get_band_result() or {}
        log.info(f"band_output has keys: {list(band_output.keys())}")
        worker.write_result(band_output, store_dir / "band_out.json")

    except Exception as e:
        log.error(f"Structure {name} failed: {e}")

    log.info(f"Structure {name} done")


def run_batch(force: bool = False):
    """Run all structures locally."""
    for name in structs:
        run_tick(name, force=force)


def submit_jobs(force: bool = False):
    manager = SlurmJobManager()
    for name in structs:
        config = manager.get_cpu_config(
            job_name=f"{name}-nscf",
            output_log=f"logs/{name}-nscf.log",
            nodelist="429pro",
            ntasks=48,
            memory="100G",
        )
        job_id = manager.submit_command(
            command=f"python {__file__} --tick {name} {'--force' if force else ''}",
            config=config,
            conda_env="htvasp",
            workdir=".",
        )
        if not job_id:
            log.error(f"Failed to submit job for {name}")
        else:
            log.info(f"Submitted job for {name} with ID {job_id}")


# %%

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Magnetic workflow")
    parser.add_argument("--tick", type=str, help="Run single structure")
    parser.add_argument("--batch", action="store_true", help="Run all locally")
    parser.add_argument("--slurm", action="store_true", help="Submit to Slurm")
    parser.add_argument("--force", action="store_true", help="Force re-run")
    args = parser.parse_args()

    if args.batch:
        run_batch(force=args.force)
    elif args.slurm:
        submit_jobs(force=args.force)
    elif args.tick:
        run_tick(args.tick, force=args.force)
    else:
        parser.print_help()
