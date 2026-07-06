import json
import time
import shutil
import logging
import argparse
from pathlib import Path
from pymatgen.core import Structure
from htvasp.model import Endmember
from htvasp.workflows import QhaWorker
from htvasp.slurm import SlurmJobManager

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

VASP_ARGS = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-cpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-cpu && srun vasp_gam'",
}

GLOBAL_INCAR = {
    "ENCUT": 400,
    # Electronic
    "ISMEAR": 1,
    "SIGMA": 0.2,
    "ALGO": "Normal",
    "NELM": 120,
    "NELMIN": 6,
    "NELMDL": -6,
    # Ionic
    "IBRION": 2,
    "ISIF": 3,
    "NSW": 80,
    "POTIM": 0.2,
    "EDIFF": 1e-6,
    "EDIFFG": -0.02,
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
    "LORBIT": None,
    "LOPTICS": False,
    "LVTOT": False,
    "GGA": "PE",
    "AMIX": 0.1,
    "BMIX": 1e-4,
    "AMIX_MAG": 0.4,
    "BMIX_MAG": 1e-4,
}

RELAX_INCAR = {"KPAR": 1, "NCORE": 1}
PHONON_INCAR = {"KPAR": 1, "NCORE": 1}

STRUCTURE_NAMES = [
    # "SER-Al",
    # "SER-Nb",
    # "SER-Ti",
    # "SER-V",
    # "BCC-Al-Al",
    # "BCC-Al-Nb",
    # "BCC-Al-Ti",
    # "BCC-Al-V",
    # "BCC-Nb-Nb",
    # "BCC-Nb-Ti",
    # "BCC-Nb-V",
    # "BCC-Ti-Ti",
    "BCC-Ti-V",
    # "BCC-V-V",
]


root_dir = Path(__file__).parent.resolve()
flow_dir = Path("/nfs_ssd/tmp")
conda_env = root_dir.parent / ".conda"


def run_qha(name: str, rerun: bool = False):
    endmember = Endmember()
    store_dir = root_dir / "data" / name
    json_path = store_dir / f"{name}-qha.json"

    # # Skip if already done
    # if not rerun and json_path.exists():
    #     log.info(f"StaticWorker for Structure {name} already done")
    #     return 0

    log.info(f"QhaWorker for Structure {name} started")

    structure = endmember.get_poscar(name, Path(root_dir) / "data" / "poscars")
    
    try:
        qha_worker = QhaWorker(
            vasp_args=VASP_ARGS,
            global_incar=GLOBAL_INCAR,
            relax_incar=RELAX_INCAR,
            eos_incar=RELAX_INCAR,
            phonon_incar=PHONON_INCAR,
            temperature_range=(0, 3100, 50),
            supercell_matrix=((2, 0, 0), (0, 2, 0), (0, 0, 2)),
        )
        qha_worker.run_flow(
            name=name,
            structure=structure,
            flow_dir=flow_dir,
            store_dir=store_dir,
            resume=not rerun,
        )
        qha_output = qha_worker.get_result() or {}
        qha_worker.write_result(data=qha_output, json_path=json_path)
    except Exception as e:
        log.error(f"QhaWorker for Structure {name} failed")
        log.error(e)
        return 1
    
    log.info(f"QhaWorker for Structure {name} done")
    return 0


def submit_jobs(force: bool = False):
    manager = SlurmJobManager()
    for name in STRUCTURE_NAMES:
        jid = manager.submit_command(
            command=f"python {__file__} --tick {name} {'--rerun' if force else ''}",
            config=manager.get_cpu_config(
                job_name=f"{name}-qha",
                output_log=f"logs/{name}-qha.log",
                ntasks=32,
                memory="64G",
            ),
            conda_env=str(conda_env),
            work_dir=root_dir,
        )
        time.sleep(1)
        if not jid:
            log.error(f"Failed to submit qha job for Structure {name}")
            continue
        log.info(f"Submitted qha job for Structure {name} with ID {jid}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Qha workflow")
    parser.add_argument("--slurm", action="store_true", help="Submit to Slurm")
    parser.add_argument("--tick", type=str, help="Run qha for structure name")
    parser.add_argument("--rerun", action="store_true", help="Re-run calculations")
    args = parser.parse_args()

    if args.slurm:
        submit_jobs(force=args.rerun)
    elif args.tick:
        run_qha(args.tick, rerun=args.rerun)
    else:
        parser.print_help()
