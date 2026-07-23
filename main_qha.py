import json
import time
import shutil
import logging
import argparse
from pathlib import Path
from pymatgen.core import Structure
from htvasp.model import Endmember
from htvasp.workflows import QhaWorker, StaticWorker
from htvasp.slurm import SlurmJobManager

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

VASP_ARGS = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-cpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-cpu && srun vasp_gam'",
}
# {"Al": 1.0, "Co": 3.0, "Cr": 3.0, "Cu": 5.0, "Fe": 5.0, "Mn": 5.0, "Nb": 3.0, "Ni": 2.0, "Ta": 3.0, "Ti": 3.0, "V": 3.0, "W": 3.0, "Zr": 1.0}
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
    "MAGMOM": {
        "Al": 1.0,
        "Co": 3.0,
        "Fe": 5.0,
        "Mn": 5.0,
        "Nb": 3.0,
        "Ni": 2.0,
        "Ti": 3.0,
        "V": 3.0,
    },
    # Precision
    "ISYM": 0,
    "LREAL": "Auto",
    "PREC": "Normal",
    "SYMPREC": 1e-7,
    # Output
    "LWAVE": False,
    "LCHARG": False,
    "LOPTICS": False,
    "LVTOT": False,
    "GGA": "PE",
    "AMIX": 0.1,
    "BMIX": 1e-4,
    "AMIX_MAG": 0.4,
    "BMIX_MAG": 1e-4,
}

RELAX_INCAR = {"KPAR": 4, "NCORE": 2, "LORBIT": 10}
PHONON_INCAR = {"KPAR": 2, "NCORE": 4}
STATIC_INCAR = {"LORBIT": 11}

STRUCTURE_NAMES = [
    "SER-Al",
    "SER-Nb",
    "SER-Ti",
    "SER-V",
    "BCC-Al-Al",
    "BCC-Al-Nb",
    "BCC-Al-Ti",
    "BCC-Al-V",
    "BCC-Nb-Nb",
    "BCC-Nb-Ti",
    "BCC-Nb-V",
    "BCC-Ti-Ti",
    "BCC-Ti-V",
    "BCC-V-V",
    #
    # "SER-Co",
    # "SER-Fe",
    # "SER-Mn",
    # "SER-Ni",
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


root_dir = Path(__file__).parent.resolve()
flow_dir = Path("/nfs_ssd/tmp")
conda_env = root_dir / ".conda"


def run_static(name: str, rerun: bool = False):
    static_store = root_dir / "data" / "endmembers" / name / "staticflow"
    static_json = static_store.parent / f"{name}-static.json"

    if not rerun and static_json.exists():
        log.info(f"StaticWorker for Structure {name} already done")
        return

    log.info(f"StaticWorker for Structure {name} started")

    endmember = Endmember()
    structure = endmember.get_poscar(name, Path(root_dir) / "data" / "poscars")
    try:
        static_worker = StaticWorker(
            vasp_args=VASP_ARGS,
            global_incar=GLOBAL_INCAR,
            relax_incar=RELAX_INCAR,
            static_incar=STATIC_INCAR,
        )
        static_worker.run_flow(
            name=name,
            structure=structure,
            flow_dir=flow_dir,
            store_dir=static_store,
            resume=not rerun,
        )
        static_output = static_worker.get_result() or {}
        static_worker.write_result(data=static_output, json_path=static_json)
    except Exception as e:
        log.error(f"StaticWorker for Structure {name} failed")
        log.error(e)

    log.info(f"StaticWorker for Structure {name} done")
    return


def run_qha(name: str, rerun: bool = False):
    qha_store = root_dir / "data" / "endmembers" / name / "qhaflow"
    qha_json = qha_store.parent / f"{name}-qha.json"

    if not rerun and qha_json.exists():
        log.info(f"QhaWorker for Structure {name} already done")
        return

    log.info(f"QhaWorker for Structure {name} started")

    endmember = Endmember()
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
            store_dir=qha_store,
            resume=not rerun,
        )
        qha_output = qha_worker.get_result() or {}
        qha_worker.write_result(data=qha_output, json_path=qha_json)
    except Exception as e:
        log.error(f"QhaWorker for Structure {name} failed")
        log.error(e)

    log.info(f"QhaWorker for Structure {name} done")
    return


def check_jobs(job: str):

    for name in STRUCTURE_NAMES:
        store_dir = root_dir / "data" / "endmembers" / name / f"{job}flow"
        json_path = store_dir.parent / f"{name}-{job}.json"
        if not json_path.exists:
            continue
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        total_magnetization = data["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"]
        magnetization = data["calcs_reversed"][0]["output"]["outcar"]["magnetization"]
        log.info(f"{name}: {total_magnetization}")  # 总磁矩
        # log.info(f"{magnetization}\n")              # 轨道分解 [s,p,d,tot]



def submit_job(name: str, job: str, manager: SlurmJobManager, rerun: bool = False):
    jid = manager.submit_command(
        command=f"python {__file__} -n {name} -j {job} {'-r' if rerun else ''}",
        config=manager.get_cpu_config(
            job_name=f"{name}-{job}",
            output_log=f"logs/{name}-{job}.log",
            ntasks=8,
            memory="16G",
            nodelist="429pro",
        ),
        conda_env=str(conda_env),
        work_dir=root_dir,
    )
    time.sleep(1)
    if not jid:
        log.error(f"Failed to submit qha job for Structure {name}")
        return -1
    log.info(f"Submitted qha job for Structure {name} with ID {jid}")
    return jid


def submit_jobs(job: str, rerun: bool = False):
    manager = SlurmJobManager()
    for name in STRUCTURE_NAMES:
        submit_job(name=name, job=job, manager=manager, rerun=rerun)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Qha workflow")
    parser.add_argument("--name", "-n", type=str, help="Run for structure name")
    parser.add_argument("--job", "-j", default="static", type=str, help="Type of job name")
    parser.add_argument("--check", "-c", action="store_true", help="Check job json")
    parser.add_argument("--slurm", "-s", action="store_true", help="Submit to Slurm")
    parser.add_argument("--rerun", "-r", action="store_true", help="Re-run calculations")
    args = parser.parse_args()

    if args.slurm:
        submit_jobs(args.job, rerun=args.rerun)
    elif args.check:
        check_jobs(args.job)
    elif args.name:
        job = args.job.lower()
        if job.startswith("q"):
            run_qha(args.name, rerun=args.rerun)
        elif job.startswith("s"):
            run_static(args.name, rerun=args.rerun)
    else:
        parser.print_help()
