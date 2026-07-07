import json
import time
import shutil
import logging
import argparse
from pathlib import Path
from pymatgen.core import Structure
from htvasp.workflows import StaticWorker, NscfWorker
from htvasp.slurm import SlurmJobManager, SlurmConfig

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

# GPU VASP for relax + static (fast)
VASP_ARGS_GPU = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-gpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-gpu && srun vasp_gam'",
}

# CPU VASP for NSCF (band/DOS) — avoid GPU memory overflow with dense k-points
VASP_ARGS_CPU = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-cpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-cpu && srun vasp_gam'",
}

GLOBAL_INCAR = {
    "ENCUT": 520,
    "PREC": "Accurate",
    "ALGO": "Fast",
    "NELM": 100,
    "EDIFF": 1e-6,
    "ISPIN": 2,
    "MAGMOM": {"Co": 3.0, "Mn": 5.0, "Ni": 2.0, "C": 0.6, "H": 0.6, "O": 0.6, "S": 0.6},
    "AMIX": 0.1,
    "BMIX": 1e-4,
    "AMIX_MAG": 0.4,
    "BMIX_MAG": 1e-4,
    "LREAL": "Auto",
    # "KPAR": 4,
    # "NCORE": 4,
    "GGA": "PE",
    "IVDW": 12,
    "LDAU": True,
    "LDAUTYPE": 2,
    "LDAUPRINT": 1,
    "LASPH": True,
    "LMAXMIX": 4,
    "LDAUL": {"Co": 2, "Mn": 2, "Ni": 2, "C": -1, "H": -1, "O": -1, "S": -1},
    "LDAUU": {"Co": 3.32, "Mn": 3.90, "Ni": 6.20},
    "LDAUJ": {"Co": 0.0, "Mn": 0.0, "Ni": 0.0},
}

RELAX_INCAR = {
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "IBRION": 2,
    "ISIF": 3,
    "NELM": 100,
    "NSW": 150,
    "EDIFF": 1e-5,
    "EDIFFG": -0.05,
}

STATIC_INCAR = {
    "ISMEAR": -5,
    "ISTART": 1,
    "IBRION": -1,
    "ALGO": "Normal",
    "NSW": 0,
    "NELM": 200,
    "LWAVE": True,
    "LCHARG": True,
    "LORBIT": 11,
}

NSCF_INCAR = {
    "ISMEAR": -5,
    "IBRION": -1,
    "NSW": 0,
    "ICHARG": 11,
    "LORBIT": 11,
    "KPAR": 2,
    "NCORE": 4,
}


root_dir = Path(__file__).parent.resolve()
flow_dir = Path("/nfs_ssd/tmp")
conda_env = str(root_dir.parent / ".conda")
gpu_config = SlurmConfig(
    job_name="gpu_job",
    output_log="gpu_job.log",
    nodes=1,
    ntasks=1,
    memory="10G",
    partition="partGPU",
    gpus_per_task=1,
    conda_env=conda_env,
    module_name="vasp-gpu",
)
cpu_config = SlurmConfig(
    job_name="cpu_job",
    output_log="cpu_job.log",
    nodes=1,
    ntasks=32,
    memory="150G",
    partition="partCPU",
    conda_env=conda_env,
    module_name="vasp-cpu",
)

poscar_path = root_dir / "data" / "poscars"
structs = {
    "CMCH_bulk": {
        "old_store": root_dir / "data" / "CoMnH2CO5",
        "poscar": poscar_path / "CoMnH2CO5-311.vasp",
        "incar": {
            "MAGMOM": {"Co": 3.0, "Mn": 5.0, "C": 0.6, "H": 0.6, "O": 0.6},
            # "NELECT": 540 , # 12*9+12*7+24*1+12*4+60*6= 540
            "LDAUL": {"Co": 2, "Mn": 2, "C": -1, "H": -1, "O": -1},
            "LDAUU": {"Co": 3.32, "Mn": 3.90},
            "LDAUJ": {"Co": 0.0, "Mn": 0.0},
            "NBANDS": 600,
        },
    },
    # "LDH_bulk": {
    #     "old_store": root_dir / "data" / "CoNiOH2",
    #     "poscar": poscar_path / "CoNiOH2-322.vasp",
    #     "incar": {
    #         "MAGMOM": {"Co": 3.0, "Ni": 2.0, "H": 0.6, "O": 0.6},
    #         "LDAUL": {"Co": 2, "Ni": 2, "H": -1, "O": -1},
    #         "LDAUU": {"Co": 3.32, "Ni": 6.20},
    #         "LDAUJ": {"Co": 0.0, "Ni": 0.0},
    #     },
    # },
    # "LDH_S_bulk": {
    #     "old_store": root_dir / "data" / "CoNiOH2S-noH",
    #     "poscar": poscar_path / "CoNiOH2-322-S-noH.vasp",
    #     "incar": {
    #         "MAGMOM": {"Co": 3.0, "Ni": 2.0, "H": 0.6, "O": 0.6, "S": 0.6},
    #         # "NELECT": 422, # 六方 # 9*9 + 9*10 + 35*6 + 35*1 + 1*6
    #         "NELECT": 563,  # 正交 # 12*9+12*10+47*1+47*6+1*6 = 12*(9+10)+47+48*6 = 228+47+288=563
    #         "LDAUL": {"Co": 2, "Ni": 2, "H": -1, "O": -1, "S": -1},
    #         "LDAUU": {"Co": 3.32, "Ni": 6.20},
    #         "LDAUJ": {"Co": 0.0, "Ni": 0.0},
    #     },
    # },
}


def run_static(name: str, device: str = "gpu", rerun: bool = False) -> int:
    store_dir = root_dir / "data" / name
    json_path = store_dir / "static_out.json"

    # # Skip if already done
    # if not rerun and json_path.exists():
    #     log.info(f"StaticWorker for Structure {name} already done")
    #     return 0

    log.info(f"StaticWorker for Structure {name} start")

    # Use relaxed structure from static_out.json if available (much faster convergence)
    json_path = store_dir / "static_out.json"
    if json_path.exists():
        with open(json_path) as f:
            static_data = json.load(f)
        structure = Structure.from_dict(static_data["output"]["structure"])
        log.info(f"Using relaxed structure from {json_path}")
    else:
        structure = Structure.from_file(structs[name]["poscar"])
        log.info(f"Using initial structure from POSCAR")
    try:
        static_worker = StaticWorker(
            vasp_args=VASP_ARGS_CPU if device == "cpu" else VASP_ARGS_GPU,
            global_incar={**GLOBAL_INCAR, **structs[name]["incar"]},
            relax_incar=RELAX_INCAR,
            static_incar=STATIC_INCAR,
        )
        static_worker.run_flow(
            name=name,
            structure=structure,
            # prev_dir=store_dir / "3-static",
            flow_dir=flow_dir,
            store_dir=store_dir,
            resume=not rerun,
        )
        static_output = static_worker.get_result("static") or {}
        static_worker.write_result(static_output, store_dir / "static_out.json")
    except Exception as e:
        log.error(f"StaticWorker for Structure {name} failed")
        log.error(e)
        return 1

    log.info(f"StaticWorker for Structure {name} done")
    return 0


def run_nscf(name: str, device: str = "cpu", rerun: bool = False) -> int:
    store_dir = root_dir / "data" / name
    json_path = store_dir / "band_out.json"

    # Skip if already done
    if not rerun and json_path.exists():
        log.info(f"NscfWorker for Structure {name} already done")
        return 0

    log.info(f"NscfWorker for Structure {name} start")

    with open(store_dir / "static_out.json", "r", encoding="utf-8") as f:
        static_output = json.load(f)

    structure = Structure.from_dict(static_output["output"]["structure"])

    try:
        worker = NscfWorker(
            vasp_args=VASP_ARGS_CPU if device == "cpu" else VASP_ARGS_GPU,
            global_incar={**GLOBAL_INCAR, **structs[name]["incar"]},
            nscf_band_incar=NSCF_INCAR,
            nscf_dos_incar=NSCF_INCAR,
        )
        worker.run_flow(
            name=name,
            structure=structure,
            prev_dir=store_dir / "3-static",
            flow_dir=flow_dir,
            store_dir=store_dir,
            resume=not rerun,
        )
        dos_output = worker.get_dos_result() or {}
        worker.write_result(dos_output, store_dir / "dos_out.json")

        band_output = worker.get_band_result() or {}
        worker.write_result(band_output, store_dir / "band_out.json")

    except Exception as e:
        log.error(f"NscfWorker for Structure {name} failed")
        log.error(e)
        return 1

    log.info(f"NscfWorker for Structure {name} done")
    return 0


def submit_jobs(rerun: bool = False):
    manager = SlurmJobManager()

    for name in structs:
        # submit static job
        static_job_name = f"{name}-static"
        static_cmd = f"python {__file__} --tick {name} --static --device gpu {'--rerun' if rerun else ''}"
        static_log = str(root_dir / "logs" / f"{static_job_name}.log")
        gpu_config.job_name = static_job_name
        gpu_config.output_log = static_log
        static_jid = manager.submit_command(
            command=static_cmd,
            config=gpu_config,
            work_dir=root_dir,
        )
        time.sleep(1)
        if not static_jid:
            log.error(f"Failed to submit static job for {name}")
            continue
        log.info(f"Submitted static job for {name} with ID {static_jid}")

        # submit NSCF job dependent on static job
        nscf_job_name = f"{name}-nscf"
        nscf_cmd = f"python {__file__} --tick {name} --nscf --device cpu {'--rerun' if rerun else ''}"
        nscf_log = str(root_dir / "logs" / f"{nscf_job_name}.log")
        cpu_config.job_name = nscf_job_name
        cpu_config.output_log = nscf_log
        cpu_config.dependency = f"afterok:{static_jid}"
        nscf_jid = manager.submit_command(
            command=nscf_cmd,
            config=cpu_config,
            work_dir=root_dir,
        )
        time.sleep(1)
        if not nscf_jid:
            log.error(f"Failed to submit NSCF job for {name}")
            continue
        log.info(f"Submitted NSCF job for {name} with ID {nscf_jid} depending on {static_jid}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Nscf workflow")
    parser.add_argument("--slurm", action="store_true", help="Submit to Slurm")
    parser.add_argument("--static", action="store_true", help="Run static calculation")
    parser.add_argument("--nscf", action="store_true", help="Run NSCF calculation")
    parser.add_argument("--tick", type=str, help="Run tick for structure name")
    parser.add_argument("--device", type=str, default="cpu", help="Use cpu/gpu")
    parser.add_argument("--rerun", action="store_true", help="Re-run calculations")
    args = parser.parse_args()

    if args.slurm:
        submit_jobs(rerun=args.rerun)
    if args.static:
        run_static(args.tick, device=args.device, rerun=args.rerun)
    elif args.nscf:
        run_nscf(args.tick, device=args.device, rerun=args.rerun)
    else:
        parser.print_help()
