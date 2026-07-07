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
    "NELM": 200,
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
    # VASPsol: 隐式溶剂（水）
    "LSOL": True,
    "EB_K": 78.4,          # 水的介电常数 (25°C)
    "LAMBDA_D_K": 3.0,     # Debye 屏蔽长度 (Å)，碱性电解液典型值
    "TAU": 0,              # Debye 模型
}

RELAX_INCAR = {
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "IBRION": 2,
    "ISIF": 2,
    "NELM": 100,
    "NSW": 150,
    "EDIFF": 1e-5,
    "EDIFFG": -0.05,
    "LDIPOL": True,
    "IDIPOL": 3,
}

STATIC_INCAR = {
    "ISMEAR": -5,
    "ISTART": 1,
    "IBRION": -1,
    "ALGO": "Normal",
    "NSW": 0,
    "NELM": 200,
    "LVTOT": True,
    "LWAVE": True,
    "LCHARG": True,
    "LORBIT": 11,
    "LDIPOL": True,
    "IDIPOL": 3,
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
    memory="20G",
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
    memory="100G",
    partition="partCPU",
    conda_env=conda_env,
    module_name="vasp-cpu",
)

poscar_path = root_dir / "data" / "poscars"
structs = {
    # "OH_Co_pristine": {
    #     "poscar": poscar_path / "OH_Co3_pristine.vasp",
    #     "incar": {
    #         "MAGMOM": {"Co": 3.0, "Ni": 2.0, "H": 0.6, "O": 0.6},
    #         "NELECT": 289,
    #         "LDAUL": {"Co": 2, "Ni": 2, "H": -1, "O": -1},
    #         "LDAUU": {"Co": 3.32, "Ni": 6.20},
    #         "LDAUJ": {"Co": 0.0, "Ni": 0.0},
    #     },
    # },
    # "OH_Ni_pristine": {
    #     "poscar": poscar_path / "OH_Ni6_pristine.vasp",
    #     "incar": {
    #         "MAGMOM": {"Co": 3.0, "Ni": 2.0, "H": 0.6, "O": 0.6},
    #         "NELECT": 289,
    #         "LDAUL": {"Co": 2, "Ni": 2, "H": -1, "O": -1},
    #         "LDAUU": {"Co": 3.32, "Ni": 6.20},
    #         "LDAUJ": {"Co": 0.0, "Ni": 0.0},
    #     },
    # },
    # "OH_Co_S_doped_near": {
    #     "poscar": poscar_path / "OH_Co3_S_doped.vasp",
    #     "incar": {
    #         "MAGMOM": {"Co": 3.0, "Ni": 2.0, "H": 0.6, "O": 0.6, "S": 0.6},
    #         "NELECT": 288,  # 正交 # 6*9+6*10+24*1+24*6+1*6 = 6*(9+10)+24+25*6 = 114+24+150 = 288
    #         "LDAUL": {"Co": 2, "Ni": 2, "H": -1, "O": -1, "S": -1},
    #         "LDAUU": {"Co": 3.32, "Ni": 6.20},
    #         "LDAUJ": {"Co": 0.0, "Ni": 0.0},
    #     },
    # },
    # "OH_Ni_S_doped_near": {
    #     "poscar": poscar_path / "OH_Ni6_S_doped.vasp",
    #     "incar": {
    #         "MAGMOM": {"Co": 3.0, "Ni": 2.0, "H": 0.6, "O": 0.6, "S": 0.6},
    #         "NELECT": 288,
    #         "LDAUL": {"Co": 2, "Ni": 2, "H": -1, "O": -1, "S": -1},
    #         "LDAUU": {"Co": 3.32, "Ni": 6.20},
    #         "LDAUJ": {"Co": 0.0, "Ni": 0.0},
    #     },
    # },
    # "OH_Co_S_doped_far": {
    #     "poscar": poscar_path / "OH_Co5_S_doped.vasp",
    #     "incar": {
    #         "MAGMOM": {"Co": 3.0, "Ni": 2.0, "H": 0.6, "O": 0.6, "S": 0.6},
    #         "NELECT": 288,
    #         "LDAUL": {"Co": 2, "Ni": 2, "H": -1, "O": -1, "S": -1},
    #         "LDAUU": {"Co": 3.32, "Ni": 6.20},
    #         "LDAUJ": {"Co": 0.0, "Ni": 0.0},
    #     },
    # },
    # "OH_Ni_S_doped_far": {
    #     "poscar": poscar_path / "OH_Ni1_S_doped.vasp",
    #     "incar": {
    #         "MAGMOM": {"Co": 3.0, "Ni": 2.0, "H": 0.6, "O": 0.6, "S": 0.6},
    #         "NELECT": 288,
    #         "LDAUL": {"Co": 2, "Ni": 2, "H": -1, "O": -1, "S": -1},
    #         "LDAUU": {"Co": 3.32, "Ni": 6.20},
    #         "LDAUJ": {"Co": 0.0, "Ni": 0.0},
    #     },
    # },
    # S_flip: S在顶面（与OH吸附同侧）
    "OH_Co_S_flip_near": {
        "poscar": poscar_path / "OH_Co3_S_flip.vasp",
        "incar": {
            "MAGMOM": {"Co": 3.0, "Ni": 2.0, "H": 0.6, "O": 0.6, "S": 0.6},
            "NELECT": 288,
            "LDAUL": {"Co": 2, "Ni": 2, "H": -1, "O": -1, "S": -1},
            "LDAUU": {"Co": 3.32, "Ni": 6.20},
            "LDAUJ": {"Co": 0.0, "Ni": 0.0},
        },
    },
    "OH_Ni_S_flip_near": {
        "poscar": poscar_path / "OH_Ni6_S_flip.vasp",
        "incar": {
            "MAGMOM": {"Co": 3.0, "Ni": 2.0, "H": 0.6, "O": 0.6, "S": 0.6},
            "NELECT": 288,
            "LDAUL": {"Co": 2, "Ni": 2, "H": -1, "O": -1, "S": -1},
            "LDAUU": {"Co": 3.32, "Ni": 6.20},
            "LDAUJ": {"Co": 0.0, "Ni": 0.0},
        },
    },
}


def run_static(name: str, device: str = "gpu", rerun: bool = False) -> int:
    store_dir = root_dir / "data" / name
    json_path = store_dir / "static_out.json"

    # Skip if already done
    if not rerun and json_path.exists():
        log.info(f"StaticWorker for Structure {name} already done")
        return 0

    log.info(f"StaticWorker for Structure {name} start")

    # Use relaxed structure from static_out.json if available (much faster convergence)
    old_json_path = structs[name].get("old_store", Path()) / "static_out.json"
    if old_json_path.exists():
        with open(old_json_path) as f:
            static_data = json.load(f)
        structure = Structure.from_dict(static_data["output"]["structure"])
        log.info(f"Using relaxed structure from {old_json_path}")
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
        prev_dir = None if not old_json_path.exists() else old_json_path.parent / "3-static"
        static_worker.run_flow(
            name=name,
            structure=structure,
            prev_dir=prev_dir,
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



def submit_jobs(rerun: bool = False):
    manager = SlurmJobManager()

    for name in structs:
        # submit static job — VASPsol 用 CPU（GPU 内存不足）
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="OH absorption workflow")
    parser.add_argument("--slurm", action="store_true", help="Submit to Slurm")
    parser.add_argument("--static", action="store_true", help="Run static calculation")
    parser.add_argument("--tick", type=str, help="Run tick for structure name")
    parser.add_argument("--device", type=str, default="cpu", help="Use cpu/gpu")
    parser.add_argument("--rerun", action="store_true", help="Re-run calculations")
    args = parser.parse_args()

    if args.slurm:
        submit_jobs(rerun=args.rerun)
    if args.static:
        run_static(args.tick, device=args.device, rerun=args.rerun)
    else:
        parser.print_help()
