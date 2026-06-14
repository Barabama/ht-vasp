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

# GPU VASP for relax + static (fast)
VASP_ARGS_GPU = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-gpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module load vasp-gpu && srun vasp_gam'",
}

# CPU VASP for NSCF (band/DOS) — avoid GPU memory overflow with dense k-points
VASP_ARGS_CPU = {
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
    "MAGMOM": {"Co": 3.0, "Mn": 5.0, "Ni": 2.0, "C": 0.6, "H": 0.6, "O": 0.6, "S": 0.6},
    "AMIX": 0.1,
    "BMIX": 1e-4,
    "AMIX_MAG": 0.4,
    "BMIX_MAG": 1e-4,
    "LREAL": "Auto",
    "KPAR": 4,
    "NCORE": 2,
    "GGA": "PE",
    "IVDW": 12,
    "LDAU": True,
    "LDAUTYPE": 2,
    "LDAUPRINT": 1,
    "LASPH": True,
    "LMAXMIX": 4,
    "LDAUL": {"Co": 2, "Mn": 2, "Ni": 2, "C": -1, "H": -1, "O": -1, "S": -1},
    "LDAUU": {"Co": 3.32, "Mn": 5.00, "Ni": 6.20},
    "LDAUJ": {"Co": 0.0, "Mn": 0.0, "Ni": 0.0},
}

RELAX_INCAR = {
    "ISMEAR": 0,
    "SIGMA": 0.05,
    "IBRION": 2,
    "ISIF": 2,
    "NELM": 100,
    "NSW": 50,
    "EDIFFG": -0.02,
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
}

NSCF_INCAR = {
    "ISMEAR": -5,
    "IBRION": -1,
    "NSW": 0,
    "ICHARG": 11,
    "LORBIT": 11,
}


root_dir = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S")
hetero_path = root_dir / "data" / "heterostructures"

structs = {
    "hetero_intrinsic": {
        "poscar": str(hetero_path / "hetero_intrinsic.vasp"),
        "incar": {
            "MAGMOM": {"Co": 3.0, "Mn": 5.0, "Ni": 2.0, "C": 0.6, "H": 0.6, "O": 0.6},
            "LDAUL": {"Co": 2, "Mn": 2, "Ni": 2, "C": -1, "H": -1, "O": -1},
            "LDAUU": {"Co": 3.32, "Mn": 5.00, "Ni": 6.20},
            "LDAUJ": {"Co": 0.0, "Mn": 0.0, "Ni": 0.0},
        },
    },
    "hetero_s_doped": {
        "poscar": str(hetero_path / "hetero_S_doped.vasp"),
        "incar": {
            "MAGMOM": {"Co": 3.0, "Mn": 5.0, "Ni": 2.0, "C": 0.6, "H": 0.6, "O": 0.6, "S": 0.6},
            "NELECT": 629,  # Mn_pv(13): 6*13+12*9+6*10+35*1+6*4+1*6+53*6 = 78+108+60+35+24+6+318
            "LDAUL": {"Co": 2, "Mn": 2, "Ni": 2, "C": -1, "H": -1, "O": -1, "S": -1},
            "LDAUU": {"Co": 3.32, "Mn": 5.00, "Ni": 6.20},
            "LDAUJ": {"Co": 0.0, "Mn": 0.0, "Ni": 0.0},
        },
    },
    "hetero_s_exposed": {
        "poscar": str(hetero_path / "hetero_S_exposed.vasp"),
        "incar": {
            "MAGMOM": {"Co": 3.0, "Mn": 5.0, "Ni": 2.0, "C": 0.6, "H": 0.6, "O": 0.6, "S": 0.6},
            "NELECT": 629,
            "LDAUL": {"Co": 2, "Mn": 2, "Ni": 2, "C": -1, "H": -1, "O": -1, "S": -1},
            "LDAUU": {"Co": 3.32, "Mn": 5.00, "Ni": 6.20},
            "LDAUJ": {"Co": 0.0, "Mn": 0.0, "Ni": 0.0},
        },
    },
}

# %%

from pymatgen.core import Structure


def run_tick(name: str, device: str = "cpu", rerun: bool = False):
    flow_dir = Path("/tmp")
    store_dir = root_dir / "data" / name
    json_path = store_dir / "band_out.json"

    # Skip if already done
    if not rerun and json_path.exists():
        log.info(f"Structure {name} already done")
        return

    log.info(f"Structure {name} start")

    if rerun and store_dir.exists():
        shutil.rmtree(store_dir)

    structure = Structure.from_file(structs[name]["poscar"])

    try:
        worker = NscfWorker(
            vasp_args=VASP_ARGS_CPU if device == "cpu" else VASP_ARGS_GPU,
            global_incar={**GLOBAL_INCAR, **structs[name]["incar"]},
            relax_incar=RELAX_INCAR,
            static_incar=STATIC_INCAR,
            nscf_band_incar=NSCF_INCAR,
            nscf_dos_incar=NSCF_INCAR,
            band_kpath_kwargs={"path_type": "hinuma"},
            relax_reciprocal_density=120,
            static_reciprocal_density=150,
        )
        worker.run_flow(
            name=name,
            structure=structure,
            flow_dir=flow_dir,
            store_dir=store_dir,
            resume=not rerun,
        )
        static_output = worker.get_result("static") or {}
        worker.write_result(static_output, store_dir / "static_out.json")

        dos_output = worker.get_dos_result() or {}
        worker.write_result(dos_output, store_dir / "dos_out.json")

        band_output = worker.get_band_result() or {}
        worker.write_result(band_output, store_dir / "band_out.json")

    except Exception as e:
        log.error(f"Structure {name} failed")
        log.error(e)

    log.info(f"Structure {name} done")


def submit_jobs(device: str = "cpu", rerun: bool = False):
    manager = SlurmJobManager()
    conda_env = str(root_dir.parent / ".conda")
    workdir = str(root_dir)

    for name in structs:
        job_name = f"{name}-nscf"
        command = f"python {__file__} --tick {name} --device {device} {'--rerun' if rerun else ''}"
        output_log = str(root_dir / "logs" / f"{job_name}.log")

        config = manager.get_cpu_config(
            job_name=job_name,
            output_log=output_log,
            nodes=1,
            ntasks=32,
            memory="100G",
            conda_env=conda_env,
            module_name="vasp-cpu",
        ) if device == "cpu" else manager.get_gpu_config(
            job_name=job_name,
            output_log=output_log,
            nodes=1,
            ntasks=1,
            memory="20G",
            conda_env=conda_env,
            module_name="vasp-gpu",
        )

        jid = manager.submit_command(
            command=command,
            config=config,
            workdir=workdir,
        )
        if not jid:
            log.error(f"Failed to submit nscf job for {name}")
        else:
            log.info(f"Submitted nscf job for {name} with ID {jid}")


# %%

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Heterostructure Nscf workflow")
    parser.add_argument("--slurm", action="store_true", help="Submit to Slurm")
    parser.add_argument("--tick", type=str, help="Run tick for structure name")
    parser.add_argument("--device", type=str, default="cpu", help="Use cpu/gpu")
    parser.add_argument("--rerun", action="store_true", help="Re-run calculations")
    args = parser.parse_args()

    if args.slurm:
        submit_jobs(device=args.device, rerun=args.rerun)
    elif args.tick:
        run_tick(args.tick, device=args.device, rerun=args.rerun)
    else:
        parser.print_help()
