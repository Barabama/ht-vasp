"""
Deformation charge density
"""

import logging
from pathlib import Path
from typing import Any

from jobflow import Flow
from pymatgen.core import Structure
from atomate2.vasp.jobs.core import StaticMaker
from atomate2.vasp.sets.core import StaticSetGenerator

from htvasp.workflows.base import Worker

log = logging.getLogger(__name__)


class StaticWorker(Worker):
    """
    Worker for structural static (ISIF=2).
    """

    def __init__(
        self,
        worker_name: str = "Static-worker",
        vasp_args: dict[str, Any] | None = None,
        potcar_functional="PBE_64",
        global_incar: dict[str, Any] | None = None,
        static_incar: dict[str, Any] | None = None,
        **kwargs,
    ):
        # Initialize base Worker
        super().__init__(
            worker_name=worker_name,
            vasp_args=vasp_args,
            potcar_functional=potcar_functional,
            global_incar=global_incar,
        )
        static_incar = static_incar or {}

        # Static structural relaxation
        static_maker = StaticMaker(
            name="static",
            run_vasp_kwargs=self.run_vasp_kwargs,
            stop_children_kwargs={"handle_unsuccessful": False},
            copy_vasp_kwargs={"additional_vasp_files": ("WAVECAR",)},
            input_set_generator=StaticSetGenerator(
                user_potcar_functional=self.potcar_functional,
                user_incar_settings={
                    **self.global_incar,
                    "ISTART": 1,
                    "IBRION": -1,
                    "ISIF": 2,
                    "NSW": 0,
                    "LWAVE": False,
                    "LCHARG": True,
                    **static_incar,
                },
            ),
        )
        self.flow_maker = static_maker

    def _make_flow(self, structure: Structure, prev_dir: Path | str | None = None) -> Flow:
        return self.flow_maker.make(structure, prev_dir)

    def get_result(self, output_job_name: str = "static") -> dict[str, Any] | None:
        return super().get_result(output_job_name)


# %%
import json
import logging
import argparse
from pymatgen.core import Structure
from htvasp.slurm import SlurmJobManager, SlurmConfig

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")


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
    memory="120G",
    partition="partCPU",
    conda_env=conda_env,
    module_name="vasp-cpu",
)

poscar_path = root_dir / "data" / "poscars"
structs = {
    # "CMCH_intrinsic": {
    #     "poscar": poscar_path / "CMCH_strained.vasp",
    #     "incar": {
    #         "MAGMOM": {"Co": 3.0, "Mn": 5.0, "C": 0.6, "H": 0.6, "O": 0.6},
    #         "LDAUL": {"Co": 2, "Mn": 2, "C": -1, "H": -1, "O": -1},
    #         "LDAUU": {"Co": 3.32, "Mn": 3.90},
    #         "LDAUJ": {"Co": 0.0, "Mn": 0.0},
    #     },
    # },
    # "CMCH_S_doped": {
    #     "poscar": poscar_path / "CMCH_S_doped.vasp",
    #     "incar": {
    #         "MAGMOM": {"Co": 3.0, "Mn": 5.0, "C": 0.6, "H": 0.6, "O": 0.6},
    #         "LDAUL": {"Co": 2, "Mn": 2, "C": -1, "H": -1, "O": -1},
    #         "LDAUU": {"Co": 3.32, "Mn": 3.90},
    #         "LDAUJ": {"Co": 0.0, "Mn": 0.0},
    #     },
    # },
    "CMCH_S_exposed": {
        "poscar": poscar_path / "CMCH_S_exposed.vasp",
        "incar": {
            "MAGMOM": {"Co": 3.0, "Mn": 5.0, "C": 0.6, "H": 0.6, "O": 0.6},
            "LDAUL": {"Co": 2, "Mn": 2, "C": -1, "H": -1, "O": -1},
            "LDAUU": {"Co": 3.32, "Mn": 3.90},
            "LDAUJ": {"Co": 0.0, "Mn": 0.0},
        },
    },
    "LDH_intrinsic": {
        "poscar": poscar_path / "LDH_intrinsic.vasp",
        "incar": {
            "MAGMOM": {"Co": 3.0, "Ni": 2.0, "H": 0.6, "O": 0.6},
            "LDAUL": {"Co": 2, "Ni": 2, "H": -1, "O": -1},
            "LDAUU": {"Co": 3.32, "Ni": 6.20},
            "LDAUJ": {"Co": 0.0, "Ni": 0.0},
        },
    },
    "LDH_S_doped": {
        "poscar": poscar_path / "LDH_S_doped.vasp",
        "incar": {
            "MAGMOM": {"Co": 3.0, "Ni": 2.0, "H": 0.6, "O": 0.6, "S": 0.6},
            "NELECT": 281,  # 正交 # 6*9+6*10+23*1+23*6+1*6 = 6*(9+10)+23+24*6 = 114+23+144=281
            "LDAUL": {"Co": 2, "Ni": 2, "H": -1, "O": -1, "S": -1},
            "LDAUU": {"Co": 3.32, "Ni": 6.20},
            "LDAUJ": {"Co": 0.0, "Ni": 0.0},
        },
    },
    "LDH_S_exposed": {
        "poscar": poscar_path / "LDH_S_exposed.vasp",
        "incar": {
            "MAGMOM": {"Co": 3.0, "Ni": 2.0, "H": 0.6, "O": 0.6, "S": 0.6},
            "NELECT": 281,
            "LDAUL": {"Co": 2, "Ni": 2, "H": -1, "O": -1, "S": -1},
            "LDAUU": {"Co": 3.32, "Ni": 6.20},
            "LDAUJ": {"Co": 0.0, "Ni": 0.0},
        },
    },
}


def run_static(name: str, device: str = "gpu", rerun: bool = False):
    store_dir = root_dir / "data" / name
    json_path = store_dir / "static_out.json"

    # Skip if already done
    if not rerun and json_path.exists():
        log.info(f"StaticWorker for Structure {name} already done")
    else:
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

        log.info(f"StaticWorker for Structure {name} done")


def submit_jobs(rerun: bool = False):
    manager = SlurmJobManager()

    for name in structs:
        # submit static job
        static_job_name = f"{name}-static"
        static_cmd = (
            f"python {__file__} --static {name} --device gpu {'--rerun' if rerun else ''}"
        )
        static_log = str(root_dir / "logs" / f"{static_job_name}.log")
        gpu_config.job_name = static_job_name
        gpu_config.output_log = static_log
        static_jid = manager.submit_command(
            command=static_cmd,
            config=gpu_config,
            work_dir=root_dir,
        )
        if not static_jid:
            log.error(f"Failed to submit static job for {name}")
            continue
        log.info(f"Submitted static job for {name} with ID {static_jid}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Slab Nscf workflow")
    parser.add_argument("--slurm", action="store_true", help="Submit to Slurm")
    parser.add_argument("--static", type=str, help="Run static calculation")
    parser.add_argument("--device", type=str, default="cpu", help="Use cpu/gpu")
    parser.add_argument("--rerun", action="store_true", help="Re-run calculations")
    args = parser.parse_args()

    if args.slurm:
        submit_jobs(rerun=args.rerun)
    elif args.static:
        run_static(args.static, device=args.device, rerun=args.rerun)
    else:
        parser.print_help()

