import json
import shutil
import logging
import argparse
from pathlib import Path
from datetime import datetime

from htvasp.model import Endmember
from htvasp.workflows import QhaWorker
from htvasp.slurm import SlurmJobManager

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

global_incar = {
            "ENCUT": 450,
            "ISTART": 0,
            "ICHARG": 2,
            # Electronic
            "ISMEAR": 1,
            "SIGMA": 0.1,
            "ALGO": "Normal",
            "NELM": 200,
            "NELMIN": 6,
            "NELMDL": -6,
            # Ionic
            "IBRION": 2,
            "ISIF": 2,
            "NSW": 100,
            "POTIM": 0.2,
            "EDIFF": 1e-6,
            "EDIFFG": -0.01,
            # Magnetic
            "ISPIN": 2,
            # Precision
            "ISYM": 0,
            "LREAL": "Auto",
            "PREC": "Accurate",
            "SYMPREC": 1e-5,
            # Output
            "LWAVE": False,
            "LCHARG": False,
            "LORBIT": None,
            "LOPTICS": False,
            "LVTOT": False,
    "KPAR": 4,
    "NCORE": 2,
    "GGA": "PE",
    "AMIX": 0.2,
    "BMIX": 1e-4,
    "AMIX_MAG": 0.8,
    "BMIX_MAG": 1e-4,
}
relax_incar = {
    "KPAR": 4,
    "NCORE": 2,
}
phonon_incar = {
    "KPAR": 2,
    "NCORE": 4,
}

struct_names = [
    "SER-Co",
    "SER-Fe",
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


def run_locally(force=False):
    """Run all structures locally (legacy function, kept for compatibility)."""
    run_qha(force=force, execution_mode="local")


def run_qha(
    force: bool = False,
    execution_mode: str = "local",
    fw_config_dir: Path | None = None,
    wait_for_completion: bool = False,
):
    """Run all structures in specified mode.

    Args:
        force: Force re-run even if results exist
        execution_mode: "local" or "fireworks"
        fw_config_dir: Directory containing FireWorks config files (fireworks mode only)
        wait_for_completion: Wait for FireWorks workflow to complete (fireworks mode only)
    """
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
            if result.get("state", "failed") == "successful":
                log.info(f"Structure {name} already done")
                continue

        log.info(f"Structure {name} start (mode: {execution_mode})")

        # Run QhaWorker
        if force and flowdir.exists():
            shutil.rmtree(flowdir)
        flowdir.mkdir(parents=True, exist_ok=True)
        struct = endmember.get_poscar(name, posdir)
        try:
            worker = QhaWorker(
                worker_name=f"{name}-qha",
                vasp_args=vasp_args,
                global_incar=global_incar,
                relax_incar=relax_incar,
                eos_incar=relax_incar,
                phonon_incar=phonon_incar,
                temperature_range=(0, 3100, 50),
                supercell_matrix=((2, 0, 0), (0, 2, 0), (0, 0, 2)),
                execution_mode=execution_mode,
                fireworks_config_dir=fw_config_dir,
            )
            worker.run_flow(
                name,
                struct,
                flowdir,
                resume=not force,
                wait_for_completion=wait_for_completion,
            )
            qha_data = worker.get_result()
            if not qha_data:
                result = {"name": name, "state": "failed", "struct": struct.as_dict()}
            else:
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
    config = manager.get_cpu_config(
        job_name="gml-qha",
        output_log="em-qha.log",
        nodes=1,
        ntasks=48,
    )
    job_id = manager.submit_command(
        command=f"python {__file__} --local {'--force' if force else ''}",
        config=config,
        conda_env="htvasp",
        workdir=".",
    )
    log.info(f"Submitted job with ID: {job_id}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run QHA workflow")
    parser.add_argument("--local", action="store_true", help="Run locally (legacy)")
    parser.add_argument("--slurm", action="store_true", help="Submit as single Slurm job")
    parser.add_argument(
        "--mode",
        choices=["local", "fireworks"],
        default="local",
        help="Execution mode: local or fireworks",
    )
    parser.add_argument(
        "--fw-config-dir",
        type=Path,
        default=None,
        help="FireWorks configuration directory",
    )
    parser.add_argument(
        "--wait",
        action="store_true",
        help="Wait for FireWorks workflow completion",
    )
    parser.add_argument("--force", action="store_true", help="Force re-run even if already done")
    args = parser.parse_args()

    # Determine execution mode
    if args.local or args.slurm:
        execution_mode = "local"
    else:
        execution_mode = args.mode

    if execution_mode == "local":
        if args.slurm:
            submit_job(args.force)
        else:
            run_qha(args.force, execution_mode="local")
    elif execution_mode == "fireworks":
        if not args.fw_config_dir:
            log.error("FireWorks mode requires --fw-config-dir")
            exit(1)
        run_qha(
            args.force,
            execution_mode="fireworks",
            fw_config_dir=args.fw_config_dir,
            wait_for_completion=args.wait,
        )
