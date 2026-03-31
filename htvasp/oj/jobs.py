"""
OJ Jobs

Job-decorated functions for OstravaJ workflow steps.
"""

import json
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Any, Literal

from jobflow import job
from pymatgen.core import Structure

from htvasp.oj.config import OJConfig
from htvasp.oj.input_set import write_oj_inputs

log = logging.getLogger(__name__)

OJ_SCRIPT = Path(__file__).parent.parent.parent.joinpath("ostravaj", "ostravaj.sh")


def _run_vasp(flip_dir: Path, vasp_cmd: str) -> dict[str, Any]:
    """Run VASP in a single directory"""
    result = subprocess.run(
        vasp_cmd,
        shell=True,
        cwd=flip_dir,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:
        return {"flip_dir": str(flip_dir), "success": False, "error": result.stderr}

    energy = None
    outcar = flip_dir.joinpath("OUTCAR")
    if outcar.exists():
        try:
            from pymatgen.io.vasp import Outcar

            energy = Outcar(outcar).final_energy
        except Exception:
            pass

    return {"flip_dir": str(flip_dir), "success": True, "energy": energy}


@job
def oj_generate(
    structure: Structure,
    config: OJConfig,
    workdir: str = "./oj_work",
    potcar_functional: Literal["PBE", "PBE_54", "PBE_64"] = "PBE_64",
) -> dict[str, Any]:
    """Generate magnetic configurations using OstravaJ"""
    workdir = Path(workdir)
    input_dir = workdir.joinpath("input")
    input_dir.mkdir(parents=True, exist_ok=True)

    write_oj_inputs(structure, input_dir, config, potcar_functional, kppa=config.kppa)

    run_dir = workdir.joinpath("run_oj")
    cmd = [str(OJ_SCRIPT), "generate", "-i", str(input_dir), "-r", str(run_dir)]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"OJ generate failed: {result.stderr}")

    # Copy input directory to run_dir for solve command
    run_input_dir = run_dir.joinpath("input")
    if not run_input_dir.exists():
        shutil.copytree(input_dir, run_input_dir)
        log.info(f"Copied input directory to {run_input_dir}")

    flip_dirs = sorted(run_dir.glob("flip*"))
    log.info(f"Generated {len(flip_dirs)} magnetic configurations")

    return {
        "workdir": str(workdir.resolve()),
        "run_dir": str(run_dir.resolve()),
        "flip_dirs": [str(d.resolve()) for d in flip_dirs],
        "num_configs": len(flip_dirs),
    }


@job
def oj_vasp(flip_dirs: list[str], vasp_cmd: str = "vasp_std") -> dict[str, Any]:
    """Run VASP for all magnetic configurations (serial)"""
    results = [_run_vasp(Path(d), vasp_cmd) for d in flip_dirs]
    successful = sum(1 for r in results if r["success"])

    log.info(f"VASP batch: {successful}/{len(flip_dirs)} successful")

    return {
        "flip_dirs": flip_dirs,
        "results": results,
        "num_successful": successful,
        "num_total": len(flip_dirs),
    }


@job
def oj_solve(run_dir: str) -> dict[str, Any]:
    """Solve for exchange parameters using OstravaJ"""
    run_path = Path(run_dir)
    cmd = [str(OJ_SCRIPT), "solve", "-r", str(run_path)]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return {"error": "Solve failed", "stderr": result.stderr}

    solution_file = run_path.joinpath("OJ_solution.json")
    if solution_file.exists():
        with open(solution_file) as f:
            return json.load(f)

    return {"error": "No solution file"}


@job
def oj_workflow(
    structure: Structure,
    config: OJConfig,
    workdir: str = "./oj_work",
    vasp_cmd: str = "vasp_std",
) -> dict[str, Any]:
    """Run complete OJ workflow in a single job"""
    gen = oj_generate.function(structure, config, workdir)
    vasp = oj_vasp.function(gen["flip_dirs"], vasp_cmd)
    solution = oj_solve.function(gen["run_dir"])

    solution["workdir"] = workdir
    solution["num_configs"] = gen["num_configs"]
    solution["vasp_success_rate"] = vasp["num_successful"] / vasp["num_total"]

    return solution
