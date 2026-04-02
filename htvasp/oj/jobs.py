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
from pymatgen.io.vasp import Outcar

from htvasp.oj.input_set import OJInputSetGenerator, write_oj_input_set

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
        return {
            "flip_dir": str(flip_dir),
            "success": False,
            "error": result.stderr,
            "returncode": result.returncode,
        }

    energy = None
    outcar = flip_dir.joinpath("OUTCAR")
    if outcar.exists():
        try:
            energy = Outcar(outcar).final_energy
        except Exception as e:
            log.warning(f"Failed to parse OUTCAR in {flip_dir}: {e}")

    return {
        "flip_dir": str(flip_dir),
        "success": True,
        "energy": energy,
        "returncode": 0,
    }


@job
def oj_generate(
    structure: Structure, input_set_generator: OJInputSetGenerator
) -> dict[str, Any]:
    """Generate magnetic configurations using OstravaJ

    This job creates input files and runs OstravaJ generate command.
    Uses OJInputSetGenerator (inherits from VaspInputGenerator) for
    atomate2 compatibility.

    Args:
        structure: Input structure
        input_set_generator: OJInputSetGenerator instance

    Returns:
        Dictionary with run_dir, flip_dirs, and num_configs
    """
    job_dir = Path.cwd()

    temp_input_dir = job_dir.parent.joinpath(f"temp_input_{job_dir.name}")
    temp_input_dir.mkdir(parents=True, exist_ok=True)

    write_oj_input_set(
        structure,
        temp_input_dir,
        input_set_generator,
        include_oj_conf=True,
    )

    run_dir = job_dir
    cmd = [str(OJ_SCRIPT), "generate", "-i", str(temp_input_dir), "-r", str(run_dir)]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if temp_input_dir.exists():
        shutil.rmtree(temp_input_dir)

    if result.returncode != 0:
        raise RuntimeError(f"OJ generate failed: {result.stderr}")

    flip_dirs = sorted(run_dir.glob("flip*"))
    log.info(f"Generated {len(flip_dirs)} magnetic configurations in {run_dir}")

    return {
        "run_dir": str(run_dir.resolve()),
        "flip_dirs": [str(d.resolve()) for d in flip_dirs],
        "num_configs": len(flip_dirs),
    }


@job
def oj_vasp(flip_dirs: list[str], vasp_cmd: str) -> dict[str, Any]:
    """Run VASP for all magnetic configurations in serial

    Args:
        flip_dirs: List of flip directory paths
        vasp_cmd: VASP command to run

    Returns:
        Dictionary with results summary
    """
    results = []
    for flip_dir in flip_dirs:
        flip_path = Path(flip_dir)
        log.info(f"Running VASP in {flip_dir}")
        result = _run_vasp(flip_path, vasp_cmd)
        results.append(result)

        if result["success"]:
            log.info(
                f"VASP completed successfully in {flip_dir}, energy={result.get('energy')}"
            )
        else:
            log.error(f"VASP failed in {flip_dir}: {result.get('error')}")

    successful = sum(1 for r in results if r["success"])
    log.info(f"VASP batch: {successful}/{len(flip_dirs)} successful")

    return {
        "flip_dirs": flip_dirs,
        "results": results,
        "num_successful": successful,
        "num_total": len(flip_dirs),
    }


@job
def oj_solve(run_dir: str, vasp_results: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    """Solve for exchange parameters using OstravaJ

    Args:
        run_dir: Directory containing flip configurations
        vasp_results: Optional list of VASP results (for validation)

    Returns:
        Dictionary with J parameters and Tc values
    """
    run_path = Path(run_dir)

    if vasp_results:
        failed = [r for r in vasp_results if not r.get("success")]
        if failed:
            log.warning(f"{len(failed)} VASP calculations failed, solution may be inaccurate")

    cmd = [str(OJ_SCRIPT), "solve", "-r", str(run_path)]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return {
            "error": "Solve failed",
            "stderr": result.stderr,
            "run_dir": run_dir,
        }

    solution_file = run_path.joinpath("OJ_solution.json")
    if solution_file.exists():
        with open(solution_file) as f:
            solution = json.load(f)
        solution["run_dir"] = run_dir
        return solution

    return {
        "error": "No solution file",
        "run_dir": run_dir,
    }
