"""
OJ Jobs

Job-decorated functions for OstravaJ workflow steps.
"""

import json
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Any

from jobflow import job
from pymatgen.core import Structure

from htvasp.oj.input_set import OJInputSetGenerator, write_oj_input_set

log = logging.getLogger(__name__)

OJ_SCRIPT = Path(__file__).parent.parent.parent.joinpath("ostravaj", "ostravaj.sh")


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
def oj_solve(run_dir: str, flip_task_docs: list[Any] | None = None) -> dict[str, Any]:
    """Solve for exchange parameters using OstravaJ

    Args:
        run_dir: Directory containing flip configurations
        flip_task_docs: List of TaskDoc objects or dicts from FlipMaker jobs.
            Each TaskDoc contains dir_name pointing to the flip directory
            with vasprun.xml.

    Returns:
        Dictionary with J parameters and Tc values
    """
    from emmet.core.tasks import TaskDoc

    run_path = Path(run_dir)

    # Filter successful flips based on TaskDoc state
    successful_flips = []
    if flip_task_docs:
        for task_doc in flip_task_docs:
            # Handle both TaskDoc objects and dicts
            if isinstance(task_doc, TaskDoc):
                # TaskDoc is a Pydantic model - access attributes directly
                state = getattr(task_doc, "state", None)
                flip_dir = getattr(task_doc, "dir_name", None)
            elif isinstance(task_doc, dict):
                # Fallback for dict format
                state = task_doc.get("state", "failed")
                flip_dir = task_doc.get("dir_name")
            else:
                log.warning(f"Unexpected task_doc type: {type(task_doc)}, skipping")
                continue

            # Check if calculation was successful
            if state == "successful" and flip_dir and Path(flip_dir).exists():
                successful_flips.append(flip_dir)

    if not successful_flips:
        log.error("No successful flip calculations found")
        return {
            "error": "No successful flips",
            "run_dir": run_dir,
        }

    log.info(f"Found {len(successful_flips)} successful flip calculations")

    # Call ostravaj solve - it will scan run_dir for flip directories
    # Note: Currently ostravaj CLI auto-scans run_dir, so we rely on all
    # successful flips being present there. Future enhancement could pass
    # explicit flip dirs via --flip-dirs argument.
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
        solution["num_successful_flips"] = len(successful_flips)
        solution["successful_flip_dirs"] = successful_flips
        return solution

    return {
        "error": "No solution file",
        "run_dir": run_dir,
    }
