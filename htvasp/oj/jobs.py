"""OJ Jobs - Job-decorated functions for OstravaJ workflow steps."""

import gzip
import json
import logging
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path
from typing import Any

from atomate2.vasp.jobs.core import StaticMaker
from jobflow import Flow, Response, job
from pymatgen.core import Structure

from htvasp.oj.input_set import OJInputSetGenerator, write_oj_input_set

log = logging.getLogger(__name__)


@lru_cache(maxsize=1)
def _get_oj_script_path() -> Path:
    """Get cached OJ script path."""
    return Path(__file__).parent.parent.parent / "ostravaj" / "ostravaj.sh"


@job
def oj_generate(
    structure: Structure, input_set_generator: OJInputSetGenerator
) -> dict[str, Any]:
    """Generate magnetic configurations using OstravaJ.

    Args:
        structure: Input structure
        input_set_generator: OJInputSetGenerator instance

    Returns:
        Dictionary with run_dir, flip_dirs, and num_configs
    """
    job_dir = Path.cwd()
    temp_input_dir = job_dir.parent / f"temp_input_{job_dir.name}"
    temp_input_dir.mkdir(parents=True, exist_ok=True)

    write_oj_input_set(structure, temp_input_dir, input_set_generator, include_oj_conf=True)

    cmd = [str(_get_oj_script_path()), "generate", "-i", str(temp_input_dir), "-r", str(job_dir)]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if temp_input_dir.exists():
        shutil.rmtree(temp_input_dir)

    if result.returncode != 0:
        raise RuntimeError(f"OJ generate failed: {result.stderr}")

    flip_dirs = sorted(job_dir.glob("flip*"))
    log.info(f"Generated {len(flip_dirs)} magnetic configurations in {job_dir}")

    return {
        "run_dir": str(job_dir.resolve()),
        "flip_dirs": [str(d.resolve()) for d in flip_dirs],
        "num_configs": len(flip_dirs),
    }


@job
def create_flip_jobs(
    flip_dirs: list[str],
    run_vasp_kwargs: dict[str, Any],
    input_set_generator: OJInputSetGenerator,
) -> Response:
    """Create StaticMaker jobs in new oj_flip_{idx} directories.

    Args:
        flip_dirs: List of original flip directory paths from oj_generate.
        run_vasp_kwargs: VASP execution parameters.
        input_set_generator: Input set generator for VASP calculations.

    Returns:
        Response with detour Flow containing StaticMaker jobs.
    """
    job_dir = Path.cwd()
    flip_jobs = []
    new_flip_dirs = []
    temp_flip_dirs = []

    try:
        for idx, flip_dir in enumerate(flip_dirs):
            flip_path = Path(flip_dir)
            new_flip_dir = job_dir / f"oj_flip_{idx}"
            new_flip_dir.mkdir(parents=True, exist_ok=True)
            temp_flip_dirs.append(new_flip_dir)

            # Copy required input files
            for filename in ["POSCAR", "INCAR", "KPOINTS", "POTCAR"]:
                src_file = flip_path / filename
                if src_file.exists():
                    shutil.copy(src_file, new_flip_dir / filename)
                else:
                    log.warning(f"File not found: {src_file}")

            # Create StaticMaker job
            structure = Structure.from_file(new_flip_dir / "POSCAR")
            flip_maker = StaticMaker(
                name=f"oj_flip_{idx}",
                run_vasp_kwargs=run_vasp_kwargs,
                input_set_generator=input_set_generator,
            )
            flip_jobs.append(flip_maker.make(structure))
            new_flip_dirs.append(str(new_flip_dir.resolve()))
    except Exception as e:
        log.error(f"Failed to create flip jobs: {e}")
        raise

    # Output: [job0.output, job1.output, ..., dir0, dir1, ...]
    output_list = [job.output for job in flip_jobs] + new_flip_dirs
    return Response(detour=Flow(jobs=flip_jobs, output=output_list), output=output_list)


def _decompress_gz_files(flip_dir: Path) -> None:
    """Decompress .gz files in flip directory if needed."""
    for name in ["vasprun.xml", "POSCAR", "INCAR"]:
        gz_path = flip_dir / f"{name}.gz"
        target_path = flip_dir / name
        if gz_path.exists() and not target_path.exists():
            log.debug(f"Decompressing {gz_path}")
            with gzip.open(gz_path, "rb") as f_in, open(target_path, "wb") as f_out:
                f_out.write(f_in.read())


@job
def oj_solve(flow_output: list[Any], num_flips: int) -> dict[str, Any]:
    """Solve for exchange parameters using OstravaJ.

    Args:
        flow_output: Flattened output from create_flip_jobs.
        num_flips: Number of flip jobs.

    Returns:
        Dictionary with J parameters and Tc values.
    """
    if not flow_output:
        log.error("Empty flow_output received")
        return {"error": "Empty flow_output", "successful_flip_dirs": []}

    log.debug(f"flow_output length: {len(flow_output)}, num_flips: {num_flips}")

    # Split into task docs and flip dirs
    flip_task_docs = flow_output[:num_flips]
    flip_dirs = flow_output[num_flips:]

    # Filter successful calculations
    successful_flip_dirs = []
    for idx, task_doc in enumerate(flip_task_docs):
        if isinstance(task_doc, dict):
            state = task_doc.get("state", "failed")
            flip_dir = task_doc.get("dir_name")
        else:
            state = getattr(task_doc, "state", None)
            flip_dir = getattr(task_doc, "dir_name", None)

        log.debug(f"Task {idx}: state={state}, dir={flip_dir}")

        if state == "successful" and flip_dir:
            if ":" in flip_dir:
                flip_dir = flip_dir.split(":", 1)[1]
            flip_path = Path(flip_dir)
            if flip_path.exists():
                successful_flip_dirs.append(flip_dir)
            else:
                log.warning(f"Flip directory not found: {flip_dir}")

    if not successful_flip_dirs:
        log.error("No successful flip calculations found")
        return {"error": "No successful flips", "successful_flip_dirs": []}

    log.info(f"Found {len(successful_flip_dirs)} successful flip calculations")

    # Setup run_dir and input directory
    run_dir_path = Path(successful_flip_dirs[0]).parent

    # Find and copy input directory
    input_dir = None
    for sibling in run_dir_path.iterdir():
        if sibling.is_dir() and (sibling / "input").exists():
            input_dir = sibling / "input"
            break

    if input_dir and not (run_dir_path / "input").exists():
        shutil.copytree(input_dir, run_dir_path / "input")
        log.debug(f"Copied input directory to: {run_dir_path / 'input'}")

    # Decompress files in flip directories
    for flip_dir in successful_flip_dirs:
        _decompress_gz_files(Path(flip_dir))

    # Run solve command
    cmd = [
        str(_get_oj_script_path()),
        "solve",
        "-r", str(run_dir_path),
        "--flip-dirs",
        *successful_flip_dirs,
    ]

    log.debug(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        log.error(f"Solve failed: {result.stderr}")
        return {
            "error": "Solve failed",
            "stderr": result.stderr,
            "stdout": result.stdout,
            "successful_flip_dirs": successful_flip_dirs,
        }

    # Read solution file
    solution_path = run_dir_path / "OJ_solution.json"
    if solution_path.exists():
        solution = json.load(open(solution_path))
        solution["num_successful_flips"] = len(successful_flip_dirs)
        solution["successful_flip_dirs"] = successful_flip_dirs
        return solution

    return {"error": "No solution file", "successful_flip_dirs": successful_flip_dirs}
