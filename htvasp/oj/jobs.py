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
    structure: Structure,
    input_set_generator: OJInputSetGenerator | None = None,
    j_count: int | None = 2,
    dist_cutoff: float | None = None,
    magnetic_ion_types: list[str] | None = None,
    noncollinear: bool = False,
    base_spin: float | list[float] = 1.0,
    extend_poscar: tuple[int, int, int] = (2, 2, 2),
    reciprocal_density: float = 100,
    reciprocal_density_metal: float = 400,
    auto_metal_kpoints: bool = True,
    force_gamma: bool = True,
    user_incar_settings: dict[str, Any] | None = None,
    potcar_functional: Literal["PBE", "PBE_54", "PBE_64"] = "PBE_64",
) -> dict[str, Any]:
    """Generate magnetic configurations using OstravaJ

    This job creates input files and runs OstravaJ generate command.
    Uses OJInputSetGenerator (inherits from VaspInputGenerator) for
    atomate2 compatibility.

    Args:
        structure: Input structure
        input_set_generator: OJInputSetGenerator instance (atomate2 style)
        j_count: Number of J pairs to consider
        dist_cutoff: Distance cutoff for J pairs (mutually exclusive with j_count)
        magnetic_ion_types: List of magnetic ion types
        noncollinear: Whether to use noncollinear magnetism
        base_spin: Base spin value
        extend_poscar: Supercell extension factors (nx, ny, nz)
        reciprocal_density: K-point density for insulators
        reciprocal_density_metal: K-point density for metals
        auto_metal_kpoints: Automatically use higher density for metals
        force_gamma: Force gamma-centered k-point mesh
        user_incar_settings: Custom INCAR settings
        potcar_functional: POTCAR functional type

    Returns:
        Dictionary with run_dir, flip_dirs, and num_configs
    """
    job_dir = Path.cwd()

    temp_input_dir = job_dir.parent.joinpath(f"temp_input_{job_dir.name}")
    temp_input_dir.mkdir(parents=True, exist_ok=True)

    # Use provided input_set_generator or create one from parameters
    if input_set_generator is None:
        generator = OJInputSetGenerator(
            j_count=j_count,
            dist_cutoff=dist_cutoff,
            magnetic_ion_types=magnetic_ion_types or [],
            noncollinear=noncollinear,
            base_spin=base_spin,
            extend_poscar=extend_poscar,
            reciprocal_density=reciprocal_density,
            reciprocal_density_metal=reciprocal_density_metal,
            auto_metal_kpoints=auto_metal_kpoints,
            force_gamma=force_gamma,
            user_incar_settings=user_incar_settings or {},
            user_potcar_functional=potcar_functional,
        )
    else:
        generator = input_set_generator

    write_oj_input_set(
        structure,
        temp_input_dir,
        generator,
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
def oj_vasp(
    flip_dirs: list[str],
    vasp_cmd: str = "vasp_std",
) -> dict[str, Any]:
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
def oj_solve(
    run_dir: str,
    vasp_results: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
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


@job
def oj_workflow(
    structure: Structure,
    j_count: int | None = 2,
    dist_cutoff: float | None = None,
    magnetic_ion_types: list[str] | None = None,
    noncollinear: bool = False,
    base_spin: float | list[float] = 1.0,
    extend_poscar: tuple[int, int, int] = (2, 2, 2),
    reciprocal_density: float = 100,
    reciprocal_density_metal: float = 400,
    auto_metal_kpoints: bool = True,
    force_gamma: bool = True,
    user_incar_settings: dict[str, Any] | None = None,
    vasp_cmd: str = "vasp_std",
    potcar_functional: Literal["PBE", "PBE_54", "PBE_64"] = "PBE_64",
) -> dict[str, Any]:
    """Run complete OJ workflow in a single job (simplified mode)

    This runs all steps in one job without Flow composition.
    Use for simple serial execution.

    Args:
        structure: Input structure
        j_count: Number of J pairs to consider
        dist_cutoff: Distance cutoff for J pairs
        magnetic_ion_types: List of magnetic ion types
        noncollinear: Whether to use noncollinear magnetism
        base_spin: Base spin value
        extend_poscar: Supercell extension factors
        reciprocal_density: K-point density for insulators
        reciprocal_density_metal: K-point density for metals
        auto_metal_kpoints: Automatically use higher density for metals
        force_gamma: Force gamma-centered k-point mesh
        user_incar_settings: Custom INCAR settings
        vasp_cmd: VASP command
        potcar_functional: POTCAR functional type

    Returns:
        Dictionary with solution results
    """
    gen = oj_generate.function(
        structure,
        j_count=j_count,
        dist_cutoff=dist_cutoff,
        magnetic_ion_types=magnetic_ion_types,
        noncollinear=noncollinear,
        base_spin=base_spin,
        extend_poscar=extend_poscar,
        reciprocal_density=reciprocal_density,
        reciprocal_density_metal=reciprocal_density_metal,
        auto_metal_kpoints=auto_metal_kpoints,
        force_gamma=force_gamma,
        user_incar_settings=user_incar_settings,
        potcar_functional=potcar_functional,
    )
    vasp = oj_vasp.function(gen["flip_dirs"], vasp_cmd)
    solution = oj_solve.function(gen["run_dir"], vasp["results"])

    solution["workdir"] = str(Path.cwd())
    solution["num_configs"] = gen["num_configs"]
    solution["vasp_success_rate"] = (
        vasp["num_successful"] / vasp["num_total"] if vasp["num_total"] > 0 else 0
    )

    return solution
