#!/usr/bin/env python3
"""DFPT-QHA portal: one endpoint per invocation, local or Slurm.

Mirrors ``main_qha.py``: ``VASP_ARGS`` + ``GLOBAL_INCAR`` + ``run_tick`` +
argparse (``--slurm`` / ``--dry-run`` / ``--rerun`` / ``--run-id``).

The scientific recipe is NOT taken from this file.  It comes from the pool's
``pool_recipes.json`` for the named endpoint, so the settings a run used are the
settings the recipe records -- there is no second place for them to drift.

Usage::

    python main_dfpt_qha.py -n BCC-Al-Ta-2 --dry-run     # print the plan, no VASP
    python main_dfpt_qha.py -n BCC-Al-Ta-2               # run locally
    python main_dfpt_qha.py -n BCC-Al-Ta-2 --slurm       # submit to Slurm

The flow is::

    EOS (relax at each volume) -> DFPT per volume -> validate -> thermal -> QHA

``--slurm`` submits one job that runs this script without ``--slurm``, so the
jobflow flow (including its dynamic per-volume expansion) executes inside the
allocation.
"""

import argparse
import json
import logging
import os
import shlex
from dataclasses import asdict, is_dataclass
from pathlib import Path

from pymatgen.core import Structure

from htvasp.model.dfpt_qha_recipe import Recipe
from htvasp.slurm import (
    SlurmJobManager,
    artifact_marker_matches,
    current_artifact_matches,
    mark_artifact_from_ledger,
    new_run_id,
    payload_fingerprint,
    stable_run_id,
)
from htvasp.workflows import DfptQhaWorker

DFPT_QHA_RECIPE_VERSION = "main_dfpt_qha-v1"

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

root_dir = Path(__file__).parent.resolve()
flow_dir = Path("/nfs_ssd/tmp")
conda_env = root_dir / ".conda"

POOL_DIR = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0")
POOL_RECIPES = POOL_DIR / "pool_recipes.json"
POOL_STRUCTURES = POOL_DIR / "pool_structures.json"

# Same verified launcher as main_qha.py / main_w90.py.  Without it atomate2
# calls a bare `vasp_std`, which is not on PATH here and is not run under srun.
VASP_ARGS = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-cpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-cpu && srun vasp_gam'",
}

# EOS deformation window and supercell.  The DFPT k-mesh, POTCAR, smearing,
# magnetic moments and displacement all come from the recipe.
LINEAR_STRAIN = (-0.05, 0.05)
NUMBER_OF_FRAMES = 6
SUPERCELL_MATRIX = ((2, 0, 0), (0, 2, 0), (0, 0, 2))
MESH = (12, 12, 12)
TEMPERATURES = list(range(0, 1001, 100))
PRESSURE_GPA = 0.0

DEFAULT_ENDPOINT = "BCC-Al-Ta-2"


def load_endpoint(name: str) -> tuple[Structure, Recipe]:
    """The endpoint's cell and recipe, both from the pool."""

    if not POOL_RECIPES.is_file():
        raise FileNotFoundError(f"missing pool recipes: {POOL_RECIPES}")
    if not POOL_STRUCTURES.is_file():
        raise FileNotFoundError(f"missing pool structures: {POOL_STRUCTURES}")
    recipes = json.loads(POOL_RECIPES.read_text())["endpoints"]
    structures = json.loads(POOL_STRUCTURES.read_text())["endpoints"]
    if name not in recipes:
        raise KeyError(f"endpoint {name!r} is absent from {POOL_RECIPES}")
    if name not in structures:
        raise KeyError(f"endpoint {name!r} is absent from {POOL_STRUCTURES}")
    recipe = Recipe.from_dict(recipes[name]["recipe"])
    structure = Structure.from_dict(structures[name]["structure"])
    if len(structure) != recipe.kpoints.cell_atoms:
        raise ValueError(
            f"{name}: recipe cell_atoms={recipe.kpoints.cell_atoms} but the pool "
            f"structure has {len(structure)} atoms"
        )
    return structure, recipe


def _stage_paths(name: str, stage: str, run_id: str | None = None):
    base = root_dir / "data" / "endpoints" / name / f"{stage}flow"
    if run_id:
        if Path(run_id).name != run_id or not run_id.strip():
            raise ValueError(f"Invalid run_id: {run_id!r}")
        run_root = base / "runs" / run_id
        return run_root, run_root / f"{name}-{stage}.json", base / "slurm_submissions.json"
    return base, base.parent / f"{name}-{stage}.json", base / "slurm_submissions.json"


def run_tick(
    name: str,
    rerun: bool = False,
    dry_run: bool = False,
    run_id: str | None = None,
    number_of_frames: int | None = None,
    mesh: tuple[int, int, int] | None = None,
):
    """Run one endpoint's DFPT-QHA flow locally (or print the plan)."""

    structure, recipe = load_endpoint(name)
    frames = NUMBER_OF_FRAMES if number_of_frames is None else number_of_frames
    dfpt_mesh = MESH if mesh is None else mesh
    if rerun:
        run_id = run_id or new_run_id("dfpt-qha")
    store_dir, json_path, ledger_path = _stage_paths(name, "dfpt-qha", run_id)

    if dry_run:
        log.info("endpoint=%s formula=%s atoms=%d", name, structure.formula, len(structure))
        log.info(
            "recipe: class=%s supercell=%s kpoints=%s/%s supercell_kpoints=%s/%s "
            "policy=%s",
            recipe.structure_class, list(recipe.supercell),
            list(recipe.kpoints.grid), list(recipe.kpoints.shift),
            list(recipe.kpoints.supercell_grid), list(recipe.kpoints.supercell_shift),
            recipe.magnetic_policy,
        )
        log.info("eos frames=%d strain=%s mesh=%s", frames, LINEAR_STRAIN, dfpt_mesh)
        worker = DfptQhaWorker(
            recipe=recipe,
            vasp_args=VASP_ARGS,
            potcar_functional=recipe.potcar.functional,
            supercell_matrix=SUPERCELL_MATRIX,
            linear_strain=LINEAR_STRAIN,
            number_of_frames=frames,
            mesh=dfpt_mesh,
            temperatures=TEMPERATURES,
            pressure_gpa=PRESSURE_GPA,
        )
        flow = worker._make_flow(structure)
        log.info("flow=%s jobs=%d", flow.name, len(flow.jobs))
        for index, job in enumerate(flow.jobs, 1):
            log.info("  %d. %s", index, job.name)
        log.info("dry-run OK: no VASP or Slurm launched")
        return

    if not rerun and json_path.exists() and (
        (run_id is None and not os.environ.get("HTVASP_INTENT_KEY"))
        or current_artifact_matches(ledger_path, json_path)
    ):
        log.info("DFPT-QHA for %s already done (use --rerun to re-run)", name)
        return

    log.info("DFPT-QHA for %s started", name)
    try:
        worker = DfptQhaWorker(
            recipe=recipe,
            vasp_args=VASP_ARGS,
            potcar_functional=recipe.potcar.functional,
            supercell_matrix=SUPERCELL_MATRIX,
            linear_strain=LINEAR_STRAIN,
            number_of_frames=frames,
            mesh=dfpt_mesh,
            temperatures=TEMPERATURES,
            pressure_gpa=PRESSURE_GPA,
        )
        worker.run_flow(
            name=name,
            structure=structure,
            flow_dir=flow_dir,
            store_dir=store_dir,
            resume=not rerun,
        )
        output = worker.get_result() or {}
        worker.write_result(data=output, json_path=json_path)
        mark_artifact_from_ledger(ledger_path, json_path)
        if output.get("ok"):
            eos = output.get("eos", {})
            log.info(
                "  V0=%.4f A^3  B0=%.2f GPa  T=%s..%s K",
                eos.get("V0_angstrom3", float("nan")),
                eos.get("B0_GPa", float("nan")),
                output["temperatures"][0], output["temperatures"][-1],
            )
        else:
            log.warning("  analysis did not succeed: %s", output.get("error"))
    except Exception as exc:  # noqa: BLE001 - the Worker already logs; keep the tick alive
        log.error("DFPT-QHA for %s failed", name)
        log.error(exc)

    log.info("DFPT-QHA for %s done", name)


def check_jobs(stage: str = "dfpt-qha") -> None:
    for name in [DEFAULT_ENDPOINT]:
        json_path = _stage_paths(name, stage)[1]
        if not json_path.exists():
            log.info("%s: not run", name)
            continue
        data = json.loads(json_path.read_text(encoding="utf-8"))
        if not data.get("ok"):
            log.info("%s: analysis failed: %s", name, data.get("error"))
            continue
        log.info(
            "%s: V0=%.4f A^3 B0=%.2f GPa, %d temperatures",
            name, data["eos"]["V0_angstrom3"], data["eos"]["B0_GPa"],
            len(data["temperatures"]),
        )


def _submission_plan(
    name: str,
    manager: SlurmJobManager,
    rerun: bool = False,
    number_of_frames: int | None = None,
    mesh: tuple[int, int, int] | None = None,
) -> dict:
    structure, recipe = load_endpoint(name)
    frames = NUMBER_OF_FRAMES if number_of_frames is None else number_of_frames
    dfpt_mesh = MESH if mesh is None else mesh
    run_payload = {
        "system": name,
        "stage": "dfpt-qha",
        "recipe": DFPT_QHA_RECIPE_VERSION,
        "endpoint_recipe": recipe.model_dump(mode="json"),
        "settings": {
            "vasp_args": VASP_ARGS,
            "linear_strain": LINEAR_STRAIN,
            "number_of_frames": frames,
            "supercell_matrix": SUPERCELL_MATRIX,
            "mesh": dfpt_mesh,
            "temperatures": TEMPERATURES,
            "pressure_gpa": PRESSURE_GPA,
        },
    }
    run_id = new_run_id("dfpt-qha") if rerun else stable_run_id("dfpt-qha", run_payload)
    run_root, result_path, ledger_path = _stage_paths(name, "dfpt-qha", run_id)
    command = (
        f"{shlex.quote(str(root_dir / '.conda' / 'bin' / 'python'))} "
        f"{shlex.quote(str(Path(__file__).resolve()))} -n {shlex.quote(name)}"
    )
    command += f" --run-id {shlex.quote(run_id)}"
    if number_of_frames is not None:
        command += f" --frames {int(frames)}"
    if mesh is not None:
        command += f" --mesh {int(dfpt_mesh[0])} {int(dfpt_mesh[1])} {int(dfpt_mesh[2])}"
    if rerun:
        command += " --rerun"
    # DFPT on a 16-atom supercell at ENCUT 520 with the recipe k-mesh: the
    # verified e1-dfpt run used 48 tasks / 96G / 6h on partCPU.
    config = manager.get_cpu_config(
        job_name=f"{name}-dfpt-qha",
        output_log=f"logs/{name}-dfpt-qha.log",
        ntasks=48,
        memory="96G",
        time_limit="24:00:00",
        module_name="vasp-cpu",
    )
    return {
        "name": name,
        "stage": "dfpt-qha",
        "recipe_version": DFPT_QHA_RECIPE_VERSION,
        "run_id": run_id,
        "command": command,
        "config": config,
        "config_dict": asdict(config) if is_dataclass(config) else dict(vars(config)),
        "artifact_path": result_path,
        "run_root": run_root,
        "ledger_path": ledger_path,
        "settings_fingerprint": payload_fingerprint(run_payload["settings"]),
        "manifest_fingerprint": payload_fingerprint({
            "system": name,
            "stage": "dfpt-qha",
            "recipe": DFPT_QHA_RECIPE_VERSION,
        }),
        "endpoint_recipe": recipe.model_dump(mode="json"),
    }


def _print_submission_plan(plan: dict, *, conflict: str | None = None) -> None:
    payload = {
        "target": plan["name"],
        "stage": plan["stage"],
        "recipe_version": plan["recipe_version"],
        "run_id": plan["run_id"],
        "command": plan["command"],
        "artifact_path": str(Path(plan["artifact_path"]).resolve()),
        "run_root": str(Path(plan["run_root"]).resolve()),
        "ledger_path": str(Path(plan["ledger_path"]).resolve()),
        "resources": plan["config_dict"],
    }
    if conflict:
        payload["conflict"] = conflict
    print(json.dumps(payload, ensure_ascii=False, indent=2, default=str))


def submit_job(
    name: str,
    manager: SlurmJobManager,
    rerun: bool = False,
    dry_run: bool = False,
    number_of_frames: int | None = None,
    mesh: tuple[int, int, int] | None = None,
):
    plan = _submission_plan(
        name, manager, rerun=rerun,
        number_of_frames=number_of_frames, mesh=mesh,
    )
    if dry_run:
        _print_submission_plan(plan)
        return plan

    result_path = Path(plan["artifact_path"])
    result = manager.submit_intent(
        system_name=plan["name"],
        stage=plan["stage"],
        recipe_version=plan["recipe_version"],
        settings_fingerprint=plan["settings_fingerprint"],
        manifest_fingerprint=plan["manifest_fingerprint"],
        command=plan["command"],
        config=plan["config"],
        ledger_path=plan["ledger_path"],
        conda_env=str(conda_env),
        work_dir=root_dir,
        run_id=plan["run_id"],
        retry=rerun,
        artifact_path=result_path,
        artifact_checker=lambda entry: artifact_marker_matches(result_path, entry),
        artifact_shortcut=True,
    )
    if not result.ok:
        log.error("Refused DFPT-QHA job for %s: %s", name, result.state)
        return -1
    log.info("DFPT-QHA job for %s: %s (%s)", name, result.state, result.job_id)
    return result.job_id


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DFPT-QHA workflow")
    parser.add_argument("--name", "-n", type=str, default=None, help="Endpoint name")
    parser.add_argument("--slurm", "-s", action="store_true", help="Submit to Slurm")
    parser.add_argument("--dry-run", action="store_true", help="Print the plan; no VASP")
    parser.add_argument("--flow-graph", action="store_true", help="Alias for --dry-run")
    parser.add_argument("--rerun", "-r", action="store_true", help="Force a re-run")
    parser.add_argument("--check", "-c", action="store_true", help="Check the result JSON")
    parser.add_argument("--run-id", help="Operate on an explicit immutable run")
    parser.add_argument(
        "--frames", type=int, default=None,
        help=f"EOS deformation frames (default {NUMBER_OF_FRAMES}); the QHA "
             f"analysis needs at least 4",
    )
    parser.add_argument(
        "--mesh", type=int, nargs=3, default=None, metavar=("N1", "N2", "N3"),
        help=f"DFPT thermal-integration mesh (default {list(MESH)}); it does NOT "
             f"change the VASP k-mesh, which comes from the recipe",
    )
    args = parser.parse_args()

    if args.slurm:
        if not args.name:
            parser.error("--slurm needs -n/--name")
        submit_job(
            args.name, SlurmJobManager(), rerun=args.rerun, dry_run=args.dry_run,
            number_of_frames=args.frames,
            mesh=tuple(args.mesh) if args.mesh else None,
        )
    elif args.dry_run or args.flow_graph:
        run_tick(
            args.name or DEFAULT_ENDPOINT, dry_run=True, run_id=args.run_id,
            number_of_frames=args.frames, mesh=tuple(args.mesh) if args.mesh else None,
        )
    elif args.check:
        check_jobs()
    elif args.name:
        run_tick(
            args.name, rerun=args.rerun, run_id=args.run_id,
            number_of_frames=args.frames, mesh=tuple(args.mesh) if args.mesh else None,
        )
    else:
        parser.print_help()
