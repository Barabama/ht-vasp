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
import hashlib
import json
import logging
import os
import shlex
from dataclasses import asdict, is_dataclass
from pathlib import Path

import numpy as np
from pymatgen.core import Structure

from htvasp.dfpt_qha.qha import continuous_thermal_eos
from htvasp.model.dfpt_qha_recipe import Recipe
from htvasp.utils.vasp_psp import psp_dir_exports
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

# Initial magnetic moments per element (muB).  Copied verbatim from
# main_qha.py (itself from the verified submit_tb2j.sh) so the DFPT and QHA/TB2J
# entries start from the SAME spin guess: a zero start falls into the
# non-magnetic solution for Cr-containing cells, which is 0.28 eV/atom above the
# magnetic ground state.  Elements outside the 7-element pool are carried for
# the shared pool; the values for the pool elements are identical to the
# package default in htvasp.dfpt_qha.input_set.
MAGMOM = {
    "Al": 1.0,
    "Co": 3.0,
    "Cr": 3.0,
    "Cu": 5.0,
    "Fe": 5.0,
    "Hf": 1.0,
    "Mn": 5.0,
    "Mo": 3.0,
    "Nb": 3.0,
    "Ni": 2.0,
    "Re": 1.0,
    "Ta": 3.0,
    "Ti": 3.0,
    "V": 3.0,
    "W": 3.0,
    "Zr": 1.0,
}

#: Endpoints whose spin screen chose the antiferromagnetic branch.  The screen
#: ran free nm / fm / afm relaxations and compared energies; ``afm`` won for
#: these two, and both are the same 2-atom bcc Cr cell, so the pattern is one
#: sign per sublattice.  Kept as a table rather than derived from
#: ``magnetic_policy``: that field says "a magnetic branch was confirmed", not
#: *which* one, and a ferromagnetic start on these two converges back to the
#: non-magnetic solution (the screen measured a 0.0 meV/atom gain for fm).
AFM_ENDPOINTS = {"BCC-Cr-Cr-2", "SER-Cr-1"}


def site_moments(name: str, structure: Structure) -> tuple[dict[str, float] | list[float], Structure | None]:
    """Initial moments for an endpoint, and the cell a per-site list refers to.

    Returns the per-element table for every ordinary endpoint.  For the two
    antiferromagnetic ones it returns a per-site list -- ``+3`` on one
    sublattice and ``-3`` on the other -- together with the cell that list was
    written for, because the DFPT stage runs on a supercell and the pattern has
    to be re-expanded there.
    """
    if name not in AFM_ENDPOINTS:
        return MAGMOM, None
    species = [str(site.specie) for site in structure]
    if species.count("Cr") < 2:
        raise ValueError(
            f"{name} is listed as antiferromagnetic but its cell has "
            f"{species.count('Cr')} Cr atoms; the staggered pattern needs two"
        )
    # Stagger by site order, the same convention the spin screen used, so the
    # DFPT run starts from the state that was actually measured to be lowest.
    moments = [3.0 if index % 2 == 0 else -3.0 for index in range(len(species))]
    return moments, structure


# EOS deformation window and supercell.  The DFPT k-mesh, POTCAR, smearing and
# displacement come from the recipe; MAGMOM comes from the table above.
LINEAR_STRAIN = (-0.05, 0.05)
NUMBER_OF_FRAMES = 6
SUPERCELL_MATRIX = ((2, 0, 0), (0, 2, 0), (0, 0, 2))
MESH = (12, 12, 12)
TEMPERATURES = list(range(0, 1001, 100))
PRESSURE_GPA = 0.0

#: MPI ranks per job.  The cluster has two nodes (48 and 64 cores) and both are
#: shared with partGPU, so a 48-rank job monopolises a node and the pool runs
#: strictly serially.  32 fits two jobs at once and was measured not to slow the
#: per-volume DFPT noticeably.
DEFAULT_NTASKS = 32

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
    ncore: int | None = None,
    kpar: int | None = None,
    resume_from: Path | None = None,
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
        moments, moments_cell = site_moments(name, structure)
        worker = DfptQhaWorker(
            recipe=recipe,
            vasp_args=VASP_ARGS,
            potcar_functional=recipe.potcar.functional,
            magmom=moments,
            magmom_cell=moments_cell,
            ncore=ncore,
            kpar=kpar,
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
        moments, moments_cell = site_moments(name, structure)
        worker = DfptQhaWorker(
            recipe=recipe,
            vasp_args=VASP_ARGS,
            potcar_functional=recipe.potcar.functional,
            magmom=moments,
            magmom_cell=moments_cell,
            ncore=ncore,
            kpar=kpar,
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
            # The run id is stable for a given endpoint+settings, so it names
            # the scratch directory.  ``flow.uuid`` is regenerated on every
            # invocation and would orphan the completed volumes of a job that
            # was killed at its time limit.
            flow_dir_name=run_id,
            resume_from=resume_from,
        )
        output = worker.get_result() or {}
        if output:
            worker.write_result(data=output, json_path=json_path)
            mark_artifact_from_ledger(ledger_path, json_path)
        else:
            # Nothing was stored: the flow died without a verdict to report.
            # Marking the artifact here would make the next invocation skip the
            # endpoint as "already done".
            log.error("  no result was stored; the endpoint is not marked done")
        if output.get("ok"):
            eos = output.get("eos", {})
            log.info(
                "  V0=%.4f A^3  B0=%.2f GPa  T=%s..%s K  applicability=%s",
                eos.get("V0_angstrom3", float("nan")),
                eos.get("B0_GPa", float("nan")),
                output["temperatures"][0], output["temperatures"][-1],
                output.get("qha_applicability"),
            )
            if output.get("dynamically_unstable_conditional"):
                log.warning(
                    "  conditional product: volumes %s are not harmonic "
                    "(branches excluded from the thermal sums).  F(V,T) is "
                    "not a harmonic free energy for those volumes; do not "
                    "promote this endpoint to a production TDB.",
                    output.get("non_harmonic_volume_indices"),
                )
        elif output.get("gate") == "dfpt_validation":
            log.warning(
                "  stopped at the DFPT validation gate: %d volume(s) rejected",
                output.get("n_failed", 0),
            )
            for failure in output.get("failures") or []:
                log.warning("    %s: %s", failure["job_name"], failure["error"])
        elif output.get("gate") == "dfpt_stability":
            # Only reached when the analysis job itself did not run; the gate no
            # longer stops the flow.  The verdicts are still the best record.
            log.warning(
                "  no analysis was stored; stability verdict: %s, "
                "%d of %d volumes non-harmonic",
                output.get("qha_applicability"),
                output.get("n_non_harmonic", 0), output.get("n_volumes", 0),
            )
            for report in output.get("reports") or []:
                if not report.get("qha_allowed"):
                    log.warning(
                        "    volume %s: %s (%s)",
                        report.get("index"), report.get("status"),
                        report.get("reason"),
                    )
            if output.get("confirmation_volume_indices"):
                log.warning(
                    "    borderline volumes needing a larger supercell: %s",
                    output["confirmation_volume_indices"],
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
        if data.get("gate") == "dfpt_validation":
            log.info(
                "%s: stopped at the DFPT validation gate, %d volume(s) rejected",
                name, data.get("n_failed", 0),
            )
            for failure in data.get("failures") or []:
                log.info("  %s: %s", failure["job_name"], failure["error"])
            continue
        if data.get("gate") == "dfpt_stability":
            log.info(
                "%s: no analysis stored; %d/%d volumes non-harmonic (%s)",
                name, data.get("n_non_harmonic", 0), data.get("n_volumes", 0),
                data.get("qha_applicability"),
            )
            for report in data.get("reports") or []:
                if not report.get("qha_allowed"):
                    log.info(
                        "  volume %s: %s (%s)",
                        report.get("index"), report.get("status"),
                        report.get("reason"),
                    )
            if data.get("confirmation_volume_indices"):
                log.info(
                    "  borderline volumes needing a larger supercell: %s",
                    data["confirmation_volume_indices"],
                )
            continue
        if not data.get("ok"):
            log.info("%s: analysis failed: %s", name, data.get("error"))
            continue
        log.info(
            "%s: V0=%.4f A^3 B0=%.2f GPa, %d temperatures, applicability=%s",
            name, data["eos"]["V0_angstrom3"], data["eos"]["B0_GPa"],
            len(data["temperatures"]), data.get("qha_applicability"),
        )
        if data.get("dynamically_unstable_conditional"):
            log.warning(
                "  conditionally unstable: volumes %s; not a production TDB "
                "endmember without an explicit decision",
                data.get("non_harmonic_volume_indices"),
            )


def _submission_plan(
    name: str,
    manager: SlurmJobManager,
    rerun: bool = False,
    number_of_frames: int | None = None,
    mesh: tuple[int, int, int] | None = None,
    ntasks: int = DEFAULT_NTASKS,
    ncore: int | None = None,
    kpar: int | None = None,
    resume_from: Path | None = None,
    nodelist: str | None = None,
    memory: str = "96G",
    run_id: str | None = None,
) -> dict:
    structure, recipe = load_endpoint(name)
    moments, _ = site_moments(name, structure)
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
            "ntasks": int(ntasks),
            "temperatures": TEMPERATURES,
            "pressure_gpa": PRESSURE_GPA,
        },
    }
    # Only recorded when the caller actually sets them.  Adding the keys
    # unconditionally would change ``stable_run_id`` for all 84 endpoints and
    # throw away every validated store; this way the default run ids stay put
    # and an explicit parallelisation change still mints a new one.
    if ncore is not None:
        run_payload["settings"]["ncore"] = int(ncore)
    if kpar is not None:
        run_payload["settings"]["kpar"] = int(kpar)
    # The node is a resource choice, not a scientific one, so it is recorded in
    # the intent identity but not in ``settings`` -- putting it there would mint
    # a new run id and discard the store every time a job moves between nodes.
    if nodelist:
        run_payload["nodelist"] = nodelist
    if run_id is not None:
        # An explicit id exists to point a job back at a store that was archived
        # under a different id -- e.g. before a pool rewrite changed
        # ``endpoint_recipe`` and therefore ``stable_run_id``, while the volumes
        # already on disk stayed valid.  ``_prepare_resume`` matches stored jobs
        # by (name, index), so such a run really does skip its finished volumes.
        # It never overrides the settings that produced those volumes: the job
        # still runs with today's recipe, which is why this is a bookkeeping
        # escape hatch and not a way to change physics.
        if rerun:
            raise ValueError("--run-id and --rerun are mutually exclusive")
        if Path(run_id).name != run_id or not run_id.strip():
            raise ValueError(f"Invalid run_id: {run_id!r}")
    else:
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
    # The job re-runs this script without --slurm, so the parallelisation has to
    # travel in the command line or the child would regenerate the defaults.
    if ncore is not None:
        command += f" --ncore {int(ncore)}"
    if kpar is not None:
        command += f" --kpar {int(kpar)}"
    # The child re-runs this script, so a resume source has to travel with it.
    if resume_from is not None:
        command += f" --resume-from {shlex.quote(str(resume_from))}"
    if rerun:
        command += " --rerun"
    # 16-atom supercell, ENCUT 520, recipe k-mesh.  ntasks is a throughput
    # choice, not a scientific one: it changes only how the ranks are laid out,
    # and the INCAR's KPAR=2 is independent of it.  32 lets two jobs share the
    # two nodes instead of monopolising one.
    config = manager.get_cpu_config(
        job_name=f"{name}-dfpt-qha",
        output_log=f"logs/{name}-dfpt-qha.log",
        ntasks=int(ntasks),
        memory=memory,
        time_limit="100:00:00",
        module_name="vasp-cpu",
        # the compute node's own pymatgen config may point at a POTCAR
        # directory that does not exist on it (B5071 does); pin the value the
        # submitting node resolved.
        extra_commands=psp_dir_exports(),
    )
    if nodelist:
        config.nodelist = nodelist
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
        # The MAGMOM table lives in the MANIFEST fingerprint, not in
        # ``settings``: the settings dict feeds ``stable_run_id``, and the table
        # is a constant whose effective values for the pool elements did not
        # change when it was made explicit.  Putting it in ``settings`` would
        # mint a new run_id for every already-validated endpoint and throw away
        # its store.  Here it is still part of the intent identity, so a real
        # change to the table still produces a new ledger entry.
        "manifest_fingerprint": payload_fingerprint({
            "system": name,
            "stage": "dfpt-qha",
            "recipe": DFPT_QHA_RECIPE_VERSION,
            # The *effective* initial moments, not the per-element table: the
            # antiferromagnetic endpoints carry a per-site pattern, and a run
            # started from +3/+3 is a different calculation from one started
            # from +3/-3 even though both use the same element table.
            "magmom": (
                list(moments)
                if isinstance(moments, (list, tuple))
                else dict(MAGMOM)
            ),
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


def _latest_analysis_json(name: str) -> Path | None:
    """The most recent analysis result for an endpoint, if any."""

    base = root_dir / "data" / "endpoints" / name / "dfpt-qhaflow" / "runs"
    files = sorted(
        base.glob(f"*/{name}-dfpt-qha.json"),
        key=lambda path: path.stat().st_mtime,
    )
    return files[-1] if files else None


def refit_vt(names: list[str]) -> list[dict]:
    """Re-fit each endpoint's F(V,T) per temperature and write a derived artifact.

    Pure post-processing: it reads the analysis JSON that already holds the
    per-volume free energies, fits the Birch-Murnaghan form at every
    temperature, and writes ``<name>-dfpt-qha-vt.json`` beside it.  The
    analysis JSON is not modified -- runs are immutable -- and the derived file
    records the source path and its SHA-256.
    """

    results: list[dict] = []
    for name in names:
        source = _latest_analysis_json(name)
        if source is None:
            log.warning("%s: no analysis JSON", name)
            results.append({"endpoint": name, "ok": False, "error": "no analysis JSON"})
            continue
        document = json.loads(source.read_text())
        if not document.get("ok") or "free_energies_eV" not in document:
            log.warning("%s: analysis did not succeed; nothing to refit", name)
            results.append(
                {"endpoint": name, "ok": False, "error": "analysis did not succeed"}
            )
            continue
        try:
            output = continuous_thermal_eos(
                document["volumes"],
                document["static_energies_eV"],
                document["free_energies_eV"],
                document["temperatures"],
                pressure_gpa=float(document.get("pressure_gpa", 0.0)),
                n_atoms_cell=int(document.get("n_atoms_cell", 1)),
            )
        except Exception as exc:  # noqa: BLE001 - one endpoint must not stop a batch
            log.error("%s: refit failed: %s", name, exc)
            results.append(
                {"endpoint": name, "ok": False, "error": f"{type(exc).__name__}: {exc}"}
            )
            continue
        output.update(
            {
                "ok": True,
                "endpoint": name,
                "source_json": str(source),
                "source_json_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                "qha_applicability": document.get("qha_applicability"),
                "non_harmonic_volume_indices": document.get(
                    "non_harmonic_volume_indices", []
                ),
            }
        )
        target = source.with_name(f"{name}-dfpt-qha-vt.json")
        target.write_text(json.dumps(output, indent=2, ensure_ascii=False, default=str))
        checks = output["checks"]
        index = min(
            range(len(output["temperatures"])),
            key=lambda i: abs(output["temperatures"][i] - 300.0),
        )
        log.info(
            "  %-14s res=%.1f meV/atom  V(%.0fK)=%.3f A^3  B=%.1f GPa  "
            "mono=%s inside=%s B>0=%s",
            name,
            checks["max_residual_meV_per_atom"],
            output["temperatures"][index],
            output["equilibrium_volume_angstrom3"][index],
            output["bulk_modulus_GPa"][index],
            checks["volume_monotonic_in_temperature"],
            checks["all_inside_sampled_window"],
            checks["all_bulk_moduli_positive"],
        )
        if not checks["residual_within_tolerance"]:
            log.warning(
                "    F(V,T) is not EOS-shaped at some temperatures "
                "(max residual %.1f meV/atom); V(T) and G(T) are then an "
                "artefact of the fit, not of the data",
                checks["max_residual_meV_per_atom"],
            )
        results.append(output)
    return results


def read_batch(path: Path) -> list[str]:
    """Endpoint names from a batch file: one per line, '#' comments."""

    names: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        name = line.split("#", 1)[0].strip()
        if name:
            names.append(name)
    if not names:
        raise ValueError(f"batch file {path} lists no endpoints")
    return names


def submit_batch(
    path: Path,
    manager: SlurmJobManager,
    rerun: bool = False,
    resume: bool = False,
    dry_run: bool = False,
    number_of_frames: int | None = None,
    mesh: tuple[int, int, int] | None = None,
    ntasks: int = DEFAULT_NTASKS,
    ncore: int | None = None,
    kpar: int | None = None,
    resume_from: Path | None = None,
    nodelist: str | None = None,
    memory: str = "96G",
) -> list[dict]:
    """Submit every endpoint in a batch file, one Slurm job each.

    The ledger's artifact shortcut makes this safe to re-run: an endpoint whose
    JSON is already present is reported as completed and not resubmitted.
    """

    names = read_batch(path)
    results: list[dict] = []
    for index, name in enumerate(names, 1):
        log.info("[%d/%d] %s", index, len(names), name)
        try:
            job_id = submit_job(
                name, manager, rerun=rerun, resume=resume, dry_run=dry_run,
                number_of_frames=number_of_frames, mesh=mesh, ntasks=ntasks,
                ncore=ncore, kpar=kpar, resume_from=resume_from,
                nodelist=nodelist,
            )
        except Exception as exc:  # noqa: BLE001 - one endpoint must not stop a batch
            log.error("  %s: submission failed: %s", name, exc)
            results.append({"name": name, "job_id": None, "state": "error",
                            "error": str(exc)})
            continue
        state = "dry-run" if dry_run else ("refused" if job_id == -1 else "submitted")
        results.append({"name": name, "job_id": job_id, "state": state})
        if not dry_run:
            log.info("  -> %s (%s)", job_id, state)
    submitted = [r for r in results if r["state"] == "submitted"]
    log.info(
        "batch %s: %d endpoints, %d submitted, %d refused, %d errors",
        path.name, len(results), len(submitted),
        sum(1 for r in results if r["state"] == "refused"),
        sum(1 for r in results if r["state"] == "error"),
    )
    return results


def submit_job(
    name: str,
    manager: SlurmJobManager,
    rerun: bool = False,
    resume: bool = False,
    dry_run: bool = False,
    number_of_frames: int | None = None,
    mesh: tuple[int, int, int] | None = None,
    ntasks: int = DEFAULT_NTASKS,
    ncore: int | None = None,
    kpar: int | None = None,
    resume_from: Path | None = None,
    nodelist: str | None = None,
    memory: str = "96G",
    run_id: str | None = None,
):
    plan = _submission_plan(
        name, manager, rerun=rerun,
        number_of_frames=number_of_frames, mesh=mesh, ntasks=ntasks,
        ncore=ncore, kpar=kpar, resume_from=resume_from, nodelist=nodelist,
        memory=memory, run_id=run_id,
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
        # ``rerun`` starts a NEW run id (fresh store); ``resume`` retries a
        # failed/timed-out entry with the SAME run id so run_locally_custom can
        # reuse the DFPT volumes already on disk.
        retry=rerun or resume,
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
    parser.add_argument(
        "--batch", metavar="FILE", default=None,
        help="Submit every endpoint listed in FILE (one per line) as its own "
             "Slurm job; already-finished endpoints are skipped by the ledger",
    )
    parser.add_argument("--slurm", "-s", action="store_true", help="Submit to Slurm")
    parser.add_argument("--dry-run", action="store_true", help="Print the plan; no VASP")
    parser.add_argument("--flow-graph", action="store_true", help="Alias for --dry-run")
    parser.add_argument("--rerun", "-r", action="store_true",
                        help="Force a re-run under a NEW run id (fresh store)")
    parser.add_argument(
        "--resume", action="store_true",
        help="Retry a failed or timed-out Slurm entry with the SAME run id, so "
             "completed DFPT volumes are reused",
    )
    parser.add_argument("--check", "-c", action="store_true", help="Check the result JSON")
    parser.add_argument(
        "--refit-vt", action="store_true",
        help="Re-fit F(V,T) per temperature and write <name>-dfpt-qha-vt.json "
             "beside each analysis JSON; accepts -n or --batch, no VASP",
    )
    parser.add_argument("--run-id", help="Operate on an explicit immutable run")
    parser.add_argument(
        "--memory", default="96G", metavar="SIZE",
        help="Slurm --mem for the job (default 96G).  Not a scientific input; "
             "it is recorded in the config but, unlike ntasks/ncore/kpar, it is "
             "not part of the run id.",
    )
    parser.add_argument(
        "--nodelist", default=None, metavar="NODES",
        help="Pin the job to specific node(s), e.g. --nodelist 429e.  Resource "
             "only: it does not change the run id.  Needed here because the two "
             "compute nodes differ in memory (429pro 128 GB, 429e 192 GB) and a "
             "heavy FCC cell is OOM-killed on the smaller one.",
    )
    parser.add_argument(
        "--resume-from", type=Path, default=None, metavar="DIR",
        help="Directory holding the store.json to resume from, for a run whose "
             "scratch directory does not sit at the deterministic location "
             "(e.g. one left by an older submission).  Only meaningful with "
             "--resume.",
    )
    parser.add_argument(
        "--ntasks", type=int, default=DEFAULT_NTASKS,
        help=f"MPI ranks for the Slurm job (default {DEFAULT_NTASKS}); a "
             f"throughput choice, not a scientific one",
    )
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
    parser.add_argument(
        "--ncore", type=int, default=None,
        help="NCORE for the relaxation stages (default from RELAX_ELECTRONIC). "
             "Parallelisation only; it does not change the physics.  It also "
             "sets how many ranks share a band group, so NBANDS must be "
             "divisible by ntasks/(1*NCORE) for VASP not to adjust it.",
    )
    parser.add_argument(
        "--kpar", type=int, default=None,
        help="KPAR for the DFPT stages (default from DFPT_ELECTRONIC).  "
             "Parallelisation only.  NCORE and KPAR are mutually exclusive, so "
             "the two knobs are applied to different stages and never both to "
             "one INCAR.",
    )
    args = parser.parse_args()

    if args.resume_from is not None and not args.resume:
        parser.error("--resume-from only applies to --resume")

    if args.refit_vt:
        if args.slurm:
            parser.error(
                "--refit-vt is local post-processing; do not combine it with --slurm"
            )
        if args.name and args.batch:
            parser.error("--refit-vt takes either -n or --batch, not both")
        if args.batch:
            refit_vt(read_batch(Path(args.batch)))
        elif args.name:
            refit_vt([args.name])
        else:
            parser.error("--refit-vt needs -n/--name or --batch FILE")
        raise SystemExit(0)

    if args.batch:
        if args.name:
            parser.error("--batch and -n/--name are mutually exclusive")
        if not args.slurm:
            parser.error("--batch needs --slurm")
        if args.rerun and args.resume:
            parser.error("--rerun and --resume are mutually exclusive")
        submit_batch(
            Path(args.batch), SlurmJobManager(), rerun=args.rerun,
            resume=args.resume, dry_run=args.dry_run,
            number_of_frames=args.frames,
            mesh=tuple(args.mesh) if args.mesh else None,
            ntasks=args.ntasks, ncore=args.ncore, kpar=args.kpar,
            resume_from=args.resume_from, nodelist=args.nodelist,
            memory=args.memory,
        )
    elif args.slurm:
        if not args.name:
            parser.error("--slurm needs -n/--name")
        if args.rerun and args.resume:
            parser.error("--rerun and --resume are mutually exclusive")
        if args.run_id and (args.rerun or not args.resume):
            parser.error("--run-id with --slurm needs --resume: it exists to point "
                         "a resume back at an archived store")
        submit_job(
            args.name, SlurmJobManager(), rerun=args.rerun, resume=args.resume,
            dry_run=args.dry_run,
            number_of_frames=args.frames,
            mesh=tuple(args.mesh) if args.mesh else None,
            ntasks=args.ntasks, ncore=args.ncore, kpar=args.kpar,
            resume_from=args.resume_from, nodelist=args.nodelist,
            memory=args.memory, run_id=args.run_id,
        )
    elif args.dry_run or args.flow_graph:
        run_tick(
            args.name or DEFAULT_ENDPOINT, dry_run=True, run_id=args.run_id,
            number_of_frames=args.frames, mesh=tuple(args.mesh) if args.mesh else None,
            ncore=args.ncore, kpar=args.kpar,
        )
    elif args.check:
        check_jobs()
    elif args.name:
        run_tick(
            args.name, rerun=args.rerun, run_id=args.run_id,
            number_of_frames=args.frames, mesh=tuple(args.mesh) if args.mesh else None,
            ncore=args.ncore, kpar=args.kpar, resume_from=args.resume_from,
        )
    else:
        parser.print_help()
