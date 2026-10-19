#!/usr/bin/env python3
"""Spin screening for endpoints whose magnetic policy is unresolved.

The pool refuses to generate inputs for a Cr-bearing endpoint whose spin state
was never measured (``magnetic_policy = 'unresolved'``).  That guard is correct:
a zero or untested start converges to the non-magnetic solution, which for
B2 AlCr is 283.7 meV/atom above the magnetic ground state, and the whole QHA
would then describe the wrong electronic state.

This portal measures the branch instead of guessing it.  For each endpoint it
runs the same short relaxation from two or three spin starts -- non-magnetic,
ferromagnetic, and (where the cell holds at least two Cr) a staggered
antiferromagnetic pattern -- and writes the comparison to
``data/endpoints/<name>/spin-screen.json``.  Updating ``pool_recipes.json``
from that file is a separate, explicit step.

The recipe's policy is forced to ``confirmed_branch`` for the screen itself:
this script exists precisely to bypass the guard while measuring it.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Any

from jobflow import Flow, JobStore, job
from monty.json import jsanitize
from maggma.stores import JSONStore, MemoryStore
from pymatgen.core import Structure

from htvasp.dfpt_qha.jobs import DfptRelaxMaker
from htvasp.model.dfpt_qha_recipe import Recipe
from htvasp.slurm import SlurmJobManager
from htvasp.utils.vasp_psp import psp_dir_exports
from htvasp.utils.local import run_locally_custom
from htvasp.utils.vasp_handlers import with_handlers

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

root_dir = Path(__file__).parent.resolve()
flow_dir = Path("/nfs_ssd/tmp")
conda_env = root_dir / ".conda"

POOL_DIR = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0")
POOL_RECIPES = POOL_DIR / "pool_recipes.json"
POOL_STRUCTURES = POOL_DIR / "pool_structures.json"

VASP_ARGS = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-cpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-cpu && srun vasp_gam'",
}

#: Same 16-element table as main_qha.py / main_dfpt_qha.py.
MAGMOM = {
    "Al": 1.0, "Co": 3.0, "Cr": 3.0, "Cu": 5.0, "Fe": 5.0, "Hf": 1.0,
    "Mn": 5.0, "Mo": 3.0, "Nb": 3.0, "Ni": 2.0, "Re": 1.0, "Ta": 3.0,
    "Ti": 3.0, "V": 3.0, "W": 3.0, "Zr": 1.0,
}

#: A branch must beat the non-magnetic one by this much at the relaxed
#: geometry.  Below it the two are indistinguishable at the accuracy of the
#: recipe (k-mesh density contributes up to ~6 meV/atom, Stage 0A).
DECISION_THRESHOLD_MEV_PER_ATOM = 5.0

#: A state whose relaxed lattice deviates from the non-magnetic one by more
#: than this is a different structure, not a spin arrangement of the same one.
STRUCTURE_TOLERANCE = 0.02

# The screen uses the pipeline's own relax INCAR (IBRION=2, ISIF=3, NSW=80,
# EDIFFG=-0.02).  NSW is a cap, not a cost: the relaxation stops at EDIFFG, and
# reusing the production settings means a branch that wins here was measured
# with the same electronic parameters the DFPT run will use.

DEFAULT_NTASKS = 32


def load_endpoint(name: str) -> tuple[Structure, Recipe]:
    """The endpoint's cell and recipe from the pool (see main_dfpt_qha.load_endpoint)."""

    if not POOL_RECIPES.is_file():
        raise FileNotFoundError(f"missing pool recipes: {POOL_RECIPES}")
    recipes = json.loads(POOL_RECIPES.read_text())["endpoints"]
    structures = json.loads(POOL_STRUCTURES.read_text())["endpoints"]
    if name not in recipes or name not in structures:
        raise KeyError(f"endpoint {name!r} is absent from the pool")
    return (
        Structure.from_dict(structures[name]["structure"]),
        Recipe.from_dict(recipes[name]["recipe"]),
    )


def unresolved_endpoints() -> list[str]:
    recipes = json.loads(POOL_RECIPES.read_text())["endpoints"]
    return sorted(
        name for name, entry in recipes.items()
        if entry["recipe"]["magnetic_policy"] == "unresolved"
    )


def spin_states(structure: Structure) -> dict[str, Any]:
    """Initial moments for each state that makes sense for this cell."""

    species = [str(site.specie) for site in structure]
    states: dict[str, Any] = {
        # non-magnetic: an explicit zero, not an absent key
        "nm": {element: 0.0 for element in dict.fromkeys(species)},
        "fm": {element: MAGMOM[element] for element in dict.fromkeys(species)},
    }
    cr_sites = [index for index, element in enumerate(species) if element == "Cr"]
    if len(cr_sites) >= 2:
        # A staggered pattern by site order.  It is a screen, not a search over
        # antiferromagnetic orderings; the sign pattern is recorded so a
        # decision that rests on it can be re-examined.
        states["afm"] = [
            (3.0 if position % 2 == 0 else -3.0) if element == "Cr" else 0.0
            for position, element in enumerate(species)
        ]
    return states


def _taskdoc_fields(document: Any) -> tuple[
    float | None, int, list[float], list[float]
]:
    """Energy, atom count and site moments from an atomate2 TaskDoc.

    The site moments are not in ``TaskDoc.output``; emmet keeps them under
    ``calcs_reversed[-1].output.outcar.magnetization`` as per-site s/p/d/tot
    blocks.  A screening that read ``output.magnetic_moments`` would silently
    report ``None`` for every state.
    """

    # jobflow stores a job's output as JSON and monty's decoder rebuilds the
    # objects it recognises, so what arrives here can be a live emmet TaskDoc
    # rather than a dict.  ``.get`` on that raises AttributeError.
    doc = document or {}
    if not isinstance(doc, dict):
        doc = jsanitize(doc, strict=True, allow_bson=True, recursive_msonable=True)
    if not isinstance(doc, dict):
        raise ValueError(f"expected a TaskDoc mapping, got {type(document).__name__}")
    task_output = doc.get("output") or {}
    if not isinstance(task_output, dict):
        task_output = jsanitize(
            task_output, strict=True, allow_bson=True, recursive_msonable=True
        )
    output = task_output.get("output") or {}
    energy = output.get("energy", task_output.get("energy"))
    structure = output.get("structure") or doc.get("structure") or {}
    n_atoms = len(structure.get("sites") or []) or int(doc.get("nsites") or 1)
    lattice = structure.get("lattice") or {}
    # pymatgen serialises the cell as a/b/c, not as "abc"
    abc = list(lattice.get("abc") or [])
    if not abc:
        abc = [lattice[key] for key in ("a", "b", "c") if lattice.get(key) is not None]
    moments: list[float] = []
    # calcs_reversed sits at the TaskDoc level, not inside its output block
    calcs = doc.get("calcs_reversed") or task_output.get("calcs_reversed") or []
    if calcs:
        outcar = (calcs[-1].get("output") or {}).get("outcar") or {}
        moments = [
            float(entry.get("tot", 0.0))
            for entry in outcar.get("magnetization") or []
        ]
    return energy, n_atoms, moments, abc


@job
def summarize_states(
    endpoint: str, documents: dict[str, Any], *, threshold_mev_per_atom: float
) -> dict[str, Any]:
    rows: dict[str, Any] = {}
    for label, document in documents.items():
        energy, n_atoms, moments, abc = _taskdoc_fields(document)
        rows[label] = {
            "energy_eV": energy,
            "energy_eV_per_atom": None if energy is None else energy / n_atoms,
            "total_moment_muB": None if not moments else float(sum(moments)),
            "max_site_moment_muB": None if not moments else float(
                max(abs(value) for value in moments)
            ),
            "lattice_abc": abc,
        }

    # A starter spin state can drive a structural collapse instead of a
    # magnetic solution.  The staggered pattern in an L1_2 cell does exactly
    # that: the cell goes cubic -> tetragonal by ~20%, the moments vanish, and
    # the energy drops by hundreds of meV/atom.  That is a different phase, not
    # a magnetic branch, so it must not win the comparison.
    reference = rows.get("nm", {}).get("lattice_abc") or []
    for row in rows.values():
        abc = row.get("lattice_abc") or []
        if reference and abc and len(abc) == len(reference):
            deviation = max(
                abs(value - base) / base
                for value, base in zip(abc, reference) if base
            )
            row["lattice_deviation_vs_nm"] = float(deviation)
            row["structure_changed"] = bool(deviation > STRUCTURE_TOLERANCE)
        else:
            row["lattice_deviation_vs_nm"] = None
            row["structure_changed"] = False
    admissible = {label: row for label, row in rows.items() if not row["structure_changed"]}
    excluded = sorted(set(rows) - set(admissible))
    if rows.get("nm", {}).get("energy_eV_per_atom") is None or "nm" in excluded:
        return {
            "schema_version": "htvasp-spin-screen/v1",
            "endpoint": endpoint,
            "ok": False,
            "error": "the non-magnetic reference relaxation is missing or changed structure",
            "states": rows,
            "excluded_states": excluded,
        }
    nm = rows["nm"]["energy_eV_per_atom"]
    deltas = {
        label: 1e3 * (row["energy_eV_per_atom"] - nm)
        for label, row in admissible.items()
        if label != "nm" and row["energy_eV_per_atom"] is not None
    }
    best = min(admissible, key=lambda label: admissible[label]["energy_eV_per_atom"])
    magnetic_gain = -min(deltas.values()) if deltas else 0.0
    decision = (
        "confirmed_branch"
        if best != "nm" and magnetic_gain >= threshold_mev_per_atom
        else "nonmagnetic"
    )
    return {
        "schema_version": "htvasp-spin-screen/v1",
        "endpoint": endpoint,
        "ok": True,
        "states": rows,
        "delta_vs_nm_meV_per_atom": deltas,
        "best_state": best,
        "excluded_states": excluded,
        "structure_tolerance": STRUCTURE_TOLERANCE,
        "magnetic_gain_meV_per_atom": magnetic_gain,
        "threshold_meV_per_atom": threshold_mev_per_atom,
        "decision": decision,
        "note": (
            "decision is a screening outcome; propagating it into "
            "pool_recipes.json is an explicit separate step"
        ),
    }


def build_flow(name: str, structure: Structure, recipe: Recipe) -> Flow:
    screen_recipe = Recipe.from_dict(
        {**recipe.model_dump(mode="json"), "magnetic_policy": "confirmed_branch"}
    )
    jobs = []
    references: dict[str, Any] = {}
    for label, magmom in spin_states(structure).items():
        maker = DfptRelaxMaker(
            name=f"spin {label}",
            recipe=screen_recipe,
            isif=3,
            magmom=magmom,
            run_vasp_kwargs=with_handlers(run_vasp_kwargs=VASP_ARGS),
        )
        job = maker.make(structure)
        job.name = f"spin {label}"
        jobs.append(job)
        references[label] = job.output
    summary = summarize_states(
        name, references, threshold_mev_per_atom=DECISION_THRESHOLD_MEV_PER_ATOM
    )
    summary.name = f"summarize spin states: {name}"
    jobs.append(summary)
    return Flow(jobs, output=summary.output, name=f"spin screen {name}")


def run_screen(names: list[str], rerun: bool = False, dry_run: bool = False) -> None:
    for index, name in enumerate(names, 1):
        log.info("[%d/%d] spin screen %s", index, len(names), name)
        structure, recipe = load_endpoint(name)
        flow = build_flow(name, structure, recipe)
        if dry_run:
            log.info("  states: %s", ", ".join(
                job.name for job in flow.jobs if job.name.startswith("spin ")
            ))
            continue
        run_root = root_dir / "data" / "endpoints" / name / "spin-screenflow"
        # read_only=False lets maggma create an empty store.json; the default
        # read_only=True raises FileNotFoundError on a first run.
        store = JobStore(
            JSONStore(str(run_root / "store.json"), read_only=False),
            additional_stores={"data": MemoryStore()},
        )
        try:
            run_locally_custom(
                flow, store=store, root_dir=run_root, resume=not rerun
            )
            output = store.get_output(
                uuid=flow.output.uuid, which=flow.output.index, load=True
            )
            # jobflow returns a one-element list for a single-index lookup
            if isinstance(output, list):
                output = output[0] if len(output) == 1 else output
            if not isinstance(output, dict):
                raise ValueError(
                    f"the summariser returned {type(output).__name__}, expected a dict"
                )
        except Exception as exc:  # noqa: BLE001 - one endpoint must not stop the screen
            log.error("  %s: screen failed: %s", name, exc)
            output = {
                "schema_version": "htvasp-spin-screen/v1",
                "endpoint": name,
                "ok": False,
                "error": f"{type(exc).__name__}: {exc}",
            }
        finally:
            store.close()
        path = root_dir / "data" / "endpoints" / name / "spin-screen.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(output, indent=2, ensure_ascii=False, default=str))
        log.info("  -> %s  decision=%s  gain=%.1f meV/atom",
                 path.name, output.get("decision"), output.get("magnetic_gain_meV_per_atom", 0.0))


def reanalyse(names: list[str]) -> None:
    """Recompute the verdict from stored relaxations.  No VASP, no new jobs.

    Used when the decision logic changes: the TaskDocs are already in the
    stores, so a corrected rule costs seconds instead of a re-run.
    """

    store = JobStore(MemoryStore())
    store.connect()
    for name in names:
        path = root_dir / "data" / "endpoints" / name / "spin-screenflow" / "store.json"
        if not path.is_file():
            log.warning("%s: no screen store at %s", name, path)
            continue
        documents: dict[str, Any] = {}
        for document in json.loads(path.read_text()):
            job_name = document.get("name") or ""
            if job_name.startswith("spin "):
                documents[job_name.split()[1]] = document.get("output")
        output = summarize_states(
            name, documents, threshold_mev_per_atom=DECISION_THRESHOLD_MEV_PER_ATOM
        ).run(store).output
        target = root_dir / "data" / "endpoints" / name / "spin-screen.json"
        target.write_text(json.dumps(output, indent=2, ensure_ascii=False, default=str))
        log.info(
            "  %-14s %-17s gain=%6.1f meV/atom excluded=%s",
            name, output.get("decision"),
            output.get("magnetic_gain_meV_per_atom", 0.0),
            output.get("excluded_states") or [],
        )


def submit_screen(names: list[str], ntasks: int = DEFAULT_NTASKS, dry_run: bool = False) -> None:
    manager = SlurmJobManager()
    config = manager.get_cpu_config(
        job_name="spin-screen",
        output_log="logs/spin-screen.log",
        ntasks=int(ntasks),
        memory="96G",
        time_limit="24:00:00",
        module_name="vasp-cpu",
        extra_commands=psp_dir_exports(),
    )
    command = (
        f"{conda_env / 'bin' / 'python'} {Path(__file__).resolve()} "
        f"--batch {' '.join(names)}"
    )
    if dry_run:
        print(json.dumps({"command": command, "names": names}, indent=2))
        return
    job_id = manager.submit_command(
        command=command, config=config, conda_env=str(conda_env),
        module_name="vasp-cpu", work_dir=str(root_dir),
    )
    log.info("spin screen submitted: %s", job_id)


def main() -> None:
    parser = argparse.ArgumentParser(description="Spin screening for unresolved endpoints")
    parser.add_argument("--batch", nargs="*", default=None,
                        help="endpoint names (default: every unresolved endpoint)")
    parser.add_argument("--slurm", action="store_true", help="submit the whole screen as one job")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--rerun", action="store_true")
    parser.add_argument(
        "--reanalyse", action="store_true",
        help="recompute the verdict from the stored relaxations (no VASP)",
    )
    parser.add_argument("--ntasks", type=int, default=DEFAULT_NTASKS)
    args = parser.parse_args()

    names = args.batch or unresolved_endpoints()
    if not names:
        log.info("no unresolved endpoints; nothing to do")
        return
    log.info("%d endpoint(s): %s", len(names), ", ".join(names))
    if args.reanalyse:
        reanalyse(names)
    elif args.slurm:
        submit_screen(names, ntasks=args.ntasks, dry_run=args.dry_run)
    else:
        run_screen(names, rerun=args.rerun, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
