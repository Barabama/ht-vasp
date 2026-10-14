#!/usr/bin/env python3
"""DFPT-QHA portal: dry-run/flow-graph only at this milestone.

This entry point deliberately refuses real execution and Slurm submission until
the DFPT volume flow has a reviewed terminal QHA consumer.  It is still useful
now: it validates one endpoint recipe and builds the complete deferred Flow
without launching VASP.
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

from pymatgen.core import Structure

from htvasp.dfpt_qha.maker import DfptQhaMaker
from htvasp.model.dfpt_qha_recipe import Recipe

log = logging.getLogger("main_dfpt_qha")
ROOT = Path(__file__).resolve().parent
POOL = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_recipes.json")
STRUCTURES = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_structures.json")


def load_endpoint(name: str) -> tuple[Structure, Recipe]:
    recipes = json.loads(POOL.read_text())["endpoints"]
    structures = json.loads(STRUCTURES.read_text())["endpoints"]
    if name not in recipes or name not in structures:
        raise KeyError(f"endpoint is absent from pool: {name}")
    recipe = Recipe.from_dict(recipes[name]["recipe"])
    structure = Structure.from_dict(structures[name]["structure"])
    return structure, recipe


def dry_run(name: str) -> None:
    structure, recipe = load_endpoint(name)
    problems = []
    # The structure/recipe cross-check remains a model gate before flow build.
    # The compact formal model currently exposes the recipe shape contract;
    # geometry-specific checking stays in the legacy audit until migrated.
    if len(structure) != recipe.kpoints.cell_atoms:
        problems.append("recipe cell_atoms disagrees with the pool structure")
    if problems:
        raise SystemExit(f"refusing flow build for {name}: {'; '.join(problems)}")
    flow = DfptQhaMaker(recipe=recipe).make(
        structure,
        supercell_matrix=[
            [recipe.supercell[0], 0, 0],
            [0, recipe.supercell[1], 0],
            [0, 0, recipe.supercell[2]],
        ],
    )
    log.info("endpoint=%s formula=%s atoms=%d", name, structure.formula, len(structure))
    log.info("recipe=%s", recipe.schema_version)
    log.info("flow=%s jobs=%d", flow.name, len(flow.jobs))
    for index, job in enumerate(flow.jobs, 1):
        log.info("  %d. %s", index, job.name)
    log.info("dry-run OK: no VASP or Slurm launched")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--name", "-n", required=True, help="pool endpoint name")
    parser.add_argument("--dry-run", action="store_true", help="build and print the flow")
    parser.add_argument("--flow-graph", action="store_true", help="alias for --dry-run")
    parser.add_argument("--slurm", action="store_true", help="not enabled at this milestone")
    parser.add_argument("--run", action="store_true", help="not enabled at this milestone")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
    if args.slurm or args.run:
        parser.error("real DFPT-QHA execution is intentionally disabled until the terminal QHA consumer is reviewed")
    if not (args.dry_run or args.flow_graph):
        parser.error("use --dry-run or --flow-graph; real execution is not enabled")
    dry_run(args.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
