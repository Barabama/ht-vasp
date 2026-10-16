"""DFPT jobflow boundaries.

The collector consumes the output of an atomate2 ``BaseVaspMaker`` job.  It does
not manufacture a finite-displacement dataset: DFPT's ``<dynmat>`` is read
from that job's own ``vasprun.xml``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from atomate2.common.flows.eos import CommonEosMaker
from jobflow import Flow, Response, job

from htvasp.dfpt_qha.jobs import DfptVaspMaker, collect_dfpt_task_doc, validate_dfpt_task_doc
from htvasp.model.dfpt_qha_recipe import Recipe


@job

def spawn_dfpt_volume_jobs(
    eos_output: dict[str, Any],
    *,
    recipe: dict[str, Any],
    run_vasp_kwargs: dict[str, Any] | None = None,
) -> Response:
    """Replace one deferred job with one DFPT flow per EOS volume.

    The EOS contract is the same dict-of-lists emitted by atomate2's
    ``CommonEosMaker``.  Each volume is passed to a fresh ``BaseVaspMaker``
    job; the collector then reads that job's TaskDoc directory.  No phonopy
    displacement structures are generated here.
    """

    structures = eos_output.get("relax", {}).get("structure", [])
    directories = eos_output.get("relax", {}).get("dir_name", [])
    if len(structures) != len(directories):
        raise ValueError("EOS structures and dir_name arrays have different lengths")
    parsed = Recipe.from_dict(recipe)
    supercell = tuple(int(x) for x in parsed.supercell)
    jobs = []
    outputs = []
    for index, (structure, directory) in enumerate(zip(structures, directories)):
        # DFPT must run on the recipe's phonon supercell.  Running it on the
        # 2-atom cell yields a Gamma-only dynamical matrix, which cannot see
        # the finite-q instability the recipe exists to characterise.
        supercell_structure = (
            structure if supercell == (1, 1, 1) else structure.make_supercell(supercell)
        )
        maker = DfptVaspMaker(
            name=f"DFPT volume {index + 1}/{len(structures)}",
            recipe=parsed,
            run_vasp_kwargs=dict(run_vasp_kwargs or {}),
        )
        vasp_job = maker.make(supercell_structure, prev_dir=directory)
        collector = collect_dfpt_task_doc(
            vasp_job.output,
            recipe=recipe,
        )
        collector.name = f"collect DFPT volume {index + 1}/{len(structures)}"
        validation = validate_dfpt_task_doc(collector.output, recipe=recipe)
        validation.name = f"validate DFPT volume {index + 1}/{len(structures)}"
        jobs.extend([vasp_job, collector, validation])
        outputs.append(validation.output)
    return Response(replace=Flow(jobs, output=outputs, name="DFPT volume jobs"))


__all__ = ["spawn_dfpt_volume_jobs"]
