"""DFPT jobflow boundaries.

The collector consumes the output of an atomate2 ``BaseVaspMaker`` job.  It does
not manufacture a finite-displacement dataset: DFPT's ``<dynmat>`` is read from
that job's own ``vasprun.xml``.

Flow shape, per endpoint::

    CommonEosMaker                 E(V) and the relaxed cell at each volume
      -> spawn_dfpt_volume_jobs    one DFPT run per volume, on the recipe supercell
           -> collect_dfpt_task_doc    vasprun <dynmat>/<hessian> -> DfptTaskDoc
           -> validate_dfpt_task_doc   fail-closed: cell, shape, provenance
      -> collect_validated_volumes  gather the survivors
      -> classify_dfpt_stability    one verdict per volume; gates the endpoint
      -> dfpt_thermal_properties    F(V,T) from the dynamical matrices
      -> analyze_dfpt_qha           G(T), V(T) and the EOS fit

The validation job is a pass-through, so a stopped volume stops the thermal and
QHA jobs behind it.  Stability is one job after collection, not one per
volume: the decision belongs to the endpoint (a soft volume fails the endpoint
rather than being dropped from the fit), and the VASP time is already spent by
then, so a per-volume fail-fast would buy nothing.  The classification's own
output is the tier-B artifact a halted flow reports through
``DfptQhaWorker.get_result``.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from jobflow import Flow, Response, job
from pymatgen.core import Structure

from htvasp.dfpt_qha.jobs import (
    DfptVaspMaker,
    collect_dfpt_task_doc,
    validate_dfpt_task_doc,
)
from htvasp.dfpt_qha.phonons import (
    as_structure,
    assert_same_ordering,
    phonopy_for,
    supercell_structure,
)
from htvasp.dfpt_qha.stability import (
    DfptStabilityReport,
    diagnose_stability,
    unresolved_report,
)
from htvasp.dfpt_qha.qha import (
    DEFAULT_MESH,
    DEFAULT_TEMPERATURES,
    analyze_dfpt_qha,
    dfpt_thermal_properties,
)
from htvasp.model.dfpt_qha_recipe import Recipe


@job
def spawn_dfpt_volume_jobs(
    eos_output: dict[str, Any],
    *,
    recipe: dict[str, Any],
    run_vasp_kwargs: dict[str, Any] | None = None,
    supercell_matrix: Any = None,
    mesh: tuple[int, int, int] = DEFAULT_MESH,
    temperatures: list[float] | None = None,
    pressure_gpa: float = 0.0,
) -> Response:
    """Replace one deferred job with the per-volume DFPT chain plus its consumers.

    Each volume gets its own ``BaseVaspMaker`` job; nothing is relaxed here,
    because the EOS stage already produced the relaxed cell at that volume.
    """

    relax = eos_output.get("relax", {})
    structures = relax.get("structure", [])
    directories = relax.get("dir_name", [])
    energies = relax.get("energy", [])
    volumes = relax.get("volume", [])
    if not (len(structures) == len(directories) == len(energies) == len(volumes)):
        raise ValueError(
            "EOS output arrays disagree in length: "
            f"structure={len(structures)} dir_name={len(directories)} "
            f"energy={len(energies)} volume={len(volumes)}"
        )

    parsed = Recipe.from_dict(recipe)
    supercell = tuple(int(x) for x in parsed.supercell)
    if supercell_matrix is None:
        supercell_matrix = [
            [supercell[0], 0, 0], [0, supercell[1], 0], [0, 0, supercell[2]]
        ]
    temperatures = list(temperatures or DEFAULT_TEMPERATURES)
    vasp_kwargs = dict(run_vasp_kwargs or {})

    jobs: list[Any] = []
    validations: list[Any] = []
    cells: list[dict[str, Any]] = []
    for index, (structure, directory) in enumerate(zip(structures, directories)):
        # ``structure`` may be a dict or a reconstructed Structure depending on
        # how the EOS job's output was stored; the rest of this loop needs the
        # object either way.
        structure = as_structure(structure)
        # The DFPT job must run on phonopy's OWN supercell.  pymatgen's
        # make_supercell produces the same crystal in a different atom order,
        # and the force constants that come back are indexed by whatever order
        # VASP read -- attaching them across the two orderings still gives
        # plausible frequencies, so the mismatch has to be designed out rather
        # than checked for later.  See htvasp.dfpt_qha.phonons.
        cell = supercell_structure(structure, supercell_matrix)
        maker = DfptVaspMaker(
            name=f"DFPT volume {index + 1}/{len(structures)}",
            recipe=parsed,
            run_vasp_kwargs=vasp_kwargs,
        )
        vasp_job = maker.make(cell, prev_dir=directory)
        collector = collect_dfpt_task_doc(
            vasp_job.output,
            recipe=recipe,
            cell=structure.as_dict(),
            supercell_matrix=supercell_matrix,
            # The DFPT run recomputes the ground state on the supercell, so its
            # energy per atom must match the static energy of this same volume.
            # A diverged SCF that still stopped early would otherwise pass the
            # TaskDoc state check.
            expected_static_energy_per_atom=float(energies[index]) / len(structure),
        )
        collector.name = f"collect DFPT volume {index + 1}/{len(structures)}"
        validation = validate_dfpt_task_doc(collector.output, recipe=recipe)
        validation.name = f"validate DFPT volume {index + 1}/{len(structures)}"
        jobs.extend([vasp_job, collector, validation])
        validations.append(validation.output)
        cells.append(
            {
                "index": index,
                "cell": structure.as_dict(),
                "supercell": cell.as_dict(),
                "volume_angstrom3": float(volumes[index]),
                "static_energy_eV": float(energies[index]),
            }
        )

    collected = collect_validated_volumes(validations, cells=cells)
    collected.name = "collect validated volumes"
    stability = classify_dfpt_stability(
        collected.output,
        recipe=recipe,
        supercell_matrix=supercell_matrix,
        mesh=tuple(int(x) for x in mesh),
    )
    stability.name = "classify DFPT stability"
    thermal = dfpt_thermal_properties(
        stability.output,
        supercell_matrix=supercell_matrix,
        mesh=mesh,
        temperatures=temperatures,
    )
    thermal.name = "DFPT thermal properties"
    analysis = analyze_dfpt_qha(thermal.output, pressure_gpa=pressure_gpa)
    analysis.name = "analyze DFPT-QHA"

    jobs.extend([collected, stability, thermal, analysis])
    return Response(replace=Flow(jobs, output=analysis.output, name="DFPT volume jobs"))


@job
def classify_dfpt_stability(
    collected: dict[str, Any],
    *,
    recipe: dict[str, Any],
    supercell_matrix: Any,
    mesh: tuple[int, int, int] = DEFAULT_MESH,
) -> Response:
    """One stability verdict per validated volume, and the endpoint decision.

    The gate is here rather than per volume for two reasons.  A soft volume
    fails the *endpoint*: dropping it and fitting the rest leaves a hole in
    F(V,T) that the analysis cannot see, and the result still looks smooth.
    And the VASP time is already spent by this point, so a per-volume fail-fast
    would save nothing while multiplying the serialization boundaries.

    Each verdict is a stored artifact.  When the endpoint fails, the analysis
    job never exists, so this output is what ``DfptQhaWorker.get_result``
    returns -- including the exact-q frequencies and eigenvector subspace.
    """

    mesh = tuple(int(x) for x in mesh)
    reports: list[dict[str, Any]] = []
    volumes: list[dict[str, Any]] = []
    for record in collected.get("volumes") or []:
        index = int(record["index"])
        volume = float(record["volume_angstrom3"])
        force_constants = np.array(record["force_constants"], dtype=float)
        try:
            cell = as_structure(record["cell"])
            supercell = as_structure(record["supercell"])
            phonon = phonopy_for(cell, supercell_matrix)
            # The matrix is indexed for whatever order VASP read; attaching it
            # to a differently-ordered supercell still gives plausible
            # frequencies at Gamma, so this is asserted, not assumed.
            assert_same_ordering(
                supercell, supercell_structure(cell, supercell_matrix)
            )
            phonon.force_constants = force_constants
            report = diagnose_stability(
                phonon,
                force_constants,
                structure=supercell,
                recipe=recipe,
                supercell_matrix=supercell_matrix,
                mesh=mesh,
                ordering_ok=True,
                index=index,
                volume_angstrom3=volume,
            )
        except Exception as exc:  # noqa: BLE001 - the report is the failure record
            report = unresolved_report(
                f"stability diagnosis could not be completed: {exc}",
                force_constants=force_constants,
                structure_hash=record.get("ordering_fingerprint", ""),
                recipe_hash="",
                index=index,
                volume_angstrom3=volume,
            )
        verdict = DfptStabilityReport.model_validate(report).model_dump(mode="json")
        reports.append(verdict)
        volumes.append({**record, "stability": verdict})

    blocked = [r for r in reports if not r["qha_allowed"]]
    confirmation = [r for r in reports if r.get("confirmation_required")]
    output = {
        **{key: value for key, value in collected.items() if key != "volumes"},
        "schema_version": "htvasp-dfpt-stability/v1",
        "gate": "dfpt_stability",
        "volumes": volumes,
        "reports": reports,
        "n_volumes": len(reports),
        "n_blocked": len(blocked),
        "blocked_volume_indices": [r["index"] for r in blocked],
        "confirmation_volume_indices": [r["index"] for r in confirmation],
        "qha_allowed": not blocked,
        "mesh": list(mesh),
        "supercell_matrix": [
            [int(x) for x in row] for row in np.asarray(supercell_matrix, dtype=int)
        ],
        "supercell_convergence": "not_checked",
        "warnings": [
            warning for r in reports for warning in r.get("warnings") or []
        ],
    }
    if blocked:
        message = (
            f"{len(blocked)} of {len(reports)} volumes are not "
            f"harmonic-QHA-applicable: "
            + "; ".join(
                f"volume {r['index']} {r['status']} ({r['reason']})"
                for r in blocked
            )
        )
        output["error"] = message
        output["ok"] = False
    else:
        output["ok"] = True
    return Response(output=output, stop_children=bool(blocked))


@job
def collect_validated_volumes(
    validations: list[dict[str, Any]],
    *,
    cells: list[dict[str, Any]],
) -> dict[str, Any]:
    """Gather the volumes that passed validation, with their static energies.

    A volume whose validation returned ``ok=False`` is dropped here rather than
    carrying an unvalidated dynamical matrix into the fit.  In a real run the
    stopped job halts the flow before this point; this check is what keeps the
    contract true when a flow is assembled or replayed by hand.
    """

    kept: list[dict[str, Any]] = []
    rejected: list[int] = []
    for index, validation in enumerate(validations):
        if not validation or not validation.get("ok"):
            rejected.append(index)
            continue
        document = validation.get("dfpt_doc") or {}
        kept.append(
            {
                **cells[index],
                "force_constants": document["force_constants"],
                "n_atoms_supercell": document["force_constants_provenance"]["n_atoms"],
                "ordering_fingerprint": document["force_constants_provenance"][
                    "ordering_fingerprint"
                ],
                "force_constants_source": document["force_constants_provenance"][
                    "source"
                ],
                # custodian's corrections are recorded by the collector, and
                # they have to survive into the artifact: the per-volume
                # collector output is not part of the final JSON.
                "input_audit": (document.get("provenance") or {}).get("input_audit"),
            }
        )
    return {
        "volumes": kept,
        "rejected_volume_indices": rejected,
        "n_requested": len(validations),
        "n_validated": len(kept),
    }


__all__ = [
    "classify_dfpt_stability",
    "collect_validated_volumes",
    "spawn_dfpt_volume_jobs",
]
