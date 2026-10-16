"""DFPT jobflow boundaries.

The collector consumes the output of an atomate2 ``BaseVaspMaker`` job.  It does
not manufacture a finite-displacement dataset: the dynamical matrix is read from
that job's own ``DYNMAT``.

Flow shape, per endpoint::

    CommonEosMaker                 E(V) and the relaxed cell at each volume
      -> spawn_dfpt_volume_jobs    one DFPT run per volume, on the recipe supercell
           -> collect_dfpt_task_doc    DYNMAT -> force constants -> DfptTaskDoc
           -> validate_dfpt_task_doc   fail-closed: cell, shape, provenance
      -> collect_validated_volumes  gather the survivors
      -> dfpt_thermal_properties    F(V,T) from the dynamical matrices
      -> analyze_dfpt_qha           G(T), V(T) and the EOS fit

The validation job is a pass-through, so a stopped volume stops the thermal and
QHA jobs behind it -- which is the whole point of the gate.

What VASP actually writes
-------------------------
``IBRION=8`` writes **DYNMAT**, not ``FORCE_CONSTANTS``.  The verified reference
directory has both only because ``FORCE_CONSTANTS`` was produced 90 minutes
after the run by an offline conversion; assuming VASP writes it means the
collector finds nothing after the allocation has already been spent.

``DYNMAT`` holds, for each perturbation ``(atom, direction)``, the force
response on every atom, in eV/Angstrom^2 because the perturbation is a unit
displacement.  That is the same quantity phonopy calls force constants, so the
matrix is recovered by replaying the perturbations through phonopy with
``distance=1.0`` and attaching the responses.  No reordering is needed here
because the DFPT job runs on phonopy's OWN supercell (see
``htvasp.dfpt_qha.phonons``), which is what makes VASP's atom order and
phonopy's index order the same list.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from atomate2.vasp.jobs.base import BaseVaspMaker
from atomate2.utils.path import strip_hostname
from jobflow import Response, job
from monty.io import zopen
from monty.json import jsanitize
from pymatgen.core import Structure

from htvasp.dfpt_qha.input_set import DfptInputSetGenerator, DfptRelaxSetGenerator
from htvasp.dfpt_qha.phonons import (
    assert_same_ordering,
    ordering_fingerprint,
    phonopy_for,
    supercell_structure,
)
from htvasp.dfpt_qha.task_doc import DfptForceConstantsProvenance, DfptTaskDoc
from htvasp.model.dfpt_qha_recipe import Recipe

#: 1 eV/Angstrom^2 expressed in GPa (for the bulk modulus, if needed).
DYNMAT_UNIT = "eV/Angstrom^2"


def _find_file(directory: Path, name: str) -> Path | None:
    """``name`` or ``name.gz`` inside ``directory``, whichever exists.

    BaseVaspMaker gzips a fixed file list at the end of every job, so a
    collector that only looks for the plain name breaks the moment that list
    changes.  Both spellings are accepted.
    """

    for candidate in (directory / name, directory / f"{name}.gz"):
        if candidate.is_file():
            return candidate
    return None


def parse_dynmat(path: Path | str) -> tuple[list[dict[str, Any]], int]:
    """Read VASP's ``DYNMAT``: the DFPT force response per perturbation.

    Layout::

        line 0   n_species  n_atoms  n_perturbations
        line 1   masses
        then, per perturbation:
            atom direction  displacement(3)
            n_atoms lines of (force_x force_y force_z)

    Returns ``(perturbations, n_atoms)`` with each perturbation as
    ``{"atom": 1-based index, "direction": 1-based, "displacement": [3],
    "forces": (n_atoms, 3) array}``.
    """

    with zopen(path, mode="rt") as handle:
        lines = handle.read().splitlines()
    if len(lines) < 3:
        raise ValueError(f"DYNMAT is too short to parse: {path}")

    header = lines[0].split()
    n_atoms = int(header[1])
    n_perturbations = int(header[2])
    if n_atoms < 1 or n_perturbations < 1:
        raise ValueError(f"DYNMAT header is not usable: {lines[0]!r}")

    perturbations: list[dict[str, Any]] = []
    index = 2
    for _ in range(n_perturbations):
        if index >= len(lines):
            raise ValueError(f"DYNMAT ended after {len(perturbations)} perturbations")
        head = lines[index].split()
        atom, direction = int(head[0]), int(head[1])
        displacement = [float(x) for x in head[2:5]] if len(head) >= 5 else [0.0, 0.0, 0.0]
        index += 1
        block = lines[index:index + n_atoms]
        if len(block) != n_atoms:
            raise ValueError(
                f"DYNMAT perturbation {len(perturbations) + 1} has {len(block)} "
                f"force lines, expected {n_atoms}"
            )
        forces = np.array([[float(x) for x in row.split()[:3]] for row in block], dtype=float)
        index += n_atoms
        perturbations.append(
            {
                "atom": atom,
                "direction": direction,
                "displacement": displacement,
                "forces": forces,
            }
        )
    return perturbations, n_atoms


def force_constants_from_dynmat(
    cell: Structure,
    supercell: Structure,
    perturbations: list[dict[str, Any]],
    *,
    supercell_matrix: Any,
) -> np.ndarray:
    """Force constants from a DFPT ``DYNMAT``, in the supercell's atom order.

    phonopy is asked to generate its own perturbation list at unit displacement
    and then handed the DFPT responses in that same order.  The VASP order and
    phonopy's order are the same because the job ran on
    ``phonopy.supercell``; that is asserted rather than assumed.
    """

    phonon = phonopy_for(cell, supercell_matrix)
    phonopy_supercell = supercell_structure(cell, supercell_matrix)
    assert_same_ordering(supercell, phonopy_supercell)

    n_supercell = len(phonopy_supercell)
    phonon.generate_displacements(distance=1.0)
    expected = len(phonon.supercells_with_displacements)
    if expected != len(perturbations):
        raise ValueError(
            f"DYNMAT has {len(perturbations)} perturbations but phonopy expects "
            f"{expected} at unit displacement; the two are not the same "
            f"calculation and the responses cannot be attached"
        )
    forces = []
    # phonopy reports displacements as [atom_index, dx, dy, dz] with a
    # zero-based atom index; VASP's DYNMAT numbers atoms and directions from 1.
    for index, (record, phonopy_disp) in enumerate(zip(perturbations, phonon.displacements)):
        atom, vector = int(phonopy_disp[0]), [float(x) for x in phonopy_disp[1:4]]
        if record["atom"] != atom + 1:
            raise ValueError(
                f"perturbation {index + 1} displaces atom {record['atom']} but "
                f"phonopy expects atom {atom + 1}"
            )
        if not np.allclose(record["displacement"], vector, atol=1e-6):
            raise ValueError(
                f"perturbation {index + 1} displaces along "
                f"{record['displacement']} but phonopy expects {vector}; the two "
                f"are not the same calculation"
            )
        array = np.asarray(record["forces"], dtype=float)
        if array.shape != (n_supercell, 3):
            raise ValueError(
                f"perturbation {index + 1} carries forces of shape {array.shape}, "
                f"expected {(n_supercell, 3)}"
            )
        forces.append(array)
    phonon.forces = forces
    phonon.produce_force_constants()
    return np.array(phonon.force_constants, dtype=float)


@dataclass
class DfptRelaxMaker(BaseVaspMaker):
    """Ionic relaxation for the EOS stage: ISIF=3 (free) or ISIF=2 (fixed volume).

    Deliberately NOT the DFPT maker.  Wiring ``IBRION=8`` into the EOS stage
    produces an E(V) curve of unrelaxed strained geometries, because VASP
    ignores ISIF/NSW for a dynamical-matrix run.
    """

    name: str = "DFPT-QHA relax"
    recipe: Recipe | None = None
    isif: int = 2
    use_supercell_kpoints: bool = False
    input_set_generator: Any = None
    run_vasp_kwargs: dict[str, Any] = field(default_factory=dict)
    stop_children_kwargs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.input_set_generator is None:
            if self.recipe is None:
                raise ValueError("DfptRelaxMaker requires an explicit Recipe")
            self.input_set_generator = DfptRelaxSetGenerator(
                recipe=self.recipe,
                isif=self.isif,
                use_supercell_kpoints=self.use_supercell_kpoints,
            )
        # BaseVaspMaker is a plain dataclass in this atomate2 release and
        # defines no post-init hook.


@dataclass
class DfptVaspMaker(BaseVaspMaker):
    """A standard atomate2 BaseVaspMaker configured for IBRION=8.

    A dataclass (not a hand-written ``__init__``) because jobflow serializes
    makers through monty, which requires every constructor argument to exist
    as a dataclass field on the instance.
    """

    name: str = "DFPT job"
    recipe: Recipe | None = None
    use_supercell_kpoints: bool = True
    input_set_generator: Any = None
    run_vasp_kwargs: dict[str, Any] = field(default_factory=dict)
    stop_children_kwargs: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.input_set_generator is None:
            if self.recipe is None:
                raise ValueError("DfptVaspMaker requires an explicit Recipe")
            self.input_set_generator = DfptInputSetGenerator(
                recipe=self.recipe, use_supercell_kpoints=self.use_supercell_kpoints
            )


@job
def collect_dfpt_task_doc(
    task_doc: Any,
    *,
    recipe: dict[str, Any],
    cell: dict[str, Any],
    supercell_matrix: Any = None,
    vasprun_path: str | Path | None = None,
) -> dict[str, Any]:
    """Convert one completed atomate2 TaskDoc into a DFPT TaskDoc.

    The collector is a pure jobflow analysis boundary: it reads the job's own
    ``DYNMAT`` and replays it through phonopy.  It never launches VASP and never
    scans a calculation tree it was not handed.

    ``cell`` is the primitive cell the supercell was built from.  It is passed
    explicitly rather than recovered from the supercell so that the ordering
    assertion compares two independently-known structures.
    """

    if hasattr(task_doc, "model_dump"):
        task_payload = task_doc.model_dump(mode="json")
        directory = task_doc.dir_name
    else:
        task_payload = dict(task_doc)
        directory = task_payload.get("dir_name")
    # jobflow stores a job's output as JSON and monty's decoder rebuilds the
    # objects it recognises, so what arrives here can be a plain dict OR a
    # reconstructed TaskDoc whose ``input.kpoints`` is a pymatgen Kpoints.  The
    # latter survives until something tries to serialize it again, and then
    # fails with "Unable to serialize unknown type: Kpoints" -- three jobs
    # downstream, in the validation gate.  Sanitizing once, at the boundary,
    # keeps the contract "this document is JSON" true from here on.
    task_payload = jsanitize(
        task_payload, strict=False, allow_bson=True, recursive_msonable=True
    )
    if not isinstance(task_payload, dict):
        raise ValueError(
            f"DFPT TaskDoc did not reduce to a mapping: {type(task_payload).__name__}"
        )
    # TaskDoc.dir_name is an emmet URI ("host:/abs/path"); atomate2 strips the
    # host for the same value in copy_vasp_outputs, so the collector must too.
    job_dir = Path(strip_hostname(str(vasprun_path or directory)))

    path = _find_file(job_dir, "DYNMAT")
    if path is None:
        raise ValueError(
            f"DFPT job directory has no DYNMAT: {job_dir}.  VASP writes it for "
            f"IBRION=8; its absence means the run did not reach the "
            f"dynamical-matrix step."
        )

    structure = task_payload.get("structure")
    if structure is None:
        raise ValueError("DFPT TaskDoc carries no final structure to fingerprint")
    supercell = Structure.from_dict(structure)

    parsed = Recipe.from_dict(recipe)
    factor = int(parsed.kpoints.supercell_factor)
    if supercell_matrix is None:
        supercell_matrix = [[factor, 0, 0], [0, factor, 0], [0, 0, factor]]

    perturbations, n_atoms_dynmat = parse_dynmat(path)
    if n_atoms_dynmat != len(supercell):
        raise ValueError(
            f"DYNMAT describes {n_atoms_dynmat} atoms but the TaskDoc structure "
            f"has {len(supercell)}"
        )
    force_constants = force_constants_from_dynmat(
        Structure.from_dict(cell), supercell, perturbations,
        supercell_matrix=supercell_matrix,
    )
    shape = tuple(int(x) for x in force_constants.shape)
    if len(shape) != 4 or shape[0] != shape[1] or shape[2:] != (3, 3):
        raise ValueError(f"unexpected DFPT force_constants shape: {shape}")
    # Acoustic sum rule: sum_j Phi_ij must vanish for every i.  A matrix that
    # fails it is not a force-constant matrix, whatever its shape.
    residual = float(np.abs(force_constants.sum(axis=1)).max())
    scale = float(np.abs(force_constants).max()) or 1.0
    if residual / scale > 1e-2:
        raise ValueError(
            f"DFPT force constants violate the acoustic sum rule "
            f"(max|sum_j Phi_ij| = {residual:.4g} vs max|Phi| = {scale:.4g})"
        )
    doc = DfptTaskDoc(
        task_doc=task_payload,
        recipe=recipe,
        force_constants=force_constants.tolist(),
        force_constants_provenance=DfptForceConstantsProvenance(
            source="DYNMAT",
            parser="phonopy replay at unit displacement",
            n_atoms=shape[0],
            shape=shape,
            unit=DYNMAT_UNIT,
            ordering_fingerprint=ordering_fingerprint(supercell),
            path=str(path),
            max_acoustic_sum_residual=residual,
        ),
        provenance={
            "task_doc_dir": str(job_dir),
            "task_state": task_payload.get("state"),
            "n_perturbations": len(perturbations),
        },
    )
    return doc.model_dump(mode="json")


@job
def validate_dfpt_task_doc(
    dfpt_doc: dict[str, Any], *, recipe: dict[str, Any] | None = None
) -> Response:
    """Stop descendants unless the DFPT result matches the recipe it claims.

    The shape check is not decoration.  A DFPT run on the wrong cell still
    returns a square ``(n, n, 3, 3)`` force-constant array, so a Gamma-only
    result from the 2-atom cell would otherwise be recorded as a valid volume
    and silently lose the finite-q instability the recipe exists to find.
    """

    try:
        parsed = DfptTaskDoc.model_validate(dfpt_doc)
        if not parsed.force_constants:
            raise ValueError("empty force constants")
        if recipe is not None:
            expected = Recipe.from_dict(recipe)
            expected_atoms = expected.kpoints.cell_atoms * int(
                expected.kpoints.supercell_factor
            ) ** 3
            if parsed.force_constants_provenance.n_atoms != expected_atoms:
                raise ValueError(
                    f"DFPT ran on {parsed.force_constants_provenance.n_atoms} atoms "
                    f"but the recipe's phonon supercell is {expected_atoms} atoms "
                    f"(cell {expected.kpoints.cell_atoms} x factor "
                    f"{expected.kpoints.supercell_factor})"
                )
    except Exception as exc:
        return Response(output={"ok": False, "error": str(exc)}, stop_children=True)
    try:
        document = parsed.model_dump(mode="json")
    except Exception as exc:  # noqa: BLE001 - an unserializable document is a finding
        return Response(
            output={
                "ok": False,
                "error": f"validated DFPT document could not be serialized: {exc}",
            },
            stop_children=True,
        )
    return Response(
        output={
            "ok": True,
            "schema_version": parsed.schema_version,
            "n_atoms": parsed.force_constants_provenance.n_atoms,
            "dfpt_doc": document,
        }
    )


__all__ = [
    "DYNMAT_UNIT",
    "DfptRelaxMaker",
    "DfptVaspMaker",
    "collect_dfpt_task_doc",
    "force_constants_from_dynmat",
    "parse_dynmat",
    "validate_dfpt_task_doc",
]
