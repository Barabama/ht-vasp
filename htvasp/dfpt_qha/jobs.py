"""One VASP DFPT job and its dynamical-matrix collector."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from dataclasses import dataclass, field

import numpy as np
from atomate2.vasp.jobs.base import BaseVaspMaker
from atomate2.utils.path import strip_hostname
from jobflow import Response, job
from monty.io import zopen
from phonopy.file_IO import parse_FORCE_CONSTANTS
from pymatgen.core import Structure

from htvasp.dfpt_qha.phonons import ordering_fingerprint
from htvasp.dfpt_qha.task_doc import DfptForceConstantsProvenance, DfptTaskDoc
from htvasp.dfpt_qha.input_set import DfptInputSetGenerator, DfptRelaxSetGenerator
from htvasp.model.dfpt_qha_recipe import Recipe


def _find_file(directory: Path, name: str) -> Path | None:
    """``name`` or ``name.gz`` inside ``directory``, whichever exists.

    BaseVaspMaker gzips a fixed file list at the end of every job (vasprun.xml
    among them), so a collector that only looks for the plain name breaks the
    moment that list changes.  Both spellings are accepted.
    """

    for candidate in (directory / name, directory / f"{name}.gz"):
        if candidate.is_file():
            return candidate
    return None


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
    vasprun_path: str | Path | None = None,
) -> dict[str, Any]:
    """Convert one completed atomate2 TaskDoc into a DFPT TaskDoc.

    ``vasprun_path`` is optional because the normal BaseVaspMaker output keeps
    the calculation directory in ``TaskDoc.dir_name``.  The collector remains
    a pure jobflow analysis boundary: it parses ``vasprun.xml`` with pymatgen,
    never launches VASP and never scans an arbitrary calculation tree.
    """

    if hasattr(task_doc, "model_dump"):
        task_payload = task_doc.model_dump(mode="json")
        directory = task_doc.dir_name
    else:
        task_payload = dict(task_doc)
        directory = task_payload.get("dir_name")
    # TaskDoc.dir_name is an emmet URI ("host:/abs/path"); atomate2 strips the
    # host for the same value in copy_vasp_outputs, so the collector must too.
    job_dir = Path(strip_hostname(str(vasprun_path or directory)))

    # FORCE_CONSTANTS, not the vasprun <dynmat>.  Both exist and both are
    # square (n, n, 3, 3), but they are different objects: the vasprun hessian
    # violates the acoustic sum rule by 65% of its own scale, so attaching it
    # to phonopy gives frequencies that are plausible and wrong.  VASP writes
    # FORCE_CONSTANTS in phonopy's own format for IBRION=8, and it satisfies
    # the sum rule to 3e-5.
    path = _find_file(job_dir, "FORCE_CONSTANTS")
    if path is None:
        raise ValueError(
            f"DFPT job directory has no FORCE_CONSTANTS: {job_dir}.  VASP writes "
            f"it for IBRION=8; its absence means the run did not reach the "
            f"dynamical-matrix step."
        )
    force_constants = np.array(parse_FORCE_CONSTANTS(str(path)), dtype=float)
    shape = tuple(int(x) for x in force_constants.shape)
    if len(shape) != 4 or shape[0] != shape[1] or shape[2:] != (3, 3):
        raise ValueError(f"unexpected DFPT force_constants shape: {shape}")
    # Acoustic sum rule: sum_j Phi_ij must vanish for every i.  A matrix that
    # fails it is not a force-constant matrix, whatever its shape.
    residual = float(np.abs(force_constants.sum(axis=1)).max())
    scale = float(np.abs(force_constants).max()) or 1.0
    if residual / scale > 1e-3:
        raise ValueError(
            f"DFPT force constants violate the acoustic sum rule "
            f"(max|sum_j Phi_ij| = {residual:.4g} vs max|Phi| = {scale:.4g}); "
            f"this is not a usable force-constant matrix"
        )
    structure = task_payload.get("structure")
    if structure is None:
        raise ValueError("DFPT TaskDoc carries no final structure to fingerprint")
    fingerprint = ordering_fingerprint(Structure.from_dict(structure))
    doc = DfptTaskDoc(
        task_doc=task_payload,
        recipe=recipe,
        force_constants=force_constants.tolist(),
        force_constants_provenance=DfptForceConstantsProvenance(
            source="FORCE_CONSTANTS",
            parser="phonopy.file_IO.parse_FORCE_CONSTANTS",
            n_atoms=shape[0],
            shape=shape,
            ordering_fingerprint=fingerprint,
            path=str(path),
            max_acoustic_sum_residual=residual,
        ),
        provenance={
            "task_doc_dir": str(job_dir),
            "task_state": task_payload.get("state"),
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
    return Response(
        output={
            "ok": True,
            "schema_version": parsed.schema_version,
            "n_atoms": parsed.force_constants_provenance.n_atoms,
            "dfpt_doc": parsed.model_dump(mode="json"),
        }
    )


__all__ = [
    "DfptRelaxMaker",
    "DfptVaspMaker",
    "collect_dfpt_task_doc",
    "validate_dfpt_task_doc",
]
