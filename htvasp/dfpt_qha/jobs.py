"""One VASP DFPT job and its dynamical-matrix collector."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from atomate2.vasp.jobs.base import BaseVaspMaker, get_vasp_task_document
from jobflow import Response, job
from pymatgen.io.vasp.outputs import Vasprun

from htvasp.dfpt_qha.task_doc import DfptForceConstantsProvenance, DfptTaskDoc
from htvasp.dfpt_qha.input_set import DfptInputSetGenerator
from htvasp.model.dfpt_qha_recipe import Recipe


class DfptVaspMaker(BaseVaspMaker):
    """A standard atomate2 BaseVaspMaker configured for IBRION=8."""

    def __init__(self, recipe: Recipe, **kwargs: Any) -> None:
        kwargs.setdefault("input_set_generator", DfptInputSetGenerator(recipe=recipe))
        kwargs.setdefault("stop_children_kwargs", {})
        super().__init__(**kwargs)
        self.recipe = recipe

    def make(self, structure: Any, prev_dir: str | Path | None = None):
        return super().make(structure, prev_dir=prev_dir)


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
    path = Path(vasprun_path or directory) / "vasprun.xml"
    if not path.is_file():
        raise ValueError(f"DFPT TaskDoc directory has no vasprun.xml: {path}")
    vasprun = Vasprun(str(path), parse_potcar_file=False)
    force_constants = getattr(vasprun, "force_constants", None)
    if force_constants is None:
        raise ValueError(f"DFPT vasprun has no <dynmat> force constants: {path}")
    shape = tuple(int(x) for x in force_constants.shape)
    if len(shape) != 4 or shape[0] != shape[1] or shape[2:] != (3, 3):
        raise ValueError(f"unexpected DFPT force_constants shape: {shape}")
    doc = DfptTaskDoc(
        task_doc=task_payload,
        recipe=recipe,
        force_constants=force_constants.tolist(),
        force_constants_provenance=DfptForceConstantsProvenance(
            source="vasprun.xml",
            parser="pymatgen.Vasprun.force_constants",
            n_atoms=shape[0],
            shape=shape,
            vasprun_path=str(path),
        ),
        provenance={"task_doc_dir": str(directory), "task_state": task_payload.get("state")},
    )
    return doc.model_dump(mode="json")


@job
def validate_dfpt_task_doc(dfpt_doc: dict[str, Any]) -> Response:
    """Stop descendants unless the DFPT force constants and provenance exist."""

    try:
        parsed = DfptTaskDoc.model_validate(dfpt_doc)
        if not parsed.force_constants:
            raise ValueError("empty force constants")
    except Exception as exc:
        return Response(output={"ok": False, "error": str(exc)}, stop_children=True)
    return Response(output={"ok": True, "schema_version": parsed.schema_version})


__all__ = ["DfptVaspMaker", "collect_dfpt_task_doc", "validate_dfpt_task_doc"]
