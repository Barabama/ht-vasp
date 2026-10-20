"""DFPT TaskDoc provenance envelope for VASP dynamical matrices."""

from __future__ import annotations

import json
from typing import Any, Literal

from monty.json import jsanitize
from pydantic import BaseModel, ConfigDict, Field, field_validator


def as_json_mapping(value: Any) -> dict[str, Any]:
    """Reduce whatever a jobflow boundary handed us to a JSON-safe mapping.

    Three shapes arrive, and which one depends on the runner rather than on the
    code:

    * a plain ``dict`` -- a job that wrote one;
    * a ``TaskDoc`` pydantic model -- ``BaseVaspMaker``'s output schema;
    * a dict monty's decoder rebuilt from JSON, whose nested ``input.kpoints``,
      ``structure`` and ``orig_inputs.incar`` are live pymatgen objects again.

    The third is the one that bites: it survives until something serializes the
    document a second time, and then raises ``Unable to serialize unknown type:
    Kpoints`` -- three jobs downstream, inside the validation gate, where the
    cause is hardest to read.  Sanitizing at this boundary makes "this document
    is JSON" true for everything after it.

    ``jsanitize`` is asked for ``strict=True`` because ``strict=False`` falls
    back to ``str(obj)`` for anything it does not recognise, which silently
    turns a pydantic model into its repr.  ``strict=True`` raises instead, and
    the repr of a whole TaskDoc is exactly the kind of quiet corruption this
    project keeps finding.

    The final JSON round-trip renders the datetime objects jsanitize keeps and
    guarantees the result survives the next job boundary.
    """

    if isinstance(value, BaseModel):
        value = value.model_dump(mode="json")
    payload = jsanitize(value, strict=True, allow_bson=True, recursive_msonable=True)
    if not isinstance(payload, dict):
        raise ValueError(f"expected a mapping, got {type(value).__name__}")
    return json.loads(json.dumps(payload, default=str))


class DfptForceConstantsProvenance(BaseModel):
    """Where a DFPT force-constant matrix came from, and in what indexing.

    ``ordering_fingerprint`` is not decoration.  VASP indexes the matrix by its
    POSCAR atom order while phonopy indexes its own supercell differently, so a
    matrix that is attached across the two orderings is still nonzero, symmetric
    and gives plausible frequencies -- it is simply the wrong assignment.  The
    fingerprint is what lets a consumer prove the two agree before using it.

    ``source`` is DYNMAT, not FORCE_CONSTANTS: VASP writes the latter only after
    an offline conversion, so a collector that waits for it finds nothing after
    the allocation has been spent.
    """

    model_config = ConfigDict(extra="forbid")

    #: VASP's own Hessian from ``vasprun.xml`` is the authority; the text
    #: ``DYNMAT`` replay is the fallback and always recorded separately.
    source: Literal["DYNMAT", "vasprun.xml dynmat hessian"]
    parser: Literal[
        "phonopy replay at unit displacement",
        "phonopy VASP vasprun parser",
    ]
    n_atoms: int = Field(gt=0)
    shape: tuple[int, int, int, int]
    unit: str = "eV/Angstrom^2"
    ordering_fingerprint: str = Field(min_length=1)
    path: str | None = None
    #: the acoustic sum rule is the one cheap check that a matrix is a
    #: force-constant matrix at all
    max_acoustic_sum_residual: float | None = None


class DfptTaskDoc(BaseModel):
    """Serializable output contract layered on atomate2's TaskDoc."""

    model_config = ConfigDict(extra="forbid")

    schema_version: Literal["htvasp-dfpt-task/v1"] = "htvasp-dfpt-task/v1"
    task_doc: dict[str, Any]
    recipe: dict[str, Any]
    force_constants: list[list[list[list[float]]]]
    force_constants_provenance: DfptForceConstantsProvenance
    provenance: dict[str, Any] = Field(default_factory=dict)

    @field_validator("task_doc", "recipe", "provenance", mode="before")
    @classmethod
    def _coerce_to_mapping(cls, value: Any) -> dict[str, Any]:
        return as_json_mapping(value)


__all__ = ["DfptForceConstantsProvenance", "DfptTaskDoc", "as_json_mapping"]
