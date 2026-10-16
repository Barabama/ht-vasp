"""DFPT TaskDoc provenance envelope for VASP dynamical matrices."""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field


class DfptForceConstantsProvenance(BaseModel):
    """Where a DFPT force-constant matrix came from, and in what indexing.

    ``ordering_fingerprint`` is not decoration.  VASP indexes the matrix by its
    POSCAR atom order while phonopy indexes its own supercell differently, so a
    matrix that is attached across the two orderings is still nonzero, symmetric
    and gives plausible frequencies -- it is simply the wrong assignment.  The
    fingerprint is what lets a consumer prove the two agree before using it.
    """

    model_config = ConfigDict(extra="forbid")

    source: Literal["FORCE_CONSTANTS"]
    parser: Literal["phonopy.file_IO.parse_FORCE_CONSTANTS"]
    n_atoms: int = Field(gt=0)
    shape: tuple[int, int, int, int]
    unit: str = "eV/Angstrom^2"
    ordering_fingerprint: str = Field(min_length=1)
    path: str | None = None
    #: the acoustic sum rule is the one cheap physical check that a matrix is a
    #: force-constant matrix at all; the vasprun <dynmat> fails it by 65%
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


__all__ = ["DfptForceConstantsProvenance", "DfptTaskDoc"]
