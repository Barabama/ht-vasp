"""DFPT TaskDoc provenance envelope for VASP dynamical matrices."""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field


class DfptForceConstantsProvenance(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source: Literal["vasprun.xml"]
    parser: Literal["pymatgen.Vasprun.force_constants"]
    n_atoms: int = Field(gt=0)
    shape: tuple[int, int, int, int]
    hessian_unit: str = "eV/Angstrom^2/amu"
    vasprun_path: str | None = None


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
