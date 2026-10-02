"""Output model for one SQS Relax -> Static calculation."""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator
from pymatgen.core import Structure


class SqsStaticObservation(BaseModel):
    """Portable observation emitted after a successful Static task."""

    model_config = ConfigDict(extra="forbid")

    calculation_id: str = Field(min_length=1)
    batch_fingerprint: str | None = None
    intent_key: str | None = None
    state: Literal["successful"]
    energy_total_eV: float
    energy_per_atom_eV: float
    n_sites: int = Field(gt=0)
    final_structure: dict[str, Any]
    task_document_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")

    @field_validator("energy_total_eV", "energy_per_atom_eV")
    @classmethod
    def _validate_finite_energy(cls, value: float) -> float:
        import math

        value = float(value)
        if not math.isfinite(value):
            raise ValueError("energy values must be finite")
        return value

    @model_validator(mode="after")
    def _validate_consistency(self) -> "SqsStaticObservation":
        import math

        if not math.isclose(
            self.energy_total_eV / self.n_sites,
            self.energy_per_atom_eV,
            rel_tol=0.0,
            abs_tol=1e-10,
        ):
            raise ValueError("energy_per_atom_eV is inconsistent with total energy and n_sites")
        try:
            structure = Structure.from_dict(self.final_structure)
        except Exception as exc:
            raise ValueError("final_structure is not a valid pymatgen structure") from exc
        if len(structure) != self.n_sites:
            raise ValueError("final_structure site count is inconsistent with n_sites")
        return self


# Short name used by the calculation layer; retain the explicit historical
# name as an import-compatible alias.
SqsObservation = SqsStaticObservation


__all__ = ["SqsObservation", "SqsStaticObservation"]
