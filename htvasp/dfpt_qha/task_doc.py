"""Standard output models for DFPT-QHA jobs.

VASP jobs return atomate2/emmet TaskDoc instances.  These models describe the
small, stable production envelope around those TaskDocs and phonon/QHA analysis
outputs; they are not replacements for TaskDoc.
"""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field


class DfptValidationDoc(BaseModel):
    """Fail-closed validation result attached before thermal/QHA analysis."""

    model_config = ConfigDict(extra="forbid")

    ok: bool
    recipe_schema: str = Field(min_length=1)
    problems: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)


class DfptPhononDoc(BaseModel):
    """Stable DFPT phonon envelope around the common phonon document."""

    model_config = ConfigDict(extra="forbid")

    schema_version: Literal["htvasp-dfpt-phonon/v1"] = "htvasp-dfpt-phonon/v1"
    phonon_document: dict[str, Any]
    validation: DfptValidationDoc
    recipe: dict[str, Any]
    provenance: dict[str, Any] = Field(default_factory=dict)


class DfptQhaDoc(BaseModel):
    """Stable production envelope around atomate2's QHA analysis output."""

    model_config = ConfigDict(extra="forbid")

    schema_version: Literal["htvasp-dfpt-qha/v1"] = "htvasp-dfpt-qha/v1"
    qha_document: dict[str, Any]
    recipe: dict[str, Any]
    validation: DfptValidationDoc
    provenance: dict[str, Any] = Field(default_factory=dict)


__all__ = ["DfptPhononDoc", "DfptQhaDoc", "DfptValidationDoc"]
