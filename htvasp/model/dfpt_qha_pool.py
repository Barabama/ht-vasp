"""Typed endpoint pool values; structure construction remains in the CLI layer."""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field

POOL_SCHEMA = "legacy-qha-pool/pool-structures/v1"


class PoolEndpoint(BaseModel):
    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    name: str
    structure_class: str = Field(alias="class")
    elements: list[str]
    phase: str
    structure: dict[str, Any]
    n_atoms: int = Field(gt=0)
    lattice: str
    a: float
    volume_per_atom: float
    source: str
    vpa_source: dict[str, str] | str | None = None
    is_proxy: bool
    physical_state_id: str | None = None
    site_role: dict[str, Any] | None = None
    elem_mpid_note: str | None = None
    detour: str | None = None


class PoolManifest(BaseModel):
    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    schema_version: Literal[POOL_SCHEMA] = Field(POOL_SCHEMA, alias="schema")
    n_endpoints: int = Field(gt=0)
    by_class: dict[str, int]
    decisions: dict[str, Any]
    physical_states: dict[str, list[str]]
    duplicate_states: dict[str, list[str]]
    endpoints: dict[str, PoolEndpoint]

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "PoolManifest":
        result = cls.model_validate(data)
        if result.n_endpoints != len(result.endpoints):
            raise ValueError("n_endpoints does not match endpoints mapping")
        return result


__all__ = ["POOL_SCHEMA", "PoolEndpoint", "PoolManifest"]
