"""Versioned SQS mixing result contract for downstream CALPHAD tooling."""

from __future__ import annotations

import math
from datetime import datetime, timezone
from typing import Any, Literal, Mapping

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from htvasp.model.sqs import SqsPlan
from htvasp.sqs.task_doc import SqsStaticObservation


SCHEMA_VERSION = "htvasp-sqs-mixing/v1"
EV_PER_ATOM_TO_J_PER_MOL = 96485.33212331002


def _finite(value: float, field_name: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{field_name} must be finite")
    return result


class SqsMixingError(ValueError):
    """Raised when an SQS mixing dataset is incomplete or inconsistent."""


class Sublattice(BaseModel):
    model_config = ConfigDict(extra="forbid")

    site_ratio: float
    constituents: list[str]

    @field_validator("site_ratio")
    @classmethod
    def validate_ratio(cls, value: float) -> float:
        value = _finite(value, "site_ratio")
        if value <= 0:
            raise ValueError("site_ratio must be positive")
        return value

    @field_validator("constituents")
    @classmethod
    def validate_constituents(cls, value: list[str]) -> list[str]:
        result = [str(item).strip() for item in value]
        if not result or any(not item for item in result) or len(set(result)) != len(result):
            raise ValueError("constituents must be distinct and non-empty")
        return result


class CalphadMapping(BaseModel):
    """Explicit semantic mapping from SQS energies to a CALPHAD G model."""

    model_config = ConfigDict(extra="forbid")

    phase_name: Literal["BCC_A2"]
    phase_model: Literal["substitutional_single_sublattice"]
    sublattices: list[Sublattice]
    components: tuple[str, str]
    composition_coordinate: Literal["x_A"]
    energy_basis: Literal["J_per_mol_atoms"]
    excess_model: Literal["redlich_kister"]
    parameter_target: Literal["G"]
    maximum_order: Literal[1]

    @field_validator("components")
    @classmethod
    def validate_components(cls, value: tuple[str, str]) -> tuple[str, str]:
        result = tuple(str(item).strip() for item in value)
        if len(result) != 2 or not all(result) or result[0] == result[1]:
            raise ValueError("components must be two different symbols")
        return result

    @model_validator(mode="after")
    def validate_mapping(self) -> "CalphadMapping":
        if len(self.sublattices) != 1 or self.sublattices[0].site_ratio != 1:
            raise ValueError("BCC_A2 requires one sublattice with site_ratio=1")
        if tuple(self.sublattices[0].constituents) != self.components:
            raise ValueError("sublattice constituents must match components")
        return self


class CalculationProvenance(BaseModel):
    model_config = ConfigDict(extra="forbid")

    settings_fingerprint: str = Field(min_length=1)
    batch_fingerprint: str | None = None
    run_id: str | None = None
    # Keep the recipe identity in the downstream artifact.  The batch
    # manifest has always carried this field, and collect must not discard it
    # when converting the stage manifest into the formal SQS contract.
    recipe_version: str = Field(default="unknown", min_length=1)
    workflow: dict[str, Any]
    global_incar: dict[str, Any] = Field(default_factory=dict)
    relax_incar: dict[str, Any]
    static_incar: dict[str, Any]
    relax_incar_overrides: dict[str, Any] = Field(default_factory=dict)
    static_incar_overrides: dict[str, Any] = Field(default_factory=dict)
    kpoints: dict[str, Any] = Field(default_factory=dict)
    vasp: dict[str, Any] = Field(default_factory=dict)
    potcar: dict[str, Any] = Field(default_factory=dict)


class SqsProvenance(BaseModel):
    model_config = ConfigDict(extra="forbid")

    engine: str = Field(min_length=1)
    engine_version: str = "unknown"
    iterations: int = Field(gt=0)
    seed: int | None = None
    seed_applied: bool | None = None
    # Some sqsgenerator versions cannot apply a seed through a documented
    # configuration key.  Preserve the producer's explanatory note instead
    # of dropping it during collect.
    seed_note: str | None = None
    config: dict[str, Any] = Field(default_factory=dict)
    objective: float | None = None
    structure_checksum: str = Field(pattern=r"^[0-9a-f]{64}$")


class StaticObservation(BaseModel):
    model_config = ConfigDict(extra="forbid")

    run_id: str = Field(min_length=1)
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
    def validate_energy(cls, value: float, info) -> float:
        return _finite(value, info.field_name)

    @model_validator(mode="after")
    def validate_consistency(self) -> "StaticObservation":
        if not math.isclose(
            self.energy_total_eV / self.n_sites,
            self.energy_per_atom_eV,
            rel_tol=0.0,
            abs_tol=1e-10,
        ):
            raise ValueError("energy_per_atom_eV is inconsistent with total energy")
        return self

    @classmethod
    def from_sqs_observation(
        cls, observation: SqsStaticObservation, *, run_id: str | None = None
    ) -> "StaticObservation":
        return cls(
            run_id=run_id or observation.calculation_id,
            batch_fingerprint=observation.batch_fingerprint,
            intent_key=observation.intent_key,
            state=observation.state,
            energy_total_eV=observation.energy_total_eV,
            energy_per_atom_eV=observation.energy_per_atom_eV,
            n_sites=observation.n_sites,
            final_structure=observation.final_structure,
            task_document_sha256=observation.task_document_sha256,
        )


class PureReference(BaseModel):
    model_config = ConfigDict(extra="forbid")

    component: str
    final_structure: dict[str, Any]
    static: StaticObservation


class SqsMixingPoint(BaseModel):
    model_config = ConfigDict(extra="forbid")

    point_id: str = Field(min_length=1)
    x_A: float
    site_counts: dict[str, int]
    sqs: SqsProvenance
    initial_structure: dict[str, Any]
    final_structure: dict[str, Any]
    static: StaticObservation
    h_mix_eV_per_atom: float
    h_mix_J_per_mol: float

    @field_validator("x_A", "h_mix_eV_per_atom", "h_mix_J_per_mol")
    @classmethod
    def validate_finite(cls, value: float, info) -> float:
        return _finite(value, info.field_name)


class DatasetValidation(BaseModel):
    model_config = ConfigDict(extra="forbid")

    settings_fingerprint: str = Field(min_length=1)
    batch_fingerprint: str | None = None
    all_static_successful: Literal[True]
    warnings: list[str] = Field(default_factory=list)


class SqsMixingData(BaseModel):
    """The only JSON dataset shape intended for em-tdb consumption."""

    model_config = ConfigDict(extra="forbid")

    schema_version: Literal[SCHEMA_VERSION]
    dataset_id: str = Field(min_length=1)
    created_at: datetime
    calphad_mapping: CalphadMapping
    calculation: CalculationProvenance
    references: dict[str, PureReference]
    points: list[SqsMixingPoint] = Field(min_length=3)
    validation: DatasetValidation

    @model_validator(mode="after")
    def validate_dataset(self) -> "SqsMixingData":
        components = self.calphad_mapping.components
        if set(self.references) != set(components):
            raise ValueError("references must contain exactly the mapped components")
        if self.validation.settings_fingerprint != self.calculation.settings_fingerprint:
            raise ValueError("validation and calculation fingerprints differ")
        if (
            self.validation.batch_fingerprint is not None
            and self.calculation.batch_fingerprint is not None
            and self.validation.batch_fingerprint != self.calculation.batch_fingerprint
        ):
            raise ValueError("validation and calculation batch fingerprints differ")
        endpoint_energies = {
            element: self.references[element].static.energy_per_atom_eV for element in components
        }
        point_ids: set[str] = set()
        x_values: set[float] = set()
        for point in self.points:
            if point.point_id in point_ids or point.x_A in x_values:
                raise ValueError("point IDs and x_A values must be unique")
            point_ids.add(point.point_id)
            x_values.add(point.x_A)
            counts = {str(key): value for key, value in point.site_counts.items()}
            if set(counts) != set(components) or any(value <= 0 for value in counts.values()):
                raise ValueError("site_counts must match positive component counts")
            if sum(counts.values()) != point.static.n_sites:
                raise ValueError("site_counts must sum to static.n_sites")
            expected_x = counts[components[0]] / point.static.n_sites
            if not math.isclose(point.x_A, expected_x, rel_tol=0.0, abs_tol=1e-12):
                raise ValueError("x_A does not match site_counts")
            hmix = (
                point.static.energy_total_eV
                - sum(counts[element] * endpoint_energies[element] for element in components)
            ) / point.static.n_sites
            if not math.isclose(point.h_mix_eV_per_atom, hmix, rel_tol=0.0, abs_tol=1e-10):
                raise ValueError("h_mix_eV_per_atom cannot be reproduced")
            expected_j = hmix * EV_PER_ATOM_TO_J_PER_MOL
            if not math.isclose(point.h_mix_J_per_mol, expected_j, rel_tol=0.0, abs_tol=1e-7):
                raise ValueError("h_mix_J_per_mol cannot be reproduced")
        return self


def build_bcc_a2_mapping(components: tuple[str, str]) -> CalphadMapping:
    return CalphadMapping(
        phase_name="BCC_A2",
        phase_model="substitutional_single_sublattice",
        sublattices=[{"site_ratio": 1, "constituents": list(components)}],
        components=components,
        composition_coordinate="x_A",
        energy_basis="J_per_mol_atoms",
        excess_model="redlich_kister",
        parameter_target="G",
        maximum_order=1,
    )


def build_sqs_mixing_data(
    plan: SqsPlan,
    *,
    reference_observations: Mapping[str, SqsStaticObservation],
    point_observations: Mapping[str, SqsStaticObservation],
    point_structures: Mapping[str, Mapping[str, Any]],
    point_provenance: Mapping[str, Mapping[str, Any]],
    calculation: CalculationProvenance | Mapping[str, Any],
    dataset_id: str | None = None,
) -> SqsMixingData:
    """Recompute Hmix from validated observations and build the v1 dataset."""

    if plan.system.phase != "BCC":
        raise SqsMixingError("htvasp-sqs-mixing/v1 currently supports BCC_A2 only")
    components = plan.system.canonical_elements
    calculation = CalculationProvenance.model_validate(calculation)
    references: dict[str, PureReference] = {}
    for element in components:
        observation = reference_observations.get(element)
        if observation is None:
            raise SqsMixingError(f"Missing endpoint observation for {element}")
        references[element] = PureReference(
            component=element,
            final_structure=observation.final_structure,
            static=StaticObservation.from_sqs_observation(observation),
        )
    points: list[SqsMixingPoint] = []
    endpoint_per_atom = {
        element: references[element].static.energy_per_atom_eV for element in components
    }
    for composition_point in plan.points:
        point_id = f"{plan.stable_id}-{composition_point.label}"
        observation = point_observations.get(point_id)
        if observation is None:
            raise SqsMixingError(f"Missing point observation for {point_id}")
        counts = {
            components[0]: composition_point.element_a_count,
            components[1]: composition_point.element_b_count,
        }
        hmix = (
            observation.energy_total_eV
            - sum(counts[element] * endpoint_per_atom[element] for element in components)
        ) / observation.n_sites
        provenance = SqsProvenance.model_validate(point_provenance.get(point_id, {}))
        initial = point_structures.get(point_id)
        if initial is None:
            raise SqsMixingError(f"Missing initial structure for {point_id}")
        points.append(
            SqsMixingPoint(
                point_id=point_id,
                x_A=composition_point.fraction_a,
                site_counts=counts,
                sqs=provenance,
                initial_structure=dict(initial),
                final_structure=observation.final_structure,
                static=StaticObservation.from_sqs_observation(observation),
                h_mix_eV_per_atom=hmix,
                h_mix_J_per_mol=hmix * EV_PER_ATOM_TO_J_PER_MOL,
            )
        )
    return SqsMixingData(
        schema_version=SCHEMA_VERSION,
        dataset_id=dataset_id or plan.stable_id,
        created_at=datetime.now(timezone.utc),
        calphad_mapping=build_bcc_a2_mapping(components),
        calculation=calculation,
        references=references,
        points=points,
        validation=DatasetValidation(
            settings_fingerprint=calculation.settings_fingerprint,
            batch_fingerprint=calculation.batch_fingerprint,
            all_static_successful=True,
        ),
    )


__all__ = [
    "CalculationProvenance",
    "CalphadMapping",
    "DatasetValidation",
    "EV_PER_ATOM_TO_J_PER_MOL",
    "PureReference",
    "SCHEMA_VERSION",
    "SqsMixingData",
    "SqsMixingError",
    "SqsMixingPoint",
    "SqsProvenance",
    "StaticObservation",
    "build_bcc_a2_mapping",
    "build_sqs_mixing_data",
]
