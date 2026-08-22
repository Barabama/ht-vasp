"""Data model and engine adapter for binary SQS mixing calculations."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import logging
import math
from datetime import datetime, timezone
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Any, Literal, Mapping, Sequence

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator
from pymatgen.core import Element, Structure
from pymatgen.io.vasp import Poscar

from htvasp.model.templates import TemplateRegistry, structure_checksum

log = logging.getLogger(__name__)

SCHEMA_VERSION = "htvasp-sqs-mixing/v1"
EV_PER_ATOM_TO_J_PER_MOL = 96485.33212331002

# The prototype snapshots contain different numbers of sites. Keep the
# default geometric supercells explicit instead of implying one fixed atom
# count for every phase.
DEFAULT_SUPERCELLS: dict[str, tuple[int, int, int]] = {
    "BCC": (2, 2, 2),
    "FCC": (2, 2, 2),
    "HCP": (2, 2, 1),
}


class SqsMixingError(ValueError):
    """Raised when an SQS dataset is incomplete or inconsistent."""


def _json_checksum(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True, default=str).encode()
    ).hexdigest()


def _finite(value: float, field_name: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{field_name} must be finite")
    return number


def _structure_dict(value: Any) -> dict[str, Any]:
    if isinstance(value, Structure):
        return value.as_dict()
    if isinstance(value, Mapping) and "sites" in value:
        return dict(value)
    if hasattr(value, "as_dict"):
        result = value.as_dict()
        if isinstance(result, dict):
            return result
    raise SqsMixingError("Expected a pymatgen Structure or Structure.as_dict()")


@dataclass(frozen=True)
class SqsSystem:
    """A binary disordered system; unlike Endmember names, A/B are symmetric."""

    phase: str
    element_a: str
    element_b: str

    def __post_init__(self) -> None:
        phase = str(self.phase).strip().upper()
        if phase not in {"BCC", "FCC", "HCP"}:
            raise ValueError(f"Unsupported SQS phase: {self.phase!r}")
        raw_a = str(self.element_a).strip()
        raw_b = str(self.element_b).strip()
        if not raw_a or not raw_b:
            raise ValueError("SQS element symbols must be non-empty")
        element_a = raw_a[0].upper() + raw_a[1:].lower()
        element_b = raw_b[0].upper() + raw_b[1:].lower()
        try:
            Element(element_a)
            Element(element_b)
        except ValueError as exc:
            raise ValueError(f"Invalid SQS element pair: {self.element_a!r}, {self.element_b!r}") from exc
        if element_a == element_b:
            raise ValueError("SQS requires two distinct elements")
        element_a, element_b = sorted((element_a, element_b))
        object.__setattr__(self, "phase", phase)
        object.__setattr__(self, "element_a", element_a)
        object.__setattr__(self, "element_b", element_b)

    @classmethod
    def from_name(cls, name: str) -> "SqsSystem":
        parts = str(name).split("-")
        if len(parts) != 3:
            raise ValueError("SQS name must be phase-A-B, e.g. FCC-Co-Fe")
        return cls(parts[0], parts[1], parts[2])

    @property
    def canonical_elements(self) -> tuple[str, str]:
        return tuple(sorted((self.element_a, self.element_b)))

    @property
    def canonical_name(self) -> str:
        a, b = self.canonical_elements
        return f"{self.phase}-{a}-{b}"

    @property
    def key(self) -> tuple[str, str, str]:
        return (self.phase, *self.canonical_elements)

    def reference_name(self, element: str) -> str:
        if element not in {self.element_a, self.element_b}:
            raise ValueError(f"{element!r} is not part of {self.canonical_name}")
        return f"{self.phase}-{element}-{element}"


@dataclass(frozen=True)
class SqsCompositionPoint:
    element_a_count: int
    element_b_count: int
    n_sites: int

    def __post_init__(self) -> None:
        if self.element_a_count <= 0 or self.element_b_count <= 0:
            raise ValueError("SQS composition points must contain both elements")
        if self.n_sites <= 0 or self.element_a_count + self.element_b_count != self.n_sites:
            raise ValueError("Composition counts must sum to a positive n_sites")

    @classmethod
    def from_counts(cls, element_a_count: int, element_b_count: int) -> "SqsCompositionPoint":
        counts = (element_a_count, element_b_count)
        if any(isinstance(value, bool) or not isinstance(value, int) for value in counts):
            raise ValueError("Composition counts must be integers")
        return cls(element_a_count, element_b_count, element_a_count + element_b_count)

    @classmethod
    def from_fraction(cls, fraction_a: float | str | Fraction, n_sites: int) -> "SqsCompositionPoint":
        fraction = Fraction(str(fraction_a)) if not isinstance(fraction_a, Fraction) else fraction_a
        count = fraction * int(n_sites)
        if count.denominator != 1:
            raise ValueError(f"Fraction {fraction_a!r} cannot be represented with {n_sites} sites")
        count_a = int(count)
        return cls.from_counts(count_a, int(n_sites) - count_a)

    @property
    def fraction_a(self) -> float:
        return self.element_a_count / self.n_sites

    @property
    def label(self) -> str:
        return f"x{self.fraction_a:.6f}".rstrip("0").rstrip(".")

    def as_dict(self) -> dict[str, Any]:
        return {
            "element_a_count": self.element_a_count,
            "element_b_count": self.element_b_count,
            "n_sites": self.n_sites,
            "fraction_a": self.fraction_a,
            "label": self.label,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "SqsCompositionPoint":
        return cls(int(data["element_a_count"]), int(data["element_b_count"]), int(data["n_sites"]))


def default_supercell(phase: str) -> tuple[int, int, int]:
    """Return the versioned default supercell for a crystal prototype."""

    key = str(phase).strip().upper()
    try:
        return DEFAULT_SUPERCELLS[key]
    except KeyError as exc:
        raise ValueError(f"Unsupported SQS phase: {phase!r}") from exc


@dataclass
class SqsPlan:
    system: SqsSystem
    template_structure: Structure
    points: tuple[SqsCompositionPoint, ...]
    supercell: tuple[int, int, int]
    iterations: int = 10000
    seed: int = 0
    engine: str = "sqsgenerator"
    engine_config: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.iterations <= 0:
            raise ValueError("SQS iterations must be positive")
        if len(self.supercell) != 3 or any(
            not isinstance(x, int) or isinstance(x, bool) or x <= 0 for x in self.supercell
        ):
            raise ValueError("supercell must contain three positive integers")
        if not self.points:
            raise ValueError("SQS plan must contain at least one composition point")
        expanded_sites = len(self.expanded_template())
        if any(point.n_sites != expanded_sites for point in self.points):
            raise ValueError("Every SQS composition point must match expanded template site count")

    @classmethod
    def default(
        cls,
        system: SqsSystem,
        *,
        registry: TemplateRegistry | None = None,
        supercell: tuple[int, int, int] | None = None,
        iterations: int = 10000,
        seed: int = 0,
    ) -> "SqsPlan":
        registry = registry or TemplateRegistry()
        template = registry.get(system.phase)
        cell = tuple(supercell) if supercell is not None else default_supercell(system.phase)
        expanded = template.copy()
        expanded.make_supercell(cell)
        n_sites = len(expanded)
        # Three interior points are sufficient for the current BCC_A2
        # Redlich-Kister (L0/L1) production mapping.
        fractions = (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4))
        points = tuple(SqsCompositionPoint.from_fraction(x, n_sites) for x in fractions)
        return cls(system, template, points, cell, iterations=iterations, seed=seed)

    @property
    def n_sites(self) -> int:
        return self.points[0].n_sites

    def expanded_template(self) -> Structure:
        structure = self.template_structure.copy()
        structure.make_supercell(self.supercell)
        return structure

    def targets(self) -> tuple[dict[str, Any], ...]:
        return tuple(point.as_dict() for point in self.points)

    def reference_structures(self) -> dict[str, Structure]:
        result: dict[str, Structure] = {}
        for element in self.system.canonical_elements:
            structure = self.expanded_template()
            structure.replace_species({str(site.specie): Element(element) for site in structure})
            result[element] = structure
        return result

    def write_manifest(self, path: Path | str) -> Path:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "schema_version": 1,
            "kind": "htvasp.sqs_plan",
            "system": self.system.canonical_name,
            "template": {
                "structure": self.template_structure.as_dict(),
                "structure_checksum": structure_checksum(self.template_structure),
            },
            "supercell": list(self.supercell),
            "n_sites": self.n_sites,
            "points": [point.as_dict() for point in self.points],
            "iterations": self.iterations,
            "seed": self.seed,
            "engine": self.engine,
            "engine_config": self.engine_config,
        }
        destination.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
        return destination

    @classmethod
    def read_manifest(cls, path: Path | str) -> "SqsPlan":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        if payload.get("kind") != "htvasp.sqs_plan":
            raise ValueError("Not an htvasp SQS plan manifest")
        template = Structure.from_dict(payload["template"]["structure"])
        if structure_checksum(template) != payload["template"]["structure_checksum"]:
            raise ValueError("SQS plan template checksum mismatch")
        return cls(
            system=SqsSystem.from_name(payload["system"]),
            template_structure=template,
            points=tuple(SqsCompositionPoint.from_dict(x) for x in payload["points"]),
            supercell=tuple(payload["supercell"]),
            iterations=int(payload["iterations"]),
            seed=int(payload["seed"]),
            engine=str(payload["engine"]),
            engine_config=dict(payload.get("engine_config", {})),
        )


class SqsGenerator:
    """Strict adapter around ``sqsgenerator``; no random fallback is provided."""

    def __init__(self, engine: str = "sqsgenerator"):
        self.engine = engine

    def generate(
        self,
        template: Structure,
        point: SqsCompositionPoint,
        *,
        element_a: str,
        element_b: str,
        seed: int = 0,
        iterations: int = 10000,
    ) -> tuple[Structure, dict[str, Any]]:
        if self.engine != "sqsgenerator":
            raise ValueError(f"Unsupported SQS engine: {self.engine}")
        try:
            from sqsgenerator import StructureFormat, optimize
        except ImportError as exc:
            raise RuntimeError(
                "sqsgenerator is required for SQS generation; install it explicitly"
            ) from exc

        config = {
            "iterations": int(iterations),
            "iteration_mode": "random",
            "sublattice_mode": "interact",
            "composition": {element_a: point.element_a_count, element_b: point.element_b_count},
            "structure": {
                "lattice": template.lattice.matrix.tolist(),
                "coords": template.frac_coords.tolist(),
                # The prototype species describe MP role sites, while the
                # SQS engine expects structure species to belong to the target
                # composition.  Keep the geometry and initialize all sites as
                # A; sqsgenerator then decorates exactly n_A/n_B sites.
                "species": [element_a] * len(template),
            },
        }
        try:
            # Both supported sqsgenerator releases accept a configuration
            # mapping; older releases additionally require parse_config().
            try:
                pack = optimize(config)
            except TypeError:
                from sqsgenerator import parse_config

                pack = optimize(parse_config(config))
            best = pack.best()
            structure = Poscar.from_str(best.structure().dump(StructureFormat.poscar)).structure
        except Exception as exc:
            raise RuntimeError(f"sqsgenerator failed for {point.label}: {exc}") from exc
        counts = {element_a: 0, element_b: 0}
        for site in structure:
            symbol = str(site.specie)
            if symbol in counts:
                counts[symbol] += 1
        if len(structure) != point.n_sites or counts != {
            element_a: point.element_a_count,
            element_b: point.element_b_count,
        }:
            raise RuntimeError(
                f"sqsgenerator returned composition {counts} for {point.label}; "
                f"expected {{{element_a}: {point.element_a_count}, {element_b}: {point.element_b_count}}}"
            )
        objective = getattr(best, "objective", None)
        if not isinstance(objective, (int, float)):
            objective = None
        try:
            engine_version = importlib.metadata.version("sqsgenerator")
        except importlib.metadata.PackageNotFoundError:
            engine_version = "unknown"
        metadata = {
            "engine": self.engine,
            "engine_version": engine_version,
            "seed": int(seed),
            "seed_applied": False,
            "seed_note": "sqsgenerator has no documented seed config key; seed is recorded for provenance",
            "iterations": int(iterations),
            "config": config,
            # Keep the on-disk generation metadata aligned with the formal
            # htvasp-sqs-mixing/v1 provenance field consumed by em-tdb.
            "structure_sha256": structure_checksum(structure),
            "objective": float(objective) if objective is not None else None,
        }
        return structure, metadata


def _extract_energy(result: Any, *, n_sites: int | None = None) -> float:
    """Extract a scalar total energy from common atomate2 result shapes."""

    if isinstance(result, (int, float)):
        return float(result)
    if isinstance(result, Mapping):
        for key in ("energy", "final_energy", "e_0_energy"):
            value = result.get(key)
            if isinstance(value, (int, float)):
                return float(value)
        for key in ("final_energy_per_atom", "energy_per_atom"):
            value = result.get(key)
            if isinstance(value, (int, float)):
                if n_sites is None:
                    raise ValueError(f"{key} requires n_sites to recover total energy")
                return float(value) * n_sites
        for key in ("output", "result", "taskdoc", "calcs_reversed"):
            if key in result:
                try:
                    return _extract_energy(result[key], n_sites=n_sites)
                except (KeyError, TypeError, ValueError):
                    continue
    raise ValueError("Static result contains no recognized total or per-atom energy")


def _result_structure(result: Any, fallback: Any) -> dict[str, Any]:
    if isinstance(result, Mapping):
        candidate = result.get("structure")
        if candidate is not None:
            try:
                return _structure_dict(candidate)
            except SqsMixingError:
                pass
        for key in ("output", "result", "taskdoc", "task_document"):
            nested = result.get(key)
            if isinstance(nested, Mapping):
                try:
                    return _result_structure(nested, fallback)
                except SqsMixingError:
                    pass
    return _structure_dict(fallback)


class Sublattice(BaseModel):
    model_config = ConfigDict(extra="forbid")
    site_ratio: float
    constituents: list[str]

    @field_validator("site_ratio")
    @classmethod
    def positive_ratio(cls, value: float) -> float:
        value = _finite(value, "site_ratio")
        if value <= 0:
            raise ValueError("site_ratio must be positive")
        return value

    @field_validator("constituents")
    @classmethod
    def normalize_constituents(cls, value: list[str]) -> list[str]:
        result = [str(item).strip().upper() for item in value]
        if not result or any(not item for item in result) or len(result) != len(set(result)):
            raise ValueError("constituents must be distinct, non-empty component symbols")
        return result


class CalphadMapping(BaseModel):
    """Explicit BCC_A2/RK semantic mapping carried by every dataset."""

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
    def normalize_components(cls, value: tuple[str, str]) -> tuple[str, str]:
        result = tuple(str(item).strip().upper() for item in value)
        if len(result) != 2 or not all(result) or result[0] == result[1]:
            raise ValueError("components must be two different non-empty symbols")
        return result

    @model_validator(mode="after")
    def validate_bcc_a2(self) -> "CalphadMapping":
        if len(self.sublattices) != 1 or self.sublattices[0].site_ratio != 1:
            raise ValueError("P2 requires one BCC_A2 sublattice with site_ratio=1")
        if tuple(self.sublattices[0].constituents) != self.components:
            raise ValueError("sublattice constituents must match components in order")
        return self


class CalculationProvenance(BaseModel):
    model_config = ConfigDict(extra="forbid")
    settings_fingerprint: str = Field(min_length=8)
    workflow: dict[str, Any]
    relax_incar: dict[str, Any]
    static_incar: dict[str, Any]
    vasp: dict[str, Any] = Field(default_factory=dict)
    potcar: dict[str, Any] = Field(default_factory=dict)


class SqsProvenance(BaseModel):
    model_config = ConfigDict(extra="forbid")
    engine: str = Field(min_length=1)
    engine_version: str = "unknown"
    iterations: int = Field(gt=0)
    seed: int | None = None
    config: dict[str, Any] = Field(default_factory=dict)
    objective: float | None = None
    structure_sha256: str = Field(min_length=8)


class StaticObservation(BaseModel):
    model_config = ConfigDict(extra="forbid")
    run_id: str = Field(min_length=1)
    state: Literal["successful"]
    energy_total_eV: float
    energy_per_atom_eV: float
    n_sites: int = Field(gt=0)
    task_document_sha256: str = Field(min_length=8)

    @model_validator(mode="after")
    def validate_energy(self) -> "StaticObservation":
        _finite(self.energy_total_eV, "energy_total_eV")
        _finite(self.energy_per_atom_eV, "energy_per_atom_eV")
        if not math.isclose(self.energy_total_eV / self.n_sites, self.energy_per_atom_eV, abs_tol=1e-10):
            raise ValueError("energy_per_atom_eV is inconsistent with total energy")
        return self


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

    @field_validator("x_A")
    @classmethod
    def interior_fraction(cls, value: float) -> float:
        value = _finite(value, "x_A")
        if not 0.0 < value < 1.0:
            raise ValueError("x_A must be strictly between zero and one")
        return value


class DatasetValidation(BaseModel):
    model_config = ConfigDict(extra="forbid")
    settings_fingerprint: str = Field(min_length=8)
    all_static_successful: Literal[True]
    warnings: list[str] = Field(default_factory=list)


class SqsMixingData(BaseModel):
    """Versioned JSON input contract consumed by em-tdb."""

    model_config = ConfigDict(extra="forbid")
    schema_version: Literal[SCHEMA_VERSION]
    dataset_id: str = Field(min_length=1)
    created_at: datetime
    calphad_mapping: CalphadMapping
    calculation: CalculationProvenance
    references: dict[str, PureReference]
    points: list[SqsMixingPoint] = Field(min_length=3, max_length=5)
    validation: DatasetValidation

    @model_validator(mode="after")
    def validate_dataset(self) -> "SqsMixingData":
        a, b = self.calphad_mapping.components
        if set(self.references) != {a, b}:
            raise ValueError("references must contain exactly the mapped components")
        if self.validation.settings_fingerprint != self.calculation.settings_fingerprint:
            raise ValueError("validation and calculation settings fingerprints differ")
        ea = self.references[a].static.energy_per_atom_eV
        eb = self.references[b].static.energy_per_atom_eV
        ids: set[str] = set()
        xs: set[float] = set()
        for point in self.points:
            if point.point_id in ids or point.x_A in xs or not 0 < point.x_A < 1:
                raise ValueError("SQS point IDs/x_A must be unique interior compositions")
            ids.add(point.point_id)
            xs.add(point.x_A)
            counts = {str(key).upper(): value for key, value in point.site_counts.items()}
            if set(counts) != {a, b} or any(value < 0 for value in counts.values()):
                raise ValueError("site_counts must match components and be non-negative")
            if sum(counts.values()) != point.static.n_sites:
                raise ValueError("site_counts must sum to static.n_sites")
            if not math.isclose(point.x_A, counts[a] / point.static.n_sites, abs_tol=1e-12):
                raise ValueError("x_A does not match site_counts")
            hmix = (
                point.static.energy_total_eV
                - counts[a] * ea
                - counts[b] * eb
            ) / point.static.n_sites
            if not math.isclose(point.h_mix_eV_per_atom, hmix, abs_tol=1e-10):
                raise ValueError("h_mix_eV_per_atom cannot be reproduced")
            if not math.isclose(point.h_mix_J_per_mol, hmix * EV_PER_ATOM_TO_J_PER_MOL, abs_tol=1e-7):
                raise ValueError("h_mix_J_per_mol cannot be reproduced")
        return self

    @classmethod
    def from_runs(
        cls,
        plan: SqsPlan,
        reference_static_results: Mapping[str, Any],
        point_static_results: Mapping[str, Any],
        point_structures: Mapping[str, Structure] | None = None,
        *,
        reference_structures: Mapping[str, Structure] | None = None,
        point_provenance: Mapping[str, Mapping[str, Any]] | None = None,
        dataset_id: str | None = None,
        calculation: CalculationProvenance | Mapping[str, Any] | None = None,
    ) -> "SqsMixingData":
        if plan.system.phase != "BCC":
            raise SqsMixingError("P2-A em-tdb mapping currently supports only BCC_A2 SQS datasets")
        components = tuple(item.upper() for item in plan.system.canonical_elements)
        mapping = CalphadMapping(
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
        if calculation is None:
            calculation = CalculationProvenance(
                settings_fingerprint="unspecified-settings",
                workflow={"worker": "SqsMixingWorker", "flow": "single_relax_then_static"},
                relax_incar={}, static_incar={},
            )
        calculation = CalculationProvenance.model_validate(calculation)
        reference_structures = reference_structures or plan.reference_structures()
        references: dict[str, PureReference] = {}
        for element in plan.system.canonical_elements:
            if element not in reference_static_results:
                raise ValueError(f"Missing static reference result for {element}")
            fallback = reference_structures[element]
            observation, final = cls._observation(
                reference_static_results[element], f"ref-{element}", fallback
            )
            references[element.upper()] = PureReference(
                component=element.upper(), final_structure=final, static=observation
            )
        point_structures = point_structures or {}
        point_provenance = point_provenance or {}
        points: list[SqsMixingPoint] = []
        a, b = plan.system.canonical_elements
        for point in plan.points:
            if point.label not in point_static_results or point.label not in point_structures:
                raise ValueError(f"Missing static result or structure for {point.label}")
            initial = _structure_dict(point_structures[point.label])
            observation, final = cls._observation(
                point_static_results[point.label], point.label, point_structures[point.label]
            )
            counts = {a.upper(): point.element_a_count, b.upper(): point.element_b_count}
            hmix = (
                observation.energy_total_eV
                - counts[a.upper()] * references[a.upper()].static.energy_per_atom_eV
                - counts[b.upper()] * references[b.upper()].static.energy_per_atom_eV
            ) / observation.n_sites
            provenance = dict(point_provenance.get(point.label, {}))
            provenance.setdefault("engine", plan.engine)
            provenance.setdefault("iterations", plan.iterations)
            provenance.setdefault("seed", plan.seed + plan.points.index(point))
            provenance.setdefault("config", {})
            provenance.setdefault("objective", None)
            provenance.setdefault("structure_sha256", structure_checksum(point_structures[point.label]))
            points.append(
                SqsMixingPoint(
                    point_id=point.label,
                    x_A=point.fraction_a,
                    site_counts=counts,
                    sqs=SqsProvenance.model_validate(provenance),
                    initial_structure=initial,
                    final_structure=final,
                    static=observation,
                    h_mix_eV_per_atom=hmix,
                    h_mix_J_per_mol=hmix * EV_PER_ATOM_TO_J_PER_MOL,
                )
            )
        return cls(
            schema_version=SCHEMA_VERSION,
            dataset_id=dataset_id or plan.system.canonical_name,
            created_at=datetime.now(timezone.utc),
            calphad_mapping=mapping,
            calculation=calculation,
            references=references,
            points=points,
            validation=DatasetValidation(
                settings_fingerprint=calculation.settings_fingerprint,
                all_static_successful=True,
            ),
        )

    @staticmethod
    def _observation(result: Any, run_id: str, fallback: Any) -> tuple[StaticObservation, dict[str, Any]]:
        final = _result_structure(result, fallback)
        n_sites = len(final.get("sites", []))
        if n_sites <= 0:
            raise SqsMixingError(f"Static result {run_id} has no final sites")
        energy = _finite(_extract_energy(result, n_sites=n_sites), "static energy")
        state = result.get("state", "successful") if isinstance(result, Mapping) else "successful"
        if state != "successful":
            raise SqsMixingError(f"Static task {run_id} did not succeed")
        observation = StaticObservation(
            run_id=run_id,
            state="successful",
            energy_total_eV=energy,
            energy_per_atom_eV=energy / n_sites,
            n_sites=n_sites,
            task_document_sha256=_json_checksum(result),
        )
        return observation, final

    def write_json(self, path: Path | str) -> Path:
        destination = Path(path).resolve()
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(self.model_dump_json(indent=2) + "\n", encoding="utf-8")
        return destination

    @classmethod
    def read_json(cls, path: Path | str) -> "SqsMixingData":
        return cls.model_validate_json(Path(path).read_text(encoding="utf-8"))


__all__ = [
    "SCHEMA_VERSION",
    "EV_PER_ATOM_TO_J_PER_MOL",
    "DEFAULT_SUPERCELLS",
    "CalphadMapping",
    "CalculationProvenance",
    "DatasetValidation",
    "PureReference",
    "SqsMixingError",
    "SqsCompositionPoint",
    "SqsGenerator",
    "SqsMixingData",
    "SqsMixingPoint",
    "SqsPlan",
    "SqsProvenance",
    "StaticObservation",
    "Sublattice",
    "SqsSystem",
    "default_supercell",
]
