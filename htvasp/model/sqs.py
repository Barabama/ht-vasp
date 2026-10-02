"""Pure domain models for planning SQS mixing calculations.

This module describes structures and calculation inputs.  It intentionally
does not import an SQS engine, atomate2, jobflow, Slurm, or VASP result
parsers.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Any, Literal, Mapping

from pymatgen.core import Element, Structure

from htvasp.model.templates import TemplateRecord, TemplateSnapshot, structure_checksum


DEFAULT_SUPERCELLS: dict[str, tuple[int, int, int]] = {
    "BCC": (2, 2, 2),
    "FCC": (2, 2, 2),
    "HCP": (2, 2, 1),
}


def _canonical_json(payload: Any) -> bytes:
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        default=str,
    ).encode("utf-8")


def payload_checksum(payload: Any) -> str:
    """Return a stable SHA-256 checksum for JSON-compatible data."""

    return hashlib.sha256(_canonical_json(payload)).hexdigest()


def default_supercell(phase: str) -> tuple[int, int, int]:
    key = str(phase).strip().upper()
    try:
        return DEFAULT_SUPERCELLS[key]
    except KeyError as exc:
        raise ValueError(f"Unsupported SQS phase: {phase!r}") from exc


def _normalize_element(value: str) -> str:
    text = str(value).strip()
    if not text:
        raise ValueError("SQS element symbols must be non-empty")
    symbol = text[0].upper() + text[1:].lower()
    try:
        Element(symbol)
    except ValueError as exc:
        raise ValueError(f"Invalid SQS element symbol: {value!r}") from exc
    return symbol


@dataclass(frozen=True)
class SqsSystem:
    """A binary disordered system whose component order is canonicalized."""

    phase: str
    element_a: str
    element_b: str

    def __post_init__(self) -> None:
        phase = str(self.phase).strip().upper()
        if phase not in DEFAULT_SUPERCELLS:
            raise ValueError(f"Unsupported SQS phase: {self.phase!r}")
        elements = sorted((_normalize_element(self.element_a), _normalize_element(self.element_b)))
        if elements[0] == elements[1]:
            raise ValueError("SQS requires two distinct elements")
        object.__setattr__(self, "phase", phase)
        object.__setattr__(self, "element_a", elements[0])
        object.__setattr__(self, "element_b", elements[1])

    @classmethod
    def from_name(cls, name: str) -> "SqsSystem":
        parts = str(name).strip().split("-")
        if len(parts) != 3:
            raise ValueError("SQS name must be phase-A-B, e.g. BCC-Al-Nb")
        return cls(*parts)

    @property
    def canonical_elements(self) -> tuple[str, str]:
        return self.element_a, self.element_b

    @property
    def canonical_name(self) -> str:
        return f"{self.phase}-{self.element_a}-{self.element_b}"

    @property
    def key(self) -> tuple[str, str, str]:
        return self.phase, self.element_a, self.element_b

    def as_dict(self) -> dict[str, str]:
        return {
            "phase": self.phase,
            "element_a": self.element_a,
            "element_b": self.element_b,
            "canonical_name": self.canonical_name,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "SqsSystem":
        return cls(str(data["phase"]), str(data["element_a"]), str(data["element_b"]))


@dataclass(frozen=True)
class SqsCompositionPoint:
    """One interior composition represented by integer site counts."""

    element_a_count: int
    element_b_count: int
    n_sites: int

    def __post_init__(self) -> None:
        values = (self.element_a_count, self.element_b_count, self.n_sites)
        if any(isinstance(value, bool) or not isinstance(value, int) for value in values):
            raise ValueError("SQS composition counts and n_sites must be integers")
        if self.element_a_count <= 0 or self.element_b_count <= 0:
            raise ValueError("SQS composition points must contain both elements")
        if self.n_sites <= 0 or self.element_a_count + self.element_b_count != self.n_sites:
            raise ValueError("Composition counts must sum to a positive n_sites")

    @classmethod
    def from_counts(cls, element_a_count: int, element_b_count: int) -> "SqsCompositionPoint":
        return cls(element_a_count, element_b_count, element_a_count + element_b_count)

    @classmethod
    def from_fraction(
        cls, fraction_a: float | str | Fraction, n_sites: int
    ) -> "SqsCompositionPoint":
        if isinstance(n_sites, bool) or not isinstance(n_sites, int) or n_sites <= 0:
            raise ValueError("n_sites must be a positive integer")
        fraction = fraction_a if isinstance(fraction_a, Fraction) else Fraction(str(fraction_a))
        if not 0 < fraction < 1:
            raise ValueError("fraction_a must be strictly between zero and one")
        count = fraction * n_sites
        if count.denominator != 1:
            raise ValueError(f"Fraction {fraction_a!r} cannot be represented with {n_sites} sites")
        return cls.from_counts(int(count), n_sites - int(count))

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


@dataclass(frozen=True)
class SqsStructureArtifact:
    """A generated or materialized structure reference with provenance."""

    path: str
    structure_checksum: str
    n_sites: int
    composition: Mapping[str, int]
    provenance: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not str(self.path).strip():
            raise ValueError("SQS structure artifact path must be non-empty")
        if len(self.structure_checksum) != 64:
            raise ValueError("structure_checksum must be a SHA-256 hex digest")
        if isinstance(self.n_sites, bool) or not isinstance(self.n_sites, int) or self.n_sites <= 0:
            raise ValueError("n_sites must be a positive integer")
        raw_counts = dict(self.composition)
        if any(isinstance(value, bool) or not isinstance(value, int) for value in raw_counts.values()):
            raise ValueError("composition counts must be integers")
        counts = {str(key): value for key, value in raw_counts.items()}
        if any(value < 0 for value in counts.values()) or sum(counts.values()) != self.n_sites:
            raise ValueError("composition counts must be non-negative and sum to n_sites")
        object.__setattr__(self, "composition", counts)
        object.__setattr__(self, "provenance", dict(self.provenance))

    def as_dict(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "structure_checksum": self.structure_checksum,
            "n_sites": self.n_sites,
            "composition": dict(self.composition),
            "provenance": dict(self.provenance),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "SqsStructureArtifact":
        return cls(
            path=str(data["path"]),
            structure_checksum=str(data["structure_checksum"]),
            n_sites=int(data["n_sites"]),
            composition=dict(data["composition"]),
            provenance=dict(data.get("provenance", {})),
        )


CalculationKind = Literal["endpoint", "point"]


@dataclass(frozen=True)
class SqsCalculationSpec:
    """Stable description of one endpoint or SQS-point calculation."""

    kind: CalculationKind
    calculation_id: str
    structure: SqsStructureArtifact
    phase: str
    settings_fingerprint: str
    system: SqsSystem | None = None
    point: SqsCompositionPoint | None = None
    component: str | None = None
    supercell: tuple[int, int, int] | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        kind = str(self.kind).strip().lower()
        if kind not in {"endpoint", "point"}:
            raise ValueError("calculation kind must be 'endpoint' or 'point'")
        if not str(self.calculation_id).strip():
            raise ValueError("calculation_id must be non-empty")
        phase = str(self.phase).strip().upper()
        if not phase:
            raise ValueError("calculation phase must be non-empty")
        if not str(self.settings_fingerprint).strip():
            raise ValueError("settings_fingerprint must be non-empty")
        if self.supercell is not None:
            if len(self.supercell) != 3 or any(
                isinstance(value, bool) or not isinstance(value, int) or value <= 0
                for value in self.supercell
            ):
                raise ValueError("supercell must contain three positive integers")
        if kind == "endpoint":
            if self.component is None or self.system is not None or self.point is not None:
                raise ValueError("endpoint specs require component and forbid system/point")
            component = _normalize_element(self.component)
            if self.structure.composition != {component: self.structure.n_sites}:
                raise ValueError("endpoint structure composition must contain only its component")
            expected_id = f"{phase}-{component}"
            object.__setattr__(self, "component", component)
        else:
            if self.system is None or self.point is None or self.component is not None:
                raise ValueError("point specs require system and point and forbid component")
            if self.system.phase != phase:
                raise ValueError("spec phase must match point system phase")
            if self.point.n_sites != self.structure.n_sites:
                raise ValueError("point composition must match structure site count")
            expected_composition = {
                self.system.element_a: self.point.element_a_count,
                self.system.element_b: self.point.element_b_count,
            }
            if self.structure.composition != expected_composition:
                raise ValueError("point structure composition must match the composition point")
            expected_id = f"{self.system.canonical_name}-{self.point.label}"
        if self.calculation_id != expected_id:
            raise ValueError(f"calculation_id must be {expected_id!r} for this spec")
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "phase", phase)
        object.__setattr__(self, "metadata", dict(self.metadata))

    @property
    def stable_id(self) -> str:
        return self.calculation_id

    @property
    def fingerprint(self) -> str:
        return payload_checksum(self.as_dict(include_id=False))

    def as_dict(self, *, include_id: bool = True) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "kind": self.kind,
            "phase": self.phase,
            "settings_fingerprint": self.settings_fingerprint,
            "structure": self.structure.as_dict(),
            "system": self.system.as_dict() if self.system is not None else None,
            "point": self.point.as_dict() if self.point is not None else None,
            "component": self.component,
            "supercell": list(self.supercell) if self.supercell is not None else None,
            "metadata": dict(self.metadata),
        }
        if include_id:
            payload["calculation_id"] = self.calculation_id
        return payload

    def write_json(self, path: Path | str) -> Path:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(self.as_dict(), indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
        return destination

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "SqsCalculationSpec":
        return cls(
            kind=data["kind"],
            calculation_id=str(data["calculation_id"]),
            structure=SqsStructureArtifact.from_dict(data["structure"]),
            phase=str(data["phase"]),
            settings_fingerprint=str(data["settings_fingerprint"]),
            system=SqsSystem.from_dict(data["system"]) if data.get("system") else None,
            point=SqsCompositionPoint.from_dict(data["point"]) if data.get("point") else None,
            component=data.get("component"),
            supercell=tuple(data["supercell"]) if data.get("supercell") else None,
            metadata=dict(data.get("metadata", {})),
        )

    @classmethod
    def endpoint(
        cls,
        *,
        phase: str,
        component: str,
        structure: SqsStructureArtifact,
        settings_fingerprint: str,
        supercell: tuple[int, int, int] | None = None,
        metadata: Mapping[str, Any] | None = None,
    ) -> "SqsCalculationSpec":
        normalized_phase = str(phase).strip().upper()
        normalized_component = _normalize_element(component)
        return cls(
            kind="endpoint",
            calculation_id=f"{normalized_phase}-{normalized_component}",
            structure=structure,
            phase=normalized_phase,
            settings_fingerprint=settings_fingerprint,
            component=normalized_component,
            supercell=supercell,
            metadata=metadata or {},
        )

    @classmethod
    def point_spec(
        cls,
        *,
        system: SqsSystem,
        point: SqsCompositionPoint,
        structure: SqsStructureArtifact,
        settings_fingerprint: str,
        supercell: tuple[int, int, int] | None = None,
        metadata: Mapping[str, Any] | None = None,
    ) -> "SqsCalculationSpec":
        return cls(
            kind="point",
            calculation_id=f"{system.canonical_name}-{point.label}",
            structure=structure,
            phase=system.phase,
            settings_fingerprint=settings_fingerprint,
            system=system,
            point=point,
            supercell=supercell,
            metadata=metadata or {},
        )


@dataclass(frozen=True)
class SqsPlan:
    """Serializable calculation plan for one binary SQS system."""

    system: SqsSystem
    template: TemplateSnapshot
    points: tuple[SqsCompositionPoint, ...]
    supercell: tuple[int, int, int]
    iterations: int = 10000
    seed: int = 0
    engine: str = "unspecified"
    engine_config: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if len(self.supercell) != 3 or any(
            isinstance(value, bool) or not isinstance(value, int) or value <= 0
            for value in self.supercell
        ):
            raise ValueError("supercell must contain three positive integers")
        if isinstance(self.iterations, bool) or not isinstance(self.iterations, int) or self.iterations <= 0:
            raise ValueError("iterations must be a positive integer")
        if not self.points:
            raise ValueError("SQS plan must contain at least one composition point")
        if self.template.phase != self.system.phase:
            raise ValueError("template phase must match SQS system phase")
        expanded_sites = len(self.expanded_template)
        if any(point.n_sites != expanded_sites for point in self.points):
            raise ValueError("Every SQS composition point must match expanded template site count")
        object.__setattr__(self, "engine_config", dict(self.engine_config))

    @classmethod
    def default(
        cls,
        system: SqsSystem,
        template: TemplateSnapshot,
        *,
        supercell: tuple[int, int, int] | None = None,
        iterations: int = 10000,
        seed: int = 0,
        engine: str = "unspecified",
        engine_config: Mapping[str, Any] | None = None,
    ) -> "SqsPlan":
        cell = tuple(supercell) if supercell is not None else default_supercell(system.phase)
        expanded = template.structure.copy()
        expanded.make_supercell(cell)
        points = tuple(SqsCompositionPoint.from_fraction(fraction, len(expanded)) for fraction in (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4)))
        return cls(
            system=system,
            template=template,
            points=points,
            supercell=cell,
            iterations=iterations,
            seed=seed,
            engine=engine,
            engine_config=engine_config or {},
        )

    @property
    def expanded_template(self) -> Structure:
        structure = self.template.structure.copy()
        structure.make_supercell(self.supercell)
        return structure

    @property
    def n_sites(self) -> int:
        return len(self.expanded_template)

    @property
    def stable_id(self) -> str:
        return self.system.canonical_name

    @property
    def fingerprint(self) -> str:
        return payload_checksum(self.as_dict(include_fingerprint=False))

    def as_dict(self, *, include_fingerprint: bool = True) -> dict[str, Any]:
        payload = {
            "schema_version": 1,
            "kind": "htvasp.sqs_plan",
            "system": self.system.as_dict(),
            "template": self.template.as_dict(),
            "supercell": list(self.supercell),
            "n_sites": self.n_sites,
            "points": [point.as_dict() for point in self.points],
            "iterations": self.iterations,
            "seed": self.seed,
            "engine": self.engine,
            "engine_config": dict(self.engine_config),
        }
        if include_fingerprint:
            payload["fingerprint"] = self.fingerprint
        return payload

    def write_manifest(self, path: Path | str) -> Path:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(self.as_dict(), indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
        return destination

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "SqsPlan":
        if data.get("kind") != "htvasp.sqs_plan":
            raise ValueError("Not an htvasp SQS plan manifest")
        template_data = data["template"]
        template = TemplateSnapshot(
            record=TemplateRecord.from_dict(template_data["record"]),
            structure=Structure.from_dict(template_data["structure"]),
        )
        plan = cls(
            system=SqsSystem.from_dict(data["system"]),
            template=template,
            points=tuple(SqsCompositionPoint.from_dict(point) for point in data["points"]),
            supercell=tuple(data["supercell"]),
            iterations=int(data["iterations"]),
            seed=int(data["seed"]),
            engine=str(data["engine"]),
            engine_config=dict(data.get("engine_config", {})),
        )
        fingerprint = data.get("fingerprint")
        if fingerprint is not None and fingerprint != plan.fingerprint:
            raise ValueError("SQS plan fingerprint mismatch")
        return plan

    @classmethod
    def read_manifest(cls, path: Path | str) -> "SqsPlan":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


__all__ = [
    "CalculationKind",
    "DEFAULT_SUPERCELLS",
    "SqsCalculationSpec",
    "SqsCompositionPoint",
    "SqsPlan",
    "SqsStructureArtifact",
    "SqsSystem",
    "default_supercell",
    "payload_checksum",
    "structure_checksum",
]
