"""File-backed SQS stage orchestration.

This module prepares and validates SQS batches without owning VASP execution.
The actual calculation remains delegated to :class:`SqsMixingWorker` in a
later execution step.  The persisted manifest is the contract between stages.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Callable, Literal, Mapping, Sequence

from pymatgen.core import Structure

from htvasp.model.sqs import (
    SqsCalculationSpec,
    SqsPlan,
    SqsStructureArtifact,
    payload_checksum,
)
from htvasp.sqs.generator import GeneratedSqs, SqsGenerator, SqsgeneratorAdapter
from htvasp.sqs.task_doc import SqsStaticObservation
from htvasp.slurm.manager import atomic_json_write, file_lock


BATCH_SCHEMA_VERSION = 1
BATCH_KIND = "htvasp.sqs_batch"
TaskState = Literal["pending", "successful", "failed"]


class SqsOrchestrationError(RuntimeError):
    """Raised when a persisted batch cannot safely advance stages."""


def endpoint_cache_key(
    *,
    phase: str,
    component: str,
    supercell: Sequence[int],
    template_checksum: str,
    settings_fingerprint: str,
) -> str:
    """Return a cache key containing every input that affects an endpoint."""

    payload = {
        "phase": str(phase).upper(),
        "component": str(component),
        "supercell": list(supercell),
        "template_checksum": str(template_checksum),
        "settings_fingerprint": str(settings_fingerprint),
    }
    digest = payload_checksum(payload)[:16]
    return f"{payload['phase']}-{payload['component']}-{digest}"


def _json_dump(path: Path, payload: Mapping[str, Any]) -> None:
    atomic_json_write(path, payload)


def _artifact_path(root: Path, artifact: SqsStructureArtifact) -> Path:
    path = Path(artifact.path)
    if path.is_absolute():
        resolved = path.resolve()
    else:
        resolved = (root / path).resolve()
    try:
        resolved.relative_to(root.resolve())
    except ValueError as exc:
        raise SqsOrchestrationError(f"Artifact path escapes batch root: {artifact.path}") from exc
    return resolved


def _safe_path(root: Path, relative_path: str) -> Path:
    """Resolve a persisted relative path without allowing root escape."""

    resolved = (root / Path(relative_path)).resolve()
    try:
        resolved.relative_to(root.resolve())
    except ValueError as exc:
        raise SqsOrchestrationError(f"Persisted path escapes batch root: {relative_path}") from exc
    return resolved


def _structure_composition(structure: Structure) -> dict[str, int]:
    return {
        str(element): int(amount)
        for element, amount in structure.composition.get_el_amt_dict().items()
    }


def _write_structure(
    root: Path,
    artifact: SqsStructureArtifact,
    structure: Structure,
    *,
    overwrite: bool = False,
) -> None:
    destination = _artifact_path(root, artifact)
    if destination.exists():
        try:
            existing = Structure.from_file(destination)
        except Exception as exc:
            raise SqsOrchestrationError(f"Cannot read existing artifact: {destination}") from exc
        from htvasp.model.templates import structure_checksum

        if structure_checksum(existing) != artifact.structure_checksum and not overwrite:
            raise SqsOrchestrationError(f"Existing artifact checksum mismatch: {destination}")
        if not overwrite:
            return
    destination.parent.mkdir(parents=True, exist_ok=True)
    structure.to(filename=destination, fmt="poscar")


def _decorated_endpoint(template: Structure, element: str) -> Structure:
    structure = template.copy()
    for index in range(len(structure)):
        structure.replace(index, element)
    return structure


@dataclass(frozen=True)
class SqsEndpointRecord:
    """One deduplicated endpoint and its stage state."""

    cache_key: str
    spec: SqsCalculationSpec
    consumers: tuple[str, ...] = ()
    state: TaskState = "pending"
    observation_path: str | None = None
    failure_path: str | None = None

    def __post_init__(self) -> None:
        if self.spec.kind != "endpoint":
            raise ValueError("Endpoint records require an endpoint calculation spec")
        if not self.cache_key:
            raise ValueError("Endpoint cache_key must be non-empty")
        if self.state == "successful" and not self.observation_path:
            raise ValueError("Successful endpoints require observation_path")
        if self.state == "failed" and not self.failure_path:
            raise ValueError("Failed endpoints require failure_path")
        if self.state not in {"pending", "successful", "failed"}:
            raise ValueError(f"Unsupported endpoint state: {self.state!r}")
        object.__setattr__(self, "consumers", tuple(sorted(set(self.consumers))))

    def as_dict(self) -> dict[str, Any]:
        return {
            "cache_key": self.cache_key,
            "spec": self.spec.as_dict(),
            "consumers": list(self.consumers),
            "state": self.state,
            "observation_path": self.observation_path,
            "failure_path": self.failure_path,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "SqsEndpointRecord":
        return cls(
            cache_key=str(data["cache_key"]),
            spec=SqsCalculationSpec.from_dict(data["spec"]),
            consumers=tuple(data.get("consumers", ())),
            state=data.get("state", "pending"),
            observation_path=data.get("observation_path"),
            failure_path=data.get("failure_path"),
        )


@dataclass(frozen=True)
class SqsPointRecord:
    """One SQS point calculation and its persisted result state."""

    spec: SqsCalculationSpec
    state: TaskState = "pending"
    observation_path: str | None = None
    failure_path: str | None = None

    def __post_init__(self) -> None:
        if self.spec.kind != "point":
            raise ValueError("Point records require a point calculation spec")
        if self.state == "successful" and not self.observation_path:
            raise ValueError("Successful points require observation_path")
        if self.state == "failed" and not self.failure_path:
            raise ValueError("Failed points require failure_path")
        if self.state not in {"pending", "successful", "failed"}:
            raise ValueError(f"Unsupported point state: {self.state!r}")

    @property
    def calculation_id(self) -> str:
        return self.spec.calculation_id

    def as_dict(self) -> dict[str, Any]:
        return {
            "spec": self.spec.as_dict(),
            "state": self.state,
            "observation_path": self.observation_path,
            "failure_path": self.failure_path,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "SqsPointRecord":
        return cls(
            spec=SqsCalculationSpec.from_dict(data["spec"]),
            state=data.get("state", "pending"),
            observation_path=data.get("observation_path"),
            failure_path=data.get("failure_path"),
        )


@dataclass(frozen=True)
class SqsSystemRecord:
    """Persisted plan plus generated point calculation specifications."""

    plan: SqsPlan
    settings_fingerprint: str
    points: tuple[SqsPointRecord, ...]
    endpoint_cache_keys: tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.settings_fingerprint.strip():
            raise ValueError("settings_fingerprint must be non-empty")
        if len(self.points) != len(self.plan.points):
            raise ValueError("Every SQS point must have one point record")
        object.__setattr__(self, "endpoint_cache_keys", tuple(sorted(set(self.endpoint_cache_keys))))

    @property
    def system_id(self) -> str:
        return self.plan.stable_id

    @property
    def point_specs(self) -> tuple[SqsCalculationSpec, ...]:
        return tuple(record.spec for record in self.points)

    def pending_points(self) -> tuple[SqsPointRecord, ...]:
        return tuple(record for record in self.points if record.state == "pending")

    def as_dict(self) -> dict[str, Any]:
        return {
            "system_id": self.system_id,
            "plan": self.plan.as_dict(),
            "settings_fingerprint": self.settings_fingerprint,
            "points": [point.as_dict() for point in self.points],
            "endpoint_cache_keys": list(self.endpoint_cache_keys),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "SqsSystemRecord":
        raw_points = data.get("points")
        if raw_points is None:
            raw_points = [{"spec": spec} for spec in data.get("point_specs", ())]
        return cls(
            plan=SqsPlan.from_dict(data["plan"]),
            settings_fingerprint=str(data["settings_fingerprint"]),
            points=tuple(SqsPointRecord.from_dict(x) for x in raw_points),
            endpoint_cache_keys=tuple(data["endpoint_cache_keys"]),
        )


@dataclass(frozen=True)
class SqsBatchManifest:
    """Serializable state shared by SQS generation and calculation stages."""

    systems: tuple[SqsSystemRecord, ...]
    endpoints: tuple[SqsEndpointRecord, ...]
    schema_version: int = BATCH_SCHEMA_VERSION
    calculation: Mapping[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "kind": BATCH_KIND,
            "systems": [record.as_dict() for record in self.systems],
            "endpoints": [record.as_dict() for record in self.endpoints],
            "calculation": dict(self.calculation),
        }

    @property
    def fingerprint(self) -> str:
        return payload_checksum(self.as_dict())

    @property
    def batch_fingerprint(self) -> str:
        """Fingerprint immutable batch intent, excluding runtime task state.

        ``fingerprint`` remains a complete snapshot hash for diagnostics.  A
        submission must instead use this value because endpoint/point states
        and observation paths change as a batch advances.
        """

        systems = []
        for system in self.systems:
            systems.append(
                {
                    "system_id": system.system_id,
                    "plan": system.plan.as_dict(),
                    "settings_fingerprint": system.settings_fingerprint,
                    "endpoint_cache_keys": list(system.endpoint_cache_keys),
                    "point_specs": [point.spec.as_dict() for point in system.points],
                }
            )
        endpoints = [
            {
                "cache_key": endpoint.cache_key,
                "spec": endpoint.spec.as_dict(),
                "consumers": list(endpoint.consumers),
            }
            for endpoint in self.endpoints
        ]
        return payload_checksum(
            {
                "schema_version": self.schema_version,
                "kind": BATCH_KIND,
                "systems": systems,
                "endpoints": endpoints,
                "calculation": dict(self.calculation),
            }
        )

    @property
    def intent_fingerprint(self) -> str:
        """Compatibility alias for callers that call the batch hash intent."""

        return self.batch_fingerprint

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "SqsBatchManifest":
        if data.get("kind") != BATCH_KIND or int(data.get("schema_version", -1)) != BATCH_SCHEMA_VERSION:
            raise SqsOrchestrationError("Unsupported SQS batch manifest")
        manifest = cls(
            systems=tuple(SqsSystemRecord.from_dict(x) for x in data["systems"]),
            endpoints=tuple(SqsEndpointRecord.from_dict(x) for x in data["endpoints"]),
            calculation=dict(data.get("calculation", {})),
        )
        return manifest

    @classmethod
    def read_json(cls, path: Path | str) -> "SqsBatchManifest":
        source = Path(path)
        try:
            data = json.loads(source.read_text(encoding="utf-8"))
        except FileNotFoundError as exc:
            raise SqsOrchestrationError(f"SQS batch manifest does not exist: {source}") from exc
        return cls.from_dict(data)

    def write_json(self, path: Path | str) -> Path:
        destination = Path(path)
        _json_dump(destination, self.as_dict())
        return destination

    def _record_map(self) -> dict[str, SqsEndpointRecord]:
        return {record.cache_key: record for record in self.endpoints}

    def validation_errors(self, root: Path | None = None) -> list[str]:
        errors: list[str] = []
        system_ids = [record.system_id for record in self.systems]
        if len(system_ids) != len(set(system_ids)):
            errors.append("duplicate SQS system IDs")
        endpoint_map = self._record_map()
        if len(endpoint_map) != len(self.endpoints):
            errors.append("duplicate endpoint cache keys")
        for record in self.systems:
            if len(record.endpoint_cache_keys) != 2:
                errors.append(f"{record.system_id}: expected two endpoint cache keys")
            for key in record.endpoint_cache_keys:
                endpoint = endpoint_map.get(key)
                if endpoint is None:
                    errors.append(f"{record.system_id}: missing endpoint {key}")
                elif record.system_id not in endpoint.consumers:
                    errors.append(f"{key}: missing consumer {record.system_id}")
            for spec in record.point_specs:
                if spec.settings_fingerprint != record.settings_fingerprint:
                    errors.append(f"{spec.calculation_id}: settings fingerprint mismatch")
            for point in record.points:
                if point.spec.settings_fingerprint != record.settings_fingerprint:
                    errors.append(f"{point.calculation_id}: settings fingerprint mismatch")
                if root is not None and point.state == "successful":
                    errors.extend(_validate_point_observation(root, point))
        for endpoint in self.endpoints:
            template_checksum = endpoint.spec.metadata.get("template_checksum")
            if template_checksum is None:
                errors.append(f"{endpoint.cache_key}: missing template checksum")
                continue
            expected_key = endpoint_cache_key(
                phase=endpoint.spec.phase,
                component=str(endpoint.spec.component),
                supercell=endpoint.spec.supercell or (),
                template_checksum=str(template_checksum),
                settings_fingerprint=endpoint.spec.settings_fingerprint,
            )
            if expected_key != endpoint.cache_key:
                errors.append(f"{endpoint.cache_key}: cache key does not match endpoint inputs")
            if root is not None and endpoint.state == "successful":
                errors.extend(_validate_observation(root, endpoint))
        if root is not None:
            for endpoint in self.endpoints:
                errors.extend(_validate_artifact(root, endpoint.spec.structure, endpoint.spec.calculation_id))
            for system in self.systems:
                for point in system.points:
                    errors.extend(_validate_artifact(root, point.spec.structure, point.calculation_id))
                try:
                    plan_manifest = root / "systems" / system.system_id / "manifest.json"
                    if not plan_manifest.is_file():
                        errors.append(f"{system.system_id}: missing system manifest")
                except OSError as exc:
                    errors.append(f"{system.system_id}: cannot inspect system manifest: {exc}")
        return errors

    def assert_valid(self, root: Path | None = None) -> None:
        errors = self.validation_errors(root)
        if errors:
            raise SqsOrchestrationError("Invalid SQS batch manifest:\n- " + "\n- ".join(errors))

    def assert_points_ready(
        self,
        root: Path | None = None,
        *,
        expected_settings_fingerprint: str | None = None,
    ) -> None:
        self.assert_valid(root)
        if expected_settings_fingerprint is not None:
            for system in self.systems:
                if system.settings_fingerprint != expected_settings_fingerprint:
                    raise SqsOrchestrationError(
                        f"{system.system_id}: settings fingerprint does not match requested run"
                    )
        pending = [record.cache_key for record in self.endpoints if record.state != "successful"]
        if pending:
            raise SqsOrchestrationError(
                "Cannot start SQS points before all endpoints succeed: " + ", ".join(pending)
            )

    def pending_endpoint_specs(self) -> tuple[SqsCalculationSpec, ...]:
        return tuple(record.spec for record in self.endpoints if record.state == "pending")

    def point_specs(self, root: Path | None = None) -> tuple[SqsCalculationSpec, ...]:
        self.assert_points_ready(root)
        return tuple(
            point.spec
            for system in self.systems
            for point in system.points
            if point.state == "pending"
        )

    def status(self) -> dict[str, int]:
        successful = sum(record.state == "successful" for record in self.endpoints)
        return {
            "systems": len(self.systems),
            "endpoints": len(self.endpoints),
            "endpoints_successful": successful,
            "endpoints_pending": sum(record.state == "pending" for record in self.endpoints),
            "endpoints_failed": sum(record.state == "failed" for record in self.endpoints),
            "points": sum(len(system.point_specs) for system in self.systems),
            "points_successful": sum(
                point.state == "successful" for system in self.systems for point in system.points
            ),
            "points_pending": sum(
                point.state == "pending" for system in self.systems for point in system.points
            ),
            "points_failed": sum(
                point.state == "failed" for system in self.systems for point in system.points
            ),
        }

    def with_endpoint_success(
        self,
        cache_key: str,
        observation: Mapping[str, Any] | SqsStaticObservation,
        *,
        root: Path,
    ) -> "SqsBatchManifest":
        endpoint_map = self._record_map()
        record = endpoint_map.get(cache_key)
        if record is None:
            raise SqsOrchestrationError(f"Unknown endpoint cache key: {cache_key}")
        normalized = (
            observation
            if isinstance(observation, SqsStaticObservation)
            else SqsStaticObservation.model_validate(observation)
        )
        if normalized.batch_fingerprint is None:
            normalized = normalized.model_copy(
                update={"batch_fingerprint": self.batch_fingerprint}
            )
        elif normalized.batch_fingerprint != self.batch_fingerprint:
            raise SqsOrchestrationError(
                f"{cache_key}: observation batch fingerprint mismatch"
            )
        if normalized.calculation_id != record.spec.calculation_id:
            raise SqsOrchestrationError(
                f"Observation ID {normalized.calculation_id!r} does not match {record.spec.calculation_id!r}"
            )
        if normalized.n_sites != record.spec.structure.n_sites:
            raise SqsOrchestrationError(f"{cache_key}: observation site count mismatch")
        try:
            final_structure = Structure.from_dict(normalized.final_structure)
        except Exception as exc:  # pragma: no cover - model validation normally catches this
            raise SqsOrchestrationError(f"{cache_key}: invalid final structure") from exc
        if _structure_composition(final_structure) != dict(record.spec.structure.composition):
            raise SqsOrchestrationError(f"{cache_key}: observation composition mismatch")
        observation_path = Path("endpoints") / f"{cache_key}.observation.json"
        _json_dump(root / observation_path, normalized.model_dump(mode="json"))
        updated = replace(record, state="successful", observation_path=str(observation_path))
        return replace(
            self,
            endpoints=tuple(updated if item.cache_key == cache_key else item for item in self.endpoints),
        )

    def with_endpoint_failure(
        self,
        cache_key: str,
        error: str,
        *,
        root: Path,
    ) -> "SqsBatchManifest":
        endpoint = self._record_map().get(cache_key)
        if endpoint is None:
            raise SqsOrchestrationError(f"Unknown endpoint cache key: {cache_key}")
        failure_path = Path("endpoints") / f"{cache_key}.failure.json"
        _json_dump(root / failure_path, {"calculation_id": endpoint.spec.calculation_id, "error": error})
        updated = replace(endpoint, state="failed", failure_path=str(failure_path))
        return replace(
            self,
            endpoints=tuple(updated if item.cache_key == cache_key else item for item in self.endpoints),
        )

    def with_point_success(
        self,
        calculation_id: str,
        observation: Mapping[str, Any] | SqsStaticObservation,
        *,
        root: Path,
    ) -> "SqsBatchManifest":
        normalized = (
            observation
            if isinstance(observation, SqsStaticObservation)
            else SqsStaticObservation.model_validate(observation)
        )
        if normalized.batch_fingerprint is None:
            normalized = normalized.model_copy(
                update={"batch_fingerprint": self.batch_fingerprint}
            )
        elif normalized.batch_fingerprint != self.batch_fingerprint:
            raise SqsOrchestrationError(
                f"{calculation_id}: observation batch fingerprint mismatch"
            )
        if normalized.calculation_id != calculation_id:
            raise SqsOrchestrationError(
                f"Observation ID {normalized.calculation_id!r} does not match {calculation_id!r}"
            )
        for system in self.systems:
            for point in system.points:
                if point.calculation_id != calculation_id:
                    continue
                if normalized.n_sites != point.spec.structure.n_sites:
                    raise SqsOrchestrationError(f"{calculation_id}: observation site count mismatch")
                final_structure = Structure.from_dict(normalized.final_structure)
                if _structure_composition(final_structure) != dict(point.spec.structure.composition):
                    raise SqsOrchestrationError(f"{calculation_id}: observation composition mismatch")
                observation_path = Path("points") / f"{calculation_id}.observation.json"
                _json_dump(root / observation_path, normalized.model_dump(mode="json"))
                updated = replace(point, state="successful", observation_path=str(observation_path))
                updated_system = replace(
                    system,
                    points=tuple(updated if item.calculation_id == calculation_id else item for item in system.points),
                )
                return replace(
                    self,
                    systems=tuple(updated_system if item.system_id == system.system_id else item for item in self.systems),
                )
        raise SqsOrchestrationError(f"Unknown point calculation ID: {calculation_id}")

    def with_point_failure(
        self,
        calculation_id: str,
        error: str,
        *,
        root: Path,
    ) -> "SqsBatchManifest":
        for system in self.systems:
            for point in system.points:
                if point.calculation_id != calculation_id:
                    continue
                failure_path = Path("points") / f"{calculation_id}.failure.json"
                _json_dump(root / failure_path, {"calculation_id": calculation_id, "error": error})
                updated = replace(point, state="failed", failure_path=str(failure_path))
                updated_system = replace(
                    system,
                    points=tuple(updated if item.calculation_id == calculation_id else item for item in system.points),
                )
                return replace(
                    self,
                    systems=tuple(updated_system if item.system_id == system.system_id else item for item in self.systems),
                )
        raise SqsOrchestrationError(f"Unknown point calculation ID: {calculation_id}")


def _validate_artifact(root: Path, artifact: SqsStructureArtifact, label: str) -> list[str]:
    errors: list[str] = []
    try:
        path = _artifact_path(root, artifact)
    except SqsOrchestrationError as exc:
        return [f"{label}: {exc}"]
    if not path.is_file():
        return [f"{label}: missing structure artifact {artifact.path}"]
    try:
        structure = Structure.from_file(path)
    except Exception as exc:
        return [f"{label}: unreadable structure artifact {artifact.path}: {exc}"]
    from htvasp.model.templates import structure_checksum

    if len(structure) != artifact.n_sites:
        errors.append(f"{label}: n_sites mismatch")
    if structure_checksum(structure) != artifact.structure_checksum:
        errors.append(f"{label}: structure checksum mismatch")
    if _structure_composition(structure) != dict(artifact.composition):
        errors.append(f"{label}: composition mismatch")
    return errors


def _validate_observation(
    root: Path,
    endpoint: SqsEndpointRecord,
    *,
    expected_batch_fingerprint: str | None = None,
) -> list[str]:
    if not endpoint.observation_path:
        return [f"{endpoint.cache_key}: successful endpoint has no observation path"]
    try:
        path = _safe_path(root, endpoint.observation_path)
    except SqsOrchestrationError as exc:
        return [f"{endpoint.cache_key}: {exc}"]
    if not path.is_file():
        return [f"{endpoint.cache_key}: missing observation {endpoint.observation_path}"]
    try:
        observation = SqsStaticObservation.model_validate(
            json.loads(path.read_text(encoding="utf-8"))
        )
        structure = Structure.from_dict(observation.final_structure)
    except Exception as exc:
        return [f"{endpoint.cache_key}: invalid observation: {exc}"]
    errors: list[str] = []
    if (
        expected_batch_fingerprint is not None
        and observation.batch_fingerprint != expected_batch_fingerprint
    ):
        errors.append(f"{endpoint.cache_key}: observation batch fingerprint mismatch")
    if observation.calculation_id != endpoint.spec.calculation_id:
        errors.append(f"{endpoint.cache_key}: observation calculation ID mismatch")
    if observation.n_sites != endpoint.spec.structure.n_sites:
        errors.append(f"{endpoint.cache_key}: observation n_sites mismatch")
    if _structure_composition(structure) != dict(endpoint.spec.structure.composition):
        errors.append(f"{endpoint.cache_key}: observation composition mismatch")
    return errors


def _validate_point_observation(
    root: Path,
    point: SqsPointRecord,
    *,
    expected_batch_fingerprint: str | None = None,
) -> list[str]:
    if not point.observation_path:
        return [f"{point.calculation_id}: successful point has no observation path"]
    try:
        path = _safe_path(root, point.observation_path)
    except SqsOrchestrationError as exc:
        return [f"{point.calculation_id}: {exc}"]
    if not path.is_file():
        return [f"{point.calculation_id}: missing observation {point.observation_path}"]
    try:
        observation = SqsStaticObservation.model_validate(
            json.loads(path.read_text(encoding="utf-8"))
        )
        structure = Structure.from_dict(observation.final_structure)
    except Exception as exc:
        return [f"{point.calculation_id}: invalid observation: {exc}"]
    errors: list[str] = []
    if (
        expected_batch_fingerprint is not None
        and observation.batch_fingerprint != expected_batch_fingerprint
    ):
        errors.append(f"{point.calculation_id}: observation batch fingerprint mismatch")
    if observation.calculation_id != point.calculation_id:
        errors.append(f"{point.calculation_id}: observation calculation ID mismatch")
    if observation.n_sites != point.spec.structure.n_sites:
        errors.append(f"{point.calculation_id}: observation n_sites mismatch")
    if _structure_composition(structure) != dict(point.spec.structure.composition):
        errors.append(f"{point.calculation_id}: observation composition mismatch")
    return errors


class SqsOrchestrator:
    """Create and gate a file-backed SQS batch."""

    def __init__(
        self,
        root: Path | str,
        *,
        generator: SqsGenerator | None = None,
        manifest_name: str = "batch.json",
    ) -> None:
        self.root = Path(root).resolve()
        self.generator = generator or SqsgeneratorAdapter()
        self.manifest_path = self.root / manifest_name

    def generate(
        self,
        plans: Sequence[SqsPlan],
        *,
        settings_fingerprint: str,
        calculation: Mapping[str, Any] | None = None,
        overwrite: bool = False,
    ) -> SqsBatchManifest:
        """Generate a batch under a cross-process lock."""

        with file_lock(self.root / ".generate.lock"):
            return self._generate_unlocked(
                plans,
                settings_fingerprint=settings_fingerprint,
                calculation=calculation,
                overwrite=overwrite,
            )

    def _generate_unlocked(
        self,
        plans: Sequence[SqsPlan],
        *,
        settings_fingerprint: str,
        calculation: Mapping[str, Any] | None = None,
        overwrite: bool = False,
    ) -> SqsBatchManifest:
        """Generate structures, deduplicate endpoints, and persist a manifest."""

        if not plans:
            raise ValueError("At least one SQS plan is required")
        if not settings_fingerprint.strip():
            raise ValueError("settings_fingerprint must be non-empty")
        if self.manifest_path.exists() and not overwrite:
            raise SqsOrchestrationError(
                f"SQS batch already exists: {self.manifest_path}; pass overwrite=True to rebuild"
            )
        systems: list[SqsSystemRecord] = []
        endpoints: dict[str, SqsEndpointRecord] = {}
        system_ids: set[str] = set()
        self.root.mkdir(parents=True, exist_ok=True)

        for plan in plans:
            if plan.stable_id in system_ids:
                raise SqsOrchestrationError(f"Duplicate SQS system: {plan.stable_id}")
            system_ids.add(plan.stable_id)
            system_root = self.root / "systems" / plan.stable_id
            points: list[SqsPointRecord] = []
            endpoint_keys: list[str] = []
            expanded = plan.expanded_template

            for component in plan.system.canonical_elements:
                key = endpoint_cache_key(
                    phase=plan.system.phase,
                    component=component,
                    supercell=plan.supercell,
                    template_checksum=plan.template.checksum,
                    settings_fingerprint=settings_fingerprint,
                )
                endpoint_keys.append(key)
                if key not in endpoints:
                    endpoint_structure = _decorated_endpoint(expanded, component)
                    endpoint_artifact = SqsStructureArtifact(
                        path=f"endpoints/{key}.vasp",
                        structure_checksum=_structure_checksum(endpoint_structure),
                        n_sites=len(endpoint_structure),
                        composition={component: len(endpoint_structure)},
                        provenance={"kind": "endpoint", "template_checksum": plan.template.checksum},
                    )
                    _write_structure(
                        self.root,
                        endpoint_artifact,
                        endpoint_structure,
                        overwrite=overwrite,
                    )
                    spec = SqsCalculationSpec.endpoint(
                        phase=plan.system.phase,
                        component=component,
                        structure=endpoint_artifact,
                        settings_fingerprint=settings_fingerprint,
                        supercell=plan.supercell,
                        metadata={
                            "template_checksum": plan.template.checksum,
                            "cache_key": key,
                        },
                    )
                    endpoints[key] = SqsEndpointRecord(cache_key=key, spec=spec)
                endpoints[key] = replace(
                    endpoints[key],
                    consumers=tuple((*endpoints[key].consumers, plan.stable_id)),
                )

            for point in plan.points:
                generated: GeneratedSqs = self.generator.generate(
                    expanded,
                    point,
                    element_a=plan.system.element_a,
                    element_b=plan.system.element_b,
                    iterations=plan.iterations,
                    seed=plan.seed,
                    engine_config=plan.engine_config,
                )
                artifact = SqsStructureArtifact(
                    path=f"systems/{plan.stable_id}/points/{point.label}.vasp",
                    structure_checksum=_structure_checksum(generated.structure),
                    n_sites=len(generated.structure),
                    composition={
                        plan.system.element_a: point.element_a_count,
                        plan.system.element_b: point.element_b_count,
                    },
                    provenance=dict(generated.provenance),
                )
                artifact.provenance.setdefault("structure_checksum", artifact.structure_checksum)
                artifact.provenance.setdefault("iterations", plan.iterations)
                artifact.provenance.setdefault("seed", plan.seed)
                artifact.provenance.setdefault("seed_applied", None)
                if artifact.n_sites != point.n_sites or dict(artifact.composition) != {
                    plan.system.element_a: point.element_a_count,
                    plan.system.element_b: point.element_b_count,
                }:
                    raise SqsOrchestrationError(f"Generator returned invalid point {point.label}")
                _write_structure(self.root, artifact, generated.structure, overwrite=overwrite)
                points.append(
                    SqsPointRecord(
                        spec=SqsCalculationSpec.point_spec(
                        system=plan.system,
                        point=point,
                        structure=artifact,
                        settings_fingerprint=settings_fingerprint,
                        supercell=plan.supercell,
                        metadata={
                            "template_checksum": plan.template.checksum,
                            "engine": generated.provenance.get("engine", plan.engine),
                        },
                        )
                    )
                )

            record = SqsSystemRecord(
                plan=plan,
                settings_fingerprint=settings_fingerprint,
                points=tuple(points),
                endpoint_cache_keys=tuple(endpoint_keys),
            )
            systems.append(record)
            _json_dump(system_root / "manifest.json", record.as_dict())

        manifest = SqsBatchManifest(
            tuple(systems),
            tuple(endpoints.values()),
            calculation=dict(calculation or {}),
        )
        manifest.assert_valid(self.root)
        manifest.write_json(self.manifest_path)
        return manifest

    def load(self) -> SqsBatchManifest:
        return SqsBatchManifest.read_json(self.manifest_path)

    def _save(self, manifest: SqsBatchManifest) -> None:
        manifest.write_json(self.manifest_path)
        for system in manifest.systems:
            _json_dump(
                self.root / "systems" / system.system_id / "manifest.json",
                system.as_dict(),
            )

    def validate(self) -> list[str]:
        return self.load().validation_errors(self.root)

    def pending_endpoints(self) -> tuple[SqsCalculationSpec, ...]:
        manifest = self.load()
        manifest.assert_valid(self.root)
        return manifest.pending_endpoint_specs()

    def pending_points(self) -> tuple[SqsCalculationSpec, ...]:
        return self.load().point_specs(self.root)

    def mark_endpoint_success(
        self,
        cache_key: str,
        observation: Mapping[str, Any] | SqsStaticObservation,
    ) -> SqsBatchManifest:
        with file_lock(self.root / ".batch.lock"):
            manifest = self.load()
            updated = manifest.with_endpoint_success(cache_key, observation, root=self.root)
            updated.assert_valid(self.root)
            self._save(updated)
            return updated

    def mark_endpoint_failure(self, cache_key: str, error: str) -> SqsBatchManifest:
        with file_lock(self.root / ".batch.lock"):
            manifest = self.load()
            updated = manifest.with_endpoint_failure(cache_key, error, root=self.root)
            self._save(updated)
            return updated

    def mark_point_success(
        self,
        calculation_id: str,
        observation: Mapping[str, Any] | SqsStaticObservation,
    ) -> SqsBatchManifest:
        with file_lock(self.root / ".batch.lock"):
            manifest = self.load()
            manifest.assert_points_ready(self.root)
            updated = manifest.with_point_success(calculation_id, observation, root=self.root)
            updated.assert_valid(self.root)
            self._save(updated)
            return updated

    def mark_point_failure(self, calculation_id: str, error: str) -> SqsBatchManifest:
        with file_lock(self.root / ".batch.lock"):
            manifest = self.load()
            manifest.assert_points_ready(self.root)
            updated = manifest.with_point_failure(calculation_id, error, root=self.root)
            self._save(updated)
            return updated

    def status(self) -> dict[str, int]:
        return self.load().status()

    def execute_pending_endpoints(
        self,
        executor: Callable[[SqsCalculationSpec], Mapping[str, Any] | SqsStaticObservation],
        *,
        continue_on_error: bool = True,
    ) -> SqsBatchManifest:
        """Run pending endpoint callbacks and persist each result immediately."""

        for record in self.load().endpoints:
            if record.state != "pending":
                continue
            try:
                observation = executor(record.spec)
                self.mark_endpoint_success(record.cache_key, observation)
            except Exception as exc:
                self.mark_endpoint_failure(record.cache_key, str(exc))
                if not continue_on_error:
                    raise
        return self.load()

    def execute_pending_points(
        self,
        executor: Callable[[SqsCalculationSpec], Mapping[str, Any] | SqsStaticObservation],
        *,
        continue_on_error: bool = True,
    ) -> SqsBatchManifest:
        """Run pending point callbacks after endpoint readiness is verified."""

        for system in self.load().systems:
            for point in system.points:
                if point.state != "pending":
                    continue
                try:
                    observation = executor(point.spec)
                    self.mark_point_success(point.calculation_id, observation)
                except Exception as exc:
                    self.mark_point_failure(point.calculation_id, str(exc))
                    if not continue_on_error:
                        raise
        return self.load()

    def collect(
        self,
        system_id: str,
        *,
        dataset_id: str | None = None,
        output_path: Path | str | None = None,
    ):
        """Build and persist one validated ``htvasp-sqs-mixing/v1`` dataset."""

        from htvasp.model.sqs_schema import (
            CalculationProvenance,
            SqsMixingData,
            build_sqs_mixing_data,
        )

        manifest = self.load()
        manifest.assert_valid(self.root)
        system = next((item for item in manifest.systems if item.system_id == system_id), None)
        if system is None:
            raise SqsOrchestrationError(f"Unknown SQS system: {system_id}")
        if any(point.state != "successful" for point in system.points):
            raise SqsOrchestrationError(f"Cannot collect {system_id} before all points succeed")

        endpoint_map = {record.cache_key: record for record in manifest.endpoints}

        def read_observation(path: str) -> SqsStaticObservation:
            observation = SqsStaticObservation.model_validate(
                json.loads(_safe_path(self.root, path).read_text(encoding="utf-8"))
            )
            if observation.batch_fingerprint != manifest.batch_fingerprint:
                raise SqsOrchestrationError(
                    "Cannot collect an observation from a different or legacy SQS batch"
                )
            return observation

        references: dict[str, SqsStaticObservation] = {}
        for key in system.endpoint_cache_keys:
            endpoint = endpoint_map[key]
            if endpoint.state != "successful" or endpoint.observation_path is None:
                raise SqsOrchestrationError(f"Endpoint {key} is not ready for collect")
            references[str(endpoint.spec.component)] = read_observation(endpoint.observation_path)

        point_observations: dict[str, SqsStaticObservation] = {}
        point_structures: dict[str, Mapping[str, Any]] = {}
        point_provenance: dict[str, Mapping[str, Any]] = {}
        for point in system.points:
            if point.observation_path is None:
                raise SqsOrchestrationError(f"Point {point.calculation_id} has no observation")
            point_observations[point.calculation_id] = read_observation(point.observation_path)
            structure_path = _artifact_path(self.root, point.spec.structure)
            point_structures[point.calculation_id] = Structure.from_file(structure_path).as_dict()
            point_provenance[point.calculation_id] = point.spec.structure.provenance

        calculation = dict(manifest.calculation)
        calculation.setdefault("settings_fingerprint", system.settings_fingerprint)
        calculation["batch_fingerprint"] = manifest.batch_fingerprint
        calculation.setdefault(
            "workflow", {"worker": "SqsMixingWorker", "flow": "relax_then_static"}
        )
        calculation.setdefault("relax_incar", {})
        calculation.setdefault("static_incar", {})
        data: SqsMixingData = build_sqs_mixing_data(
            system.plan,
            reference_observations=references,
            point_observations=point_observations,
            point_structures=point_structures,
            point_provenance=point_provenance,
            calculation=CalculationProvenance.model_validate(calculation),
            dataset_id=dataset_id or system_id,
        )
        destination = (
            Path(output_path)
            if output_path is not None
            else self.root.parent / f"{system_id}-sqs.json"
        )
        destination = destination if destination.is_absolute() else self.root / destination
        _json_dump(destination, data.model_dump(mode="json"))
        return data


def _structure_checksum(structure: Structure) -> str:
    from htvasp.model.templates import structure_checksum

    return structure_checksum(structure)


__all__ = [
    "BATCH_KIND",
    "BATCH_SCHEMA_VERSION",
    "SqsBatchManifest",
    "SqsEndpointRecord",
    "SqsOrchestrationError",
    "SqsOrchestrator",
    "SqsSystemRecord",
    "endpoint_cache_key",
]
