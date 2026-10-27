"""Local structure-template snapshots used by SQS planning.

This module deliberately has no Materials Project or other network client.
Templates are inputs to a plan and must be materialized locally before a
calculation is prepared.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from pymatgen.core import Structure


class TemplateError(ValueError):
    """Raised when a local template snapshot is missing or inconsistent."""


def _canonical_json(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        default=str,
    ).encode("utf-8")


def structure_checksum(structure: Structure) -> str:
    """Return a stable SHA-256 checksum for a pymatgen structure."""

    if not isinstance(structure, Structure):
        raise TypeError("structure_checksum expects a pymatgen Structure")
    return hashlib.sha256(_canonical_json(structure.as_dict())).hexdigest()


def file_checksum(path: Path | str) -> str:
    """Return a SHA-256 checksum for a local file."""

    source = Path(path)
    if not source.is_file():
        raise TemplateError(f"Template file does not exist: {source}")
    digest = hashlib.sha256()
    with source.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_checksum(value: str | None, field_name: str) -> None:
    if value is None:
        return
    if len(value) != 64 or any(character not in "0123456789abcdef" for character in value.lower()):
        raise TemplateError(f"{field_name} must be a SHA-256 hex digest")


@dataclass(frozen=True)
class TemplateRecord:
    """Metadata for one locally stored template structure."""

    phase: str
    path: str
    source_species: tuple[str, ...]
    source_material_id: str | None = None
    structure_checksum: str | None = None
    file_checksum: str | None = None

    def __post_init__(self) -> None:
        phase = str(self.phase).strip().upper()
        if not phase:
            raise TemplateError("Template phase must be non-empty")
        path = str(self.path).strip()
        if not path:
            raise TemplateError("Template path must be non-empty")
        species = tuple(str(item).strip() for item in self.source_species)
        if not species or any(not item for item in species):
            raise TemplateError("Template source_species must be non-empty")
        _validate_checksum(self.structure_checksum, "structure_checksum")
        _validate_checksum(self.file_checksum, "file_checksum")
        object.__setattr__(self, "phase", phase)
        object.__setattr__(self, "path", path)
        object.__setattr__(self, "source_species", species)

    def as_dict(self) -> dict[str, Any]:
        return {
            "phase": self.phase,
            "path": self.path,
            "source_species": list(self.source_species),
            "source_material_id": self.source_material_id,
            "structure_checksum": self.structure_checksum,
            "file_checksum": self.file_checksum,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "TemplateRecord":
        return cls(
            phase=str(data["phase"]),
            path=str(data["path"]),
            source_species=tuple(data["source_species"]),
            source_material_id=data.get("source_material_id"),
            structure_checksum=data.get("structure_checksum"),
            file_checksum=data.get("file_checksum"),
        )


@dataclass(frozen=True)
class TemplateSnapshot:
    """An immutable structure plus the checksums recorded for its source."""

    record: TemplateRecord
    structure: Structure

    def __post_init__(self) -> None:
        if not isinstance(self.structure, Structure):
            raise TypeError("TemplateSnapshot.structure must be a pymatgen Structure")
        expected = self.record.structure_checksum
        if expected is not None and expected != structure_checksum(self.structure):
            raise TemplateError(
                f"Template structure checksum mismatch for phase {self.record.phase}"
            )

    @property
    def phase(self) -> str:
        return self.record.phase

    @property
    def checksum(self) -> str:
        return structure_checksum(self.structure)

    def as_dict(self) -> dict[str, Any]:
        return {
            "record": {
                **self.record.as_dict(),
                "structure_checksum": self.checksum,
            },
            "structure": self.structure.as_dict(),
        }


class TemplateRegistry:
    """Resolve and validate templates under an injected local root directory.

    ``records`` may be supplied by a caller or loaded from a JSON registry;
    no default record performs network discovery.  The default root is only a
    conventional local path and can be replaced in tests or production setup.
    """

    DEFAULT_ROOT = Path(__file__).resolve().parents[2] / "data" / "sqs_templates"
    DEFAULT_RECORDS: dict[str, TemplateRecord] = {}

    def __init__(
        self,
        root: Path | str | None = None,
        records: Mapping[str, TemplateRecord | Mapping[str, Any]] | None = None,
    ) -> None:
        self.root = Path(root) if root is not None else self.DEFAULT_ROOT
        raw_records = records if records is not None else self.DEFAULT_RECORDS
        self._records = {
            str(key).strip().upper(): (
                value if isinstance(value, TemplateRecord) else TemplateRecord.from_dict(value)
            )
            for key, value in raw_records.items()
        }

    @classmethod
    def from_json(cls, path: Path | str, *, root: Path | str | None = None) -> "TemplateRegistry":
        source = Path(path)
        try:
            payload = json.loads(source.read_text(encoding="utf-8"))
        except FileNotFoundError as exc:
            raise TemplateError(f"Template registry does not exist: {source}") from exc
        records = payload.get("records", payload)
        if not isinstance(records, Mapping):
            raise TemplateError("Template registry JSON must contain a records mapping")
        return cls(root=root, records=records)

    def get_record(self, phase: str) -> TemplateRecord:
        key = str(phase).strip().upper()
        try:
            return self._records[key]
        except KeyError as exc:
            raise TemplateError(f"No local template record registered for phase {key!r}") from exc

    def get(self, phase: str, *, validate: bool = True) -> Structure:
        return self.snapshot(phase, validate=validate).structure.copy()

    def snapshot(self, phase: str, *, validate: bool = True) -> TemplateSnapshot:
        record = self.get_record(phase)
        path = Path(record.path)
        if not path.is_absolute():
            path = self.root / path
        if not path.is_file():
            raise TemplateError(f"Template file does not exist: {path}")
        if validate and record.file_checksum is not None:
            actual_file_checksum = file_checksum(path)
            if actual_file_checksum != record.file_checksum:
                raise TemplateError(f"Template file checksum mismatch: {path}")
        try:
            structure = Structure.from_file(path)
        except Exception as exc:
            raise TemplateError(f"Cannot read template structure: {path}") from exc
        snapshot = TemplateSnapshot(record, structure)
        if validate and record.structure_checksum is not None and snapshot.checksum != record.structure_checksum:
            raise TemplateError(f"Template structure checksum mismatch: {path}")
        return snapshot

    def as_dict(self) -> dict[str, Any]:
        return {"root": str(self.root), "records": {key: value.as_dict() for key, value in self._records.items()}}


__all__ = [
    "TemplateError",
    "TemplateRecord",
    "TemplateRegistry",
    "TemplateSnapshot",
    "file_checksum",
    "structure_checksum",
]
