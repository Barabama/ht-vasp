"""Jobflow jobs for standardizing SQS calculation outputs."""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping
from numbers import Real
from typing import Any

from jobflow import job
from pymatgen.core import Structure

from htvasp.sqs.task_doc import SqsStaticObservation


def _jsonable(value: Any) -> Any:
    """Convert common task-document models into deterministic JSON data."""

    if isinstance(value, Structure):
        return value.as_dict()
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if hasattr(value, "model_dump"):
        try:
            dumped = value.model_dump(mode="python")
        except TypeError:
            dumped = value.model_dump()
        return _jsonable(dumped)
    if hasattr(value, "as_dict"):
        return _jsonable(value.as_dict())
    if hasattr(value, "tolist"):
        return _jsonable(value.tolist())
    return value


def _digest(value: Any) -> str:
    payload = json.dumps(
        _jsonable(value), sort_keys=True, separators=(",", ":"), ensure_ascii=True, default=str
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _mapping(value: Any) -> Mapping[str, Any] | None:
    normalized = _jsonable(value)
    return normalized if isinstance(normalized, Mapping) else None


def _find_value(value: Any, keys: tuple[str, ...]) -> Any:
    mapping = _mapping(value)
    if mapping is None:
        return None
    for key in keys:
        if key in mapping:
            return mapping[key]
    for nested_key in ("output", "result", "taskdoc", "task_doc", "task_document"):
        if nested_key in mapping:
            found = _find_value(mapping[nested_key], keys)
            if found is not None:
                return found
    for nested_key in ("calcs_reversed",):
        nested = mapping.get(nested_key)
        if isinstance(nested, list):
            for item in reversed(nested):
                found = _find_value(item, keys)
                if found is not None:
                    return found
    return None


def _structure_from_task_doc(value: Any) -> Structure:
    candidate = _find_value(value, ("structure", "final_structure"))
    if candidate is None:
        raise ValueError("Static task document has no final structure")
    if isinstance(candidate, Structure):
        return candidate
    try:
        return Structure.from_dict(candidate)
    except Exception as exc:  # pragma: no cover - error text is the contract
        raise ValueError("Static task document contains an invalid final structure") from exc


def _energy(value: Any, n_sites: int) -> float:
    total = _find_value(value, ("energy", "final_energy", "e_0_energy", "energy_total_eV"))
    if isinstance(total, Real) and not isinstance(total, bool):
        return float(total)
    per_atom = _find_value(value, ("energy_per_atom", "final_energy_per_atom", "energy_per_atom_eV"))
    if isinstance(per_atom, Real) and not isinstance(per_atom, bool):
        return float(per_atom) * n_sites
    raise ValueError("Static task document has no total or per-atom energy")


@job
def make_static_observation(
    static_output: Any,
    calculation_id: str,
    batch_fingerprint: str | None = None,
    intent_key: str | None = None,
) -> dict[str, Any]:
    """Convert one successful Static task document into a portable observation."""

    state = _find_value(static_output, ("state",))
    state = getattr(state, "value", state)
    if state != "successful":
        raise ValueError(f"Static task did not succeed: state={state!r}")

    structure = _structure_from_task_doc(static_output)
    n_sites = len(structure)
    total_energy = _energy(static_output, n_sites)
    if not math.isfinite(total_energy):
        raise ValueError("Static task energy must be finite")

    observation = SqsStaticObservation(
        calculation_id=calculation_id,
        batch_fingerprint=batch_fingerprint,
        intent_key=intent_key,
        state="successful",
        energy_total_eV=total_energy,
        energy_per_atom_eV=total_energy / n_sites,
        n_sites=n_sites,
        final_structure=structure.as_dict(),
        task_document_sha256=_digest(static_output),
    )
    return observation.model_dump(mode="json")


def normalize_static_task_document(
    static_output: Any,
    *,
    calculation_id: str,
    batch_fingerprint: str | None = None,
    intent_key: str | None = None,
) -> dict[str, Any]:
    """Normalize a Static task document without constructing a Job.

    This helper keeps the validation logic directly testable while
    ``make_static_observation`` remains the jobflow boundary.
    """

    return make_static_observation.original(
        static_output,
        calculation_id,
        batch_fingerprint,
        intent_key,
    )


# Newer callers can use the shorter name without breaking the historical
# ``make_static_observation`` import.
sqs_observation = make_static_observation


__all__ = [
    "make_static_observation",
    "normalize_static_task_document",
    "sqs_observation",
]
