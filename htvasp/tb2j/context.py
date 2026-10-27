"""Serializable producer-side semantic context for TB2J artifacts."""

from __future__ import annotations

import copy
import hashlib
import inspect
import json
from collections.abc import Callable, Mapping
from typing import Any

from pymatgen.core import Structure

CONTEXT_SOURCE = "htvasp_workflow_context"
CONTEXT_VERSION = "htvasp-magnetic-context/v1"
PRODUCER_RECIPE_VERSION = "main_w90-v2"
COORDINATE_TOLERANCE = 1.0e-6


def _rounded(value: float) -> float:
    value = round(float(value), 12)
    return 0.0 if abs(value) <= COORDINATE_TOLERANCE else value


def _fractional_coordinates(site: Any) -> list[float]:
    result = []
    for value in site.frac_coords:
        normalized = float(value) % 1.0
        if normalized <= COORDINATE_TOLERANCE or 1.0 - normalized <= COORDINATE_TOLERANCE:
            normalized = 0.0
        result.append(_rounded(normalized))
    return result


def _site_element(site: Any) -> str:
    if len(site.species) != 1:
        raise ValueError(
            "magnetic context requires ordered sites; disordered occupancies need explicit input"
        )
    return next(iter(site.species)).symbol


def _site_snapshot(structure: Structure) -> list[dict[str, Any]]:
    counts: dict[str, int] = {}
    result: list[dict[str, Any]] = []
    for index, site in enumerate(structure):
        element = _site_element(site)
        counts[element] = counts.get(element, 0) + 1
        result.append(
            {
                "id": f"{element}{counts[element]}",
                "index": index,
                "element": element,
                "frac_coords": _fractional_coordinates(site),
            }
        )
    return result


def _structure_payload(structure: Structure, sites: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "lattice_matrix_angstrom": [
            [_rounded(value) for value in row] for row in structure.lattice.matrix
        ],
        "sites": [
            {
                "id": site["id"],
                "element": site["element"],
                "frac_coords": list(site["frac_coords"]),
            }
            for site in sites
        ],
    }


def structure_checksum(structure: Structure) -> str:
    """Return a deterministic SHA-256 for lattice and ordered site data."""

    sites = _site_snapshot(structure)
    encoded = json.dumps(
        _structure_payload(structure, sites),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _near(coords: list[float], target: tuple[float, float, float]) -> bool:
    return all(abs(value - expected) <= COORDINATE_TOLERANCE for value, expected in zip(coords, target))


def _find_unique_site(
    sites: list[dict[str, Any]], predicate: Callable[[list[float]], bool], label: str
) -> dict[str, Any]:
    matches = [site for site in sites if predicate(site["frac_coords"])]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one {label} site, found {len(matches)}")
    return matches[0]


def _base_context(
    structure: Structure,
    structure_role: str,
    producer_recipe_version: str,
) -> dict[str, Any]:
    sites = _site_snapshot(structure)
    payload = _structure_payload(structure, sites)
    checksum = structure_checksum(structure)
    return {
        "source": CONTEXT_SOURCE,
        "context_version": CONTEXT_VERSION,
        "producer_recipe_version": producer_recipe_version,
        "structure_checksum": checksum,
        "structure_checksum_algorithm": "sha256(canonical lattice+site snapshot)",
        "structure_role": structure_role,
        "sites": copy.deepcopy(sites),
        "structure": {
            "lattice_matrix_angstrom": payload["lattice_matrix_angstrom"],
            "sites": copy.deepcopy(sites),
            "checksum": checksum,
        },
        "magnetic_moment_source": (
            "TB2J atom_moments from the solve result; initial VASP MAGMOM is not an atom moment"
        ),
        "producer": {
            "source": CONTEXT_SOURCE,
            "recipe_version": producer_recipe_version,
        },
    }


def _mapping_rationale(text: str, **details: Any) -> dict[str, Any]:
    return {
        "method": "explicit fractional-coordinate geometry",
        "tolerance": COORDINATE_TOLERANCE,
        "explanation": text,
        **details,
    }


def _attach_calphad(context: dict[str, Any], mapping: dict[str, Any]) -> dict[str, Any]:
    normalized = copy.deepcopy(mapping)
    # ``ratios`` is the concise workflow spelling; the canonical JSON contract
    # uses the explicit ``sublattice_ratios`` name.  Emit both at the producer
    # boundary so downstream consumers do not need producer-specific aliases.
    if "sublattice_ratios" not in normalized and "ratios" in normalized:
        normalized["sublattice_ratios"] = copy.deepcopy(normalized["ratios"])
    normalized.setdefault("structure_role", context.get("structure_role"))
    context["calphad"] = normalized
    context["calphad_mapping"] = copy.deepcopy(normalized)
    context["status"] = mapping["status"]
    context["mapping_status"] = mapping["status"]
    context["production_ready"] = bool(mapping.get("production_ready", False))
    for key in (
        "phase",
        "type_character",
        "ratios",
        "site_to_sublattice",
        "endmember_components",
        "order_num",
    ):
        context[key] = copy.deepcopy(normalized.get(key))
    context["mapping_rationale"] = copy.deepcopy(normalized.get("mapping_rationale"))
    return context


def _not_applicable_context(context: dict[str, Any]) -> dict[str, Any]:
    return _attach_calphad(
        context,
        {
            "status": "not_applicable",
            "production_ready": False,
            "phase": None,
            "type_character": None,
            "ratios": None,
            "site_to_sublattice": {},
            "endmember_components": None,
            "order_num": None,
            "mapping_rationale": _mapping_rationale(
                "A unary SER structure is a reference state, not an ordered binary CALPHAD endmember; no B2 mapping is fabricated."
            ),
        },
    )


def _bcc_context(context: dict[str, Any]) -> dict[str, Any]:
    sites = context["sites"]
    if len(sites) != 2:
        raise ValueError("ordered BCC endmember mapping requires exactly two sites")
    corner = _find_unique_site(
        sites, lambda coords: _near(coords, (0.0, 0.0, 0.0)), "BCC corner"
    )
    body = _find_unique_site(
        sites, lambda coords: _near(coords, (0.5, 0.5, 0.5)), "BCC body-center"
    )
    rationale = _mapping_rationale(
        "The site at fractional (0,0,0) is mapped to corner_1a and the site at (1/2,1/2,1/2) to body_center_1b; endmember order follows their actual occupations.",
        corner_site_id=corner["id"],
        body_center_site_id=body["id"],
        site_occupations={corner["id"]: corner["element"], body["id"]: body["element"]},
    )
    return _attach_calphad(
        context,
        {
            "status": "ready",
            "production_ready": True,
            "phase": "B2_BCC",
            "type_character": "&",
            "ratios": [0.5, 0.5],
            "site_to_sublattice": {corner["id"]: 0, body["id"]: 1},
            "sublattice_labels": ["corner_1a", "body_center_1b"],
            "endmember_components": [corner["element"].upper(), body["element"].upper()],
            "order_num": 0,
            "mapping_rationale": rationale,
        },
    )


def _fcc_context(context: dict[str, Any]) -> dict[str, Any]:
    sites = context["sites"]
    rationale = _mapping_rationale(
        "FCC/L12 production mapping requires one 1a corner site and three 3c face sites; it is not inferred from the directory token alone."
    )
    mapping: dict[str, Any] = {
        "status": "needs_review",
        "production_ready": False,
        "phase": "L12_FCC",
        "type_character": "&",
        "ratios": [0.25, 0.75],
        "site_to_sublattice": {},
        "sublattice_labels": ["corner_1a", "face_3c"],
        "endmember_components": None,
        "order_num": 0,
        "mapping_rationale": rationale,
    }
    if len(sites) != 4:
        rationale["explanation"] += " The structure does not contain the four sites required for an explicit 1a/3c mapping."
        return _attach_calphad(context, mapping)
    try:
        corner = _find_unique_site(
            sites, lambda coords: _near(coords, (0.0, 0.0, 0.0)), "FCC 1a corner"
        )
        face = [
            site
            for site in sites
            if sum(abs(value - 0.5) <= COORDINATE_TOLERANCE for value in site["frac_coords"]) == 2
            and sum(abs(value) <= COORDINATE_TOLERANCE for value in site["frac_coords"]) == 1
        ]
        if len(face) != 3:
            raise ValueError(f"expected three FCC 3c face sites, found {len(face)}")
    except ValueError as exc:
        rationale["explanation"] += f" Explicit 1a/3c detection failed: {exc}."
        return _attach_calphad(context, mapping)
    rationale.update(
        {
            "corner_1a_site_id": corner["id"],
            "face_3c_site_ids": [site["id"] for site in face],
            "site_occupations": {
                corner["id"]: corner["element"],
                **{site["id"]: site["element"] for site in face},
            },
        }
    )
    mapping.update(
        {
            "status": "ready",
            "production_ready": True,
            "site_to_sublattice": {corner["id"]: 0, **{site["id"]: 1 for site in face}},
            "endmember_components": [corner["element"].upper(), face[0]["element"].upper()],
        }
    )
    return _attach_calphad(context, mapping)


def build_magnetic_context(
    name: str,
    structure: Structure,
    *,
    producer_recipe_version: str = PRODUCER_RECIPE_VERSION,
) -> dict[str, Any]:
    """Build an explicit context; never infer a phase for an unknown role."""

    if not isinstance(structure, Structure):
        raise TypeError("structure must be a pymatgen Structure")
    prefix = str(name).split("-", 1)[0].upper()
    role = {
        "SER": "unary_reference",
        "BCC": "ordered_bcc_endmember",
        "FCC": "ordered_fcc_endmember",
    }.get(prefix, "unknown")
    context = _base_context(structure, role, producer_recipe_version)
    if role == "unary_reference":
        return _not_applicable_context(context)
    if role == "ordered_bcc_endmember":
        try:
            return _bcc_context(context)
        except ValueError as exc:
            return _attach_calphad(
                context,
                {
                    "status": "needs_review",
                    "production_ready": False,
                    "phase": "B2_BCC",
                    "type_character": "&",
                    "ratios": [0.5, 0.5],
                    "site_to_sublattice": {},
                    "endmember_components": None,
                    "order_num": 0,
                    "mapping_rationale": _mapping_rationale(
                        f"BCC role was requested, but explicit corner/body-center mapping failed: {exc}."
                    ),
                },
            )
    if role == "ordered_fcc_endmember":
        return _fcc_context(context)
    return _attach_calphad(
        context,
        {
            "status": "needs_input",
            "production_ready": False,
            "phase": None,
            "type_character": None,
            "ratios": None,
            "site_to_sublattice": {},
            "endmember_components": None,
            "order_num": None,
            "mapping_rationale": _mapping_rationale(
                "The workflow role is unknown; explicit user input is required and no CALPHAD phase is inferred."
            ),
        },
    )


def normalize_magnetic_context(context: Mapping[str, Any] | None) -> dict[str, Any] | None:
    """Return a JSON-safe detached copy of an optional context."""

    if context is None:
        return None
    if hasattr(context, "to_dict"):
        context = context.to_dict()
    if not isinstance(context, Mapping):
        raise TypeError("magnetic_context must be a mapping or None")
    return json.loads(json.dumps(dict(context), sort_keys=True, separators=(",", ":"), allow_nan=False))


def serialize_magnetic_context(context: Mapping[str, Any]) -> str:
    normalized = normalize_magnetic_context(context)
    if normalized is None:
        raise ValueError("cannot serialize a missing magnetic_context")
    return json.dumps(normalized, sort_keys=True, separators=(",", ":"), allow_nan=False)


def deserialize_magnetic_context(value: str | bytes) -> dict[str, Any]:
    normalized = normalize_magnetic_context(json.loads(value))
    if normalized is None:
        raise ValueError("serialized magnetic_context must contain an object")
    return normalized


def context_summary(context: Mapping[str, Any]) -> dict[str, Any]:
    normalized = normalize_magnetic_context(context)
    assert normalized is not None
    calphad = normalized.get("calphad") or {}
    return {
        "role": normalized.get("structure_role"),
        "status": calphad.get("status"),
        "phase": calphad.get("phase"),
        "structure_checksum": normalized.get("structure_checksum"),
        "site_count": len(normalized.get("sites") or []),
    }


def _accepts_keyword(function: Callable[..., Any], keyword: str) -> bool:
    try:
        parameters = inspect.signature(function).parameters.values()
    except (TypeError, ValueError):
        return True
    return any(parameter.kind == inspect.Parameter.VAR_KEYWORD for parameter in parameters) or any(
        parameter.name == keyword for parameter in parameters
    )


def call_tb2j_solve_with_context(
    job_function: Callable[..., Any],
    *args: Any,
    magnetic_context: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> Any:
    """Call the new job interface, with an old-job compatibility path."""

    normalized = normalize_magnetic_context(magnetic_context)
    if normalized is None or _accepts_keyword(job_function, "magnetic_context"):
        if normalized is not None:
            kwargs["magnetic_context"] = normalized
        return job_function(*args, **kwargs)
    job_object = job_function(*args, **kwargs)
    metadata = getattr(job_object, "metadata", None)
    if isinstance(metadata, dict):
        metadata["magnetic_context"] = normalized
    return job_object


def canonicalize_with_context(
    canonicalizer: Callable[..., Any],
    data: Mapping[str, Any],
    *,
    output_path: Any = None,
    magnetic_context: Mapping[str, Any] | None = None,
    **kwargs: Any,
) -> Any:
    """Call a canonicalizer with context while retaining old signatures."""

    normalized = normalize_magnetic_context(magnetic_context)
    call_kwargs = dict(kwargs)
    if output_path is not None:
        call_kwargs["output_path"] = output_path
    if normalized is not None and _accepts_keyword(canonicalizer, "magnetic_context"):
        call_kwargs["magnetic_context"] = normalized
        return canonicalizer(dict(data), **call_kwargs)
    if normalized is not None:
        payload = copy.deepcopy(dict(data))
        payload["magnetic_context"] = normalized
        return canonicalizer(payload, **call_kwargs)
    return canonicalizer(dict(data), **call_kwargs)


__all__ = [
    "CONTEXT_SOURCE",
    "CONTEXT_VERSION",
    "PRODUCER_RECIPE_VERSION",
    "build_magnetic_context",
    "call_tb2j_solve_with_context",
    "canonicalize_with_context",
    "context_summary",
    "deserialize_magnetic_context",
    "normalize_magnetic_context",
    "serialize_magnetic_context",
    "structure_checksum",
]
