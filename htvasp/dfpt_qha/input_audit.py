"""Original-versus-actual VASP input audit.

custodian writes ``INCAR.orig``/``KPOINTS.orig``/``POSCAR.orig``/``POTCAR.orig``
before it starts correcting, then may rewrite the live files.  Those
corrections are allowed (see ``htvasp.utils.vasp_handlers``), but a DFPT
force-constant matrix is only interpretable if the inputs that produced it are
known.  This module compares the two sets and reports what moved.

It is an audit, not a gate.  A modified input does not fail the run: the
calculation happened, and the honest record is "this is what it actually used".
Only ``numerical_unresolved`` from the stability diagnosis stops a volume.

``POSCAR`` is part of the audit because ``VaspErrorHandler`` can reorder or
nudge it.  A force-constant matrix is indexed by that atom order, so a POSCAR
correction is the one that would silently mis-index everything downstream.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
from pymatgen.io.vasp import Incar, Kpoints, Poscar

log = logging.getLogger(__name__)

SCHEMA_VERSION = "htvasp-input-audit/v1"

#: Files compared structurally, original against the version actually run.
AUDITED_FILES = ("INCAR", "KPOINTS", "POSCAR")
#: Files recorded by presence and hash only: a POTCAR diff is not informative
#: without the pseudopotential library, and it is large.
HASHED_FILES = ("POTCAR",)

LATTICE_TOL = 1e-6
COORD_TOL = 1e-6
VALUE_TOL = 1e-9


def _locate(directory: Path, name: str) -> Path | None:
    """``NAME`` or ``NAME.gz``, whichever exists."""

    for candidate in (directory / name, directory / f"{name}.gz"):
        if candidate.is_file():
            return candidate
    return None


def _read_text(path: Path) -> str:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", errors="replace") as handle:
            return handle.read()
    return path.read_text(errors="replace")


def sha256_file(path: Path) -> str:
    """Hash of the decompressed content, so gzip does not change the identity."""

    if path.suffix == ".gz":
        with gzip.open(path, "rb") as handle:
            data = handle.read()
    else:
        data = path.read_bytes()
    return hashlib.sha256(data).hexdigest()


def _normalise(value: Any) -> Any:
    """A comparable form: numbers as floats, lists as tuples, the rest stripped."""

    if isinstance(value, (list, tuple)):
        return tuple(_normalise(v) for v in value)
    if isinstance(value, bool):
        return value
    try:
        return float(value)
    except (TypeError, ValueError):
        return str(value).strip()


def _same(left: Any, right: Any) -> bool:
    if isinstance(left, float) and isinstance(right, float):
        if np.isnan(left) and np.isnan(right):
            return True
        return abs(left - right) <= VALUE_TOL
    if isinstance(left, tuple) and isinstance(right, tuple):
        return len(left) == len(right) and all(
            _same(a, b) for a, b in zip(left, right)
        )
    return left == right


def _diff_mappings(original: dict[str, Any], current: dict[str, Any]) -> dict[str, Any]:
    changes: dict[str, Any] = {}
    for key in sorted(set(original) | set(current)):
        left = _normalise(original.get(key))
        right = _normalise(current.get(key))
        if not _same(left, right):
            changes[key] = {"original": original.get(key), "current": current.get(key)}
    return changes


def _incar_mapping(text: str) -> dict[str, Any]:
    return dict(Incar.from_str(text))


def _kpoints_mapping(text: str) -> dict[str, Any]:
    data = Kpoints.from_str(text).as_dict()
    return {
        key: value
        for key, value in data.items()
        if not key.startswith("@") and key != "comment"
    }


def _poscar_mapping(text: str) -> dict[str, Any]:
    structure = Poscar.from_str(text).structure
    return {
        "lattice": np.array(structure.lattice.matrix, dtype=float),
        "species": [str(site.specie) for site in structure],
        "frac_coords": np.array([site.frac_coords for site in structure], dtype=float),
    }


def _poscar_changes(original: dict[str, Any], current: dict[str, Any]) -> dict[str, Any]:
    """Lattice, species order and coordinates, each with its own tolerance.

    The species ORDER is compared as a sequence, not a multiset: an atom swap
    leaves the structure identical as a crystal and every force constant
    attached to the wrong atom.
    """

    changes: dict[str, Any] = {}
    if original["species"] != current["species"]:
        changes["species_order"] = {
            "original": original["species"],
            "current": current["species"],
        }
    if not np.allclose(original["lattice"], current["lattice"], atol=LATTICE_TOL):
        changes["lattice_max_abs_diff"] = float(
            np.abs(original["lattice"] - current["lattice"]).max()
        )
    if original["frac_coords"].shape != current["frac_coords"].shape:
        changes["frac_coords_shape"] = {
            "original": list(original["frac_coords"].shape),
            "current": list(current["frac_coords"].shape),
        }
    else:
        # wrap into (-1/2, 1/2]: 0.0 and 1.0 are the same site, and a custodian
        # rewrite that only re-wraps the coordinates must not read as a
        # full-cell displacement
        diff = original["frac_coords"] - current["frac_coords"]
        wrapped = (diff + 0.5) % 1.0 - 0.5
        delta = np.abs(wrapped).max()
        if delta > COORD_TOL:
            changes["frac_coords_max_abs_diff"] = float(delta)
    return changes


def audit_vasp_inputs(directory: Path | str) -> dict[str, Any]:
    """Compare custodian's ``.orig`` backups with the inputs that actually ran.

    Missing ``.orig`` files are reported, not treated as failures: a run whose
    job directory was trimmed still produced a usable matrix, it just cannot
    prove the inputs were untouched.  ``original_input_available`` says which
    case this is.
    """

    directory = Path(directory)
    record: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "directory": str(directory),
        "original_files": {},
        "current_files": {},
        "missing_original_files": [],
        "missing_current_files": [],
        "modified_files": [],
        "changes": {},
        "hashes": {},
        "warnings": [],
    }

    parsers = {
        "INCAR": (_incar_mapping, _diff_mappings),
        "KPOINTS": (_kpoints_mapping, _diff_mappings),
        "POSCAR": (_poscar_mapping, _poscar_changes),
    }

    for name in AUDITED_FILES:
        original_path = _locate(directory, f"{name}.orig")
        current_path = _locate(directory, name)
        if original_path is None:
            record["missing_original_files"].append(name)
            continue
        record["original_files"][name] = original_path.name
        if current_path is None:
            record["missing_current_files"].append(name)
            continue
        record["current_files"][name] = current_path.name

        parse, compare = parsers[name]
        try:
            record["hashes"].setdefault(name, {})["original"] = sha256_file(
                original_path
            )
            record["hashes"].setdefault(name, {})["current"] = sha256_file(
                current_path
            )
            changes = compare(
                parse(_read_text(original_path)), parse(_read_text(current_path))
            )
        except Exception as exc:  # noqa: BLE001 - an unreadable input is a finding
            # A truncated .gz (an interrupted atomate2 gzip_output_folder, a
            # partial copy) must be recorded, not allowed to kill the collector
            # three jobs downstream of the cause.
            changes = {"error": f"{type(exc).__name__}: {exc}"}
        if changes:
            record["modified_files"].append(name)
            record["changes"][name] = changes

    for name in HASHED_FILES:
        original_path = _locate(directory, f"{name}.orig")
        current_path = _locate(directory, name)
        if original_path is None:
            record["missing_original_files"].append(name)
            continue
        record["original_files"][name] = original_path.name
        if current_path is None:
            record["missing_current_files"].append(name)
            continue
        record["current_files"][name] = current_path.name
        try:
            original_hash = sha256_file(original_path)
            current_hash = sha256_file(current_path)
        except Exception as exc:  # noqa: BLE001 - an unreadable input is a finding
            record["changes"][name] = {"error": f"{type(exc).__name__}: {exc}"}
            record["modified_files"].append(name)
            continue
        record["hashes"].setdefault(name, {})["original"] = original_hash
        record["hashes"].setdefault(name, {})["current"] = current_hash
        if original_hash != current_hash:
            record["modified_files"].append(name)

    record["original_input_available"] = not record["missing_original_files"]
    for name in record["modified_files"]:
        warning = f"custodian modified {name}"
        record["warnings"].append(warning)
        log.warning("%s in %s", warning, directory)
    return json.loads(json.dumps(record, default=str))


__all__ = [
    "AUDITED_FILES",
    "HASHED_FILES",
    "SCHEMA_VERSION",
    "audit_vasp_inputs",
    "sha256_file",
]
