"""Tier-A stability checks and tier-B finite-q soft-mode diagnostics.

Pure analysis only: no VASP, no scheduler and no scripts imports.  The input is
a production ``Phonopy`` object with its DFPT force constants attached.  The
module deliberately distinguishes:

* numerical integrity (ordering, ASR, reciprocal convention), and
* physical stability (exact supercell-commensurate q points).

A Monkhorst-Pack mesh is diagnostic context only: its minimum is an interpolated
sample and is not allowed to stand in for the instability wavevector.
"""

from __future__ import annotations

import hashlib
import json
from typing import Any

import numpy as np
from phonopy.harmonic.dynmat_to_fc import get_commensurate_points

from htvasp.dfpt_qha.phonons import ordering_fingerprint

SOFT_FREQUENCY_TOL_THz = 1e-4
GAMMA_ACOUSTIC_TOL_THz = 0.1
EIGEN_RESIDUAL_TOL = 1e-6
PERMUTATION_DIAGNOSTIC_TOL = 1e-2


def _canonical_json(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()


def recipe_sha256(recipe: Any) -> str:
    payload = recipe.model_dump(mode="json") if hasattr(recipe, "model_dump") else recipe
    return hashlib.sha256(_canonical_json(payload)).hexdigest()


def force_constants_sha256(force_constants: Any) -> str:
    array = np.ascontiguousarray(force_constants, dtype="<f8")
    shape = np.asarray(array.shape, dtype="<i8").tobytes()
    return hashlib.sha256(shape + array.tobytes()).hexdigest()


def force_constants_checks(force_constants: Any) -> dict[str, Any]:
    fc = np.asarray(force_constants, dtype=float)
    problems: list[str] = []
    shape_ok = fc.ndim == 4 and fc.shape[0] == fc.shape[1] and tuple(fc.shape[2:]) == (3, 3)
    if not shape_ok:
        problems.append(f"unexpected force_constants shape {fc.shape}")
        return {"shape_ok": False, "problems": problems}
    scale = float(np.abs(fc).max()) or 1.0
    asr_abs = float(np.abs(fc.sum(axis=1)).max())
    permutation_abs = float(np.abs(fc - fc.transpose(1, 0, 3, 2)).max())
    asr_rel = asr_abs / scale
    permutation_rel = permutation_abs / scale
    return {
        "shape_ok": True,
        "shape": list(fc.shape),
        "max_force_constant": scale,
        "asr_abs": asr_abs,
        "asr_rel": asr_rel,
        "asr_ok": asr_rel <= 1e-2,
        "asr_report_threshold_ok": not (asr_rel > 0.10 and asr_abs > 0.01),
        "permutation_abs": permutation_abs,
        "permutation_rel": permutation_rel,
        "permutation_ok": permutation_rel <= PERMUTATION_DIAGNOSTIC_TOL,
        "force_constants_sha256": force_constants_sha256(fc),
        "problems": problems,
    }


def reciprocal_lattice(ph: Any) -> np.ndarray:
    """Return B = 2*pi*A^-T, rows are reciprocal Cartesian vectors (1/A)."""

    lattice = np.asarray(ph.primitive.cell, dtype=float)
    return 2.0 * np.pi * np.linalg.inv(lattice).T


def reciprocal_check(ph: Any, tol: float = 1e-10) -> dict[str, Any]:
    lattice = np.asarray(ph.primitive.cell, dtype=float)
    reciprocal = reciprocal_lattice(ph)
    deviation = float(np.abs(reciprocal @ lattice.T - 2.0 * np.pi * np.eye(3)).max())
    return {"ok": deviation <= tol, "max_deviation": deviation, "tol": tol}


def commensurate_qpoints(supercell_matrix: Any) -> list[tuple[float, float, float]]:
    points = get_commensurate_points(np.asarray(supercell_matrix, dtype=int))
    return sorted({tuple(float(x % 1.0) for x in q) for q in np.asarray(points)})


def _q_distance(a: np.ndarray, b: np.ndarray) -> float:
    delta = (a - b + 0.5) % 1.0 - 0.5
    return float(np.linalg.norm(delta))


def qpoint_report(ph: Any, qpoint: tuple[float, float, float], tol_thz: float = SOFT_FREQUENCY_TOL_THz) -> dict[str, Any]:
    ph.run_qpoints([list(qpoint)], with_eigenvectors=True)
    data = ph.get_qpoints_dict()
    frequencies = np.asarray(data["frequencies"])[0]
    vectors = np.asarray(data["eigenvectors"])[0]
    ph.dynamical_matrix.run(np.asarray(qpoint, dtype=float))
    dynamical = np.asarray(ph.dynamical_matrix.dynamical_matrix)
    factor = 15.633302
    bands = []
    for band, frequency in enumerate(frequencies):
        vector = vectors[:, band]
        eigenvalue = float(np.sign(frequency) * (frequency / factor) ** 2)
        column_residual = float(np.linalg.norm(dynamical @ vector - eigenvalue * vector))
        row = vectors[band, :]
        row_residual = float(np.linalg.norm(dynamical @ row - eigenvalue * row))
        bands.append({
            "band": int(band), "frequency_THz": float(frequency),
            "eigenvector": [{"re": float(z.real), "im": float(z.imag)} for z in vector],
            "column_residual": column_residual, "row_residual": row_residual,
        })
    groups: list[list[int]] = []
    for band, frequency in enumerate(frequencies):
        group = [int(i) for i, other in enumerate(frequencies) if abs(other - frequency) < 1e-3]
        if group not in groups: groups.append(group)
    soft = [int(i) for i, f in enumerate(frequencies) if f < -tol_thz]
    return {
        "q_fractional": [float(x) for x in qpoint],
        "q_cartesian_per_A": [float(x) for x in np.asarray(qpoint) @ reciprocal_lattice(ph)],
        "frequencies_THz": [float(x) for x in frequencies],
        "min_frequency_THz": float(frequencies.min()),
        "n_imaginary": len(soft), "soft_bands": soft,
        "degeneracy_groups": groups, "bands": bands,
        "eigenvector_basis": "primitive cell, 3*n_primitive components",
        "eigenvector_convention": "mass-weighted; column E[:, band]",
    }


def gamma_acoustic_report(ph: Any) -> dict[str, Any]:
    report = qpoint_report(ph, (0.0, 0.0, 0.0))
    acoustic = report["frequencies_THz"][:3]
    return {"frequencies_THz": acoustic, "max_abs_THz": max(abs(x) for x in acoustic),
            "noise_ok": max(abs(x) for x in acoustic) < GAMMA_ACOUSTIC_TOL_THz}


def mesh_diagnostic(ph: Any, mesh: tuple[int, int, int]) -> dict[str, Any]:
    ph.run_mesh(list(mesh), with_eigenvectors=False)
    data = ph.get_mesh_dict()
    frequencies = np.asarray(data["frequencies"])
    qpoints = np.asarray(data["qpoints"])
    index = np.unravel_index(int(np.argmin(frequencies)), frequencies.shape)
    q = qpoints[index[0]]
    return {
        "mesh": list(mesh),
        "grid_contains_gamma": bool(np.any(np.all(np.isclose(qpoints, 0.0), axis=1))),
        "grid_contains_half": bool(np.any(np.isclose(qpoints, 0.5))),
        "mesh_min_qpoint": [float(x) for x in q],
        "mesh_min_frequency_THz": float(frequencies[index]),
        "n_imaginary_modes": int((frequencies < -SOFT_FREQUENCY_TOL_THz).sum()),
        "n_modes": int(frequencies.size),
        "note": "mesh minimum is interpolated context, not the instability wavevector",
    }


def diagnose_stability(ph: Any, force_constants: Any, *, structure: Any, recipe: Any,
                       supercell_matrix: Any, mesh: tuple[int, int, int]) -> dict[str, Any]:
    checks = force_constants_checks(force_constants)
    reciprocal = reciprocal_check(ph)
    gamma = gamma_acoustic_report(ph)
    qpoints = commensurate_qpoints(supercell_matrix)
    q_reports = [qpoint_report(ph, q) for q in qpoints]
    finite_q = [r for r in q_reports if r["q_fractional"] != [0.0, 0.0, 0.0] and r["n_imaginary"]]
    gamma_report = next(r for r in q_reports if r["q_fractional"] == [0.0, 0.0, 0.0])
    structural_ok = checks.get("shape_ok") and checks.get("asr_ok") and checks.get("ordering_ok", True) and reciprocal["ok"]
    if not structural_ok:
        status, allowed, reason = "numerical_unresolved", False, "tier-A structural check failed"
    elif finite_q:
        status, allowed, reason = "finite_q_soft_mode_diagnostic", False, "commensurate finite-q imaginary branch"
    elif not gamma["noise_ok"] and gamma_report["n_imaginary"]:
        status, allowed, reason = "gamma_soft_mode", False, "Gamma imaginary branch exceeds acoustic noise"
    else:
        status, allowed, reason = "harmonic_QHA_applicable", True, "no commensurate imaginary branch"
    return {
        "schema_version": "htvasp-dfpt-stability/v1",
        "status": status, "qha_allowed": allowed, "reason": reason,
        "checks": {**checks, "reciprocal": reciprocal, "gamma": gamma},
        "commensurate_qpoints": q_reports,
        "mesh_diagnostic": mesh_diagnostic(ph, mesh),
        "force_constants_sha256": force_constants_sha256(force_constants),
        "structure_hash": ordering_fingerprint(structure),
        "recipe_hash": recipe_sha256(recipe),
        "units": {"frequency": "THz", "force_constants": "eV/Angstrom^2"},
    }


__all__ = ["commensurate_qpoints", "diagnose_stability", "force_constants_checks",
           "force_constants_sha256", "gamma_acoustic_report", "mesh_diagnostic",
           "qpoint_report", "reciprocal_check", "reciprocal_lattice", "recipe_sha256"]
