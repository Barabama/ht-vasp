"""Tier-A stability checks and tier-B finite-q soft-mode diagnostics.

Pure analysis only: no VASP, no scheduler and no scripts imports.  The input is
a production ``Phonopy`` object with its DFPT force constants attached.  The
module deliberately distinguishes:

* numerical integrity (ordering, ASR, reciprocal convention), and
* physical stability (exact supercell-commensurate q points).

Stability is decided in three tiers, first match wins:

1. structural / numerical failure             -> ``numerical_unresolved``
2. an imaginary branch at Gamma (optical)     -> ``gamma_soft_mode``
3. an imaginary branch at an exact
   commensurate q point                       -> ``finite_q_soft_mode``
4. a deep imaginary branch on the dense mesh  -> ``incommensurate_soft_mode_suspected``
   that no exact q point represents
5. otherwise                                  -> ``harmonic_QHA_applicable``

A Monkhorst-Pack mesh is diagnostic context only.  Its minimum is an
interpolated sample, so a shallow negative value there (BCC-Nb-Nb-2: -0.024 THz
at q = (-1/24, 1/24, 1/24)) is a truncation artifact, while a deep one is a
real instability living where the supercell cannot sample it (bcc omega[111] at
q ~ 2/3).  Only tiers 1-4 can block QHA.

A 2x2x2 supercell represents q in {0, 1/2}^3 and nothing else, so passing this
gate does not prove the crystal is stable at incommensurate q.  The report says
so explicitly (``supercell_convergence``) rather than implying otherwise.
"""

from __future__ import annotations

import hashlib
import json
from typing import Any, Literal

import numpy as np
from phonopy.harmonic.dynmat_to_fc import DynmatToForceConstants
from pydantic import BaseModel, ConfigDict, Field, model_validator

from htvasp.dfpt_qha.phonons import ordering_fingerprint

SOFT_FREQUENCY_TOL_THz = 1e-4
#: Applies to the three ACOUSTIC branches at Gamma only.  It is numerical noise
#: tolerance, not physics: an optical branch near Gamma gets no such allowance.
GAMMA_ACOUSTIC_TOL_THz = 0.1
#: A dense-mesh imaginary branch deeper than this, at a q the supercell cannot
#: represent, is treated as a suspected incommensurate instability rather than
#: an interpolation artifact.  The gap it sits in is wide: the archived
#: BCC-Nb-Nb-2 mesh artifact is -0.024 THz, the BCC-Al-Ta-2 real N-point mode is
#: -3.5 THz.
MESH_SOFT_SUSPECT_TOL_THz = 0.5
#: Below this magnitude an exact-q soft mode is "borderline": a larger supercell
#: is needed to tell a real instability from a supercell-truncation effect.
BORDERLINE_CONFIRMATION_TOL_THz = 1.0
EIGEN_RESIDUAL_TOL = 1e-6
PERMUTATION_DIAGNOSTIC_TOL = 1e-2

SCHEMA_VERSION = "htvasp-dfpt-stability/v1"


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


def commensurate_qpoints(phonon: Any) -> list[tuple[float, float, float]]:
    """Exact q points the supercell can represent, in the primitive basis.

    Taken from phonopy's own ``DynmatToForceConstants`` rather than computed
    from the recipe's supercell matrix.  ``get_commensurate_points`` documents
    its argument as "with respect to primitive cell basis vectors", and the
    recipe matrix is expressed in the *unit* cell: for an fcc A=A endpoint the
    4-atom conventional cell has a 1-atom primitive, so the 2x2x2 recipe matrix
    is 8 primitives per axis pair -- 32 commensurate points, not 8.  Passing the
    recipe matrix would silently test a third of them.

    Deriving the set from the constructed phonopy object removes the conversion
    entirely: it is whatever that object's supercell actually is.
    """

    points = DynmatToForceConstants(
        phonon.primitive, phonon.supercell
    ).commensurate_points
    return sorted(
        {tuple(float(x % 1.0) for x in q) for q in np.asarray(points, dtype=float)}
    )


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


def _translational_overlaps(ph: Any, report: dict[str, Any]) -> list[float]:
    """How much of each Gamma eigenvector is a rigid translation of the cell.

    The acoustic branches are identified by their eigenvector, not by their
    position in the sorted frequency list.  phonopy returns frequencies in
    ascending order, so an imaginary *optical* branch sorts BELOW the three
    ~0 THz acoustic modes; slicing ``frequencies[:3]`` would then hand the
    optical instability the acoustic noise allowance and report the volume as
    harmonic-QHA-applicable.

    In mass-weighted coordinates a rigid translation along alpha is
    ``sqrt(m_i) e_alpha`` for every atom i, so the projection onto those three
    vectors is 1 for an acoustic mode and 0 for an optical one.
    """

    masses = np.asarray(ph.primitive.masses, dtype=float)
    n_atoms = len(masses)
    basis = []
    for alpha in range(3):
        vector = np.zeros((n_atoms, 3), dtype=float)
        vector[:, alpha] = np.sqrt(masses)
        basis.append(vector.ravel() / np.linalg.norm(vector))

    overlaps: list[float] = []
    for band in report["bands"]:
        vector = np.array(
            [complex(x["re"], x["im"]) for x in band["eigenvector"]], dtype=complex
        )
        norm = np.linalg.norm(vector) or 1.0
        vector = vector / norm
        overlaps.append(float(max(abs(np.vdot(b, vector)) for b in basis)))
    return overlaps


def gamma_acoustic_report(ph: Any) -> dict[str, Any]:
    report = qpoint_report(ph, (0.0, 0.0, 0.0))
    frequencies = report["frequencies_THz"]
    overlaps = _translational_overlaps(ph, report)
    # the three most translation-like branches ARE the acoustic branches
    acoustic_bands = sorted(
        int(i) for i in np.argsort(overlaps)[::-1][:3]
    )
    acoustic = [frequencies[i] for i in acoustic_bands]
    optical = [
        frequencies[i] for i in range(len(frequencies)) if i not in acoustic_bands
    ]
    return {
        "frequency_THz": [float(x) for x in frequencies],
        "acoustic_bands": acoustic_bands,
        "acoustic_overlap": [float(overlaps[i]) for i in acoustic_bands],
        "acoustic_frequencies_THz": [float(x) for x in acoustic],
        "max_abs_THz": max(abs(x) for x in acoustic) if acoustic else 0.0,
        "noise_ok": (
            max(abs(x) for x in acoustic) < GAMMA_ACOUSTIC_TOL_THz
            if acoustic
            else True
        ),
        "optical_frequencies_THz": [float(x) for x in optical],
        "optical_min_THz": float(min(optical)) if optical else None,
        # an optical branch gets NO acoustic tolerance: any negative value
        # beyond print noise is a soft mode
        "optical_soft_THz": [
            float(x) for x in optical if x < -SOFT_FREQUENCY_TOL_THz
        ],
    }


def mesh_diagnostic(ph: Any, mesh: tuple[int, int, int]) -> dict[str, Any]:
    """Interpolated frequencies on a dense mesh: context, never a verdict.

    The mesh minimum is not a wavevector the force constants directly
    represent; it is what Fourier interpolation produces between them.
    """

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
        "n_mesh_imaginary_modes": int((frequencies < -SOFT_FREQUENCY_TOL_THz).sum()),
        "n_modes": int(frequencies.size),
        "note": "interpolated diagnostic only; not the instability wavevector",
    }


class DfptStabilityReport(BaseModel):
    """Serializable verdict for one DFPT volume.

    The report is a stored artifact, not a log line: when the gate halts a flow
    the analysis job is never produced, so this is what the Worker's
    ``get_result`` falls back to in order to say *why* the run stopped.

    ``qha_allowed`` is tied to ``status`` rather than carried independently.
    Two sources of truth for the same decision is how a soft-mode volume ends
    up in a free-energy fit with a clean-looking ``ok``.
    """

    model_config = ConfigDict(extra="forbid")

    schema_version: Literal["htvasp-dfpt-stability/v1"] = SCHEMA_VERSION
    status: Literal[
        "harmonic_QHA_applicable",
        "finite_q_soft_mode",
        "incommensurate_soft_mode_suspected",
        "gamma_soft_mode",
        "numerical_unresolved",
    ]
    qha_allowed: bool
    reason: str
    #: which volume this is; absent for a standalone diagnosis
    index: int | None = None
    volume_angstrom3: float | None = None
    checks: dict[str, Any]
    commensurate_qpoints: list[dict[str, Any]]
    mesh_diagnostic: dict[str, Any]
    #: identifies the matrix the verdict was computed from
    force_constants_sha256: str = Field(min_length=64, max_length=64)
    structure_hash: str = ""
    recipe_hash: str = ""
    units: dict[str, str] = Field(default_factory=dict)
    #: a borderline soft mode needs a larger supercell before it is believed
    confirmation_required: bool = False
    confirmation_reason: str = ""
    warnings: list[str] = Field(default_factory=list)
    supercell_convergence: str = "not_checked"

    @model_validator(mode="after")
    def _status_decides_qha(self) -> "DfptStabilityReport":
        expected = self.status == "harmonic_QHA_applicable"
        if self.qha_allowed != expected:
            raise ValueError(
                f"qha_allowed={self.qha_allowed} contradicts status={self.status!r}"
            )
        return self


def unresolved_report(
    reason: str,
    *,
    force_constants: Any = None,
    structure_hash: str = "",
    recipe_hash: str = "",
    checks: dict[str, Any] | None = None,
    index: int | None = None,
    volume_angstrom3: float | None = None,
) -> dict[str, Any]:
    """A ``numerical_unresolved`` verdict for a diagnosis that could not run.

    The alternative -- letting the diagnosis job raise -- loses the artifact
    exactly when it is most needed: the matrix is unusable and the only record
    of that would be a traceback.
    """

    matrix = np.asarray([] if force_constants is None else force_constants, dtype=float)
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "numerical_unresolved",
        "qha_allowed": False,
        "reason": reason,
        "index": index,
        "volume_angstrom3": volume_angstrom3,
        "checks": {"shape_ok": False, "problems": [reason], **(checks or {})},
        "commensurate_qpoints": [],
        "mesh_diagnostic": {},
        "force_constants_sha256": force_constants_sha256(matrix),
        "structure_hash": structure_hash,
        "recipe_hash": recipe_hash,
        "units": {"frequency": "THz", "force_constants": "eV/Angstrom^2"},
        "warnings": [],
    }


def _softest(q_reports: list[dict[str, Any]]) -> float:
    return min((r["min_frequency_THz"] for r in q_reports), default=0.0)


def diagnose_stability(
    ph: Any,
    force_constants: Any,
    *,
    structure: Any,
    recipe: Any,
    supercell_matrix: Any,
    mesh: tuple[int, int, int],
    ordering_ok: bool = True,
    index: int | None = None,
    volume_angstrom3: float | None = None,
) -> dict[str, Any]:
    """Three-tier stability verdict for one volume.

    ``ordering_ok`` is asserted by the caller before the matrix is attached;
    recording it here keeps it part of the structural verdict rather than an
    unstated precondition.
    """

    checks = force_constants_checks(force_constants)
    checks["ordering_ok"] = bool(ordering_ok)
    reciprocal = reciprocal_check(ph)
    gamma = gamma_acoustic_report(ph)
    qpoints = commensurate_qpoints(ph)
    q_reports = [qpoint_report(ph, q) for q in qpoints]
    gamma_report = next(r for r in q_reports if r["q_fractional"] == [0.0, 0.0, 0.0])
    finite_reports = [
        r for r in q_reports if r["q_fractional"] != [0.0, 0.0, 0.0]
    ]
    mesh_report = mesh_diagnostic(ph, mesh)

    warnings: list[str] = []
    confirmation_reason = ""
    structural_ok = bool(
        checks.get("shape_ok")
        and checks.get("asr_ok")
        and checks.get("ordering_ok", False)
        and reciprocal["ok"]
    )
    # Gamma: the first three branches are acoustic and are allowed numerical
    # noise.  Everything above them is an optical branch and gets no allowance,
    # however close to zero it sits.
    gamma_optical_soft = gamma["optical_soft_THz"]
    softest_finite = _softest(finite_reports)
    mesh_min = mesh_report["mesh_min_frequency_THz"]

    if not structural_ok:
        status = "numerical_unresolved"
        reason = "tier-A structural check failed"
    elif gamma_optical_soft:
        status = "gamma_soft_mode"
        reason = (
            f"Gamma has {len(gamma_optical_soft)} imaginary optical branch(es), "
            f"softest {min(gamma_optical_soft):.4f} THz"
        )
    elif softest_finite < -SOFT_FREQUENCY_TOL_THz:
        status = "finite_q_soft_mode"
        reason = (
            f"imaginary branch at an exact commensurate q point "
            f"({softest_finite:.4f} THz)"
        )
        if abs(softest_finite) < BORDERLINE_CONFIRMATION_TOL_THz:
            confirmation_reason = (
                f"softest exact-q branch is {softest_finite:.4f} THz, within "
                f"{BORDERLINE_CONFIRMATION_TOL_THz} THz of zero; a larger "
                f"supercell is needed to separate a real instability from "
                f"supercell truncation"
            )
    elif mesh_min < -MESH_SOFT_SUSPECT_TOL_THz:
        status = "incommensurate_soft_mode_suspected"
        reason = (
            f"dense mesh minimum {mesh_min:.4f} THz at q = "
            f"{mesh_report['mesh_min_qpoint']}, which no commensurate q point "
            f"represents; the supercell cannot decide this instability"
        )
        if abs(mesh_min) < BORDERLINE_CONFIRMATION_TOL_THz:
            confirmation_reason = (
                f"mesh minimum {mesh_min:.4f} THz is within "
                f"{BORDERLINE_CONFIRMATION_TOL_THz} THz of zero"
            )
    else:
        status = "harmonic_QHA_applicable"
        reason = "no imaginary branch at any commensurate q point"
        if mesh_min < -SOFT_FREQUENCY_TOL_THz:
            warnings.append(
                f"dense mesh minimum {mesh_min:.4f} THz at "
                f"{mesh_report['mesh_min_qpoint']} is treated as an "
                f"interpolation artifact: it is shallower than "
                f"{MESH_SOFT_SUSPECT_TOL_THz} THz and no exact commensurate q "
                f"point is imaginary"
            )

    return {
        "schema_version": SCHEMA_VERSION,
        "status": status,
        "qha_allowed": status == "harmonic_QHA_applicable",
        "reason": reason,
        "index": index,
        "volume_angstrom3": volume_angstrom3,
        "checks": {**checks, "reciprocal": reciprocal, "gamma": gamma},
        "commensurate_qpoints": q_reports,
        "mesh_diagnostic": mesh_report,
        "force_constants_sha256": force_constants_sha256(force_constants),
        "structure_hash": ordering_fingerprint(structure),
        "recipe_hash": recipe_sha256(recipe),
        "units": {"frequency": "THz", "force_constants": "eV/Angstrom^2"},
        "confirmation_required": bool(confirmation_reason),
        "confirmation_reason": confirmation_reason,
        "warnings": warnings,
        "supercell_convergence": "not_checked",
    }


__all__ = [
    "BORDERLINE_CONFIRMATION_TOL_THz",
    "DfptStabilityReport",
    "GAMMA_ACOUSTIC_TOL_THz",
    "MESH_SOFT_SUSPECT_TOL_THz",
    "SOFT_FREQUENCY_TOL_THz",
    "commensurate_qpoints",
    "diagnose_stability",
    "force_constants_checks",
    "force_constants_sha256",
    "gamma_acoustic_report",
    "mesh_diagnostic",
    "qpoint_report",
    "recipe_sha256",
    "reciprocal_check",
    "reciprocal_lattice",
    "unresolved_report",
]
