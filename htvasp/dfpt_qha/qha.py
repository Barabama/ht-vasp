"""DFPT-QHA thermal analysis from VASP force constants.

This is the terminal consumer the fail-closed gate protects.  It never builds a
finite-displacement dataset: the force constants come from each volume's DFPT
``FORCE_CONSTANTS`` (see ``jobs.collect_dfpt_task_doc``), and phonopy is used
only to integrate the thermal properties of that dynamical matrix.

Physics
-------
At each volume ``V``::

    F(V, T) = E_static(V) + F_vib(V, T)

``E_static`` is the DFT energy of the volume's relaxed cell and ``F_vib`` is
phonopy's vibrational Helmholtz free energy, which already contains the
zero-point term.  The quasi-harmonic step then minimises

    G(T, p) = min_V [ F(V, T) + p V ]

over the sampled volumes at every temperature.

Units are carried explicitly rather than implied: phonopy reports kJ/mol while
the DFT energies are eV per cell, and the two are only comparable through the
cell's atom count.
"""

from __future__ import annotations

from typing import Any

import logging

import numpy as np
from jobflow import Response, job
from pymatgen.core import Structure

from htvasp.dfpt_qha.phonons import (
    as_structure,
    assert_same_ordering,
    ordering_fingerprint,
    phonopy_for,
)
from htvasp.dfpt_qha.stability import SOFT_FREQUENCY_TOL_THz

log = logging.getLogger(__name__)

#: 1 eV/atom expressed in kJ/mol (elementary charge x Avogadro).
EV_PER_ATOM_TO_KJ_PER_MOL = 96.48533212331002
#: 1 eV/Angstrom^3 in GPa.
EV_PER_ANGSTROM3_TO_GPA = 160.21766208

DEFAULT_MESH = (12, 12, 12)
DEFAULT_TEMPERATURES = [float(t) for t in range(0, 1001, 100)]


def volume_free_energy(
    cell: Structure,
    supercell: Structure,
    force_constants: Any,
    *,
    supercell_matrix: Any,
    mesh: tuple[int, int, int],
    temperatures: list[float],
) -> dict[str, Any]:
    """Vibrational free energy of one volume, from its DFPT force constants.

    ``supercell`` is the cell the DFPT job actually ran on; the matrix is only
    attached after phonopy's own supercell is shown to have the same atom
    ordering, because a mismatch still yields plausible-looking frequencies.
    """

    phonon = phonopy_for(cell, supercell_matrix)
    phonopy_supercell = Structure(
        np.array(phonon.supercell.cell),
        list(phonon.supercell.symbols),
        [np.array(p) for p in phonon.supercell.scaled_positions],
        coords_are_cartesian=False,
    )
    assert_same_ordering(supercell, phonopy_supercell)

    fc = np.array(force_constants, dtype=float)
    n_supercell = len(phonopy_supercell)
    if fc.shape[:2] != (n_supercell, n_supercell):
        raise ValueError(
            f"force constants are indexed for {fc.shape[:2]} atoms but the "
            f"supercell has {n_supercell}"
        )
    phonon.force_constants = fc

    step = 100
    if len(temperatures) > 1:
        step = max(1, int(round(temperatures[1] - temperatures[0])))
    phonon.run_mesh(list(mesh), with_eigenvectors=False)
    phonon.run_thermal_properties(
        t_min=min(temperatures), t_max=max(temperatures), t_step=step
    )
    thermal = phonon.thermal_properties
    mesh_data = phonon.get_mesh_dict()
    mesh_frequencies = np.array(mesh_data["frequencies"])
    mesh_qpoints = np.array(mesh_data["qpoints"])
    # phonopy >= 2.4x exposes the integrated curves as a (4, n_temperatures)
    # array -- temperatures, free energy, entropy, heat capacity -- and no
    # longer as the named .free_energy / .entropy / .heat_capacity attributes.
    curve = np.array(thermal.thermal_properties, dtype=float)
    if curve.ndim != 2 or curve.shape[0] < 4:
        raise ValueError(f"unexpected phonopy thermal_properties shape: {curve.shape}")
    return {
        "free_energy_kJ_per_mol": [float(x) for x in curve[1]],
        "entropy_J_per_mol_K": [float(x) for x in curve[2]],
        "heat_capacity_J_per_mol_K": [float(x) for x in curve[3]],
        "zero_point_energy_kJ_per_mol": float(thermal.zero_point_energy),
        "thermal_temperatures": [float(x) for x in curve[0]],
        "n_atoms_supercell": int(n_supercell),
        "n_atoms_cell": len(cell),
        # phonopy's thermal properties are per mole of ITS primitive cell, which
        # is not always the recipe cell (fcc SER: 1 atom vs 4).
        "n_atoms_primitive": len(phonon.primitive),
        "supercell_ordering_fingerprint": ordering_fingerprint(supercell),
        # Interpolated context, NOT a stability verdict.  The verdict is made
        # on exact supercell-commensurate q points by
        # ``classify_dfpt_stability``; a shallow negative minimum here can be a
        # Fourier-interpolation artifact, and a real instability can sit at a q
        # the supercell cannot represent.
        "mesh_diagnostic": {
            "mesh": list(mesh),
            "mesh_min_THz": float(mesh_frequencies.min()),
            "mesh_min_qpoint": [float(x) for x in mesh_qpoints[
                int(np.unravel_index(int(np.argmin(mesh_frequencies)),
                                     mesh_frequencies.shape)[0])
            ]],
            "n_mesh_imaginary_modes": int(
                (mesh_frequencies < -SOFT_FREQUENCY_TOL_THz).sum()
            ),
            "n_modes": int(mesh_frequencies.size),
            "note": "interpolated diagnostic only; not a stability verdict",
        },
    }


def birch_murnaghan_3rd_energy(
    volume: Any, e0: float, b0: float, bp: float, v0: float
) -> Any:
    """The BM3 energy surface, in eV with ``b0`` in eV/Angstrom^3.

    Module level so a fit can be evaluated at the sampled volumes (residuals)
    or between them (the pressurised minimum), instead of only at its own
    parameters.
    """

    v = np.asarray(volume, dtype=float)
    eta = (v0 / v) ** (2.0 / 3.0) - 1.0
    return e0 + 9.0 * v0 * b0 / 16.0 * (eta ** 3 * bp + eta ** 2 * (6.0 - 4.0 * (v0 / v) ** (2.0 / 3.0)))


def birch_murnaghan_3rd(volume: np.ndarray, energy: np.ndarray) -> dict[str, float]:
    """BM3 fit of E(V), returning the minimum and the bulk modulus."""

    from scipy.optimize import curve_fit

    bm3 = birch_murnaghan_3rd_energy

    popt, _ = curve_fit(
        bm3, volume, energy,
        p0=[float(energy.min()), 1.5, 4.0, float(volume[int(np.argmin(energy))])],
        maxfev=20000,
    )
    e0, b0, bp, v0 = (float(x) for x in popt)
    return {
        "E0_eV": e0,
        "V0_angstrom3": v0,
        "B0_GPa": b0 * EV_PER_ANGSTROM3_TO_GPA,
        "B0prime": bp,
    }


def continuous_thermal_eos(
    volumes: Any,
    static_energies_eV: Any,
    free_energies_eV: Any,
    temperatures: Any,
    *,
    pressure_gpa: float = 0.0,
    n_atoms_cell: int = 1,
    residual_tol_mev_per_atom: float = 10.0,
) -> dict[str, Any]:
    """Per-temperature EOS fit of F(V,T): continuous V(T), G(T) and B(T).

    ``analyze_dfpt_qha`` reports V(T) as the argmin over the sampled volumes,
    so the curve is piecewise constant and its temperature derivative does not
    exist.  A TDB fitted to that curve inherits the jumps.  Here the
    Birch-Murnaghan form is fitted at every temperature instead, which gives a
    V(T) that can be differentiated.

    Pure post-processing: the inputs are the arrays the analysis already
    stores, so this needs no VASP run and does not rewrite the analysis JSON.

    ``volumes`` is per recipe cell (Angstrom^3), ``static_energies_eV`` is the
    DFT energy per cell, and ``free_energies_eV`` is ``[volume][temperature]``
    per cell.  At each temperature the fitted total is
    ``E_static(V) + F_vib(V,T)``; the Gibbs energy is its minimum (plus ``PV``
    when a pressure is given).

    The residual is the diagnostic that decides whether the curve is usable: a
    healthy endpoint in this pool fits to 1-6 meV/atom, while a volume whose
    F(V,T) is not EOS-shaped (a soft volume developing a second dip at large V)
    grows to tens of meV/atom, and the fitted V(T) and G(T) are then an
    artefact of the functional form rather than of the data.
    """

    volume = np.asarray(volumes, dtype=float)
    static = np.asarray(static_energies_eV, dtype=float)
    free = np.asarray(free_energies_eV, dtype=float)
    t_grid = [float(x) for x in temperatures]

    if volume.ndim != 1 or static.shape != volume.shape:
        raise ValueError("volumes and static energies must be one-dimensional and equal in length")
    if free.shape != (volume.size, len(t_grid)):
        raise ValueError(
            f"free energies have shape {free.shape}, expected "
            f"{(volume.size, len(t_grid))}"
        )
    if volume.size < 5:
        # four parameters on four points is an interpolation, not a fit: the
        # residual below could not detect a bad EOS.
        raise ValueError(
            f"a per-temperature EOS fit needs at least 5 volumes, got {volume.size}"
        )
    if np.any(np.diff(volume) <= 0):
        raise ValueError("volumes must be strictly increasing")
    if not np.all(np.isfinite(free)):
        raise ValueError("free energies contain non-finite values")

    pressure = pressure_gpa / EV_PER_ANGSTROM3_TO_GPA
    v_lo, v_hi = float(volume.min()), float(volume.max())
    grid = np.linspace(v_lo, v_hi, 2001)

    v_fit: list[float] = []
    gibbs: list[float] = []
    bulk: list[float] = []
    params: list[dict[str, float]] = []
    residuals: list[float] = []
    argmin_v: list[float] = []
    at_edge: list[bool] = []

    for i, temperature in enumerate(t_grid):
        total = static + free[:, i]
        fit = birch_murnaghan_3rd(volume, total)
        b0 = fit["B0_GPa"] / EV_PER_ANGSTROM3_TO_GPA
        model = birch_murnaghan_3rd_energy(volume, fit["E0_eV"], b0, fit["B0prime"], fit["V0_angstrom3"])
        residuals.append(float(np.abs(model - total).max() * 1e3 / max(n_atoms_cell, 1)))

        if pressure == 0.0:
            v0 = float(fit["V0_angstrom3"])
        else:
            curve = birch_murnaghan_3rd_energy(grid, fit["E0_eV"], b0, fit["B0prime"], fit["V0_angstrom3"]) + pressure * grid
            v0 = float(grid[int(np.argmin(curve))])
        at_edge.append(bool(v0 <= v_lo * 1.0005 or v0 >= v_hi * 0.9995))

        h = 1e-3 * v0
        e_plus = float(birch_murnaghan_3rd_energy(v0 + h, fit["E0_eV"], b0, fit["B0prime"], fit["V0_angstrom3"]))
        e_minus = float(birch_murnaghan_3rd_energy(v0 - h, fit["E0_eV"], b0, fit["B0prime"], fit["V0_angstrom3"]))
        second = (e_plus - 2.0 * float(birch_murnaghan_3rd_energy(v0, fit["E0_eV"], b0, fit["B0prime"], fit["V0_angstrom3"])) + e_minus) / h ** 2
        bulk.append(float(v0 * second * EV_PER_ANGSTROM3_TO_GPA))

        g = float(birch_murnaghan_3rd_energy(v0, fit["E0_eV"], b0, fit["B0prime"], fit["V0_angstrom3"])) + pressure * v0
        v_fit.append(v0)
        gibbs.append(g)
        params.append({**fit, "temperature": temperature})

        sampled = total + pressure * volume
        argmin_v.append(float(volume[int(np.argmin(sampled))]))

    v_array = np.asarray(v_fit, dtype=float)
    b_array = np.asarray(bulk, dtype=float)
    r_array = np.asarray(residuals, dtype=float)
    checks = {
        "n_volumes": int(volume.size),
        "n_temperatures": len(t_grid),
        "volumes_strictly_increasing": bool(np.all(np.diff(volume) > 0)),
        "volume_monotonic_in_temperature": bool(np.all(np.diff(v_array) >= -1e-12)),
        "all_inside_sampled_window": bool(not any(at_edge)),
        "all_bulk_moduli_positive": bool(np.all(b_array > 0.0)),
        "residual_within_tolerance": bool(np.all(r_array <= residual_tol_mev_per_atom)),
        "max_residual_meV_per_atom": float(r_array.max()),
        "residual_meV_per_atom": [float(x) for x in r_array],
        "residual_tolerance_meV_per_atom": float(residual_tol_mev_per_atom),
        "max_volume_decrease_angstrom3": float(
            max(0.0, -float(np.diff(v_array).min())) if v_array.size > 1 else 0.0
        ),
        "equilibrium_at_window_edge_temperatures": [
            t for t, edge in zip(t_grid, at_edge) if edge
        ],
    }
    return {
        "schema_version": "htvasp-dfpt-qha-vt/v1",
        "pressure_gpa": float(pressure_gpa),
        "volumes": volume.tolist(),
        "temperatures": t_grid,
        "equilibrium_volume_angstrom3": v_fit,
        "gibbs_eV": gibbs,
        "bulk_modulus_GPa": bulk,
        "eos_parameters": params,
        "argmin_volume_angstrom3": argmin_v,
        "max_abs_volume_difference_angstrom3": float(
            np.abs(v_array - np.asarray(argmin_v, dtype=float)).max()
        ),
        "checks": checks,
        "units": {
            "equilibrium_volume_angstrom3": "Angstrom^3 per cell",
            "gibbs_eV": "eV per cell",
            "bulk_modulus_GPa": "GPa",
        },
    }


@job
def dfpt_thermal_properties(
    collected: dict[str, Any],
    *,
    supercell_matrix: Any,
    mesh: tuple[int, int, int] = DEFAULT_MESH,
    temperatures: list[float] | None = None,
) -> dict[str, Any]:
    """Vibrational free energy at every validated DFPT volume.

    Each record must already have passed ``validate_dfpt_task_doc`` and
    ``classify_dfpt_stability``; a Gamma-only result cannot reach here because
    the validator compares the atom count against the recipe's phonon
    supercell.

    A verdict of "not harmonic-QHA-applicable" does NOT stop this job.  The
    free energy is still integrated and the verdict is carried forward, so one
    DFPT run yields both the diagnostic and the (labelled) G(T).  Only a
    missing verdict fails, because that means the flow was wired around the
    classification and the numbers would be unlabelled.
    """

    temperatures = list(temperatures or DEFAULT_TEMPERATURES)
    per_volume = []
    for record in collected["volumes"]:
        stability = record.get("stability") or {}
        if "status" not in stability:
            return Response(
                output={
                    "ok": False,
                    "error": (
                        f"volume {record.get('index')} reached the thermal job "
                        f"without a stability verdict (status={stability.get('status')!r}); "
                        f"the classification job was bypassed"
                    ),
                },
                stop_children=True,
            )
        result = volume_free_energy(
            as_structure(record["cell"]),
            as_structure(record["supercell"]),
            record["force_constants"],
            supercell_matrix=supercell_matrix,
            mesh=mesh,
            temperatures=temperatures,
        )
        result.update(
            {
                "index": int(record["index"]),
                "volume_angstrom3": float(record["volume_angstrom3"]),
                "static_energy_eV": float(record["static_energy_eV"]),
                "stability_status": (record.get("stability") or {}).get("status"),
            }
        )
        per_volume.append(result)
    per_volume.sort(key=lambda r: r["volume_angstrom3"])

    # phonopy builds its own temperature grid from (min, max, step), so the
    # grid it actually evaluated is not necessarily the one that was requested
    # -- asking for 0, 300, 600, 1000 yields 0, 300, 600, 900.  Reporting the
    # requested list while indexing phonopy's arrays by position would pair
    # every free energy with the wrong temperature, and the curve would still
    # look smooth.  The evaluated grid is therefore the one carried forward,
    # and it must be the same at every volume or the columns do not line up.
    grids = {tuple(r["thermal_temperatures"]) for r in per_volume}
    if len(grids) != 1:
        raise ValueError(
            f"volumes evaluated different temperature grids: "
            f"{sorted(len(g) for g in grids)} points; the free-energy columns "
            f"cannot be compared"
        )
    evaluated = list(grids.pop())
    if not np.allclose(evaluated, temperatures):
        log.warning(
            "phonopy evaluated %s K rather than the requested %s K; the "
            "evaluated grid is what the analysis will report",
            evaluated, temperatures,
        )
    return {
        "volumes": per_volume,
        "temperatures": evaluated,
        "requested_temperatures": temperatures,
        "rejected_volume_indices": collected.get("rejected_volume_indices", []),
        "n_validated": collected.get("n_validated"),
        "n_requested": collected.get("n_requested"),
        # carried through from classify_dfpt_stability: the curve is computed
        # either way, but a consumer must be able to tell what it is
        "qha_applicability": collected.get(
            "qha_applicability", "harmonic_QHA_applicable"
        ),
        "non_harmonic_volume_indices": collected.get(
            "non_harmonic_volume_indices", []
        ),
    }


@job
def analyze_dfpt_qha(
    thermal: dict[str, Any],
    *,
    pressure_gpa: float = 0.0,
) -> Response:
    """Gibbs free energy from the E(V) and F(V,T) curves.

    Fails closed on a non-physical input instead of returning a smooth-looking
    curve: too few volumes, unstable phonons, or an EOS minimum outside the
    sampled window are reported, not extrapolated over.
    """

    per_volume = thermal.get("volumes") or []
    problems: list[str] = []
    if len(per_volume) < 4:
        return Response(
            output={
                "ok": False,
                "error": f"QHA needs at least 4 validated volumes, got {len(per_volume)}",
                "rejected_volume_indices": thermal.get("rejected_volume_indices", []),
            },
            stop_children=True,
        )

    volumes = np.array([r["volume_angstrom3"] for r in per_volume], dtype=float)
    static = np.array([r["static_energy_eV"] for r in per_volume], dtype=float)
    temperatures = [float(t) for t in thermal["temperatures"]]
    n_atoms_cell = int(per_volume[0]["n_atoms_cell"])

    if np.any(np.diff(volumes) <= 0):
        problems.append("volumes are not strictly increasing after sorting")
    # Deliberately NOT a gate.  The mesh minimum is an interpolated sample: a
    # shallow negative value there is a truncation artifact, and gating on it
    # rejected runs whose exact commensurate q points are all stable.  The
    # harmonic-applicability decision was made per volume by
    # ``classify_dfpt_stability``; this only carries the numbers forward.
    mesh_diagnostics = [
        {"index": int(r["index"]), **(r.get("mesh_diagnostic") or {})}
        for r in per_volume
    ]

    # phonopy reports kJ per mole of its own PRIMITIVE cell (it divides by
    # ``sum(weights)``), while the DFT energies are per recipe cell.  The two
    # coincide for B2, where the primitive cell IS the 2-atom cell, but not for
    # an fcc SER entry, whose primitive cell is 1 atom against a 4-atom
    # conventional cell.  The ratio is therefore explicit rather than assumed.
    # Multiplying by the atom count alone -- a bug this replaced -- inflates
    # every free energy by n_atoms_cell; the zero-point energy is the check
    # that catches it (for the B2 fixture it lands near the Debye estimate of
    # the two elements, and the fcc factor-4 case is pinned by a test).
    per_volume_scale = [
        r["n_atoms_cell"] / r["n_atoms_primitive"] for r in per_volume
    ]
    if len(set(per_volume_scale)) != 1:
        problems.append(
            f"volumes disagree on the cell/primitive ratio: {per_volume_scale}"
        )
    free_energy = np.array(
        [
            [
                r["free_energy_kJ_per_mol"][i]
                / EV_PER_ATOM_TO_KJ_PER_MOL
                * (r["n_atoms_cell"] / r["n_atoms_primitive"])
                for i in range(len(temperatures))
            ]
            for r in per_volume
        ],
        dtype=float,
    )

    eos: dict[str, float] = {}
    try:
        eos = birch_murnaghan_3rd(volumes, static)
    except Exception as exc:  # noqa: BLE001 - the fit failing is the finding
        problems.append(f"BM3 fit failed: {exc}")

    if problems:
        return Response(
            output={"ok": False, "error": "; ".join(problems), "eos": eos},
            stop_children=True,
        )

    pressure_ev_per_angstrom3 = pressure_gpa / EV_PER_ANGSTROM3_TO_GPA
    gibbs: list[float] = []
    equilibrium_volume: list[float] = []
    for i in range(len(temperatures)):
        total = static + free_energy[:, i] + pressure_ev_per_angstrom3 * volumes
        best = int(np.argmin(total))
        gibbs.append(float(total[best]))
        equilibrium_volume.append(float(volumes[best]))

    # ``equilibrium_volume_angstrom3`` is the argmin over the sampled volumes,
    # so it is piecewise constant (a 7-point scan gives V(T) that jumps between
    # grid points) and its temperature derivative does not exist.  The
    # continuous curve comes from ``continuous_thermal_eos`` below, which fits
    # the Birch-Murnaghan form at each temperature from the same per-volume
    # F(V,T); that is post-processing and needs no new VASP run.
    inside = bool(volumes.min() < eos["V0_angstrom3"] < volumes.max())
    if not inside:
        problems.append(
            f"EOS minimum V0={eos['V0_angstrom3']:.4f} sits outside the sampled "
            f"window [{volumes.min():.4f}, {volumes.max():.4f}]; V0 is then an "
            f"artefact of where the scan stopped"
        )

    return Response(
        output={
            "ok": True,
            "schema_version": "htvasp-dfpt-qha-analysis/v1",
            "temperatures": temperatures,
            "volumes": volumes.tolist(),
            "static_energies_eV": static.tolist(),
            "free_energies_eV": free_energy.tolist(),
            "gibbs_eV": gibbs,
            "equilibrium_volume_angstrom3": equilibrium_volume,
            "eos": eos,
            "v0_inside_sampled_window": inside,
            "warnings": problems,
            "pressure_gpa": pressure_gpa,
            "n_atoms_cell": n_atoms_cell,
            "rejected_volume_indices": thermal.get("rejected_volume_indices", []),
            "n_validated": thermal.get("n_validated"),
            "n_requested": thermal.get("n_requested"),
            # what the product is, not just whether the fit converged
            "qha_applicability": thermal.get(
                "qha_applicability", "harmonic_QHA_applicable"
            ),
            "non_harmonic_volume_indices": thermal.get(
                "non_harmonic_volume_indices", []
            ),
            "volume_stability": [
                {
                    "index": int(r["index"]),
                    "status": r.get("stability_status"),
                }
                for r in per_volume
            ],
            "dynamically_unstable_conditional": bool(
                thermal.get("non_harmonic_volume_indices")
            ),
            # carried for the record; the mesh never decided anything
            "mesh_diagnostics": mesh_diagnostics,
            "units": {
                "free_energies_eV": "eV per cell",
                "gibbs_eV": "eV per cell",
                "equilibrium_volume_angstrom3": "Angstrom^3 per cell",
                "B0_GPa": "GPa",
            },
        }
    )


__all__ = [
    "DEFAULT_MESH",
    "DEFAULT_TEMPERATURES",
    "EV_PER_ANGSTROM3_TO_GPA",
    "EV_PER_ATOM_TO_KJ_PER_MOL",
    "analyze_dfpt_qha",
    "birch_murnaghan_3rd_energy",
    "continuous_thermal_eos",
    "birch_murnaghan_3rd",
    "dfpt_thermal_properties",
    "volume_free_energy",
]
