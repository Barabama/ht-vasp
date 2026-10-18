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
        # Fourier-interpolation artifact (BCC-Nb-Nb-2: -0.024 THz), and a real
        # instability can sit at a q the supercell cannot represent.
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


def birch_murnaghan_3rd(volume: np.ndarray, energy: np.ndarray) -> dict[str, float]:
    """BM3 fit of E(V), returning the minimum and the bulk modulus."""

    from scipy.optimize import curve_fit

    def bm3(v, e0, b0, bp, v0):
        return e0 + 9.0 * v0 * b0 / 16.0 * (
            ((v0 / v) ** (2.0 / 3.0) - 1.0) ** 3 * bp
            + ((v0 / v) ** (2.0 / 3.0) - 1.0) ** 2
            * (6.0 - 4.0 * (v0 / v) ** (2.0 / 3.0))
        )

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
    supercell, and a soft-mode volume cannot reach here because the
    classification job stops the flow.  The ``qha_allowed`` check is repeated
    so a hand-assembled flow gets the same answer.
    """

    temperatures = list(temperatures or DEFAULT_TEMPERATURES)
    per_volume = []
    for record in collected["volumes"]:
        stability = record.get("stability") or {}
        if not stability.get("qha_allowed"):
            return Response(
                output={
                    "ok": False,
                    "error": (
                        f"volume {record.get('index')} reached the thermal job "
                        f"without a harmonic-QHA verdict "
                        f"(status={stability.get('status')!r})"
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
    # that catches it (30 meV/atom here, against 32 meV/atom from the Debye
    # temperatures of Al and Ta).
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
    "birch_murnaghan_3rd",
    "dfpt_thermal_properties",
    "volume_free_energy",
]
