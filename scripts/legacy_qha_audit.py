#!/usr/bin/env python
"""Fail-closed audit of an ingested legacy DFPT-QHA dataset.

This is a READ-ONLY post-processing tool. It never runs VASP, never writes
into ``legacy-as-received/``, and never modifies any source file.

What it establishes, per volume point:

1. the VASP outputs parse and the electronic structure converged
   (the same fail-closed notion the BCC-Ti-V pilot uses, not merely
   "the job wrote a vasprun.xml");
2. the DFPT force constants reproduce the recorded ``thermal_properties``
   (num_modes / num_integrated_modes / ZPE), which is what proves the
   phonon layer we are shipping is the one that produced the curve;
3. whether the phonon spectrum has imaginary modes, how many, where, and
   with what eigenvector character.

Then, across volumes: the E(V) fit and the reported Gibbs curve.

Usage:
    legacy_qha_audit.py --root data/legacy-dfpt-qha-pilot/BCC-Co-Ti \
        [--mesh 24] [--outdir <root>/audit]
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import re
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

# ---------------------------------------------------------------- helpers ---


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def parse_vasprun(path: Path) -> dict[str, Any]:
    """Parse the few fields we need without a full VASP output library.

    Deliberately tolerant: a legacy DFPT ``vasprun.xml`` is not a valid input
    to every consumer, so we read the XML ourselves and record what we found.
    """
    out: dict[str, Any] = {"parsed": False}
    try:
        tree = ET.parse(path)
    except ET.ParseError as ex:
        out["error"] = f"XML parse error: {ex}"
        return out
    root = tree.getroot()
    out["parsed"] = True

    ver = root.find("./generator/i[@name='version']")
    out["vasp_version"] = ver.text.strip() if ver is not None and ver.text else None

    # final volume / lattice from the last <structure>
    structs = root.findall(".//structure")
    if structs:
        s = structs[-1]
        basis = s.findall("./crystal/varray[@name='basis']/v")
        lat = [[float(x) for x in v.text.split()] for v in basis]
        if len(lat) == 3:
            out["final_lattice"] = lat
            out["final_volume"] = float(abs(np.linalg.det(np.array(lat))))

    # electronic convergence: every <scstep> carries <i name="e_wo_entrp">;
    # VASP's own convergence flag is the "reached required accuracy" line,
    # which lives in OUTCAR, not here. What we can read here is the final
    # energy and whether the SCF block completed at all.
    energies = root.findall(".//calculation/energy/i[@name='e_0_energy']")
    if not energies:
        energies = root.findall(".//calculation/energy/i[@name='e_fr_energy']")
    if energies:
        out["e_0_energy_eV"] = float(energies[-1].text)

    efermi = root.findall(".//calculation/dos/i[@name='efermi']")
    if efermi:
        out["efermi_eV"] = float(efermi[-1].text)

    out["scstep_blocks"] = len(root.findall(".//calculation/scstep"))
    out["ionic_steps"] = len(root.findall(".//calculation"))
    return out


def outcar_markers(path: Path) -> dict[str, Any]:
    """The same fail-closed markers the BCC-Ti-V pilot's ``cmd_verify`` uses.

    Kept deliberately identical so that "this dataset passes verification"
    means the same thing on both sides of the pilot.

    ``IBRION=8`` does not print ``reached required accuracy`` (that line belongs
    to the ionic optimiser); the DFPT convergence marker is
    ``aborting loop because EDIFF is reached``.
    """
    text = path.read_text(errors="replace")
    degrees = [
        (int(n), int(m))
        for n, m in re.findall(r"Degree of freedom:\s*(\d+)\s*/\s*(\d+)", text)
    ]
    return {
        "job_completion_marker": (
            "General timing and accounting informations for this job" in text
        ),
        "converged_scf_marker_count": text.count("aborting loop because EDIFF is reached"),
        "unconverged_scf_marker_count": text.count(
            "aborting loop EDIFF was not reached (unconverged)"),
        "dfpt_degrees_of_freedom": [{"index": n, "total": m} for n, m in degrees],
        "dfpt_all_perturbations_solved": bool(degrees) and degrees[-1][0] == degrees[-1][1],
        "charge_density_of_overlapping_atoms": (
            "charge density of overlapping atoms calculated" in text),
    }


def parse_fc_full(path: Path) -> np.ndarray:
    """Read FORCE_CONSTANTS ourselves, in its full (n_satom, n_satom, 3, 3) form.

    ``phonopy.load`` silently compacts the matrix to (n_patom, n_satom, 3, 3),
    which makes the permutation-symmetry diagnostic below ill-defined. We want
    the matrix exactly as the file stores it, so that the ASR and permutation
    numbers we report are properties of the shipped file, not of phonopy's
    in-memory view of it.
    """
    with path.open() as f:
        head = [int(x) for x in f.readline().split()]
        n, m = (head[0], head[0]) if len(head) == 1 else (head[0], head[1])
        fc = np.zeros((n, m, 3, 3), dtype=float)
        for i in range(n):
            for j in range(m):
                f.readline()                       # "i j" index line
                for k in range(3):
                    fc[i, j, k] = [float(x) for x in f.readline().split()]
    return fc


# VASP writes both energies on one line; `energy(sigma->0)` is the one the
# legacy v-e.dat files carry (verified to sub-ueV against R3f/OUTCAR).
_ENERGY_RE = re.compile(
    r"energy  without entropy=\s*(-?[\d.]+)\s+energy\(sigma->0\)\s*=\s*(-?[\d.]+)")


def static_energy(path: Path, log_path: Path | None = None,
                  whole_log: bool = False) -> dict[str, Any]:
    """The ground-state total energy VASP printed, plus its SCF loop shape.

    For an ``IBRION=8`` run the *first* energy block is the ground state; the
    linear-response blocks that follow are never summed into a total energy.

    VASP 5.4.4 writes the per-iteration ``DAV:``/``RMM:`` table to stdout only,
    so the iteration count is taken from ``log_path`` when that file was
    captured; otherwise it is counted in the OUTCAR and will be 0. For an
    ``IBRION=8`` log the count stops at the linear-response boundary; for a
    plain static (``whole_log=True``) the whole file is the one SCF loop.
    """
    text = path.read_text(errors="replace")
    hits = _ENERGY_RE.findall(text)
    out: dict[str, Any] = {"n_energy_blocks": len(hits)}
    if hits:
        out["energy_without_entropy_eV"] = float(hits[0][0])
        out["energy_sigma0_eV"] = float(hits[0][1])

    nelm = re.search(r"NELM\s*=\s*(\d+)", text)
    out["nelm"] = int(nelm.group(1)) if nelm else None
    ediff = re.search(r"EDIFF\s*=\s*([\d.eE+-]+)", text)
    out["ediif_eV"] = float(ediff.group(1)) if ediff else None

    # The ground-state electronic loop runs up to "Linear response reoptimize"
    # (IBRION=8) or, failing that, up to the first energy block.
    cut = text.find("Linear response reoptimize wavefunctions")
    if cut < 0:
        cut = text.find("free  energy   TOTEN")
    block = text[:cut] if cut > 0 else text
    src = "outcar"
    if log_path is not None and log_path.is_file():
        log_text = log_path.read_text(errors="replace")
        if whole_log:
            block, src = log_text, str(log_path.name)
        else:
            lcut = log_text.find("Linear response reoptimize")
            if lcut < 0:
                lcut = log_text.find("Linear response DOF")
            if lcut > 0:
                block, src = log_text[:lcut], str(log_path.name)
    out["scf_iteration_source"] = src
    iters = re.findall(r"^\s*(?:DAV|RMM):\s*\d+\s+(-?[\d.]+E?[-+]?\d*)\s+(-?[\d.]+E?[-+]?\d*)",
                       block, flags=re.M)
    out["ground_state_scf_iterations"] = len(iters)
    out["ground_state_hit_nelm"] = (
        out["nelm"] is not None and out["ground_state_scf_iterations"] >= out["nelm"])
    if iters:
        # How far the tail of the loop still wanders. A run that merely ran out
        # of NELM steps while sitting still is a very different thing from one
        # that diverged, and the two must not be reported with the same words.
        tail = [float(e) for _, e in iters[-10:]]
        out["ground_state_last_dE_eV"] = float(iters[-1][1])
        out["ground_state_tail_energy_span_meV"] = float(
            (max(tail) - min(tail)) * 1000.0)
    return out


def read_v_e(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows = []
    for line in path.read_text().splitlines():
        s = line.split()
        if len(s) >= 2:
            rows.append((float(s[0]), float(s[1])))
    a = np.array(rows)
    return a[:, 0], a[:, 1]


def read_curve(path: Path) -> tuple[np.ndarray, np.ndarray]:
    return read_v_e(path)


# ------------------------------------------------------------- phonopy ------

def phonon_audit(vol_dir: Path, mesh: int, log_level: int = 0) -> dict[str, Any]:
    """Rebuild the phonon spectrum from FORCE_CONSTANTS and audit stability."""
    import phonopy

    unitcell = vol_dir / "POSCAR-unitcell"
    fc_file = vol_dir / "FORCE_CONSTANTS"
    res: dict[str, Any] = {
        "unitcell": str(unitcell.name),
        "force_constants": str(fc_file.name),
    }

    # Diagnostics on the file as shipped (full matrix, phonopy-independent).
    fc_full = parse_fc_full(fc_file)
    res["force_constants_shape_full"] = list(fc_full.shape)
    res["force_constants_max_abs"] = float(np.abs(fc_full).max())
    res["asr_residual_max"] = float(np.abs(fc_full.sum(axis=1)).max())
    res["asr_residual_relative"] = float(
        np.abs(fc_full.sum(axis=1)).max() / np.abs(fc_full).max())
    res["permutation_symmetry_max_dev"] = float(
        np.abs(fc_full - fc_full.transpose(1, 0, 3, 2)).max())

    ph = phonopy.load(
        unitcell_filename=str(unitcell),
        supercell_matrix=np.diag([2, 2, 2]),
        force_constants_filename=str(fc_file),
        produce_fc=False,
        is_nac=False,
        log_level=log_level,
    )
    res["force_constants_shape_phonopy"] = list(ph.force_constants.shape)

    ph.run_mesh([mesh, mesh, mesh])
    md = ph.get_mesh_dict()
    freqs = md["frequencies"]                          # (nq, nbands) THz
    weights = np.asarray(md["weights"], dtype=float)
    nbands = freqs.shape[1]

    neg = freqs < 0.0
    total_modes = float(weights.sum() * nbands)
    neg_modes = float((neg * weights[:, None]).sum())

    res.update({
        "mesh": [mesh, mesh, mesh],
        "mesh_shifted_mp": True,
        "nqpoints_irreducible": int(freqs.shape[0]),
        "nbands": int(nbands),
        "total_modes_weighted": int(round(total_modes)),
        "imaginary_modes_weighted": int(round(neg_modes)),
        "imaginary_fraction": float(neg_modes / total_modes) if total_modes else None,
        "qpoints_with_imaginary": int(((neg * weights[:, None]).sum(axis=1) > 0).sum()),
        "min_frequency_THz": float(freqs.min()),
        "max_frequency_THz": float(freqs.max()),
    })

    imin = int(np.unravel_index(np.argmin(freqs), freqs.shape)[0])
    res["min_frequency_qpoint"] = [round(float(x), 6) for x in md["qpoints"][imin]]

    # Gamma-point frequencies, for the acoustic-sum-rule sanity check.
    ph.run_qpoints([[0.0, 0.0, 0.0]])
    res["gamma_frequencies_THz"] = [round(float(x), 6) for x in
                                    ph.get_qpoints_dict()["frequencies"][0]]

    # High-symmetry points of the bcc primitive cell, in primitive reciprocal
    # coordinates: N = (1/2,1/2,0), H = (1/2,1/2,1/2).
    ph.run_qpoints([[0.5, 0.5, 0.0], [0.5, 0.5, 0.5]])
    qd = ph.get_qpoints_dict()
    res["N_frequencies_THz"] = [round(float(x), 6) for x in qd["frequencies"][0]]
    res["H_frequencies_THz"] = [round(float(x), 6) for x in qd["frequencies"][1]]

    # Thermal properties, exactly as the legacy driver asked for.
    # phonopy 2.48.0 spells these t_min/t_max/t_step (the legacy driver's
    # phonopy 1.x used tmin/tmax/tstep). ``cutoff_frequency`` is left at
    # phonopy's default so that the silently-dropped imaginary modes show up
    # as num_modes > num_integrated_modes, exactly as in the legacy yaml.
    ph.run_thermal_properties(t_min=0, t_max=3000, t_step=50)
    tp = ph.thermal_properties
    res["recomputed"] = {
        "num_modes": int(tp.number_of_modes),
        "num_integrated_modes": int(tp.number_of_integrated_modes),
        "cutoff_frequency_THz": (None if tp.cutoff_frequency is None
                                 else float(tp.cutoff_frequency)),
        "zero_point_energy_kJ_per_mol": float(tp.zero_point_energy),
        "n_temperatures": int(len(tp.temperatures)),
    }
    res["recomputed"]["modes_dropped"] = (
        res["recomputed"]["num_modes"] - res["recomputed"]["num_integrated_modes"])

    # Lowest-frequency mode on the mesh, with its eigenvector character.
    # (Named "soft_mode" only when it is genuinely imaginary; at a stable
    # volume it is just the acoustic branch minimum and is labelled as such.)
    worst = md["qpoints"][imin]
    fr, ev = ph.get_frequencies_with_eigenvectors(np.array(worst))
    i = int(np.argmin(fr))
    masses = np.asarray(ph.unitcell.masses, dtype=float)   # PhonopyAtoms has .masses
    vec = ev[i]
    # phonopy eigenvectors are mass-weighted; undo that to get the real
    # Cartesian displacement pattern of each primitive-cell site.
    real = np.array([vec[3 * k:3 * k + 3] / math.sqrt(masses[k])
                     for k in range(len(masses))])
    res["soft_mode"] = {
        "qpoint": [round(float(x), 6) for x in worst],
        "frequency_THz": float(fr[i]),
        "is_imaginary": bool(fr[i] < 0.0),
        "degeneracy": int(np.sum(np.abs(fr - fr[i]) < 1e-6)),
        "symbols": list(ph.unitcell.symbols),
        "real_displacement_per_site": [[round(float(c), 6) for c in v] for v in real],
        "displacement_magnitude_per_site": [round(float(np.linalg.norm(v)), 6)
                                            for v in real],
    }
    return res


# ------------------------------------------------------------------- EOS ----

def birch_murnaghan3(V: np.ndarray, E0, V0, B0, Bp) -> np.ndarray:
    eta = (V0 / V) ** (2.0 / 3.0)
    return E0 + 9.0 * V0 * B0 / 16.0 * (
        (eta - 1.0) ** 3 * Bp + (eta - 1.0) ** 2 * (6.0 - 4.0 * eta)
    )


def fit_eos(V: np.ndarray, E: np.ndarray) -> dict[str, Any]:
    from scipy.optimize import curve_fit

    i0 = int(np.argmin(E))
    p0 = [E[i0], V[i0], 1.0, 4.0]
    out: dict[str, Any] = {}
    for name, p in (("bm3", p0),):
        try:
            popt, pcov = curve_fit(birch_murnaghan3, V, E, p0=p0, maxfev=200000)
            perr = np.sqrt(np.diag(pcov))
            resid = E - birch_murnaghan3(V, *popt)
            ss_res = float((resid ** 2).sum())
            ss_tot = float(((E - E.mean()) ** 2).sum())
            out[name] = {
                "E0_eV": float(popt[0]),
                "V0_ang3": float(popt[1]),
                "B0_eV_per_ang3": float(popt[2]),
                "B0_GPa": float(popt[2] * 160.21766208),
                "Bp": float(popt[3]),
                "sigma": [float(x) for x in perr],
                "r2": float(1.0 - ss_res / ss_tot) if ss_tot else None,
                "max_abs_residual_meV": float(np.abs(resid).max() * 1000.0),
            }
        except Exception as ex:  # noqa: BLE001 - report, never crash the audit
            out[name] = {"error": f"{type(ex).__name__}: {ex}"}
    out["discrete_minimum"] = {
        "V_ang3": float(V[i0]),
        "E_eV": float(E[i0]),
        "is_interior": bool(0 < i0 < len(V) - 1),
    }
    return out


# ------------------------------------------------------------------ main ----

VOL_RE = re.compile(r"^QHA-(\d+)$")


def audit_endpoint(root: Path, mesh: int, outdir: Path) -> dict[str, Any]:
    rec = root / "legacy-as-received"
    if not rec.is_dir():
        raise SystemExit(f"missing legacy-as-received under {root}")

    report: dict[str, Any] = {
        "schema": "legacy-qha-audit/v1",
        "endpoint_root": str(root.resolve()),
        "legacy_as_received": str(rec.resolve()),
        "audited_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "audited_on_host": __import__("os").uname().nodename,
        "phonopy_version": __import__("phonopy").__version__,
        "mesh": mesh,
        "volumes": {},
        "failures": [],
        "warnings": [],
        "notes": [],
    }

    v_e_path = rec / "v-e.dat"
    if not v_e_path.is_file():
        raise SystemExit(f"missing {v_e_path}")
    V, E = read_v_e(v_e_path)
    report["eos_input"] = {
        "file": "v-e.dat",
        "sha256": sha256_file(v_e_path),
        "n_points": int(len(V)),
        "volumes_ang3": [float(x) for x in V],
        "energies_eV": [float(x) for x in E],
        "units": {"volume": "angstrom^3", "energy": "eV"},
        # Established, not assumed: every row of v-e.dat is the
        # `energy(sigma->0)` of that volume's R3f static, matched to sub-ueV.
        "provenance": (
            "row i equals R3f/OUTCAR energy(sigma->0) of the i-th volume "
            "directory; verified per row in volumes[].v_e_row_matches_r3f_static_meV"
        ),
        "normalization": "unknown_from_the_files; the cell is the 2-atom "
                         "cubic bcc cell, so a per-cell value equals a "
                         "per-2-atom value",
    }
    report["eos_fit"] = fit_eos(V, E)

    # ---- per volume --------------------------------------------------
    # Gate thresholds. Chosen so that the healthy volumes of the two
    # ingested endpoints sit an order of magnitude inside them and the one
    # known-broken volume (BCC-Co-Ti QHA-08) sits far outside.
    FC_MAX_ABS_LIMIT = 50.0        # eV/A^2; healthy volumes span 6.1 - 10.1
    ASR_RELATIVE_LIMIT = 0.10      # healthy volumes span 0.008 - 0.024
    DFPT_VS_R3F_MEV_PER_ATOM_LIMIT = 100.0   # healthy span 4.6 - 28.7 meV/atom
    V_E_VS_R3F_MEV_LIMIT = 1.0     # healthy: exact to <1 ueV

    dfpt = rec / "dfpt"
    vol_dirs = []
    if dfpt.is_dir():
        vol_dirs = sorted((d for d in dfpt.iterdir()
                           if d.is_dir() and VOL_RE.match(d.name)),
                          key=lambda d: int(VOL_RE.match(d.name).group(1)))

    for idx, vd in enumerate(vol_dirs):
        tag = VOL_RE.match(vd.name).group(1)
        entry: dict[str, Any] = {"dir": str(vd.relative_to(rec))}
        vpath = vd / "vasprun.xml"
        opath = vd / "OUTCAR"
        if not (vpath.is_file() and opath.is_file()):
            entry["status"] = "missing_vasp_output"
            report["failures"].append(f"{vd.name}: missing vasprun.xml/OUTCAR")
            report["volumes"][tag] = entry
            continue

        vr = parse_vasprun(vpath)
        mk = outcar_markers(opath)
        gs = static_energy(opath, log_path=vd / "log-DFPT")
        chgcar = vd / "CHGCAR"
        chgcar_ok = chgcar.is_file() and chgcar.stat().st_size > 0
        entry["vasp"] = {
            "vasprun_parsed": vr.get("parsed"),
            "vasp_version": vr.get("vasp_version"),
            # The vasprun <calculation>/<energy>/e_0_energy tag is NOT the
            # ground-state energy of an IBRION=8 run; it is kept only so the
            # mismatch stays visible. The OUTCAR block below is authoritative.
            "vasprun_e_0_energy_eV_UNRELIABLE_FOR_IBRION8": vr.get("e_0_energy_eV"),
            "supercell_volume_ang3": vr.get("final_volume"),
            "ionic_steps": vr.get("ionic_steps"),
            "vasprun_sha256": sha256_file(vpath),
            "chgcar_present_for_icharg1": chgcar_ok,
            "dfpt_ground_state_energy_sigma0_eV": gs.get("energy_sigma0_eV"),
            "dfpt_ground_state_energy_per_atom_eV": (
                None if gs.get("energy_sigma0_eV") is None
                else gs["energy_sigma0_eV"] / 16.0),
            "dfpt_ground_state_scf_iterations": gs.get("ground_state_scf_iterations"),
            "dfpt_ground_state_scf_iteration_source": gs.get("scf_iteration_source"),
            "dfpt_ground_state_last_dE_eV": gs.get("ground_state_last_dE_eV"),
            "dfpt_ground_state_tail_energy_span_meV": gs.get(
                "ground_state_tail_energy_span_meV"),
            "nelm": gs.get("nelm"),
            "ediif_eV": gs.get("ediif_eV"),
            "dfpt_ground_state_hit_nelm": gs.get("ground_state_hit_nelm"),
            **mk,
        }

        # ---- the R3f static that produced this v-e.dat point -----------
        r3f_dir = vd / "R3f"
        r3f_out = r3f_dir / "OUTCAR"
        if r3f_out.is_file():
            r3f_txt = r3f_out.read_text(errors="replace")
            r3f_gs = static_energy(r3f_out, log_path=r3f_dir / "R3f", whole_log=True)
            entry["r3f_static"] = {
                "outcar": str(r3f_out.relative_to(rec)),
                "energy_sigma0_eV": r3f_gs.get("energy_sigma0_eV"),
                "energy_per_atom_eV": (None if r3f_gs.get("energy_sigma0_eV") is None
                                       else r3f_gs["energy_sigma0_eV"] / 2.0),
                "ground_state_scf_iterations": r3f_gs.get("ground_state_scf_iterations"),
                "ground_state_hit_nelm": r3f_gs.get("ground_state_hit_nelm"),
                "nelm": r3f_gs.get("nelm"),
                "ediif_eV": r3f_gs.get("ediif_eV"),
                "converged_scf_marker_count": r3f_txt.count(
                    "aborting loop because EDIFF is reached"),
                "unconverged_scf_marker_count": r3f_txt.count(
                    "aborting loop EDIFF was not reached (unconverged)"),
                "job_completion_marker": (
                    "General timing and accounting informations for this job" in r3f_txt),
                "incomplete_charge_density_warning": (
                    "chargedensity file is incomplete" in r3f_txt),
                "charge_density_of_overlapping_atoms": (
                    "charge density of overlapping atoms calculated" in r3f_txt),
            }
            # does this volume's v-e.dat row really come from this R3f static?
            if idx < len(E) and r3f_gs.get("energy_sigma0_eV") is not None:
                dev = (float(E[idx]) - r3f_gs["energy_sigma0_eV"]) * 1000.0
                entry["v_e_row_matches_r3f_static_meV"] = dev
                if abs(dev) > V_E_VS_R3F_MEV_LIMIT:
                    report["failures"].append(
                        f"{vd.name}: v-e.dat row {idx} does not match its "
                        f"R3f/OUTCAR energy(sigma->0) (off by {dev:.3f} meV)")
            if not entry["r3f_static"]["job_completion_marker"]:
                report["failures"].append(f"{vd.name}: R3f/OUTCAR has no completion marker")
            if entry["r3f_static"]["converged_scf_marker_count"] == 0:
                report["failures"].append(
                    f"{vd.name}: R3f static never reached EDIFF")
            if entry["r3f_static"]["ground_state_hit_nelm"]:
                report["failures"].append(
                    f"{vd.name}: R3f static SCF ran to NELM="
                    f"{entry['r3f_static']['nelm']} without converging")
            if entry["r3f_static"]["charge_density_of_overlapping_atoms"]:
                report["failures"].append(
                    f"{vd.name}: R3f static fell back to the atomic superposition")

        # ---- the DFPT ground state vs the R3f static -------------------
        # This is the gate that catches a diverged ground-state SCF. VASP
        # 5.4.4 prints NO "unconverged" marker when the loop simply runs out
        # of NELM steps, so the marker counts alone cannot see it.
        if gs.get("ground_state_hit_nelm"):
            span = gs.get("ground_state_tail_energy_span_meV")
            tail = "" if span is None else (
                f", last dE={gs.get('ground_state_last_dE_eV')} eV, tail of the "
                f"loop still moves by {span:.1f} meV")
            report["failures"].append(
                f"{vd.name}: DFPT ground-state SCF ran to NELM="
                f"{gs.get('nelm')} without reaching EDIFF{tail}")
        if (gs.get("energy_sigma0_eV") is not None
                and entry.get("r3f_static", {}).get("energy_per_atom_eV") is not None):
            d = gs["energy_sigma0_eV"] / 16.0 - entry["r3f_static"]["energy_per_atom_eV"]
            entry["dfpt_ground_state_vs_r3f_meV_per_atom"] = d * 1000.0
            if abs(d) * 1000.0 > DFPT_VS_R3F_MEV_PER_ATOM_LIMIT:
                report["failures"].append(
                    f"{vd.name}: DFPT ground state is {d * 1000.0:.1f} meV/atom "
                    f"away from the R3f static it is supposed to share a "
                    f"geometry and charge density with")

        # ---- the pilot's own marker gates ------------------------------
        if not mk["job_completion_marker"]:
            report["failures"].append(f"{vd.name}: OUTCAR has no completion marker")
        if mk["converged_scf_marker_count"] == 0:
            report["failures"].append(f"{vd.name}: electronic SCF never reached EDIFF")
        if mk["unconverged_scf_marker_count"] > 0:
            report["failures"].append(f"{vd.name}: OUTCAR reports an unconverged SCF loop")
        if not mk["dfpt_all_perturbations_solved"]:
            report["failures"].append(
                f"{vd.name}: DFPT did not solve every perturbation "
                f"({mk['dfpt_degrees_of_freedom']})")
        if not chgcar_ok:
            report["failures"].append(
                f"{vd.name}: ICHARG=1 but no CHGCAR is present")
        if mk["charge_density_of_overlapping_atoms"]:
            report["failures"].append(
                f"{vd.name}: VASP fell back to the atomic superposition "
                f"(ICHARG=1 found no usable CHGCAR)")

        try:
            entry["phonon"] = phonon_audit(vd, mesh)
        except Exception as ex:  # noqa: BLE001
            entry["phonon"] = {"error": f"{type(ex).__name__}: {ex}"}
            report["failures"].append(f"{vd.name}: phonon audit failed")

        ph = entry.get("phonon", {})
        if ph.get("force_constants_max_abs") is not None:
            if ph["force_constants_max_abs"] > FC_MAX_ABS_LIMIT:
                report["failures"].append(
                    f"{vd.name}: FORCE_CONSTANTS max |Phi| = "
                    f"{ph['force_constants_max_abs']:.2f} eV/A^2 exceeds "
                    f"{FC_MAX_ABS_LIMIT} -- the DFPT solution is not physical")
            if ph["asr_residual_relative"] > ASR_RELATIVE_LIMIT:
                report["failures"].append(
                    f"{vd.name}: FORCE_CONSTANTS violates the acoustic sum rule "
                    f"by {ph['asr_residual_relative'] * 100:.1f}%")
            if ph.get("imaginary_modes_weighted", 0) > 0:
                report["warnings"].append(
                    f"{vd.name}: {ph['imaginary_fraction'] * 100:.3f}% of the "
                    f"weighted mesh modes are imaginary (min "
                    f"{ph['min_frequency_THz']:.4f} THz at "
                    f"{ph['min_frequency_qpoint']}); phonopy drops them from "
                    f"the thermal properties, so the recorded ZPE/entropy at "
                    f"this volume covers only "
                    f"{ph['recomputed']['num_integrated_modes']}/"
                    f"{ph['recomputed']['num_modes']} modes")

        # recorded thermal properties, for cross-checking
        tps = sorted(vd.glob("thermal_properties-*.yaml"))
        if tps:
            import yaml
            d = yaml.safe_load(tps[0].read_text())
            entry["recorded_thermal_properties"] = {
                "file": str(tps[0].relative_to(rec)),
                "num_modes": d.get("num_modes"),
                "num_integrated_modes": d.get("num_integrated_modes"),
                "cutoff_frequency": d.get("cutoff_frequency"),
                "zero_point_energy_kJ_per_mol": d.get("zero_point_energy"),
                "sha256": sha256_file(tps[0]),
            }
            rec_tp = entry["recorded_thermal_properties"]
            new_tp = entry.get("phonon", {}).get("recomputed")
            if new_tp:
                rec_tp["recomputed_matches_recorded"] = {
                    "num_modes": new_tp["num_modes"] == rec_tp["num_modes"],
                    "num_integrated_modes": (
                        new_tp["num_integrated_modes"] == rec_tp["num_integrated_modes"]),
                    "zero_point_energy": (
                        abs(new_tp["zero_point_energy_kJ_per_mol"]
                            - rec_tp["zero_point_energy_kJ_per_mol"]) < 1e-6
                        if rec_tp["zero_point_energy_kJ_per_mol"] is not None else None),
                }
        report["volumes"][tag] = entry

    # ---- R3 anchor ---------------------------------------------------
    r3 = rec / "R3" / "CONTCAR"
    if r3.is_file():
        lines = r3.read_text().splitlines()
        try:
            a = float(lines[2].split()[0])
            report["anchor"] = {"R3_CONTCAR_a_ang": a, "R3_CONTCAR_V_ang3": a ** 3}
        except (IndexError, ValueError):
            report["anchor"] = {"error": "unparseable CONTCAR"}
    else:
        report["anchor"] = {"missing": True}
        report["notes"].append(
            "no R3/CONTCAR: this endpoint has no VASP layer at all, so the "
            "anchor cell and the DFPT force constants cannot be audited here"
        )

    # ---- the phonon layer as shipped, whatever layout it arrived in ----
    # Some endpoints ship the thermal-property yaml at the top level next to
    # the curves instead of inside per-volume DFPT directories. Record what is
    # there either way, so the delivery can state whether the phonon layer is
    # present, re-derivable, or absent.
    flat = sorted(rec.glob("thermal_properties-*.yaml"))
    if flat:
        import yaml
        rows = []
        for p in flat:
            d = yaml.safe_load(p.read_text())
            nm, nim = d.get("num_modes"), d.get("num_integrated_modes")
            rows.append({
                "file": p.name,
                "sha256": sha256_file(p),
                "natom": d.get("natom"),
                "num_modes": nm,
                "num_integrated_modes": nim,
                "modes_dropped": (None if nm is None or nim is None else nm - nim),
                "zero_point_energy_kJ_per_mol": d.get("zero_point_energy"),
                "high_T_entropy_J_per_K_per_mol": d.get("high_T_entropy"),
                "temperature_grid_K": (
                    None if not d.get("thermal_properties")
                    else [d["thermal_properties"][0]["temperature"],
                          d["thermal_properties"][-1]["temperature"],
                          d["thermal_properties"][1]["temperature"]
                          - d["thermal_properties"][0]["temperature"]]),
            })
        report["top_level_thermal_properties"] = rows
        unstable = [r for r in rows if r["modes_dropped"]]
        if unstable:
            report["warnings"].append(
                f"{len(unstable)}/{len(rows)} top-level thermal_properties yaml "
                f"files record fewer integrated modes than total modes, i.e. "
                f"phonopy silently dropped imaginary modes there "
                f"({', '.join(r['file'] for r in unstable)})")
        report["notes"].append(
            "these thermal_properties yaml are shipped next to the curves, not "
            "inside per-volume DFPT directories; without the corresponding "
            "FORCE_CONSTANTS they cannot be recomputed here"
        )
    elif not vol_dirs:
        report["top_level_thermal_properties"] = []
        report["notes"].append(
            "no thermal_properties yaml anywhere under legacy-as-received"
        )

    # ---- reported curve ----------------------------------------------
    curve_candidates = [
        rec / "gibbs-temperature.dat",
        rec / "QHA-TiCo" / "gibbs-temperature.dat",
    ]
    for cp in curve_candidates:
        if cp.is_file():
            T, G = read_curve(cp)
            dG = np.diff(G)
            i = int(np.argmax(np.abs(dG)))
            report["reported_curve"] = {
                "file": str(cp.relative_to(rec)),
                "sha256": sha256_file(cp),
                "n_points": int(len(T)),
                "temperature_K": [float(T[0]), float(T[-1]), float(T[1] - T[0])],
                "gibbs_first_last_eV": [float(G[0]), float(G[-1])],
                "max_adjacent_jump_eV": float(abs(dG[i])),
                "max_adjacent_jump_at_K": [float(T[i]), float(T[i + 1])],
                "n_jumps_over_0_2_eV": int((np.abs(dG) > 0.2).sum()),
            }
            break
    else:
        report["reported_curve"] = {"missing": True}

    # A dataset with no VASP layer cannot be verified at all: reporting "pass"
    # for one would be exactly the fail-open behaviour this tool exists to
    # prevent. Such an endpoint gets its own status.
    vasp_layer = (rec / "R3").is_dir() or bool(vol_dirs)
    report["vasp_layer_present"] = vasp_layer
    if report["failures"]:
        report["status"] = "fail"
    elif not vasp_layer:
        report["status"] = "unverifiable"
        report["notes"].append(
            "no VASP layer under legacy-as-received (no R3/ and no dfpt/): the "
            "curve and the thermal properties can be read, but nothing about "
            "how they were produced can be checked from this endpoint alone"
        )
    else:
        report["status"] = "pass"
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "audit.json").write_text(
        json.dumps(report, indent=2, default=str) + "\n", encoding="utf-8")
    return report


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True, type=Path)
    ap.add_argument("--mesh", type=int, default=24)
    ap.add_argument("--outdir", type=Path, default=None)
    args = ap.parse_args()
    outdir = args.outdir or (args.root / "audit")
    rep = audit_endpoint(args.root, args.mesh, outdir)
    print(f"status: {rep['status']}")
    for f in rep["failures"]:
        print(f"  FAIL {f}")
    for w in rep["warnings"]:
        print(f"  WARN {w}")
    print(f"audit -> {outdir/'audit.json'}")


if __name__ == "__main__":
    main()
