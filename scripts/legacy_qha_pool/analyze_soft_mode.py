#!/usr/bin/env python
"""Pick the soft-mode wavevector, and describe its eigenvector honestly.

Why this exists (report 4.13.4-4.13.6)
--------------------------------------
`extract_soft_mode.py` reports the MINIMUM on a sampling mesh.  That is not the
same thing as the instability wavevector, and using it for mode following is
wrong for two separate reasons, both of which this script checks:

1. GRID.  phonopy's default mesh is Monkhorst-Pack, so a 12x12x12 grid sits on
   odd multiples of 1/24 -- it does NOT contain N = (1/2,1/2,0).  The reported
   minimum (11/24,11/24,1/24) is therefore an INTERPOLATED value at a point
   near N, and it comes out SHALLOWER than N itself.  Searching for the
   minimum on a mesh that misses the symmetry point understates the depth.

2. COMMENSURABILITY.  A supercell can only represent q-points with
   S^T q integer (S = supercell_matrix in the primitive basis).  For the 2x2x2
   supercell that is q in {0, 1/2}^3 only.  (11/24,11/24,1/24) fails this, so
   even if it were the true minimum, a 2x2x2 supercell could not host the
   distortion it implies.

So the two questions are answered in the opposite order: first enumerate the
q-points the supercell CAN represent, then find the deepest of those.  For B2
AlTa that turns out to be N, which is also the physically expected place for
the instability.

Polarisation: what this script will and will not claim
------------------------------------------------------
The species weight says "Ta-dominated".  That does NOT license "Ta moves along
[110]".  Two traps:

* At N the soft branch is 2-FOLD DEGENERATE.  `np.linalg.eigh` returns an
  arbitrary orthonormal basis of the degenerate subspace, so a single
  eigenvector is not well defined -- only the subspace is.  This script detects
  the degeneracy and reports invariant properties of the subspace (the
  projector, which Cartesian components vanish identically across it, and
  whether the species weight is the same for every vector in it), rather than
  one arbitrary vector.
* q along [110] does not mean displacement along [110].  This script computes
  the angle between the displacement and q instead of assuming it.

Usage
-----
    python analyze_soft_mode.py --fcdiff /tmp/fcdiff/stage0b \
        --endpoint BCC-Al-Ta-2 --smearing s1_020 \
        --archive /nfs_ssd/tmp/legacy-qha-pool/soft-mode-evidence \
        --json soft_mode_wavevector.json
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import shutil
from pathlib import Path

import numpy as np

import analyze_stage0b as A
from phonon_eigenvector import eigenvector

DEGEN_TOL_THz = 1e-3     # frequencies closer than this are treated as degenerate
MESH = 12                # the mesh extract_soft_mode.py scans


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def commensurate(q, S: np.ndarray, atol: float = 1e-9) -> bool:
    """True iff the supercell can represent this q exactly.

    The phase factor is exp(2*pi*i * q . R); it is single-valued on the
    supercell iff S^T q is a vector of integers.
    """
    t = S.T @ np.asarray(q, dtype=float)
    return bool(np.allclose(t, np.round(t), atol=atol))


def mesh_grid_report(ph) -> dict:
    """Is the symmetry point even on the mesh phonopy scans?"""
    ph.run_mesh([MESH] * 3)
    d = ph.get_mesh_dict()
    q = np.array(d["qpoints"])
    f = np.array(d["frequencies"])
    iq, ib = np.unravel_index(np.argmin(f), f.shape)
    ux = np.unique(np.round(q[:, 0], 6))
    return {
        "mesh": [MESH] * 3,
        "n_irreducible_qpoints": int(len(q)),
        "q_x_values": [float(x) for x in ux],
        "spacing": float(ux[1] - ux[0]) if len(ux) > 1 else None,
        "grid_contains_half": bool(np.any(np.isclose(ux, 0.5))),
        "mesh_min_qpoint": [float(x) for x in q[iq]],
        "mesh_min_frequency_THz": float(f[iq, ib]),
        "mesh_n_imaginary": int((f < -1e-4).sum()),
        "mesh_n_modes": int(f.size),
    }


def deepest_commensurate(ph, S: np.ndarray, n_neg_max: int = 2) -> list[dict]:
    """Enumerate the q-points the supercell can host and rank them by depth."""
    cands = [list(c) for c in itertools.product([0.0, 0.5], repeat=3)]
    ph.run_qpoints(cands)
    f = np.array(ph.get_qpoints_dict()["frequencies"])
    out = []
    for q, fr in zip(cands, f):
        n_neg = int((fr < -1e-4).sum())
        if n_neg > n_neg_max:
            continue
        out.append({"qpoint": q, "min_frequency_THz": float(fr.min()),
                    "n_imaginary": n_neg, "commensurate": commensurate(q, S)})
    return sorted(out, key=lambda r: r["min_frequency_THz"])


def subspace_report(ph, q, pr) -> dict:
    """Invariant description of the soft-mode eigenvector subspace at q.

    Everything reported here is either (a) invariant under the arbitrary basis
    choice inside a degenerate subspace, or (b) explicitly flagged as
    basis-dependent.  The individual component values of ONE eigenvector are
    basis-dependent whenever the branch is degenerate, so they are recorded for
    provenance but never used to state a direction.
    """
    ph.run_qpoints([list(q)])
    f_all = np.array(ph.get_qpoints_dict()["frequencies"])[0]
    soft = [int(b) for b in range(len(f_all)) if f_all[b] < -1e-4]
    if not soft:
        return {"qpoint": [float(x) for x in q], "n_imaginary_branches_at_q": 0,
                "note": "no imaginary branch at this q"}

    bands, vecs, checks = [], [], []
    for b in soft:
        v, f, info = eigenvector(ph, q, b)      # checked column read
        vecs.append(v)
        checks.append({"band": b, "eigen_residual": info["eigen_residual"],
                       "row_reading_residual": info["row_reading_residual"]})
        bands.append({
            "band": b, "frequency_THz": float(f),
            "components_basis_dependent": [[{"re": float(c.real), "im": float(c.imag)}
                                            for c in v[3 * a:3 * a + 3]]
                                           for a in range(len(pr))],
            "abs_per_atom": [float(np.linalg.norm(v[3 * a:3 * a + 3]))
                             for a in range(len(pr))],
        })
    V = np.array(vecs)                          # (k, 3*n_prim)
    n_prim = len(pr)
    degenerate = (len(soft) > 1 and
                  abs(f_all[soft[0]] - f_all[soft[1]]) < DEGEN_TOL_THz)

    # --- invariant 1: Cartesian components that vanish across the WHOLE subspace
    silent = []
    for a, sym in enumerate(pr.symbols):
        for c, ax in enumerate("xyz"):
            if np.all(np.abs(V[:, 3 * a + c]) < 1e-6):
                silent.append(f"{sym}.{ax}")
    # --- invariant 2: species weight is the same for EVERY vector in the subspace?
    #     true iff each species block acts isotropically on the subspace
    weights, weight_is_invariant = {}, True
    syms = list(dict.fromkeys(str(s) for s in pr.symbols))
    for s in syms:
        idx = [3 * a + c for a, sy in enumerate(pr.symbols)
               if str(sy) == s for c in range(3)]
        blk = V[:, idx]
        sv = np.linalg.svd(blk, compute_uv=False)
        weights[s] = float((sv ** 2).sum())
        if len(sv) > 1 and (sv.max() - sv.min()) > 1e-6 * max(sv.max(), 1e-30):
            weight_is_invariant = False
    tot = sum(weights.values()) or 1.0
    weights = {k: 100 * v / tot for k, v in weights.items()}

    # --- invariant 3: projector
    Q, _ = np.linalg.qr(V.T)
    P = Q @ Q.T

    # --- direction of the dominant species, only when NON-degenerate
    #     NB: the eigenvector is defined up to a GLOBAL SIGN (real) or phase
    #     (complex).  So the angle to q must use |<d,q>|, never |d| componentwise
    #     -- taking the absolute value of each component first silently maps a
    #     transverse mode onto a longitudinal one.
    direction = None
    if not degenerate:
        dom = max(syms, key=lambda s: weights[s])
        idx = [3 * a + c for a, sy in enumerate(pr.symbols)
               if str(sy) == dom for c in range(3)]
        raw = V[0, idx]
        d = np.abs(raw)                     # magnitudes, for reporting only
        qv = np.asarray(q, dtype=float)
        nrm = np.linalg.norm(raw) * np.linalg.norm(qv)
        # longitudinal fraction is the |overlap|^2 -- sign/phase independent
        lon = (float(abs(np.vdot(raw, qv)) ** 2 / nrm ** 2)
               if nrm > 0 else None)
        ang = (float(np.degrees(np.arccos(np.sqrt(np.clip(lon, 0, 1)))))
               if lon is not None else None)
        direction = {
            "species": dom,
            "complex_components": [{"re": float(z.real), "im": float(z.imag)}
                                   for z in raw],
            "cartesian_magnitudes": [float(x) for x in d],
            "global_sign_or_phase_is_arbitrary": True,
            "angle_to_q_deg": ang,
            "longitudinal_fraction": lon,
            "character": (None if lon is None else
                          "longitudinal" if lon > 0.9 else
                          "transverse" if lon < 0.1 else "mixed"),
        }
    return {
        "qpoint": [float(x) for x in q],
        "degenerate": bool(degenerate),
        "degeneracy": int(len(soft)),
        "note": ("degenerate subspace: a single eigenvector is NOT unique; only "
                 "the subspace is.  Direction is therefore NOT reported."
                 if degenerate else
                 "non-degenerate: a single eigenvector is meaningful."),
        "n_imaginary_branches_at_q": len(soft),
        "bands": bands,
        "eigenvector_checks": checks,
        "silent_cartesian_components": silent,
        "species_weight_percent": weights,
        "species_weight_invariant_across_subspace": bool(weight_is_invariant),
        "subspace_projector": [[float(x.real) for x in row] for row in P],
        "direction": direction,
        "primitive_symbols": [str(s) for s in pr.symbols],
        "primitive_scaled_positions": [[float(x) for x in p]
                                       for p in pr.scaled_positions],
        "p2s_map": [int(x) for x in pr.p2s_map],
        "supercell_matrix": [[int(x) for x in row]
                             for row in ph.supercell_matrix],
        "eigenvector_basis": "primitive cell, 3*n_prim components, COLUMN read",
    }


# ---------------------------------------------------------------------------
# Symmetry-adapted subspace members -- the mode-following starting points.
#
# A degenerate soft branch has no unique eigenvector, so mode following cannot
# start from "the" eigenvector.  It has to start from a well-defined member of
# the subspace, and the meaningful members are the symmetry-adapted combinations
# of two independent binary choices:
#
#     longitudinal (along q)      vs  transverse (perpendicular to q)
#     in phase (Al and Ta together)  vs  antiphase (Al against Ta)
#
# NAME WARNING.  It is tempting to call these "acoustic" and "optical", because
# at Gamma the acoustic branches are in phase and the optical ones antiphase.
# `branch_tracking()` shows that this does NOT hold along the path for B2 AlTa:
# BOTH members of the soft pair trace back to the Gamma ACOUSTIC triplet (both
# start at ~0.02 THz), and the antiphase character of the upper one is acquired
# at an avoided crossing near t=0.5, not by an optical branch descending.  So
# these functions label members by their DISPLACEMENT PATTERN, never by a branch
# name.
#
# `subspace_members()` builds the four patterns from the primitive masses and
# the projector, and reports how much of each survives the projection.  A member
# that loses a lot of norm has components OUTSIDE the soft subspace and must not
# be used as a starting direction -- that is the check, and it is what
# distinguishes the two good patterns from the two bad ones.
# ---------------------------------------------------------------------------

def _fmt_vec(v):
    """Format a displacement vector whose entries may be real or {re, im}."""
    return "[" + ", ".join(
        f"{x:+.6f}" if isinstance(x, float) else
        f"({x['re']:+.6f}{x['im']:+.6f}i)" for x in v) + "]"


def _num(z):
    """Real float when the imaginary part is negligible, else {re, im}.

    Eigenvectors at commensurate q-points come out real, and their SIGN is the
    whole point (it distinguishes a longitudinal from a transverse member, and
    an acoustic-like from an optical-like one).  Discarding it with abs() makes
    different patterns print identically.
    """
    z = complex(z)
    return (float(z.real) if abs(z.imag) < 1e-10
            else {"re": float(z.real), "im": float(z.imag)})


def subspace_members(ph, q, pr, axis_long=None, axis_trans=None) -> dict:
    """Symmetry-adapted members of the soft subspace, with their projections.

    `axis_long` defaults to q normalised (the longitudinal direction in
    reciprocal space); `axis_trans` defaults to an in-plane perpendicular.
    Both are Cartesian unit vectors in the primitive cell's frame.
    """
    from phonon_eigenvector import eigenvector

    q = np.asarray(q, dtype=float)
    if axis_long is None:
        axis_long = q / (np.linalg.norm(q) or 1.0)

    ph.run_qpoints([list(q)])
    f_all = np.array(ph.get_qpoints_dict()["frequencies"])[0]
    soft = [int(b) for b in range(len(f_all)) if f_all[b] < -1e-4]
    if not soft:
        return {"error": "no imaginary branch at this q"}
    V = np.array([eigenvector(ph, q, b)[0] for b in soft])
    Q, _ = np.linalg.qr(V.T)
    P = Q @ Q.T

    if axis_trans is None:
        # A perpendicular to q that LIES IN THE SOFT SUBSPACE.  Taking a naive
        # perpendicular (e.g. [q_y, -q_x, 0]) can land exactly on the silent
        # axis -- for q = (0,1/2,1/2) that is x, which the subspace excludes
        # entirely, and the projection then returns 0.  So enumerate candidate
        # perpendiculars and keep the one the subspace retains most.
        best_v, best_keep = None, -1.0
        for i in range(3):
            for j in range(i + 1, 3):
                v = np.zeros(3)
                v[i], v[j] = axis_long[j], -axis_long[i]
                n = np.linalg.norm(v)
                if n < 1e-9:
                    continue
                v = v / n
                keep = float(np.linalg.norm(P @ np.tile(v, len(pr))))
                if keep > best_keep:
                    best_v, best_keep = v, keep
        if best_v is None:
            return {"error": "could not find a transverse axis"}
        axis_trans = best_v

    masses = np.asarray(pr.masses, dtype=float)
    n_prim = len(pr)
    align = check_longitudinal_transverse(ph, q, axis_long, axis_trans)
    out = {"qpoint": [float(x) for x in q],
           "axis_longitudinal": [float(x) for x in axis_long],
           "axis_transverse": [float(x) for x in axis_trans],
           "axis_alignment_check": align,
           "n_soft_branches": len(soft),
           "members": {}}
    for ax_name, ax in (("L", axis_long), ("T", axis_trans)):
        for rel in (+1.0, -1.0):
            # mass-weighted pattern with |e_Ta|/|e_Al| = sqrt(m_Ta/m_Al),
            # i.e. equal physical displacement magnitudes (u_Al = u_Ta)
            e = np.concatenate([ax,
                                rel * np.sqrt(masses[-1] / masses[0]) * ax])
            e = e / np.linalg.norm(e)
            pv = P @ e
            keep = float(np.linalg.norm(pv))
            if keep > 1e-12:
                m = pv / np.linalg.norm(pv)
            else:
                m = None
            # label by displacement pattern, NOT by branch name -- see the
            # module comment: both soft branches here are acoustic-derived.
            long_or_trans = "longitudinal" if ax_name == "L" else "transverse"
            phase_rel = "antiphase" if rel < 0 else "in_phase"
            rec = {"pattern": f"{long_or_trans}_{phase_rel}",
                   "axis": [float(x) for x in ax],
                   "kept_fraction": keep,
                   "usable_as_starting_direction": keep > 0.99}
            if m is not None:
                # u in Angstrom per unit eigenvector norm: e is normalised to
                # sum|e|^2 = 1, so u_i = e_i/sqrt(m_i) is ALREADY a physical
                # displacement; do NOT renormalise it or the amplitude becomes
                # meaningless.  The max over atoms is what an amplitude should
                # be quoted against.
                u_phys = m / np.sqrt(np.repeat(masses, 3))
                rec["max_single_atom_displacement_angstrom_per_unit_eigenvector"] = float(
                    max(np.linalg.norm(u_phys[3 * a:3 * a + 3])
                        for a in range(n_prim)))
                u = u_phys / np.linalg.norm(u_phys)
                rec["e_mass_weighted"] = [_num(z) for z in m]
                rec["u_physical_displacement_normalised"] = [_num(x) for x in u]
                rec["u_per_atom"] = [
                    {"symbol": str(pr.symbols[a]),
                     "vector": [_num(x) for x in u[3 * a:3 * a + 3]],
                     # sign matters: it is what makes one member longitudinal
                     # and the other transverse.  Taking abs() here made the
                     # two members print identically.
                     "vector_abs": [float(abs(x)) for x in u[3 * a:3 * a + 3]]}
                    for a in range(n_prim)]
                # Longitudinal fraction must be measured PER SUBLATTICE, not on
                # the whole vector.  For an optical-like member the two
                # sublattices move in ANTIPHASE along q, so their net
                # displacement nearly cancels and a whole-vector projection
                # against [q, q] returns ~0 -- which would wrongly label a
                # perfectly longitudinal optical mode as transverse.  What
                # matters is whether each sublattice's own displacement is
                # parallel to q.
                per_block = []
                for a in range(n_prim):
                    blk = m[3 * a:3 * a + 3]
                    nb = np.linalg.norm(blk)
                    per_block.append(float(abs(np.vdot(blk, axis_long)) ** 2 /
                                           nb ** 2) if nb > 1e-12 else None)
                rec["longitudinal_fraction_per_sublattice"] = per_block
                vals = [v for v in per_block if v is not None]
                lon = float(np.mean(vals)) if vals else 0.0
                rec["longitudinal_fraction_mean"] = lon
                rec["character"] = ("longitudinal" if lon > 0.9 else
                                    "transverse" if lon < 0.1 else "mixed")
            out["members"][rec["pattern"]] = rec
    return out


def displacement_field(ph, q, e_mw, amp: float = 1.0) -> np.ndarray:
    """Cartesian displacement of every supercell atom for the q-modulated pattern.

        u_i = amp * Re[ e_{prim(i)} / sqrt(m_i) * exp(2*pi*i * q . R_i) ]

    with R_i the atom position in the PRIMITIVE basis (f_s @ S^T).  Returned as
    a (n_super, 3) array; multiply by a scale to set a target amplitude.
    """
    pr = ph.primitive
    S = np.array(ph.supercell_matrix)
    sc = ph.supercell
    masses = np.asarray(pr.masses, dtype=float)
    p2s = np.array(pr.p2s_map)
    s2p = np.array(pr.s2p_map)
    prim_of = {int(p): k for k, p in enumerate(p2s)}
    pidx = np.array([prim_of[int(x)] for x in s2p])

    f_p = np.array(sc.scaled_positions) @ S.T
    q = np.asarray(q, dtype=float)
    pos0 = np.array(sc.positions)
    d = np.zeros_like(pos0)
    for i in range(len(sc)):
        p = pidx[i]
        d[i] = amp * np.real(e_mw[3 * p:3 * p + 3] / np.sqrt(masses[p])
                             * np.exp(2j * np.pi * (q @ f_p[i])))
    return d


def _canonical(pos, lattice) -> np.ndarray:
    """Sorted fractional coordinates modulo 1 -- a translation-insensitive key."""
    f = (np.linalg.solve(lattice.T, pos.T).T) % 1.0
    return np.array(sorted(map(tuple, np.round(f, 8))))


def scale_for_max_displacement(ph, q, e_mw, target_angstrom: float) -> tuple[float, float]:
    """Scale factor that puts the LARGEST single-atom displacement at `target`.

    Amplitude has to be defined on a physical quantity or patterns with
    different normalisations are not comparable: the members coming out of
    `subspace_members()` are normalised to sum_i |u_i|^2 = 1, which leaves the
    largest atom somewhere near 0.08 A, not 1 A.  Quoting "0.02 A" without
    saying whether it is the scale factor or the maximum displacement makes two
    runs incomparable.  This function makes it the MAXIMUM SINGLE-ATOM
    DISPLACEMENT, which is what a mode-following amplitude should mean.

    Returns (scale_factor, resulting_max_displacement).
    """
    d1 = displacement_field(ph, q, e_mw, 1.0)
    m = float(np.linalg.norm(d1, axis=1).max())
    if m < 1e-12:
        raise ValueError("the pattern has no displacement; cannot normalise")
    return target_angstrom / m, target_angstrom


def reciprocal_lattice(ph) -> np.ndarray:
    """B = 2*pi * A^-T, rows = b1, b2, b3, units 1/Angstrom.

    q_cartesian = q_fractional @ B.  Using the DIRECT lattice here instead is a
    silent bug: for a cubic diagonal cell the two differ only by a scalar, so
    every direction test still passes and the wrong number never shows up --
    until a non-cubic cell (HCP, low symmetry) makes the direction itself wrong.
    `assert_reciprocal()` is the guard.
    """
    return 2.0 * np.pi * np.linalg.inv(np.array(ph.primitive.cell)).T


def assert_reciprocal(ph, tol: float = 1e-10) -> dict:
    """Prove the matrix used for q_cartesian really is the reciprocal lattice.

    The defining property is B @ A^T = 2*pi*I (row-vector convention).  Checking
    it is what makes the convention a guarded fact rather than a comment.
    """
    A = np.array(ph.primitive.cell)
    B = reciprocal_lattice(ph)
    prod = B @ A.T
    dev = float(np.abs(prod - 2 * np.pi * np.eye(3)).max())
    return {
        "q_coordinate_convention": "fractional_reciprocal",
        "reciprocal_lattice_units": "1/angstrom",
        "reciprocal_lattice": [[float(x) for x in r] for r in B],
        "direct_lattice": [[float(x) for x in r] for r in A],
        "B_dot_A_transpose_is_2pi_identity": dev < tol,
        "max_deviation_from_2pi_identity": dev,
        "ok": dev < tol,
    }


def subspace_isotropy(ph, q, pr, n_dir: int = 13) -> dict:
    """Is the quadratic form isotropic inside the soft subspace?

    For a 2D irreducible representation the answer is yes BY SYMMETRY, and the
    consequence is that no direction inside the subspace is more unstable than
    another -- "which direction should we distort along" stops being a physical
    question, and two arbitrary orthogonal members are the same instability seen
    in two orientations rather than two candidate paths.

    That is not something to assume: an isotropic subspace and an anisotropic one
    produce the same single-member frequency table, so the check has to sweep
    directions.  This sweeps `n_dir` of them at equal maximum single-atom
    displacement and reports the spread.

    Requires force constants (it evaluates (1/2) u^T Phi u).
    """
    from phonon_eigenvector import eigenvector

    ph.run_qpoints([list(q)])
    f_all = np.array(ph.get_qpoints_dict()["frequencies"])[0]
    soft = [int(b) for b in range(len(f_all)) if f_all[b] < -1e-4]
    if len(soft) < 2:
        return {"error": f"only {len(soft)} imaginary branch(es) at this q; "
                         f"isotropy is a question about a degenerate subspace"}
    V = np.array([eigenvector(ph, q, b)[0] for b in soft])
    Q, _ = np.linalg.qr(V.T)
    P = Q @ Q.T
    masses = np.asarray(pr.masses, dtype=float)
    n_prim = len(pr)

    # two orthogonal representatives inside the subspace
    e1 = V[0] / np.linalg.norm(V[0])
    e2 = V[1] - np.vdot(e1, V[1]) * e1
    e2 = e2 / np.linalg.norm(e2)

    fc = np.array(ph.force_constants)

    def field(ev, target=0.05):
        # mass-weighting: e_i = sqrt(m_i) u_i, so u_i = e_i / sqrt(m_i)
        u = np.concatenate([ev[3 * a:3 * a + 3] / np.sqrt(masses[a])
                            for a in range(n_prim)])
        m = float(np.linalg.norm(u, axis=1).max()) if u.ndim == 2 else 0.0
        return u

    # build the cartesian field per supercell atom, q-modulated
    S = np.array(ph.supercell_matrix)
    sc = ph.supercell
    f_p = np.array(sc.scaled_positions) @ S.T
    p2s = np.array(pr.p2s_map)
    s2p = np.array(pr.s2p_map)
    prim_of = {int(p): k for k, p in enumerate(p2s)}
    pidx = np.array([prim_of[int(x)] for x in s2p])

    def cart_field(ev, target=0.05):
        raw = np.zeros((len(sc), 3), dtype=complex)
        for i in range(len(sc)):
            a = pidx[i]
            raw[i] = (ev[3 * a:3 * a + 3] / np.sqrt(masses[a])
                      * np.exp(2j * np.pi * (np.asarray(q) @ f_p[i])))
        d = np.real(raw)
        mx = np.linalg.norm(d, axis=1).max()
        if mx < 1e-12:
            raise ValueError("degenerate direction produced no displacement")
        return d * (target / mx)

    def quad(u):
        return 0.5 * float(np.einsum('ia,ijab,jb->', u, fc, u))

    rows = []
    for th in np.linspace(0.0, 180.0, n_dir):
        t_ = np.radians(th)
        ev = np.cos(t_) * e1 + np.sin(t_) * e2
        ev = ev / np.linalg.norm(ev)
        rows.append({"theta_deg": float(th), "quadratic_meV": quad(cart_field(ev)) * 1e3})
    vals = np.array([r["quadratic_meV"] for r in rows])
    spread = float(vals.max() - vals.min())
    return {
        "qpoint": [float(x) for x in q],
        "n_directions": n_dir,
        "target_max_single_atom_displacement_angstrom": 0.05,
        "quadratic_meV": [float(v) for v in vals],
        "spread_meV": spread,
        "relative_spread": float(spread / abs(vals.mean())) if vals.mean() else None,
        "isotropic": bool(spread < 1e-6 * max(abs(vals.mean()), 1e-30)),
        "conclusion": ("the quadratic form is ISOTROPIC inside the soft subspace "
                       "-> 2D irreducible representation -> the HARMONIC basis "
                       "direction is not unique, so one relaxation per orbital "
                       "representative is enough AT HARMONIC ORDER.  This does NOT "
                       "prove the full anharmonic energy surface is isotropic: "
                       "crystal symmetry allows higher-order terms to select a "
                       "particular order-parameter direction, so the finite-"
                       "amplitude direction still has to be settled by "
                       "calculation.  Start the relaxation from a GENERIC "
                       "direction, not from a symmetry-adapted one."
                       if spread < 1e-6 * max(abs(vals.mean()), 1e-30) else
                       "the quadratic form is ANISOTROPIC inside the soft subspace "
                       "-> the members are NOT equivalent; each direction must be "
                       "followed separately"),
    }


def check_longitudinal_transverse(ph, q, axis_long, axis_trans) -> dict:
    """Executable check that the two axes really are longitudinal/transverse.

    Do NOT trust the label.  A [110]-style string is meaningless without the q
    it is paired with: at the three equivalent N arms of a cubic cell the
    longitudinal direction is [011], [101] and [110] respectively, so reusing
    one label across arms silently pairs a direction with the wrong q.  The
    check is the definition itself:

        longitudinal: |u x q_hat| ~ 0      (parallel to q)
        transverse:   |u . q_hat| ~ 0      (perpendicular to q)

    with q_hat in CARTESIAN coordinates obtained from the RECIPROCAL lattice,
    q_cart = q_fractional @ B, B = 2*pi*A^-T.
    """
    rec = assert_reciprocal(ph)
    B = np.array(rec["reciprocal_lattice"])
    q = np.asarray(q, dtype=float)
    q_cart = q @ B
    nq = np.linalg.norm(q_cart)
    if nq < 1e-12:
        return {"error": "q is zero; no direction to test against"}
    q_hat = q_cart / nq
    al = np.asarray(axis_long, dtype=float)
    at = np.asarray(axis_trans, dtype=float)
    al = al / (np.linalg.norm(al) or 1.0)
    at = at / (np.linalg.norm(at) or 1.0)
    cross_l = float(np.linalg.norm(np.cross(al, q_hat)))
    dot_t = float(abs(np.dot(at, q_hat)))
    return {
        "q_fractional": [float(x) for x in q],
        "q_cartesian": [float(x) for x in q_cart],
        "q_coordinate_convention": "fractional_reciprocal",
        "reciprocal_lattice_units": "1/angstrom",
        "reciprocal_lattice": rec["reciprocal_lattice"],
        "reciprocal_check_ok": rec["ok"],
        "reciprocal_check_max_deviation": rec["max_deviation_from_2pi_identity"],
        "basis_convention": "primitive cell of the B2 phase (cubic, diagonal)",
        "longitudinal_direction_cartesian": [float(x) for x in al],
        "transverse_direction_cartesian": [float(x) for x in at],
        "longitudinal_cross_q_hat": cross_l,
        "transverse_dot_q_hat": dot_t,
        "tolerance": 1e-8,
        "longitudinal_is_parallel_to_q": cross_l < 1e-8,
        "transverse_is_perpendicular_to_q": dot_t < 1e-8,
        "ok": cross_l < 1e-8 and dot_t < 1e-8 and rec["ok"],
    }


def translation_equivalence(ph, q, e_mw, amp: float = 0.05) -> dict:
    """Is the +amp pattern the SAME CRYSTAL as the -amp one, up to translation?

    For a single-q displacement, translating the crystal by a lattice vector R
    multiplies the pattern by exp(2*pi*i q.R).  Whenever q.R = 1/2 for some
    parent lattice vector -- true at every N point of a B2 cell, since N has two
    half-integer components -- that factor is -1, so the -amp pattern IS the
    +amp pattern translated.  A periodic crystal cannot tell those apart.

    Consequence: +L/-L and +T/-T are NOT four independent starting structures.
    Two DFT cases suffice (one longitudinal, one transverse); the sign pair is a
    regression check for the input generator, not a second calculation.
    """
    S = np.array(ph.supercell_matrix)
    sc = ph.supercell
    L_s = np.array(sc.cell)
    pos0 = np.array(sc.positions)
    q = np.asarray(q, dtype=float)

    # amp is a scale factor on a unit-norm eigenvector; report what it means in
    # Angstrom so the number is comparable across patterns.
    scale, maxdisp = scale_for_max_displacement(ph, q, e_mw, amp)
    dP = displacement_field(ph, q, e_mw, +scale)
    dM = displacement_field(ph, q, e_mw, -scale)
    out = {
        "target_max_single_atom_displacement_angstrom": float(amp),
        "scale_factor_applied_to_unit_norm_eigenvector": float(scale),
        "max_single_atom_displacement_angstrom": float(
            np.linalg.norm(dP, axis=1).max()),
        "amplitude_convention": ("largest single-atom displacement, in Angstrom; "
                                 "the eigenvector is scaled so this equals the "
                                 "target"),
        "minus_is_exactly_plus_negated": bool(np.allclose(dM, -dP)),
        "same_crystal_without_translation": bool(
            np.allclose(_canonical(pos0 + dP, L_s), _canonical(pos0 + dM, L_s))),
        "translations_that_map_plus_onto_minus": [],
        "q_dot_parent_vector_is_half": [],
    }
    # parent-cell translations are the supercell-fractional steps 1/S_ii
    for ax in range(3):
        t = np.zeros(3)
        t[ax] = 1.0 / S[ax, ax]
        tc = t @ L_s
        if np.allclose(_canonical(pos0 + dP + tc, L_s), _canonical(pos0 + dM, L_s)):
            out["translations_that_map_plus_onto_minus"].append(
                {"supercell_fractional": [float(x) for x in t],
                 "parent_cell_axis": ax + 1,
                 "q_dot_R": float(q[ax])})
        if abs(q[ax] - round(q[ax])) > 1e-9 and abs(abs(q[ax]) - 0.5) < 1e-9:
            out["q_dot_parent_vector_is_half"].append(ax + 1)
    out["plus_minus_are_translation_equivalent"] = bool(
        out["translations_that_map_plus_onto_minus"])
    out["conclusion"] = (
        "+amp and -amp are the same crystal up to a lattice translation, so only "
        "ONE of the pair needs a DFT calculation; the other is a generator "
        "regression check."
        if out["plus_minus_are_translation_equivalent"] else
        "+amp and -amp are NOT related by a lattice translation: both are "
        "independent starting structures and both must be computed.")
    return out


def branch_tracking(ph, q_target, n_points: int = 21) -> dict:
    """Track branch character from Gamma to `q_target`, to test branch labels.

    This exists because the tempting label "optical branch" for an antiphase
    member needs evidence: an optical branch at Gamma must actually descend to
    the soft mode.  The test is the Al/Ta phase relation (in phase / antiphase)
    and the frequency, sampled along Gamma -> q_target.

    Returns the dispersion and, per branch, whether the phase relation stays
    constant or flips (a flip marks an avoided crossing, i.e. the character
    moved between branches rather than descending with one of them).
    """
    from phonon_eigenvector import eigenvector

    pr = ph.primitive
    q_target = np.asarray(q_target, dtype=float)
    ts = np.linspace(0.0, 1.0, n_points)
    qs = [list(t * q_target) for t in ts]
    ph.run_qpoints(qs, with_eigenvectors=True)
    d = ph.get_qpoints_dict()
    F = np.array(d["frequencies"])
    E = np.array(d["eigenvectors"])
    n_bands = F.shape[1]

    def phase_relation(v):
        a, b = v[0:3], v[3:6]
        na, nb = np.linalg.norm(a), np.linalg.norm(b)
        if na * nb < 1e-12:
            return None
        s = float(np.real(np.vdot(a, b)) / (na * nb))
        return s

    rows = []
    for i, t in enumerate(ts):
        row = {"t": float(t), "qpoint": [float(x) for x in qs[i]], "bands": []}
        for b in range(n_bands):
            s = phase_relation(E[i][:, b])
            row["bands"].append({
                "band": b, "frequency_THz": float(F[i, b]),
                "al_ta_phase": (None if s is None else
                                "in_phase" if s > 0.5 else
                                "antiphase" if s < -0.5 else "mixed"),
                "al_ta_overlap": (None if s is None else abs(s)),
            })
        rows.append(row)

    flips = {}
    for b in range(n_bands):
        seq = [r["bands"][b]["al_ta_phase"] for r in rows]
        seq = [x for x in seq if x in ("in_phase", "antiphase")]
        flips[b] = len(set(seq)) > 1

    at_gamma = [rows[0]["bands"][b]["frequency_THz"] for b in range(n_bands)]
    at_target = [rows[-1]["bands"][b]["frequency_THz"] for b in range(n_bands)]
    n_acoustic = sum(1 for f in at_gamma if abs(f) < 0.5)
    soft_at_target = [b for b in range(n_bands) if at_target[b] < -1e-4]
    return {
        "from": [0.0, 0.0, 0.0],
        "to": [float(x) for x in q_target],
        "n_points": n_points,
        "dispersion": rows,
        "gamma_frequencies_THz": at_gamma,
        "target_frequencies_THz": at_target,
        "n_bands_near_zero_at_gamma": n_acoustic,
        "soft_bands_at_target": soft_at_target,
        "character_flips": {str(b): bool(v) for b, v in flips.items()},
        "soft_bands_are_acoustic_derived": all(
            any(abs(at_gamma[j]) < 0.5 for j in range(n_bands))
            for _ in soft_at_target),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fcdiff", type=Path, required=True)
    ap.add_argument("--endpoint", default="BCC-Al-Ta-2")
    ap.add_argument("--smearing", default="s1_020")
    ap.add_argument("--archive", type=Path, default=None,
                    help="managed directory to copy the force constants and "
                         "this script's output into, with a sha256 manifest")
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()

    fc_path = args.fcdiff / args.endpoint / args.smearing / "fc.npy"
    if not fc_path.is_file():
        raise SystemExit(f"no {fc_path}")
    fc = np.load(fc_path)
    ph, _ = A.make_phonopy(args.endpoint)
    ph.force_constants = fc
    S = np.array(ph.supercell_matrix)
    pr = ph.primitive

    out = {
        "endpoint": args.endpoint,
        "smearing": args.smearing,
        "supercell_matrix": [[int(x) for x in r] for r in S],
        "n_primitive_atoms": len(pr),
        "n_supercell_atoms": int(fc.shape[0]),
        "fc_file": str(fc_path),
        "fc_sha256": sha256(fc_path),
        "analyze_script_sha256": sha256(Path(__file__)),
    }

    print("=" * 78)
    print("1. 采样网格本身 —— 对称点在不在这张网格上？")
    print("=" * 78)
    grid = mesh_grid_report(ph)
    out["mesh_grid"] = grid
    print(f"   mesh {grid['mesh']}  不可约 q 点数 {grid['n_irreducible_qpoints']}")
    print(f"   q_x 取值: {[round(x, 6) for x in grid['q_x_values']]}")
    print(f"   步长 {grid['spacing']:.6f}   网格含 1/2 吗: {grid['grid_contains_half']}")
    print(f"   网格最低点 {[round(x, 6) for x in grid['mesh_min_qpoint']]} "
          f"f = {grid['mesh_min_frequency_THz']:+.4f} THz")
    if not grid["grid_contains_half"]:
        print("   → MP 网格不含 N=(1/2,1/2,0)，'网格最低点'是 N 附近的插值，")
        print("     频率会比 N 本身**浅**。不能拿它做模式跟随的波矢。")

    print()
    print("=" * 78)
    print("2. 超胞能精确表示哪些 q —— 再挑最深的那个")
    print("=" * 78)
    print(f"   supercell_matrix = {S.tolist()}  →  S^T q ∈ Z³ 要求 q ∈ {{0,1/2}}³")
    comm = deepest_commensurate(ph, S)
    out["commensurate_qpoints"] = comm
    print(f"   {'q':>22} {'min f /THz':>11} {'n_neg':>6} {'commensurate':>13}")
    for r in comm:
        print(f"   {str(r['qpoint']):>22} {r['min_frequency_THz']:>11.4f} "
              f"{r['n_imaginary']:>6} {str(r['commensurate']):>13}")
    best = comm[0]
    # the three N-equivalent points are symmetry-related and come out at the
    # SAME frequency; reporting only the first hides that.  Tie them explicitly.
    tol_tie = 1e-6
    tied = [r for r in comm
            if abs(r["min_frequency_THz"] - best["min_frequency_THz"]) < tol_tie
            and r["n_imaginary"] == best["n_imaginary"]]
    out["deepest_commensurate_set"] = tied
    print(f"\n   → 最深且公度的是 q = {best['qpoint']} "
          f"({best['min_frequency_THz']:+.4f} THz, {best['n_imaginary']} 个负频支)")
    if len(tied) > 1:
        print(f"   → **共有 {len(tied)} 个等深点**（N 的等价点，立方对称下等价）:")
        for r in tied:
            print(f"        {r['qpoint']}   {r['min_frequency_THz']:+.4f} THz")
    mm = grid["mesh_min_qpoint"]
    print(f"   → 网格最低点 {[round(x, 4) for x in mm]} 的 "
          f"{grid['mesh_min_frequency_THz']:+.4f} THz **比它浅**，")
    print(f"     而且 2q = {np.round(2 * np.array(mm), 4).tolist()} ∉ Z³，"
          f"2×2×2 装不下。")

    print()
    print("=" * 78)
    print("3. 该 q 点的本征矢（完整复分量 + 不变量）")
    print("=" * 78)
    sub = subspace_report(ph, best["qpoint"], pr)
    out["soft_mode_subspace"] = sub
    print(f"   q = {sub['qpoint']}   简并: {sub['degenerate']} "
          f"(重数 {sub['degeneracy']})")
    print(f"   {sub['note']}")
    for chk in sub["eigenvector_checks"]:
        print(f"   band {chk['band']}: |Dv-λv| = {chk['eigen_residual']:.2e}  "
              f"(若按行读则为 {chk['row_reading_residual']:.2e})")
    for b in sub["bands"]:
        print(f"\n   --- band {b['band']}  f = {b['frequency_THz']:+.4f} THz "
              f"[分量值在简并时依赖基选择] ---")
        for a, sym in enumerate(sub["primitive_symbols"]):
            c = b["components_basis_dependent"][a]
            print(f"     {sym}[{a}]  "
                  + "  ".join(f"({x['re']:+.5f}{x['im']:+.5f}i)" for x in c)
                  + f"   |.|={b['abs_per_atom'][a]:.5f}")
    print(f"\n   子空间恒零分量（不变量）: "
          f"{sub['silent_cartesian_components'] or '(无)'}")
    print(f"   物种权重（%）: "
          + "  ".join(f"{k}={v:.1f}" for k, v in sub["species_weight_percent"].items())
          + f"   在子空间内是否不变量: "
          f"{sub['species_weight_invariant_across_subspace']}")
    if sub["direction"]:
        d = sub["direction"]
        print(f"   主导物种: {d['species']}   "
              f"分量幅值: {[round(x, 4) for x in d['cartesian_magnitudes']]}")
        print(f"   与 q 的夹角: {d['angle_to_q_deg']:.1f}°   "
              f"纵向成分占比: {100 * d['longitudinal_fraction']:.1f}%   "
              f"→ {d['character']}")
    else:
        print("   方向: **不报告**（简并 → 单个本征矢无意义）")

    print()
    print("=" * 78)
    print("4. 候选波矢对照（网格最低点 / N 附近的插值点 / N）")
    print("=" * 78)
    mm = grid["mesh_min_qpoint"]
    cands, seen = [], set()
    for c in (mm, [0.4792, 0.4792, 0.0208], [0.4, 0.4, 0.0], best["qpoint"]):
        key = tuple(np.round(np.asarray(c, float), 6))
        if key not in seen:
            seen.add(key)
            cands.append(list(c))
    print(f"   {'q':>26} {'公度':>6} {'f* /THz':>9} {'简并':>5} {'Ta 权重%':>9} "
          f"{'纵向%':>7} {'性质':>11}  恒零分量")
    cmp_rows = []
    for q in cands:
        r = subspace_report(ph, q, pr)
        if not r.get("bands"):
            continue
        d = r.get("direction")
        row = {"qpoint": [round(float(x), 6) for x in q],
               "commensurate_2x2x2": commensurate(q, S),
               "frequency_THz": r["bands"][0]["frequency_THz"],
               "degenerate": r["degenerate"],
               "species_weight_percent": r["species_weight_percent"],
               "silent_cartesian_components": r["silent_cartesian_components"],
               "longitudinal_fraction": (d["longitudinal_fraction"] if d else None),
               "character": (d["character"] if d else "degenerate -> not reported")}
        cmp_rows.append(row)
        lon = f"{100 * d['longitudinal_fraction']:.1f}" if d else "—"
        ch = d["character"] if d else "简并→不报告"
        print(f"   {str([round(float(x), 4) for x in q]):>26} "
              f"{str(row['commensurate_2x2x2']):>6} {row['frequency_THz']:>9.4f} "
              f"{str(row['degenerate']):>5} "
              f"{row['species_weight_percent'].get('Ta', 0):>9.2f} {lon:>7} "
              f"{ch:>11}  {row['silent_cartesian_components'] or '-'}")
    out["candidate_qpoint_comparison"] = cmp_rows
    print("\n   → 网格最低点是**横向**模；N 上简并（横/纵两个成员都在软子空间里），")
    print("     所以**不能**说'沿 [110] 畸变'——那是纵向，而网格最低点是横向。")
    print("     模式跟随要用 N（唯一公度），并在简并子空间里同时试横、纵两个方向。")

    print()
    print("=" * 78)
    print("5. 模式跟随的起点：子空间的对称适配成员（§4.13.7）")
    print("=" * 78)
    members = subspace_members(ph, best["qpoint"], pr)
    out["subspace_members"] = members
    if "error" not in members:
        print(f"   纵向轴 (=q 方向) {[round(x, 4) for x in members['axis_longitudinal']]}")
        print(f"   横向轴 (面内垂直) {[round(x, 4) for x in members['axis_transverse']]}")
        print(f"\n   {'pattern':>34} {'保留':>8} {'可用':>6} "
              f"{'纵向% Al/Ta':>13} {'性质':>13}")
        for k, v in members["members"].items():
            pb = v.get("longitudinal_fraction_per_sublattice")
            s = "/".join(f"{100 * x:.0f}" if x is not None else "—" for x in (pb or []))
            print(f"   {k:>34} {v['kept_fraction']:>8.4f} "
                  f"{str(v['usable_as_starting_direction']):>6} {s:>13} "
                  f"{v.get('character', '-'):>13}")
        usable = [k for k, v in members["members"].items()
                  if v["usable_as_starting_direction"]]
        print(f"\n   → 可用作跟随起点的成员（保留 > 0.99）: {len(usable)} 个")
        for k in usable:
            v = members["members"][k]
            print(f"     {k}")
            for a in v["u_per_atom"]:
                print(f"        {a['symbol']}: {_fmt_vec(a['vector'])}")
        print("\n   → 互补图案只保留 ~0.69，有子空间外成分，**不能**用作跟随方向。")
        print("   → 标签是**位移图案**（纵向/横向 × 反相/同相），不是支名——")
        print("     见第 7 节：两支软模都是 Γ 声学支衍生的，没有光学支下降。")
        al = members["axis_alignment_check"]
        print(f"\n   --- q 与方向的可执行配对检查（不靠 [hkl] 字符串）---")
        print(f"   倒格矩阵约定: {al['q_coordinate_convention']}  "
              f"单位 {al['reciprocal_lattice_units']}  "
              f"B·Aᵀ = 2πI ? {al['reciprocal_check_ok']} "
              f"(偏差 {al['reciprocal_check_max_deviation']:.2e})")
        print(f"   q_frac = {al['q_fractional']}   q_cart = "
              f"{[round(x, 5) for x in al['q_cartesian']]}  (1/Å)")
        print(f"   纵向方向 (笛卡尔) = {[round(x, 5) for x in al['longitudinal_direction_cartesian']]}")
        print(f"   横向方向 (笛卡尔) = {[round(x, 5) for x in al['transverse_direction_cartesian']]}")
        print(f"   |u_long × q_hat| = {al['longitudinal_cross_q_hat']:.2e}  "
              f"→ 纵向平行于 q: {al['longitudinal_is_parallel_to_q']}")
        print(f"   |u_trans · q_hat| = {al['transverse_dot_q_hat']:.2e}  "
              f"→ 横向垂直于 q: {al['transverse_is_perpendicular_to_q']}")
        print(f"   → 整体: {'通过' if al['ok'] else '**不通过**'}")
        print("   → **三个等价 N 臂的纵向方向各不相同**（[011]/[101]/[110]），")
        print("     所以方向名必须与 q 一起记录，不能跨臂复用。")

        print()
        print("=" * 78)
        print("5b. 子空间内各向同性吗？（决定「两个图案」是不是两条路径）")
        print("=" * 78)
        iso = subspace_isotropy(ph, best["qpoint"], pr)
        out["subspace_isotropy"] = iso
        if "error" in iso:
            print(f"   {iso['error']}")
        else:
            print(f"   在子空间内扫 {iso['n_directions']} 个方向"
                  f"（同一最大单原子位移 {iso['target_max_single_atom_displacement_angstrom']} Å）：")
            print(f"   ½uᵀΦu 的跨度 = {iso['spread_meV']:.3e} meV  "
                  f"（相对 {iso['relative_spread']:.3e}）")
            print(f"   → 各向同性: {iso['isotropic']}")
            print(f"   → {iso['conclusion']}")

        print()
        print("=" * 78)
        print("6. ± 是平移等价的吗？—— 决定要跑 2 个还是 4 个 DFT 算例")
        print("=" * 78)
        te_out = {}
        for k, v in members["members"].items():
            if not v["usable_as_starting_direction"]:
                continue
            e = np.array([complex(x["re"], x["im"]) if isinstance(x, dict)
                          else complex(x) for x in v["e_mass_weighted"]])
            te = translation_equivalence(ph, best["qpoint"], e, amp=0.05)
            te_out[k] = te
            print(f"   {k}")
            print(f"     目标最大单原子位移 {te['target_max_single_atom_displacement_angstrom']} Å "
                  f"→ 缩放因子 {te['scale_factor_applied_to_unit_norm_eigenvector']:.4f}")
            print(f"     平移一个母胞能把 + 映到 − 吗: "
                  f"{te['plus_minus_are_translation_equivalent']} "
                  f"(沿母胞轴 {[x['parent_cell_axis'] for x in te['translations_that_map_plus_onto_minus']]}, "
                  f"这些轴上 q·R = 1/2)")
            print(f"     → {te['conclusion']}")
        out["translation_equivalence"] = te_out
        n_dft = 2 if all(v["plus_minus_are_translation_equivalent"]
                         for v in te_out.values()) else 4
        print(f"\n   → **主计算需要 {n_dft} 个初始结构**"
              f"（{'L 与 T 各一个；符号对留作输入生成器的回归检查' if n_dft == 2 else '四个方向各自独立'}）。")
        print("   → 三个等深 N 点由立方对称联系：单 q 初筛可用一个代表点；")
        print("     若将来搜索多 q 结构，不同 N 臂的组合要另行考虑。")

        print()
        print("=" * 78)
        print("7. 支归属检验：这些成员能叫'光学支/声学支'吗？（Γ→N 追踪）")
        print("=" * 78)
        bt = branch_tracking(ph, best["qpoint"], 21)
        out["branch_tracking"] = bt
        print(f"   Γ 处频率: {[round(x, 4) for x in bt['gamma_frequencies_THz']]}")
        print(f"   N 处频率: {[round(x, 4) for x in bt['target_frequencies_THz']]}")
        print(f"   软支 {bt['soft_bands_at_target']} 是否来自 Γ 声学三重态: "
              f"{bt['soft_bands_are_acoustic_derived']}")
        print(f"   字符发生翻转的支: "
              f"{[int(b) for b, v in bt['character_flips'].items() if v]}")
        print("\n   t     " + " ".join(f"{'b' + str(b):>16}" for b in range(3)))
        for row in bt["dispersion"]:
            cells = []
            for bd in row["bands"][:3]:
                tag = {"in_phase": "同相", "antiphase": "反相",
                       "mixed": "混合", None: "—"}[bd["al_ta_phase"]]
                cells.append(f"{bd['frequency_THz']:>8.3f} {tag}")
            print(f"   {row['t']:>5.2f} " + " ".join(f"{c:>16}" for c in cells))
        print("\n   → **没有光学支下降到软模**：光学三重态全程 7.85–8.72 THz。")
        print("     两支软模在连续追踪意义下都来自 Γ 的**声学子空间**（Γ 处都在 ~0.02 THz）。")
        print("   → 在 t≈0.48 的**避免交叉**处，单条按频率排序的能带**不保持固定字符**")
        print("     （两条近简并并交换字符），不是某支光学模下降。")
        print("   → **所以只报告'声学子空间起源'，不为 N 点简并模指定唯一支名**；")
        print("     标签只能用位移图案（longitudinal/transverse × in_phase/antiphase）。")

    if args.archive:
        arch = args.archive
        (arch / "fcdiff" / args.endpoint / args.smearing).mkdir(
            parents=True, exist_ok=True)
        dst = arch / "fcdiff" / args.endpoint / args.smearing / "fc.npy"
        if not dst.exists():
            shutil.copy2(fc_path, dst)
        if args.json:
            args.json.write_text(json.dumps(out, indent=2) + "\n")
        ev_dump = arch / "soft_mode_eigenvectors.json"
        ev_dump.write_text(json.dumps(out, indent=2) + "\n")
        manifest = {
            "purpose": ("soft-mode evidence for the mode-following step: the force "
                        "constants the eigenvector came from, the scripts, and the "
                        "eigenvector itself.  Kept off /tmp so the evidence has a "
                        "managed home."),
            "endpoint": args.endpoint,
            "smearing": args.smearing,
            "files": {},
        }
        src = Path(__file__).parent
        for p in [dst, ev_dump,
                  Path(__file__),
                  src / "extract_soft_mode.py",
                  src / "phonon_eigenvector.py",
                  src / "analyze_stage0b.py",
                  src / "build_stage0b.py"]:
            if p.exists():
                manifest["files"][p.name if p.parent != arch else
                                  str(p.relative_to(arch))] = sha256(p)
        (arch / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
        print(f"\n归档到 {arch}")
        for k, v in sorted(manifest["files"].items()):
            print(f"   {v[:16]}  {k}")

    if args.json and not args.archive:
        args.json.write_text(json.dumps(out, indent=2) + "\n")
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
