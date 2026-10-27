#!/usr/bin/env python
"""Read a phonopy eigenvector — with the convention CHECKED, not assumed.

The trap this module exists to close
------------------------------------
`phonopy`'s `QpointsPhonon._run` does

    eigvals, eigvecs = np.linalg.eigh(dm)
    eigenvectors[i] = eigvecs

and `np.linalg.eigh` puts the eigenvectors in the COLUMNS.  Its own YAML writer
agrees -- `self._eigenvectors[i][k * 3 + ll, j]` has the component index first
and the band index second.  So the eigenvector of band j is

    eigenvectors[iq][:, j]        # column  -- correct
    eigenvectors[iq][j, :]        # row     -- WRONG, and this is the one that
                                  #            looks like the obvious reading

Both readings give a unit vector, so the mistake is invisible in the norm.  It
is NOT invisible in the result: at B2 AlTa the row reading says the soft mode is
Ta along [011] and the column reading says Ta along [1-10] -- a different
direction, and only the column reading survives the checks below.

There is a second, independent trap in the same place: whenever a branch is
DEGENERATE, `eigh` returns an ARBITRARY orthonormal basis of the degenerate
subspace.  Then no single vector is meaningful -- only the subspace is.  Reading
"the eigenvector" of a degenerate band and reporting its direction is wrong even
with the indexing right.

So this module returns the eigenvector together with the evidence that it is
one: the residual |D v - lambda v|.  A caller that gets a large residual has
either the wrong indexing, the wrong basis, or the wrong atom order.

    from phonon_eigenvector import eigenvector
    v, f, info = eigenvector(ph, q, band)
"""
from __future__ import annotations

import numpy as np

FACTOR = 15.633302      # phonopy's sqrt(eV/A^2/AMU) -> THz


def eigenvector(ph, qpoint, band: int, check: bool = True) -> tuple[np.ndarray, float, dict]:
    """The eigenvector of `band` at `qpoint`, as a (3*natom,) complex array.

    Reads the COLUMN (see module docstring) and, unless `check=False`, asserts
    that it really is an eigenvector of the dynamical matrix.  The assert is the
    point of the module: it turns a silent convention into a loud failure.
    """
    ph.run_qpoints([list(qpoint)], with_eigenvectors=True)
    d = ph.get_qpoints_dict()
    freqs = np.array(d["frequencies"])[0]
    E = np.array(d["eigenvectors"])[0]
    n = E.shape[0]
    if not (0 <= band < n):
        raise IndexError(f"band {band} out of range for {n} bands")
    v = E[:, band]                       # <-- COLUMN
    f = float(freqs[band])
    info = {"n_bands": int(n), "shape": list(E.shape)}

    if check:
        dm = ph.dynamical_matrix
        dm.run(np.asarray(qpoint, dtype=float))
        D = np.array(dm.dynamical_matrix)
        # phonopy reports f = sign(lambda)*sqrt(|lambda|)*FACTOR, so lambda =
        # sign(f)*(f/FACTOR)^2
        lam = float(np.sign(f) * (f / FACTOR) ** 2)
        resid = float(np.linalg.norm(D @ v - lam * v))
        row_resid = float(np.linalg.norm(D @ E[band, :] - lam * E[band, :]))
        info["eigen_residual"] = resid
        info["row_reading_residual"] = row_resid
        info["column_reading_is_eigenvector"] = resid < 1e-6
        if resid >= 1e-6:
            raise AssertionError(
                f"the column of eigenvectors[{band}] is not an eigenvector "
                f"(residual {resid:.3e}); the convention assumed by this module "
                f"does not hold for this phonopy version")
    return v, f, info


def degeneracy(ph, qpoint, band: int, tol: float = 1e-3) -> list[int]:
    """Bands degenerate with `band` at `qpoint` (including itself).

    Callers MUST check this before using a single eigenvector: on a degenerate
    branch the returned vector is an arbitrary basis element, so its individual
    components carry no information -- only the subspace does.
    """
    ph.run_qpoints([list(qpoint)], with_eigenvectors=False)
    f = np.array(ph.get_qpoints_dict()["frequencies"])[0]
    return [int(b) for b in range(len(f)) if abs(f[b] - f[band]) < tol]
