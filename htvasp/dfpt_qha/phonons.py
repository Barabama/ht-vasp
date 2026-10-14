"""Read a phonopy eigenvector -- with the convention CHECKED, not assumed.

The trap this module closes
---------------------------
``phonopy``'s ``QpointsPhonon._run`` does ``eigvals, eigvecs = np.linalg.eigh(dm)``
and ``np.linalg.eigh`` puts the eigenvectors in the COLUMNS.  Its own YAML
writer agrees: ``self._eigenvectors[i][k * 3 + ll, j]`` has the component index
first and the band index second.  So the eigenvector of band ``j`` is

    eigenvectors[iq][:, j]     # column -- correct
    eigenvectors[iq][j, :]     # row    -- WRONG, and this is the reading that
                               #           looks like the obvious one

Both readings give a unit vector, so the mistake is invisible in the norm.  It
is not invisible in the result: at B2 AlTa the row reading says the soft mode
is Ta along [011] and the column reading says Ta along [1-10] -- a different
direction, and only the column reading survives the residual check that
``eigenvector`` performs.

There is a second, independent trap in the same place: whenever a branch is
DEGENERATE, ``eigh`` returns an ARBITRARY orthonormal basis of the degenerate
subspace, so no single vector is meaningful -- only the subspace is.  Reporting
"the eigenvector" of a degenerate band is wrong even with the indexing right.
``degeneracy`` exists so a caller can check before doing that.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np

# phonopy reports f = sign(lambda) * sqrt(|lambda|) * FACTOR, so the eigenvalue
# is recovered as sign(f) * (f / FACTOR) ** 2.
FACTOR = 15.633302
EIGEN_RESIDUAL_TOL = 1e-6
DEGENERACY_TOL_THZ = 1e-3


@dataclass(frozen=True)
class BandEigenvector:
    """One band's eigenvector at one q-point, with the evidence that it is one.

    The residuals are not diagnostics to log and forget: a large ``residual``
    means the wrong indexing, the wrong basis, or the wrong atom order, and a
    caller that ignores it is back to the silent-convention failure this module
    exists to prevent.
    """

    qpoint: tuple[float, float, float]
    band: int
    frequency_thz: float
    #: (3 * n_primitive_atom,) complex, in the PRIMITIVE cell basis -- not the
    #: supercell.  Map to the supercell through ``primitive.p2s_map``.
    vector: np.ndarray
    n_bands: int
    matrix_shape: tuple[int, int]
    residual: float | None = None
    row_reading_residual: float | None = None

    @property
    def is_checked_eigenvector(self) -> bool:
        """Whether the residual check actually ran and passed."""

        return self.residual is not None and self.residual < EIGEN_RESIDUAL_TOL

    def as_dict(self) -> dict[str, Any]:
        """JSON-safe view.  The vector's components are emitted real-imaginary."""

        return {
            "qpoint": list(self.qpoint),
            "band_index": self.band,
            "frequency_THz": self.frequency_thz,
            "n_bands": self.n_bands,
            "matrix_shape": list(self.matrix_shape),
            "eigen_residual": self.residual,
            "row_reading_residual": self.row_reading_residual,
            "column_reading_is_eigenvector": self.is_checked_eigenvector,
            "eigenvector": [
                float(z.real) if abs(z.imag) < 1e-10
                else {"re": float(z.real), "im": float(z.imag)}
                for z in self.vector
            ],
        }


def _dynamical_matrix(ph: Any, qpoint: Sequence[float]) -> np.ndarray:
    matrix = ph.dynamical_matrix
    matrix.run(np.asarray(qpoint, dtype=float))
    return np.array(matrix.dynamical_matrix)


def eigenvector(
    ph: Any, qpoint: Sequence[float], band: int, check: bool = True
) -> BandEigenvector:
    """The eigenvector of ``band`` at ``qpoint``, read from the COLUMN.

    Unless ``check`` is False, asserts that the result really is an eigenvector
    of the dynamical matrix.  That assertion is the point of the module: it
    turns a silent convention into a loud failure.  ``check=False`` skips it and
    leaves both residuals as None -- it does not make the answer trustworthy.
    """

    qpoint = tuple(float(x) for x in qpoint)
    ph.run_qpoints([list(qpoint)], with_eigenvectors=True)
    data = ph.get_qpoints_dict()
    frequencies = np.array(data["frequencies"])[0]
    matrix = np.array(data["eigenvectors"])[0]
    n_bands = matrix.shape[0]
    if not 0 <= band < n_bands:
        raise IndexError(f"band {band} out of range for {n_bands} bands")

    vector = matrix[:, band]  # <-- COLUMN
    frequency = float(frequencies[band])

    residual: float | None = None
    row_residual: float | None = None
    if check:
        dynamical = _dynamical_matrix(ph, qpoint)
        eigenvalue = float(np.sign(frequency) * (frequency / FACTOR) ** 2)
        residual = float(np.linalg.norm(dynamical @ vector - eigenvalue * vector))
        row_residual = float(
            np.linalg.norm(dynamical @ matrix[band, :] - eigenvalue * matrix[band, :])
        )
        if residual >= EIGEN_RESIDUAL_TOL:
            raise AssertionError(
                f"the column of eigenvectors[{band}] is not an eigenvector "
                f"(residual {residual:.3e}); the convention this module assumes "
                f"does not hold for this phonopy version"
            )

    return BandEigenvector(
        qpoint=qpoint,
        band=int(band),
        frequency_thz=frequency,
        vector=vector,
        n_bands=int(n_bands),
        matrix_shape=(int(matrix.shape[0]), int(matrix.shape[1])),
        residual=residual,
        row_reading_residual=row_residual,
    )


def degeneracy(
    ph: Any, qpoint: Sequence[float], band: int, tol: float = DEGENERACY_TOL_THZ
) -> list[int]:
    """Bands degenerate with ``band`` at ``qpoint``, including ``band`` itself.

    Callers MUST check this before using a single eigenvector: on a degenerate
    branch the returned vector is an arbitrary basis element, so its individual
    components carry no information -- only the subspace does.
    """

    ph.run_qpoints([list(qpoint)], with_eigenvectors=False)
    frequencies = np.array(ph.get_qpoints_dict()["frequencies"])[0]
    return [int(b) for b in range(len(frequencies)) if abs(frequencies[b] - frequencies[band]) < tol]


__all__ = [
    "BandEigenvector",
    "DEGENERACY_TOL_THZ",
    "EIGEN_RESIDUAL_TOL",
    "FACTOR",
    "degeneracy",
    "eigenvector",
]
