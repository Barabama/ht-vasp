"""DFPT phonon jobs built on atomate2 common phonon jobs.

The job boundaries are jobflow boundaries.  VASP execution belongs to the
``BaseVaspMaker`` passed to these jobs, not to this module.
"""

from __future__ import annotations

from typing import Any

from jobflow import Response, job

from atomate2.common.jobs.phonons import generate_frequencies_eigenvectors


@job

def collect_dfpt_phonons(
    structure: Any,
    supercell_matrix: Any,
    displacement_data: dict[str, Any],
    *,
    displacement: float = 0.01,
    sym_reduce: bool = False,
    symprec: float = 1e-3,
    use_symmetrized_structure: str | None = "conventional",
    kpath_scheme: str = "seekpath",
    code: str = "vasp",
    total_dft_energy: float = 0.0,
    **kwargs: Any,
) -> Any:
    """Collect DFPT force/displacement results with the common phonon schema."""

    return generate_frequencies_eigenvectors.original(
        structure=structure,
        supercell_matrix=supercell_matrix,
        displacement=displacement,
        sym_reduce=sym_reduce,
        symprec=symprec,
        use_symmetrized_structure=use_symmetrized_structure,
        kpath_scheme=kpath_scheme,
        code=code,
        displacement_data=displacement_data,
        total_dft_energy=total_dft_energy,
        **kwargs,
    )


@job
def validate_dfpt_output(phonon_doc: Any, *, recipe: dict[str, Any] | None = None) -> Response:
    """Fail closed before a downstream thermal/QHA analysis job.

    The common phonon document remains the output schema for the calculation;
    this job only adds the production decision and stops children when the
    document is incomplete or the recipe is absent.
    """

    problems: list[str] = []
    if phonon_doc is None:
        problems.append("phonon document is missing")
    if recipe is None:
        problems.append("recipe provenance is missing")
    if problems:
        return Response(output={"ok": False, "problems": problems}, stop_children=True)
    return Response(output={"ok": True, "problems": []})


__all__ = ["collect_dfpt_phonons", "validate_dfpt_output"]
