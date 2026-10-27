"""Custodian handler sets for the htvasp Workers.

Why this module exists
----------------------
``atomate2`` installs eleven error handlers by default.  ``htvasp.Worker`` used
to replace that list with ``[VaspErrorHandler()]`` alone, so ten of them were
silently dropped -- including ``UnconvergedErrorHandler``, which is the only one
that reacts to "the SCF did not converge".  A run that hit ``NELM`` therefore
passed through custodian untouched and relied entirely on the TaskDoc state
check at the job boundary.

Policy: custodian may correct
-----------------------------
Several handlers change *scientific inputs* to make a run proceed --
``VaspErrorHandler`` can rewrite KPOINTS, ISMEAR/SIGMA, POSCAR or POTIM;
``UnconvergedErrorHandler`` raises ALGO.  Refusing those corrections was tried
and rejected: a run that cannot proceed is worth less than a run whose actual
inputs are recorded.  Custodian saves ``INCAR.orig``, ``KPOINTS.orig``,
``POSCAR.orig`` and ``POTCAR.orig`` before it touches anything, so the DFPT
collector audits original against current and reports every difference
(``htvasp.dfpt_qha.input_audit``).  The correction is allowed; the fact that it
happened is not allowed to be invisible.

The set is therefore the official handlers, unmodified:

* ``VaspErrorHandler`` -- 45 signatures, the general input/SCF/symmetry faults.
* ``UnconvergedErrorHandler`` and ``NonConvergingErrorHandler`` -- non-convergence.
* ``FrozenJobErrorHandler`` -- a node that stops writing ``vasp.out``.
* ``StdErrHandler`` -- MPI/scheduler-level aborts.
* ``PositiveEnergyErrorHandler`` -- a diverged total energy.

``LargeSigmaHandler``, ``IncorrectSmearingHandler``, ``MeshSymmetryErrorHandler``,
``KspacingMetalHandler`` and ``PotimErrorHandler`` are still not added.  That is
not a judgement on their corrections -- it is that atomate2's defaults do not
include them, and the input audit is what makes any correction visible.
"""

from __future__ import annotations

import logging
from typing import Any, Sequence

from custodian.custodian import ErrorHandler
from custodian.vasp.handlers import (
    FrozenJobErrorHandler,
    NonConvergingErrorHandler,
    PositiveEnergyErrorHandler,
    StdErrHandler,
    UnconvergedErrorHandler,
    VaspErrorHandler,
)

log = logging.getLogger(__name__)


class DfptKpointSetChangeHandler(ErrorHandler):
    """Fail a DFPT run that asks VASP to change its k-point set mid-run.

    ``IBRION=8`` solves the Sternheimer equations at a fixed k-point set; VASP
    aborts with "VASP internal routines have requested a change of the k-point
    set" when the symmetry of the perturbed cell would require a different set.
    This project has hit it: on VASP 6.4.3 a DFPT run with ``NCORE`` set (rather
    than ``KPAR``) aborts this way.

    It is deliberately a *detector*, not a corrector.  The only corrections
    would be to change the parallelisation or the mesh, and both change what the
    run computes.  Returning no actions makes custodian raise
    ``NonRecoverableError`` -- loud, and early enough to see why.
    """

    is_monitor = True

    MESSAGES = (
        "requested a change of the k-point set",
        "change of the k-point set",
    )

    def __init__(self, output_filename: str = "vasp.out") -> None:
        self.output_filename = output_filename
        self.errors: set[str] = set()

    def check(self, directory: str = "./") -> bool:
        path = f"{directory}/{self.output_filename}"
        try:
            with open(path, errors="replace") as handle:
                text = handle.read()
        except OSError:
            return False
        hit = any(message in text for message in self.MESSAGES)
        if hit:
            self.errors = {"dfpt_kpoint_set_change"}
        return hit

    def correct(self, directory: str = "./") -> dict[str, Any]:
        return {"errors": ["dfpt_kpoint_set_change"], "actions": []}


def base_handlers() -> list[ErrorHandler]:
    """The handler set every Worker gets.

    Returns a fresh list each call -- custodian keeps per-run counters on the
    handler instances, so sharing one across jobs would carry state between
    them.
    """

    return [
        VaspErrorHandler(),
        UnconvergedErrorHandler(),
        NonConvergingErrorHandler(),
        FrozenJobErrorHandler(),
        StdErrHandler(),
        PositiveEnergyErrorHandler(),
    ]


def dfpt_handlers() -> list[ErrorHandler]:
    """Extra handlers the DFPT-QHA worker adds on top of the base set."""

    return [DfptKpointSetChangeHandler()]


def with_handlers(
    extra: Sequence[ErrorHandler] = (),
    *,
    run_vasp_kwargs: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Compose ``run_vasp_kwargs`` from the base set plus a worker's own.

    ``vasp_args`` callers may still pass a ``handlers`` key; it is appended
    rather than allowed to replace the base set, because replacing it is how the
    ten default handlers went missing in the first place.
    """

    kwargs = dict(run_vasp_kwargs or {})
    handlers = list(kwargs.pop("handlers", []))
    kwargs["handlers"] = base_handlers() + list(extra) + handlers
    # ``max_errors_per_job`` is a retry budget, not a convergence criterion.
    # Three is enough for the usual one-shot corrections, but not for the
    # NBANDS handshake: when VASP's own auto-NBANDS is not divisible by the
    # band-group count (ntasks/(KPAR*NCORE)), VASP adjusts it and says so,
    # custodian's ``too_few_bands`` fix adds one band, VASP adjusts again, and
    # the loop needs four attempts to settle.  BCC-Al-Al-2 (NELECT=6, two atoms)
    # spent its whole budget that way and aborted although every attempt had
    # converged physically -- the third already printed "reached required
    # accuracy".  Six leaves room for the handshake plus a genuine correction.
    kwargs.setdefault(
        "custodian_kwargs", {"max_errors_per_job": 6, "gzipped_output": False}
    )
    return kwargs


__all__ = [
    "DfptKpointSetChangeHandler",
    "base_handlers",
    "dfpt_handlers",
    "with_handlers",
]
