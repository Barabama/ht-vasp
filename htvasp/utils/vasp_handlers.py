"""Custodian handler sets, with the science-changing corrections guarded.

Why this module exists
----------------------
``atomate2`` installs eleven error handlers by default.  ``htvasp.Worker`` used
to replace that list with ``[VaspErrorHandler()]`` alone, so ten of them were
silently dropped -- including ``UnconvergedErrorHandler``, which is the only one
that reacts to "the SCF did not converge".  A run that hit ``NELM`` therefore
passed through custodian untouched and relied entirely on the TaskDoc state
check at the job boundary.

The reason that override existed is real, though: several default handlers
change **scientific inputs** to make a run proceed.  ``LargeSigmaHandler``
rewrites SIGMA, ``IncorrectSmearingHandler`` rewrites ISMEAR,
``MeshSymmetryErrorHandler`` can rewrite the k-mesh.  In a frozen-recipe pool
those are not fixes, they are silent recipe changes.

So the set here is the union of what is useful and what is safe:

* ``VaspErrorHandler`` -- 45 signatures, the general input/SCF/symmetry faults.
* ``UnconvergedErrorHandler`` and ``NonConvergingErrorHandler`` -- the two that
  react to non-convergence, wrapped so their ALGO changes must be *more* stable.
* ``FrozenJobErrorHandler`` -- a node that stops writing ``vasp.out``.
* ``StdErrHandler`` -- MPI/scheduler-level aborts.
* ``PositiveEnergyErrorHandler`` -- a diverged total energy.

Not included, and why: ``LargeSigmaHandler`` (rewrites the frozen SIGMA),
``IncorrectSmearingHandler`` (rewrites the frozen ISMEAR),
``MeshSymmetryErrorHandler`` (can rewrite the explicit k-mesh),
``KspacingMetalHandler`` (SCAN-only), ``PotimErrorHandler`` (rewrites the
relaxation step size, and DFPT has no ionic loop).
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
from pymatgen.io.vasp import Incar

log = logging.getLogger(__name__)

#: How much more stable an algorithm is.  ``All`` and ``Damped`` are the
#: progressively more expensive, more robust end of the ladder; ``VeryFast`` the
#: least.  A correction that moves DOWN this ordering is a downgrade.
ALGO_STABILITY: dict[str, int] = {
    "veryfast": 0,
    "fast": 1,
    "normal": 2,
    "all": 3,
    "damped": 3,
    "chi": 3,          # response-function mode, not part of the ladder
    "conjugate": 2,
    "subrot": 2,
    "eigenval": 2,
    "none": 2,
    "exact": 3,
}


class UnstableAlgoCorrection(RuntimeError):
    """A handler proposed replacing ALGO with a less stable algorithm."""


def _algo_rank(algo: Any) -> int | None:
    if not isinstance(algo, str):
        return None
    return ALGO_STABILITY.get(algo.strip().lower())


def guard_algo_actions(
    result: dict[str, Any], directory: str, *, handler: str
) -> dict[str, Any]:
    """Drop or reject corrections that make the electronic algorithm weaker.

    Custodian's own ladder is ``VeryFast -> Fast -> Normal -> All``, described in
    its source as "progressively switch to more stable but more expensive
    algorithms", so those moves are allowed.  The last-resort branch of the same
    ladder sets ``ALGO = Normal`` alongside mixing parameters, which for a run
    already at ``All`` is a downgrade -- and a downgrade applied to a DFPT run
    silently changes the dynamical matrix it is about to compute.

    Anything that lowers the rank raises ``UnstableAlgoCorrection`` rather than
    being dropped: a correction that is quietly skipped leaves the run looking
    handled when it is not.
    """

    actions = result.get("actions") or []
    current = Incar.from_file(f"{directory}/INCAR").get("ALGO")
    current_rank = _algo_rank(current)

    kept: list[dict[str, Any]] = []
    for action in actions:
        proposed = (action.get("action") or {}).get("_set") or {}
        if "ALGO" not in proposed:
            kept.append(action)
            continue
        new_rank = _algo_rank(proposed["ALGO"])
        if current_rank is None or new_rank is None:
            kept.append(action)          # unranked value: leave custodian's judgement
            continue
        if new_rank < current_rank:
            raise UnstableAlgoCorrection(
                f"{handler} proposed ALGO {current!r} -> {proposed['ALGO']!r}, "
                f"which is a less stable algorithm.  A run is only allowed to move "
                f"toward more stable settings; anything else changes the physics "
                f"rather than the convergence."
            )
        if new_rank > current_rank:
            log.info(
                "%s: ALGO %r -> %r (more stable, allowed)",
                handler, current, proposed["ALGO"],
            )
        kept.append(action)

    return {**result, "actions": kept}


class StableUnconvergedHandler(UnconvergedErrorHandler):
    """``UnconvergedErrorHandler`` that may only strengthen the algorithm."""

    def correct(self, directory: str = "./") -> dict[str, Any]:
        result = super().correct(directory=directory)
        return guard_algo_actions(result, directory, handler=type(self).__name__)


class StableNonConvergingHandler(NonConvergingErrorHandler):
    """``NonConvergingErrorHandler`` that may only strengthen the algorithm."""

    def correct(self, directory: str = "./") -> dict[str, Any]:
        result = super().correct(directory=directory)
        return guard_algo_actions(result, directory, handler=type(self).__name__)


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
    """The handler set every Worker gets: useful, and safe for a frozen recipe.

    Returns a fresh list each call -- custodian keeps per-run counters on the
    handler instances, so sharing one across jobs would carry state between
    them.
    """

    return [
        VaspErrorHandler(),
        StableUnconvergedHandler(),
        StableNonConvergingHandler(),
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
    kwargs.setdefault(
        "custodian_kwargs", {"max_errors_per_job": 3, "gzipped_output": False}
    )
    return kwargs


__all__ = [
    "ALGO_STABILITY",
    "DfptKpointSetChangeHandler",
    "StableNonConvergingHandler",
    "StableUnconvergedHandler",
    "UnstableAlgoCorrection",
    "base_handlers",
    "dfpt_handlers",
    "guard_algo_actions",
    "with_handlers",
]
