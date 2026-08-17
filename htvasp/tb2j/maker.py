"""TB2J Maker - Maker class for the Wannier90 -> TB2J exchange step."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from atomate2.vasp.jobs.core import StaticMaker
from jobflow import Flow, Maker
from pymatgen.core import Structure

from htvasp.tb2j.input_set import Tb2jInputSetGenerator
from htvasp.tb2j.jobs import tb2j_solve

log = logging.getLogger(__name__)


@dataclass
class Tb2jMaker(Maker):
    """Maker for the SCF+Wannier90 -> TB2J exchange part of the TB2J workflow.

    This maker only covers the Wannier90 SCF step and the TB2J solve step.
    The R3 structural relaxation is composed *before* this maker by the
    Worker (e.g., a RelaxMaker whose output structure is passed into
    ``make(structure, prev_dir)``).

    Args:
        name: Maker name.
        input_set_generator: Tb2jInputSetGenerator instance (writes wannier90.win
            and the Wannier90-interface INCAR keys).
        run_vasp_kwargs: VASP execution parameters for the W90 SCF step.
        elements: Magnetic elements passed to TB2J --elements. Defaults to the
            generator's ``elements`` if not given here.
        kmesh: kmesh for the TB2J real-space interpolation (default 9 9 9).
        ntasks: Number of MPI ranks for the W90 SCF. When set, forwarded to the
            ``input_set_generator`` so NBANDS is pre-rounded to a multiple of
            NPAR (= NTASKS with KPAR=NCORE=1). If None, the generator's own
            ``ntasks`` is used (default 32).
        copy_vasp_kwargs: Files copied from the previous (relax) directory into
            the W90 step. Default is no WAVECAR copy: the W90 step runs a fresh
            SCF (ISTART=0, matching the verified submit_tb2j.sh INCAR_W90) so
            the R3 WAVECAR (different NBANDS) must not be carried over.
        w90_name: Name of the W90 StaticMaker step.
        tb2j_name: Name of the tb2j_solve job.
    """

    name: str = "tb2j exchange"
    input_set_generator: Tb2jInputSetGenerator = field(default_factory=Tb2jInputSetGenerator)
    run_vasp_kwargs: dict[str, Any] = field(default_factory=dict)
    elements: list[str] = field(default_factory=list)
    kmesh: tuple[int, int, int] = (9, 9, 9)
    ntasks: int | None = None
    copy_vasp_kwargs: dict[str, Any] = field(default_factory=dict)
    w90_name: str = "w90 static"
    tb2j_name: str = "tb2j solve"

    def make(
        self,
        structure: Structure,
        prev_dir: str | Path | None = None,
    ) -> Flow:
        """Generate the W90 static SCF + TB2J solve flow.

        Args:
            structure: Structure to run the Wannier90 SCF on (typically the
                R3-relaxed structure).
            prev_dir: Optional previous calculation directory (the R3 relax
                dir). Forwarded to the W90 ``StaticMaker`` as its reference
                directory. No WAVECAR is hot-started: the default
                ``copy_vasp_kwargs`` is empty and the W90 step runs a fresh
                ISTART=0 SCF, so the R3 WAVECAR (different NBANDS) is not
                carried over.

        Returns:
            Flow([w90_static_job, tb2j_solve_job], output=tb2j_solve_job.output).
        """
        if self.ntasks is not None:
            self.input_set_generator.ntasks = self.ntasks
        w90_maker = StaticMaker(
            name=self.w90_name,
            input_set_generator=self.input_set_generator,
            run_vasp_kwargs=self.run_vasp_kwargs,
            copy_vasp_kwargs=self.copy_vasp_kwargs,
            stop_children_kwargs={"handle_unsuccessful": False},
        )
        w90_job = w90_maker.make(structure, prev_dir)

        elements = self.elements or self.input_set_generator.elements
        num_wann = self.input_set_generator.get_num_wann(structure)
        tb2j_job = tb2j_solve(
            w90_job.output.dir_name,
            elements=elements,
            kmesh=self.kmesh,
            num_wann=num_wann,
        )
        tb2j_job.name = self.tb2j_name

        return Flow([w90_job, tb2j_job], output=tb2j_job.output, name=self.name)
