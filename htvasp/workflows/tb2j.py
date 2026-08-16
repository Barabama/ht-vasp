"""
HT-VASP - TB2J Workflow

End-member magnetic exchange via Wannier90::

    R3 relax (DoubleRelax) -> SCF + Wannier90 -> TB2J (exchange.out)

The R3 relaxation uses a double relax (R3 -> R3) so the structure is fully
relaxed at ISIF=3 (with a WAVECAR hot-start between the two R3 steps). The
Wannier90 SCF runs fresh (ISTART=0, no WAVECAR copy from R3 - the R3 WAVECAR
has a different NBANDS and a hot-start conflicts with the Wannier90 interface
and custodian). NUM_WANN is auto-derived as ``n_atoms * 9`` (s+p+d orbitals
per atom) unless overridden at construction.

Parallelism (run-critical): the W90 SCF must run with NCORE=1 and KPAR=1
(the Wannier90/PEAD interface), so NPAR = NTASKS and NBANDS is pre-rounded to
a multiple of NTASKS by Tb2jInputSetGenerator.get_num_bands(). Tb2jWorker
defaults to ntasks=32, and any Slurm submission must use the same ``--ntasks``
(main_w90.py submits with ntasks=32); a mismatch would let VASP silently
change NBANDS and desync the INCAR from the fixed wannier90.win num_bands.
"""

import logging
from pathlib import Path
from typing import Any

from jobflow import Flow
from pymatgen.core import Structure
from pymatgen.io.vasp import Kpoints
from custodian.vasp.handlers import VaspErrorHandler
from atomate2.vasp.flows.core import DoubleRelaxMaker
from atomate2.vasp.jobs.core import RelaxMaker
from atomate2.vasp.sets.core import RelaxSetGenerator

from htvasp.tb2j.input_set import ORBITALS_PER_ATOM, Tb2jInputSetGenerator
from htvasp.tb2j.maker import Tb2jMaker
from htvasp.workflows.base import Worker

log = logging.getLogger(__name__)


class Tb2jWorker(Worker):
    """Worker for TB2J magnetic exchange calculations.

    Pipeline:
        1. R3 full structural relaxation (``DoubleRelaxMaker``: R3 -> R3).
        2. SCF + Wannier90 interface. Fresh SCF (ISTART=0, matching the verified
           ``tb2j_test/submit_tb2j.sh`` INCAR_W90) - a WAVECAR hot-start
           conflicts with the Wannier90 NBANDS requirement and custodian.
        3. TB2J ``wann2J.py`` solve for the exchange parameters J, including
           parsing of ``exchange.out`` into structured J pairs/tensors.

    The magnetic INCAR keys (ISPIN=2, MAGMOM as ``{element: magmom}`` dict,
    ISMEAR=1 / SIGMA=0.2, ENCUT=520, EDIFF=1e-5 for R3 / 1e-6 for W90) are
    carried in ``global_incar`` following the ``main_qha.py`` GLOBAL_INCAR
    pattern. The W90 EDIFF is forced to 1e-6 by this worker (verified value);
    R3 inherits the global EDIFF (base default 1e-5).

    Run-critical parallelism: the W90 SCF is hard-constrained to KPAR=1,
    NCORE=1 (Wannier90/PEAD interface), giving NPAR = NTASKS with NBANDS
    pre-rounded to a multiple of NTASKS. ``ntasks`` (default 32) must match
    the Slurm ``--ntasks`` the job is launched with; see ``main_w90.py``.

    Args:
        worker_name: Worker name.
        vasp_args: VASP execution arguments (vasp_cmd, vasp_gamma_cmd, ...).
        potcar_functional: POTCAR functional.
        global_incar: Global INCAR settings (magnetic keys, ENCUT, EDIFF, ...).
        relax_incar: Extra INCAR overrides for the R3 relax step.
        elements: Magnetic elements for TB2J ``--elements``. Default None =
            read from the structure (distinct species in POSCAR site order).
        num_wann: Number of Wannier orbitals. Default None = auto-derive as
            ``n_atoms * ORBITALS_PER_ATOM`` (9 per spd atom).
        kmesh: kmesh for the TB2J real-space interpolation (default 9 9 9).
        kpoints: Explicit k-mesh grid for the Wannier90 SCF (default 8x8x8,
            Gamma-centered). Wannier90 requires an explicit isotropic mesh.
        ntasks: Number of MPI ranks the W90 SCF runs on (default 32). Forwarded
            to the Tb2jInputSetGenerator so NBANDS is pre-rounded to a multiple
            of NPAR (= NTASKS since KPAR=NCORE=1), keeping the INCAR NBANDS
            equal to wannier90.win ``num_bands``. Must match the Slurm
            ``--ntasks`` used to launch the job.
    """

    def __init__(
        self,
        worker_name: str = "tb2j-worker",
        vasp_args: dict[str, Any] | None = None,
        potcar_functional="PBE",
        global_incar: dict[str, Any] | None = None,
        relax_incar: dict[str, Any] | None = None,
        elements: list[str] | None = None,
        num_wann: int | None = None,
        kmesh: tuple = (9, 9, 9),
        kpoints: tuple = (8, 8, 8),
        ntasks: int = 32,
        **kwargs,
    ):
        # Initialize base Worker
        super().__init__(
            worker_name=worker_name,
            vasp_args=vasp_args,
            potcar_functional=potcar_functional,
            global_incar=global_incar,
        )

        relax_incar = relax_incar or {}
        self.elements = list(elements) if elements else None
        self.kmesh = tuple(kmesh)
        self.kpoints = tuple(kpoints)
        self.num_wann = num_wann
        self.ntasks = ntasks

        # R3 structural relaxation (DoubleRelax: R3 -> R3)
        self.relax_maker = DoubleRelaxMaker.from_relax_maker(
            RelaxMaker(
                name="r3 relax",
                run_vasp_kwargs=self.run_vasp_kwargs,
                stop_children_kwargs={"handle_unsuccessful": False},
                copy_vasp_kwargs={"additional_vasp_files": ("WAVECAR",)},
                input_set_generator=RelaxSetGenerator(
                    user_potcar_functional=self.potcar_functional,
                    # Use the same explicit isotropic Gamma mesh as the Wannier90
                    # step (verified submit_tb2j.sh KPOINTS). The atomate2
                    # auto-kspacing mesh is anisotropic for HCP (10x10x6), which
                    # during the ISIF=3 cell relaxation biases the relaxed c/a
                    # (1.613 vs the isotropic 1.6064) and shifts the TB2J J
                    # values away from the verified benchmark.
                    user_kpoints_settings=Kpoints.gamma_automatic(kpts=self.kpoints),
                    user_incar_settings={
                        **self.global_incar,
                        "ISTART": 1,
                        "ISIF": 3,
                        "LWAVE": True,
                        # Enforce symmetry during relaxation (global_incar sets
                        # ISYM=0). ISYM=0 lets the cell drift slightly off the
                        # intended spacegroup (e.g. HCP angle 120.0001 deg).
                        # ISYM=1 matches the verified submit_tb2j.sh INCAR_R3.
                        "ISYM": 1,
                        **relax_incar,
                    },
                ),
            ),
        )

        # SCF + Wannier90 -> TB2J solve
        #
        # The Wannier90 interface in VASP always activates the PEAD routines,
        # which require NCORE=1; a larger NCORE makes VASP print "PEAD routines
        # do not work for NCORE". VASP also checks that NBANDS is divisible by
        # NPAR (= NTASKS with KPAR=NCORE=1) and silently bumps it up otherwise
        # (e.g. 38 -> 64), printing "number of bands changed".
        # Tb2jInputSetGenerator pre-rounds NBANDS to a multiple of ntasks
        # (get_num_bands), so VASP should never change it and the INCAR NBANDS
        # always matches wannier90.win num_bands. The W90 step runs a fresh SCF
        # (ISTART=0, no WAVECAR hot-start) with NCORE=1 (see the ISTART/NCORE
        # keys below). This reproduces the verified submit_tb2j.sh INCAR_W90.
        #
        # Custodian handler subset: auto_nbands and dfpt_ncore are intentionally
        # excluded from errors_subset_to_catch. The W90 step cannot legitimately
        # trigger them (NBANDS pre-rounded to a multiple of NPAR; NCORE forced
        # to 1), and custodian's automatic corrections for either would be wrong
        # for a Wannier90 run:
        # - auto_nbands ("The number of bands has been changed") would bump the
        #   INCAR NBANDS to the next NPAR multiple, re-desyncing it from the
        #   fixed wannier90.win num_bands. A real NBANDS drift here signals an
        #   ntasks mismatch (worker vs Slurm) that must fail loudly.
        # - dfpt_ncore ("PEAD routines do not work for NCORE") would UNSET
        #   NCORE and NPAR from the INCAR (custodian's fix), the opposite of
        #   the NCORE=1 the Wannier90 interface requires - silently breaking
        #   the PEAD step. Keeping it out of the catch subset surfaces such a
        #   misconfiguration (e.g. a stray global_incar NCORE override) as an
        #   explicit job failure instead of an incorrect auto-fix.
        w90_handler_subset = list(VaspErrorHandler().error_msgs)
        for _pattern in ("auto_nbands", "dfpt_ncore"):
            if _pattern in w90_handler_subset:
                w90_handler_subset.remove(_pattern)
        w90_run_vasp_kwargs = {
            **self.run_vasp_kwargs,
            "handlers": [VaspErrorHandler(errors_subset_to_catch=w90_handler_subset)],
        }
        self.tb2j_maker = Tb2jMaker(
            run_vasp_kwargs=w90_run_vasp_kwargs,
            elements=self.elements or [],
            kmesh=self.kmesh,
            ntasks=self.ntasks,
            input_set_generator=Tb2jInputSetGenerator(
                elements=self.elements or [],
                num_wann=self.num_wann,
                kpoints_grid=self.kpoints,
                ntasks=self.ntasks,
                user_potcar_functional=self.potcar_functional,
                user_incar_settings={
                    **self.global_incar,
                    # LASPH=True / ENAUG=1360 are NOT overridden here on purpose:
                    # atomate2's base VASP set injects them and they are kept (more
                    # accurate for magnetic transition metals). The verified
                    # tb2j_test baseline ran without them; the resulting small
                    # quantitative shift (e.g. BCC J +4.8%) is a known decision,
                    # not a bug (see main_w90.py GLOBAL_INCAR comment).
                    # fresh SCF, matching the verified submit_tb2j.sh INCAR_W90
                    # (a WAVECAR hot-start conflicts with NBANDS and custodian)
                    "ISTART": 0,
                    "ICHARG": 2,
                    # verified tighter SCF convergence for the Wannier90 interface
                    "EDIFF": 1e-6,
                    # Wannier90-interface requirements. These must come after
                    # global_incar: the base Worker defaults (NSW=10, IBRION=2,
                    # NELM=100, LREAL=Auto, LORBIT=None) would otherwise override
                    # Tb2jInputSetGenerator.W90_INCAR_DEFAULTS and turn this into
                    # an ionic relaxation instead of a static SCF.
                    "KPAR": 1,
                    # PEAD (used by the Wannier90 interface) requires NCORE=1
                    "NCORE": 1,
                    # match the verified INCAR_W90 (ISYM default = 1): lets the
                    # SCF use the symmetry-reduced k-mesh (~50 vs 512 k-points);
                    # wannier90 still receives the full 8x8x8 grid.
                    "ISYM": 1,
                    "NSW": 0,
                    "IBRION": -1,
                    "NELM": 200,
                    "LORBIT": 11,
                    "LREAL": False,
                    # no wave/charge output needed beyond the Wannier90 files
                    "LWAVE": False,
                    "LCHARG": False,
                },
            ),
        )

    @staticmethod
    def _structure_elements(structure: Structure) -> list[str]:
        """Distinct element symbols in POSCAR site order."""
        seen: list[str] = []
        for site in structure:
            sym = site.species_string
            if sym not in seen:
                seen.append(sym)
        return seen

    def _make_flow(self, structure: Structure, prev_dir: Path | str | None = None) -> Flow:
        # Resolve magnetic elements from the structure when not given explicitly.
        elements = self.elements if self.elements is not None else self._structure_elements(structure)
        self.tb2j_maker.elements = elements
        self.tb2j_maker.input_set_generator.elements = elements

        # Tb2jMaker.make() derives NUM_WANN from the structure at flow-build
        # time, but the relaxed structure here is an OutputReference (resolved
        # only at run time). Atom count is invariant under relaxation, so
        # resolve num_wann from the input structure.
        if self.num_wann is None:
            self.tb2j_maker.input_set_generator.num_wann = len(structure) * ORBITALS_PER_ATOM

        relax_job = self.relax_maker.make(structure, prev_dir)
        tb2j_flow = self.tb2j_maker.make(
            relax_job.output.structure,
            prev_dir=relax_job.output.dir_name,
        )
        # Nest the relax flow and the W90->TB2J flow (jobs cannot be shared
        # between flows, so the tb2j flow is added as a child flow).
        return Flow(
            [relax_job, tb2j_flow],
            output=tb2j_flow.output,
            name=f"{self.worker_name} flow",
        )

    def get_result(self, output_job_name: str = "tb2j solve") -> dict[str, Any] | None:
        return super().get_result(output_job_name)
