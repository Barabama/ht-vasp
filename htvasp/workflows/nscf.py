"""
HT-VASP - NSCF Workflows

DOS + Band electronic structure calculations.
"""

import logging
from typing import Any

from jobflow import Flow
from pymatgen.core import Structure
from atomate2.vasp.flows.core import DoubleRelaxMaker
from atomate2.vasp.jobs.core import RelaxMaker, StaticMaker, NonSCFMaker
from atomate2.vasp.sets.core import (
    RelaxSetGenerator,
    StaticSetGenerator,
    NonSCFSetGenerator,
)

from htvasp.workflows.base import Worker

log = logging.getLogger(__name__)


class NscfWorker(Worker):
    """
    Worker for DOS + Band electronic structure calculations.

    Workflow: Relax → Static → NSCF-DOS (uniform mesh)
                         └→ NSCF-Band (line mode)
    """

    def __init__(
        self,
        worker_name: str = "nscf-worker",
        vasp_args: dict[str, Any] | None = None,
        potcar_functional: str = "PBE_64",
        global_incar: dict[str, Any] | None = None,
        relax_incar: dict[str, Any] | None = None,
        static_incar: dict[str, Any] | None = None,
        nscf_dos_incar: dict[str, Any] | None = None,
        nscf_band_incar: dict[str, Any] | None = None,
        compute_band: bool = True,
        reciprocal_density: int = 200,
        dedos: float = 0.02,
        line_density: int = 20,
        **kwargs,
    ):
        """
        Initialize NscfWorker.

        Args:
            worker_name: Name of the worker
            vasp_args: VASP command and handler settings
            potcar_functional: POTCAR functional type
            global_incar: Global INCAR settings
            relax_incar: Relax-specific INCAR settings
            static_incar: Static-specific INCAR settings
            nscf_dos_incar: NSCF-DOS-specific INCAR settings
            nscf_band_incar: NSCF-Band-specific INCAR settings
            compute_band: Whether to compute band structure
            reciprocal_density: K-point density for DOS (uniform mode)
            dedos: Energy resolution for DOS (eV), used to auto-calculate NEDOS
            line_density: Line density for band structure (line mode)
            **kwargs: Additional keyword arguments
        """
        # Initialize base Worker
        super().__init__(
            worker_name=worker_name,
            vasp_args=vasp_args,
            potcar_functional=potcar_functional,
            global_incar=global_incar,
        )

        # Disable VaspErrorHandler for NSCF calculations
        # NSCF jobs should be stable and fail fast if there's a real issue
        # The "number of bands has been changed" warning is handled by VASP itself
        self.run_vasp_kwargs["handlers"] = []

        relax_incar = relax_incar or {}
        static_incar = static_incar or {}
        nscf_dos_incar = nscf_dos_incar or {}
        nscf_band_incar = nscf_band_incar or {}

        self.compute_band = compute_band
        self.reciprocal_density = reciprocal_density
        self.dedos = dedos
        self.line_density = line_density

        # Structural relaxation (ISIF=3)
        relax_maker = DoubleRelaxMaker.from_relax_maker(
            RelaxMaker(
                run_vasp_kwargs=self.run_vasp_kwargs,
                input_set_generator=RelaxSetGenerator(
                    user_potcar_functional=self.potcar_functional,
                    user_incar_settings={
                        **self.global_incar,
                        "ISIF": 3,
                        **relax_incar,
                    },
                ),
            )
        )

        # Static calculation (generate CHGCAR for NSCF)
        static_maker = StaticMaker(
            run_vasp_kwargs=self.run_vasp_kwargs,
            input_set_generator=StaticSetGenerator(
                user_potcar_functional=self.potcar_functional,
                user_incar_settings={
                    **self.global_incar,
                    "ISTART": 1,
                    "IBRION": -1,
                    "ISIF": 2,
                    "NSW": 0,
                    "LWAVE": True,
                    "LCHARG": True,
                    **static_incar,
                },
            ),
        )

        # NSCF-DOS (uniform mesh for density of states)
        nscf_dos_maker = NonSCFMaker(
            name="nscf uniform",
            run_vasp_kwargs=self.run_vasp_kwargs,
            input_set_generator=NonSCFSetGenerator(
                mode="uniform",
                reciprocal_density=reciprocal_density,
                dedos=dedos,
                user_potcar_functional=self.potcar_functional,
                user_incar_settings={
                    **self.global_incar,
                    "IBRION": -1,
                    "NSW": 0,
                    "ICHARG": 11,
                    "LORBIT": 11,
                    **nscf_dos_incar,
                },
            ),
        )

        # NSCF-Band (line mode for band structure)
        nscf_band_maker = NonSCFMaker(
            name="nscf line",
            run_vasp_kwargs=self.run_vasp_kwargs,
            input_set_generator=NonSCFSetGenerator(
                mode="line",
                line_density=line_density,
                user_potcar_functional=self.potcar_functional,
                user_incar_settings={
                    **self.global_incar,
                    "IBRION": -1,
                    "NSW": 0,
                    "ICHARG": 11,
                    "LORBIT": 11,
                    **nscf_band_incar,
                },
            ),
        )

        self.flow_makers: tuple[
            DoubleRelaxMaker, StaticMaker, NonSCFMaker, NonSCFMaker
        ] = (
            relax_maker,
            static_maker,
            nscf_dos_maker,
            nscf_band_maker,
        )

    def _make_flow(self, structure: Structure) -> Flow:
        """
        Create the NSCF flow.

        Args:
            structure: Input structure

        Returns:
            Flow object with Relax → Static → NSCF-DOS/Band jobs
        """
        relax_maker, static_maker, nscf_dos_maker, nscf_band_maker = self.flow_makers

        # Create jobs
        relax_job = relax_maker.make(structure)
        static_job = static_maker.make(
            relax_job.output.structure,
            prev_dir=relax_job.output.dir_name,
        )
        nscf_dos_job = nscf_dos_maker.make(
            static_job.output.structure,
            prev_dir=static_job.output.dir_name,
        )

        if self.compute_band:
            nscf_band_job = nscf_band_maker.make(
                static_job.output.structure,
                prev_dir=static_job.output.dir_name,
            )
            # Both NSCF jobs depend on static, but not on each other
            return Flow(
                [relax_job, static_job, nscf_dos_job, nscf_band_job],
                output={"dos": nscf_dos_job.output, "band": nscf_band_job.output},
            )
        else:
            return Flow(
                [relax_job, static_job, nscf_dos_job],
                output=nscf_dos_job.output,
            )

    def get_result(
        self, output_job_name: str = "nscf uniform"
    ) -> dict[str, Any] | None:
        """
        Get the result from the specified job.

        Args:
            output_job_name: Name of the job to get result from.
                Options: "nscf uniform" (DOS), "nscf line" (Band)

        Returns:
            Output from the specified job or None if not found
        """
        return super().get_result(output_job_name)
    
    def get_dos_result(self) -> dict[str, Any] | None:
        return super().get_result("nscf uniform")
    def get_band_result(self) -> dict[str, Any] | None:
        return super().get_result("nscf line")