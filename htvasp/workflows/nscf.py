"""
HT-VASP - NSCF Workflows

DOS + Band electronic structure calculations.
"""

import logging
from pathlib import Path
from typing import Any

from jobflow import Flow
from pymatgen.core import Structure
from atomate2.vasp.jobs.core import NonSCFMaker
from atomate2.vasp.sets.core import NonSCFSetGenerator

from htvasp.workflows.base import Worker

log = logging.getLogger(__name__)


class NscfWorker(Worker):
    """
    Worker for DOS + Band electronic structure calculations.

    Workflow: NSCF-DOS (uniform mesh) and optional NSCF-Band (line mode)
    """

    def __init__(
        self,
        worker_name: str = "nscf-worker",
        vasp_args: dict[str, Any] | None = None,
        potcar_functional: str = "PBE",
        global_incar: dict[str, Any] | None = None,
        nscf_dos_incar: dict[str, Any] | None = None,
        nscf_band_incar: dict[str, Any] | None = None,
        compute_band: bool = True,
        reciprocal_density: int = 200,
        line_density: int = 20,
        dedos: float = 0.02,
        band_kpath_kwargs: dict[str, Any] | None = None,
        **kwargs,
    ):
        """
        Initialize NscfWorker.

        Args:
            worker_name: Name of the worker
            vasp_args: VASP command and handler settings
            potcar_functional: POTCAR functional type
            global_incar: Global INCAR settings
            nscf_dos_incar: NSCF-DOS-specific INCAR settings
            nscf_band_incar: NSCF-Band-specific INCAR settings
            compute_band: Whether to compute band structure
            reciprocal_density: K-point density for DOS (uniform mode)
            line_density: Line density for band structure (line mode)
            dedos: Energy resolution for DOS (eV), used to auto-calculate NEDOS
            band_kpath_kwargs: Extra kwargs for HighSymmKpath in line mode,
                e.g. {"path_type": "hinuma"} to use SeeK-path for slabs.
            **kwargs: Additional keyword arguments
        """
        # Initialize base Worker
        super().__init__(
            worker_name=worker_name,
            vasp_args=vasp_args,
            potcar_functional=potcar_functional,
            global_incar=global_incar,
        )

        # # Disable VaspErrorHandler for NSCF calculations
        # self.run_vasp_kwargs["handlers"] = []

        nscf_dos_incar = nscf_dos_incar or {}
        nscf_band_incar = nscf_band_incar or {}

        self.compute_band = compute_band

        # NSCF-DOS (uniform mesh for density of states)
        nscf_dos_maker = NonSCFMaker(
            name="nscf uniform",
            run_vasp_kwargs=self.run_vasp_kwargs,
            stop_children_kwargs={"handle_unsuccessful": False},
            copy_vasp_kwargs={"additional_vasp_files": ("CHGCAR",)},
            input_set_generator=NonSCFSetGenerator(
                mode="uniform",
                reciprocal_density=reciprocal_density,
                dedos=dedos,
                user_potcar_functional=self.potcar_functional,
                user_incar_settings={
                    **self.global_incar,
                    "ISTART": 1,
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
            stop_children_kwargs={"handle_unsuccessful": False},
            copy_vasp_kwargs={"additional_vasp_files": ("CHGCAR",)},
            input_set_generator=NonSCFSetGenerator(
                mode="line",
                line_density=line_density,
                user_potcar_functional=self.potcar_functional,
                user_kpoints_settings={
                    "line_density": line_density,
                    "kpath_kwargs": band_kpath_kwargs or {},
                },
                user_incar_settings={
                    **self.global_incar,
                    "ISTART": 1,
                    "IBRION": -1,
                    "NSW": 0,
                    "ICHARG": 11,
                    "LORBIT": 11,
                    **nscf_band_incar,
                },
            ),
        )

        self.flow_makers: tuple[NonSCFMaker, NonSCFMaker] = (
            nscf_dos_maker,
            nscf_band_maker,
        )

    def _make_flow(self, structure: Structure, prev_dir: Path | str | None = None) -> Flow:
        nscf_dos_maker, nscf_band_maker = self.flow_makers

        # Create jobs
        nscf_dos_job = nscf_dos_maker.make(structure, prev_dir=prev_dir)

        if self.compute_band:
            nscf_band_job = nscf_band_maker.make(
                structure,
                prev_dir=prev_dir,
                mode="line",  #  avoiding default "uniform"
            )
            # Both NSCF jobs depend on static, but not on each other
            return Flow(
                [nscf_dos_job, nscf_band_job],
                output=nscf_band_job.output,
            )
        else:
            return Flow(
                [nscf_dos_job],
                output=nscf_dos_job.output,
            )

    def get_result(self, output_job_name: str = "nscf uniform") -> dict[str, Any] | None:
        return super().get_result(output_job_name)

    def get_dos_result(self) -> dict[str, Any] | None:
        return super().get_result("nscf uniform")

    def get_band_result(self) -> dict[str, Any] | None:
        return super().get_result("nscf line")
