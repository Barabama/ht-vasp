"""
OJ Maker

Maker class for OstravaJ magnetic exchange calculations.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Literal

from jobflow import Flow, Maker
from pymatgen.core import Structure

from htvasp.oj.config import OJConfig
from htvasp.oj.input_set import OJInputSetGenerator
from htvasp.oj.jobs import oj_generate, oj_vasp, oj_solve, oj_workflow

log = logging.getLogger(__name__)


@dataclass
class OJMaker(Maker):
    """Maker for OstravaJ magnetic exchange calculations

    Generates a Flow that:
    1. Generates magnetic configurations
    2. Runs VASP for each configuration (serial)
    3. Solves for exchange parameters J and Tc

    Args:
        name: Maker name
        input_set_generator: OJInputSetGenerator instance (atomate2 style)
        config: OJ configuration (for backward compatibility)
        vasp_cmd: VASP command
        potcar_functional: POTCAR functional type
        j_count: Number of J pairs to consider
        dist_cutoff: Distance cutoff for J pairs
        magnetic_ion_types: List of magnetic ion types
        noncollinear: Whether to use noncollinear magnetism
        base_spin: Base spin value
        extend_poscar: Supercell extension factors
        reciprocal_density: K-point density for insulators
        reciprocal_density_metal: K-point density for metals
        auto_metal_kpoints: Automatically use higher density for metals
        force_gamma: Force gamma-centered k-point mesh
        user_incar_settings: Custom INCAR settings (atomate2 style)

    Example:
        >>> from pymatgen.core import Structure
        >>> from htvasp.oj import OJMaker, OJInputSetGenerator
        >>>
        >>> structure = Structure.from_file("POSCAR")
        >>> input_set_gen = OJInputSetGenerator(
        ...     j_count=3,
        ...     magnetic_ion_types=["Fe", "Co"]
        ... )
        >>> maker = OJMaker(input_set_generator=input_set_gen)
        >>> flow = maker.make(structure)

        Or with atomate2 compatibility:
        >>> from atomate2.vasp.jobs.core import StaticMaker
        >>> static_maker = StaticMaker()
        >>> oj_maker = OJMaker(j_count=3, magnetic_ion_types=["Fe"])
        >>> static_job = static_maker.make(structure)
        >>> oj_flow = oj_maker.make(static_job.output.structure)
        >>> combined_flow = Flow([static_job, *oj_flow.jobs])
    """

    name: str = "oj_exchange"
    input_set_generator: OJInputSetGenerator | None = field(default=None)
    j_count: int | None = 2
    dist_cutoff: float | None = None
    magnetic_ion_types: list[str] = field(default_factory=list)
    noncollinear: bool = False
    base_spin: float | list[float] = 1.0
    extend_poscar: tuple[int, int, int] = (2, 2, 2)
    reciprocal_density: float = 100
    reciprocal_density_metal: float = 400
    auto_metal_kpoints: bool = True
    force_gamma: bool = True
    user_incar_settings: dict[str, Any] = field(default_factory=dict)
    vasp_cmd: str = "vasp_std"
    potcar_functional: Literal["PBE", "PBE_54", "PBE_64"] = "PBE_64"
    config: OJConfig | None = field(default=None)

    def __post_init__(self):
        # If input_set_generator is provided, use it directly
        if self.input_set_generator is None:
            # Create input_set_generator from other parameters
            if self.config is not None:
                self.j_count = self.config.j_count
                self.dist_cutoff = self.config.dist_cutoff
                self.magnetic_ion_types = self.config.magnetic_ion_types
                self.noncollinear = self.config.noncollinear
                self.base_spin = self.config.base_spin
                self.extend_poscar = self.config.extend_poscar
                if hasattr(self.config, 'reciprocal_density'):
                    self.reciprocal_density = self.config.reciprocal_density
                if hasattr(self.config, 'reciprocal_density_metal'):
                    self.reciprocal_density_metal = self.config.reciprocal_density_metal
                if hasattr(self.config, 'auto_metal_kpoints'):
                    self.auto_metal_kpoints = self.config.auto_metal_kpoints
                if hasattr(self.config, 'force_gamma'):
                    self.force_gamma = self.config.force_gamma
                self.user_incar_settings = {**self.config.incar, **self.user_incar_settings}
            
            # Create input_set_generator
            self.input_set_generator = OJInputSetGenerator(
                j_count=self.j_count,
                dist_cutoff=self.dist_cutoff,
                magnetic_ion_types=self.magnetic_ion_types,
                noncollinear=self.noncollinear,
                base_spin=self.base_spin,
                extend_poscar=self.extend_poscar,
                reciprocal_density=self.reciprocal_density,
                reciprocal_density_metal=self.reciprocal_density_metal,
                auto_metal_kpoints=self.auto_metal_kpoints,
                force_gamma=self.force_gamma,
                user_incar_settings=self.user_incar_settings,
                user_potcar_functional=self.potcar_functional,
            )

    def make(self, structure: Structure) -> Flow:
        """Generate the OJ workflow

        Args:
            structure: Input structure

        Returns:
            Flow containing the workflow
        """
        jobs = []

        gen_job = oj_generate(
            structure,
            input_set_generator=self.input_set_generator,
        )
        gen_job.name = f"{self.name}_generate"
        jobs.append(gen_job)

        vasp_job = oj_vasp(gen_job.output["flip_dirs"], self.vasp_cmd)
        vasp_job.name = f"{self.name}_vasp"
        jobs.append(vasp_job)

        solve_job = oj_solve(gen_job.output["run_dir"], vasp_job.output["results"])
        solve_job.name = f"{self.name}_solve"
        jobs.append(solve_job)

        return Flow(jobs, output=solve_job.output, name=self.name)


@dataclass
class OJSimpleMaker(Maker):
    """Simplified maker - runs everything in a single job

    Use for simple serial execution without Flow composition.
    All steps run in one job directory.

    Args:
        name: Maker name
        config: OJ configuration
        vasp_cmd: VASP command
        potcar_functional: POTCAR functional type
    """

    name: str = "oj_simple"
    config: OJConfig = field(default_factory=lambda: OJConfig(j_count=1))
    vasp_cmd: str = "vasp_std"
    potcar_functional: Literal["PBE", "PBE_54", "PBE_64"] = "PBE_64"

    def make(self, structure: Structure) -> Flow:
        """Generate a single-job flow

        Args:
            structure: Input structure

        Returns:
            Flow containing a single job
        """
        job_obj = oj_workflow(
            structure,
            self.config,
            self.vasp_cmd,
            self.potcar_functional,
        )
        job_obj.name = self.name
        return Flow([job_obj], output=job_obj.output, name=self.name)


@dataclass
class OJRelaxMaker(Maker):
    """Combined Maker: Relax + OJ exchange calculation

    Chains structure relaxation with OJ calculation.

    Args:
        name: Maker name
        relax_maker: Optional relax maker (atomate2 Maker)
        config: OJ configuration
        vasp_cmd: VASP command
        potcar_functional: POTCAR functional type
    """

    name: str = "oj_relax"
    relax_maker: Maker | None = None
    config: OJConfig = field(default_factory=lambda: OJConfig(j_count=1))
    vasp_cmd: str = "vasp_std"
    potcar_functional: Literal["PBE", "PBE_54", "PBE_64"] = "PBE_64"

    def make(self, structure: Structure) -> Flow:
        """Generate combined relax + OJ flow

        Args:
            structure: Input structure

        Returns:
            Flow containing the workflow
        """
        jobs = []

        if self.relax_maker is not None:
            relax_job = self.relax_maker.make(structure)
            relax_job.name = f"{self.name}_relax"
            jobs.append(relax_job)
            current_structure = relax_job.output.structure
        else:
            current_structure = structure

        gen_job = oj_generate(current_structure, self.config, self.potcar_functional)
        gen_job.name = f"{self.name}_generate"
        jobs.append(gen_job)

        vasp_job = oj_vasp(gen_job.output["flip_dirs"], self.vasp_cmd)
        vasp_job.name = f"{self.name}_vasp"
        jobs.append(vasp_job)

        solve_job = oj_solve(gen_job.output["run_dir"], vasp_job.output["results"])
        solve_job.name = f"{self.name}_solve"
        jobs.append(solve_job)

        return Flow(jobs, output=solve_job.output, name=self.name)
