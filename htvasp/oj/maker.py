"""
OJ Maker

Maker class for OstravaJ magnetic exchange calculations.
"""

import logging
from dataclasses import dataclass, field
from typing import Any

from jobflow import Flow, Maker

from pymatgen.core import Structure

from htvasp.oj.config import OJConfig
from htvasp.oj.jobs import oj_generate, oj_vasp, oj_solve, oj_workflow

log = logging.getLogger(__name__)


@dataclass
class OJMaker(Maker):
    """Maker for OstravaJ magnetic exchange calculations

    Generates a Flow that:
    1. Generates magnetic configurations
    2. Runs VASP for each configuration (serial)
    3. Solves for exchange parameters J and Tc

    Example:
        >>> from pymatgen.core import Structure
        >>> from htvasp.oj import OJMaker, OJConfig
        >>>
        >>> structure = Structure.from_file("POSCAR")
        >>> config = OJConfig(j_count=3, magnetic_ion_types=["Fe", "Co"])
        >>> maker = OJMaker(config=config)
        >>> flow = maker.make(structure)
    """

    name: str = "oj_exchange"
    config: OJConfig = field(default_factory=lambda: OJConfig(j_count=1))
    vasp_cmd: str = "vasp_std"
    workdir: str = "./oj_work"

    def make(self, structure: Structure) -> Flow:
        """Generate the OJ workflow"""
        gen_job = oj_generate(structure, self.config, self.workdir)
        gen_job.name = f"{self.name}_generate"

        vasp_job = oj_vasp(gen_job.output["flip_dirs"], self.vasp_cmd)
        vasp_job.name = f"{self.name}_vasp"

        solve_job = oj_solve(gen_job.output["run_dir"])
        solve_job.name = f"{self.name}_solve"

        return Flow([gen_job, vasp_job, solve_job], output=solve_job.output, name=self.name)


@dataclass
class OJSimpleMaker(Maker):
    """Simplified maker - runs everything in a single job

    Use for simple serial execution without Flow composition.
    """

    name: str = "oj_simple"
    config: OJConfig = field(default_factory=lambda: OJConfig(j_count=1))
    vasp_cmd: str = "vasp_std"
    workdir: str = "./oj_work"

    def make(self, structure: Structure) -> Flow:
        """Generate a single-job flow"""
        job = oj_workflow(structure, self.config, self.workdir, self.vasp_cmd)
        job.name = self.name
        return Flow([job], output=job.output, name=self.name)


@dataclass
class OJRelaxMaker(Maker):
    """Combined Maker: Relax + OJ exchange calculation

    Chains structure relaxation with OJ calculation.
    """

    name: str = "oj_relax"
    relax_maker: Maker | None = None
    config: OJConfig = field(default_factory=lambda: OJConfig(j_count=1))
    vasp_cmd: str = "vasp_std"
    workdir: str = "./oj_work"

    def make(self, structure: Structure) -> Flow:
        """Generate combined relax + OJ flow"""
        jobs = []

        if self.relax_maker is not None:
            relax_job = self.relax_maker.make(structure)
            relax_job.name = f"{self.name}_relax"
            jobs.append(relax_job)
            current_structure = relax_job.output.structure
        else:
            current_structure = structure

        gen_job = oj_generate(current_structure, self.config, f"{self.workdir}/oj")
        gen_job.name = f"{self.name}_generate"
        jobs.append(gen_job)

        vasp_job = oj_vasp(gen_job.output["flip_dirs"], self.vasp_cmd)
        vasp_job.name = f"{self.name}_vasp"
        jobs.append(vasp_job)

        solve_job = oj_solve(gen_job.output["run_dir"])
        solve_job.name = f"{self.name}_solve"
        jobs.append(solve_job)

        return Flow(jobs, output=solve_job.output, name=self.name)
