"""OJ Maker - Maker class for OstravaJ magnetic exchange calculations."""

import logging
from dataclasses import dataclass, field
from typing import Any

from jobflow import Flow, Maker
from pymatgen.core import Structure

from htvasp.oj.input_set import OJInputSetGenerator
from htvasp.oj.jobs import create_flip_jobs, oj_generate, oj_solve

log = logging.getLogger(__name__)


@dataclass
class OJMaker(Maker):
    """Maker for OstravaJ magnetic exchange calculations.

    Generates a Flow that:
    1. Generates magnetic configurations
    2. Runs VASP for each configuration as independent jobs (parallel)
    3. Solves for exchange parameters J and Tc

    Args:
        name: Maker name
        input_set_generator: OJInputSetGenerator instance
        run_vasp_kwargs: VASP execution parameters
    """

    name: str = "oj exchange"
    run_vasp_kwargs: dict[str, Any] = field(default_factory=dict)
    input_set_generator: OJInputSetGenerator = field(default_factory=OJInputSetGenerator)

    def make(self, structure: Structure) -> Flow:
        """Generate the OJ workflow with per-flip jobs.

        Args:
            structure: Input structure

        Returns:
            Flow containing: generate job -> N flip jobs (parallel) -> solve job
        """
        # Step 1: Generate magnetic configurations
        gen_job = oj_generate(structure, input_set_generator=self.input_set_generator)
        gen_job.name = f"{self.name}_generate"

        # Step 2: Create flip jobs dynamically
        flip_flow_job = create_flip_jobs(
            gen_job.output["flip_dirs"],
            self.run_vasp_kwargs,
            self.input_set_generator,
        )
        flip_flow_job.name = f"{self.name}_create_flip_jobs"

        # Step 3: Solve using successful flip results
        solve_job = oj_solve(flip_flow_job.output, num_flips=gen_job.output["num_configs"])
        solve_job.name = f"{self.name}_solve"

        return Flow(
            [gen_job, flip_flow_job, solve_job],
            output=solve_job.output,
            name=self.name,
        )
