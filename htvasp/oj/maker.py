"""
OJ Maker

Maker class for OstravaJ magnetic exchange calculations.
"""

import logging
from dataclasses import dataclass, field
from typing import Any, Literal

from jobflow import Flow, Maker
from pymatgen.core import Structure

from htvasp.oj.input_set import OJInputSetGenerator
from htvasp.oj.jobs import oj_generate, oj_vasp, oj_solve

log = logging.getLogger(__name__)


@dataclass
class OJMaker(Maker):
    """Maker for OstravaJ magnetic exchange calculations

    Generates a Flow that:
    1. Generates magnetic configurations
    2. Runs VASP for each configuration (serial)
    3. Solves for exchange parameters J and Tc

    Note: OJ does not require pre-relaxation. The input structure can be
    used directly. If relaxation is needed, users should combine with
    atomate2's relax makers separately.

    Args:
        name: Maker name
        input_set_generator: OJInputSetGenerator instance (atomate2 style)
        vasp_cmd: VASP command

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

        Combined with relaxation (if needed):
        >>> from atomate2.vasp.jobs.core import RelaxMaker
        >>> relax_maker = RelaxMaker()
        >>> relax_job = relax_maker.make(structure)
        >>> oj_flow = oj_maker.make(relax_job.output.structure)
        >>> combined_flow = Flow([relax_job, *oj_flow.jobs])
    """

    name: str = "oj exchange"
    run_vasp_kwargs: dict[str, Any] = field(default_factory=dict)   
    stop_children_kwargs: dict[str, Any] = field(default_factory=lambda: {"handle_unsuccessful": False})
    input_set_generator: OJInputSetGenerator = field(default_factory=OJInputSetGenerator)


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

        vasp_job = oj_vasp(gen_job.output["flip_dirs"], self.run_vasp_kwargs["vasp_cmd"])
        vasp_job.name = f"{self.name}_vasp"
        jobs.append(vasp_job)

        solve_job = oj_solve(gen_job.output["run_dir"], vasp_job.output["results"])
        solve_job.name = f"{self.name}_solve"
        jobs.append(solve_job)

        return Flow(jobs, output=solve_job.output, name=self.name)
