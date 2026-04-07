"""
OJ Maker

Maker class for OstravaJ magnetic exchange calculations.
"""

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

from jobflow import Flow, Maker, Response, job
from pymatgen.core import Structure

from atomate2 import SETTINGS
from atomate2.common.files import gzip_output_folder
from atomate2.vasp.jobs.base import (
    _FILES_TO_ZIP,
    get_vasp_task_document,
    vasp_job,
)
from atomate2.vasp.run import run_vasp, should_stop_children
from atomate2.vasp.sets.core import StaticSetGenerator
from htvasp.oj.input_set import OJInputSetGenerator
from htvasp.oj.jobs import oj_generate, oj_solve

if TYPE_CHECKING:
    from atomate2.vasp.sets.base import VaspInputGenerator

log = logging.getLogger(__name__)


@dataclass
class FlipMaker(Maker):
    """
    Maker for running VASP in a single OJ flip directory.

    This maker is designed to work with pre-generated flip directories
    from the OJ generate step. It does NOT generate input files - it
    assumes POSCAR, INCAR, KPOINTS, POTCAR already exist in flip_dir.

    Parameters
    ----------
    name : str
        The job name.
    flip_dir : str | Path
        Path to the flip directory containing VASP input files.
    run_vasp_kwargs : dict
        Arguments for running VASP (vasp_cmd, handlers, etc.).
    task_document_kwargs : dict
        Arguments for TaskDoc parsing.
    stop_children_kwargs : dict
        Arguments for deciding whether to stop child jobs.
    write_additional_data : dict
        Additional data to write to the current directory.
    """

    name: str = "oj flip vasp"
    flip_dir: str | Path | None = None
    run_vasp_kwargs: dict[str, Any] = field(default_factory=dict)
    task_document_kwargs: dict[str, Any] = field(default_factory=dict)
    stop_children_kwargs: dict[str, Any] = field(
        default_factory=lambda: {"handle_unsuccessful": False}
    )
    write_additional_data: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate flip_dir is provided."""
        if self.flip_dir is None:
            raise ValueError("flip_dir must be specified for FlipMaker")

    @vasp_job
    def make(
        self, structure: Structure | None = None, prev_dir: str | Path | None = None
    ) -> Response:
        """
        Run VASP in the flip directory.

        Unlike typical BaseVaspMaker implementations, this method does NOT
        write input files (they already exist from oj_generate). It directly
        runs VASP in the flip_dir.

        Parameters
        ----------
        structure : Structure | None
            Not used (kept for API compatibility with BaseVaspMaker).
        prev_dir : str | Path | None
            Not used (kept for API compatibility).

        Returns
        -------
        Response
            Contains TaskDoc with VASP results and stop_children flag.
        """
        from monty.serialization import dumpfn

        flip_path = Path(self.flip_dir).resolve()

        if not flip_path.exists():
            raise FileNotFoundError(f"Flip directory does not exist: {flip_path}")

        # Change to flip directory for VASP execution
        original_dir = Path.cwd()
        os.chdir(flip_path)

        try:
            # Write any additional data if specified
            for filename, data in self.write_additional_data.items():
                dumpfn(data, filename.replace(":", "."))

            # Run VASP (input files already exist)
            run_vasp(**self.run_vasp_kwargs)

            # Parse VASP outputs into TaskDoc
            task_doc = get_vasp_task_document(Path.cwd(), **self.task_document_kwargs)
            task_doc.task_label = self.name

            # Decide whether child jobs should proceed
            stop_children = should_stop_children(task_doc, **self.stop_children_kwargs)

            # Gzip output folder
            gzip_output_folder(
                directory=Path.cwd(),
                setting=SETTINGS.VASP_ZIP_FILES,
                files_list=_FILES_TO_ZIP,
            )

            return Response(
                stop_children=stop_children,
                stored_data={"custodian": task_doc.custodian},
                output=task_doc,
            )
        finally:
            # Always restore original directory
            os.chdir(original_dir)


@job
def create_flip_jobs(
    flip_dirs: list[str],
    run_vasp_kwargs: dict[str, Any],
    stop_children_kwargs: dict[str, Any],
    task_document_kwargs: dict[str, Any] | None = None,
) -> Flow:
    """
    Dynamically create FlipMaker jobs for each flip directory.

    This job receives the list of flip directories from oj_generate and
    creates individual FlipMaker jobs for each one. Using a job allows
    us to access the actual flip_dirs list (not OutputReference).

    Parameters
    ----------
    flip_dirs : list[str]
        List of flip directory paths from oj_generate.
    run_vasp_kwargs : dict
        VASP execution parameters.
    stop_children_kwargs : dict
        Stop children parameters.
    task_document_kwargs : dict, optional
        Task document parsing parameters.

    Returns
    -------
    Flow
        Contains one FlipMaker job per flip directory.
    """
    if task_document_kwargs is None:
        task_document_kwargs = {}

    flip_jobs = []

    for i, flip_dir in enumerate(flip_dirs):
        flip_maker = FlipMaker(
            flip_dir=flip_dir,
            run_vasp_kwargs=run_vasp_kwargs,
            stop_children_kwargs=stop_children_kwargs,
            task_document_kwargs=task_document_kwargs,
            name=f"oj_flip_{i}",
        )

        # Create job using the maker
        # Note: We call make() with None structure since FlipMaker doesn't need it
        flip_job = flip_maker.make(structure=None)
        flip_job.name = f"oj_flip_{i}"
        flip_jobs.append(flip_job)

    # Return Flow containing all flip jobs
    # The output is a list of TaskDocs from each flip job
    return Flow(jobs=flip_jobs, output=[job.output for job in flip_jobs])


@dataclass
class OJMaker(Maker):
    """Maker for OstravaJ magnetic exchange calculations

    Generates a Flow that:
    1. Generates magnetic configurations
    2. Runs VASP for each configuration as independent jobs (parallel)
    3. Solves for exchange parameters J and Tc

    Note: OJ does not require pre-relaxation. The input structure can be
    used directly. If relaxation is needed, users should combine with
    atomate2's relax makers separately.

    Args:
        name: Maker name
        input_set_generator: OJInputSetGenerator instance (atomate2 style)
        run_vasp_kwargs: VASP execution parameters (includes handlers, vasp_cmd)
        stop_children_kwargs: Error handling settings

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
    stop_children_kwargs: dict[str, Any] = field(
        default_factory=lambda: {"handle_unsuccessful": False}
    )
    input_set_generator: OJInputSetGenerator = field(default_factory=OJInputSetGenerator)
    task_document_kwargs: dict[str, Any] = field(default_factory=dict)

    def make(self, structure: Structure) -> Flow:
        """Generate the OJ workflow with per-flip jobs.

        Creates a Flow where each flip directory calculation is an independent
        job, enabling parallel execution and individual retry capability.

        Args:
            structure: Input structure

        Returns:
            Flow containing: generate job -> N flip jobs (parallel) -> solve job
        """
        jobs = []

        # Step 1: Generate magnetic configurations
        gen_job = oj_generate(
            structure,
            input_set_generator=self.input_set_generator,
        )
        gen_job.name = f"{self.name}_generate"
        jobs.append(gen_job)

        # Step 2: Create flip jobs dynamically via intermediate job
        flip_flow_job = create_flip_jobs(
            gen_job.output["flip_dirs"],
            self.run_vasp_kwargs,
            self.stop_children_kwargs,
            self.task_document_kwargs,
        )
        flip_flow_job.name = f"{self.name}_create_flip_jobs"
        jobs.append(flip_flow_job)

        # Step 3: Solve using successful flip results
        # Extract TaskDocs from flip_flow_job.output (which is a list)
        solve_job = oj_solve(
            gen_job.output["run_dir"],
            flip_flow_job.output,  # List of TaskDocs
        )
        solve_job.name = f"{self.name}_solve"
        jobs.append(solve_job)

        return Flow(jobs, output=solve_job.output, name=self.name)
