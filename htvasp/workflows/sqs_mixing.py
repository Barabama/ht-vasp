"""Single-structure SQS mixing worker: one relax followed by one static run."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from atomate2.vasp.jobs.core import RelaxMaker, StaticMaker
from atomate2.vasp.sets.core import RelaxSetGenerator, StaticSetGenerator
from jobflow import Flow
from pymatgen.core import Structure

from htvasp.workflows.base import Worker


class SqsMixingWorker(Worker):
    """Build a dedicated Relax -> Static flow for one SQS or reference structure.

    The worker deliberately has no SQS naming, composition, generator, or
    mixing-energy logic.  Those concerns belong to ``htvasp.model.sqs``.
    """

    def __init__(
        self,
        worker_name: str = "sqs-mixing-worker",
        vasp_args: dict[str, Any] | None = None,
        potcar_functional: str = "PBE",
        global_incar: dict[str, Any] | None = None,
        relax_incar: dict[str, Any] | None = None,
        static_incar: dict[str, Any] | None = None,
        relax_reciprocal_density: int | None = None,
        static_reciprocal_density: int | None = None,
        **kwargs: Any,
    ):
        super().__init__(
            worker_name=worker_name,
            vasp_args=vasp_args,
            potcar_functional=potcar_functional,
            global_incar=global_incar,
        )
        relax_incar = relax_incar or {}
        static_incar = static_incar or {}
        relax_kpoints = (
            {"reciprocal_density": relax_reciprocal_density}
            if relax_reciprocal_density is not None
            else {}
        )
        static_kpoints = (
            {"reciprocal_density": static_reciprocal_density}
            if static_reciprocal_density is not None
            else {}
        )
        self.relax_maker = RelaxMaker(
            name="sqs mixing relax",
            run_vasp_kwargs=self.run_vasp_kwargs,
            stop_children_kwargs={"handle_unsuccessful": False},
            copy_vasp_kwargs={"additional_vasp_files": ("WAVECAR",)},
            input_set_generator=RelaxSetGenerator(
                user_potcar_functional=self.potcar_functional,
                user_kpoints_settings=relax_kpoints,
                user_incar_settings={
                    **self.global_incar,
                    "ISTART": 1,
                    "ISIF": 3,
                    "LWAVE": True,
                    **relax_incar,
                },
            ),
        )
        self.static_maker = StaticMaker(
            name="sqs mixing static",
            run_vasp_kwargs=self.run_vasp_kwargs,
            stop_children_kwargs={"handle_unsuccessful": False},
            copy_vasp_kwargs={"additional_vasp_files": ("WAVECAR",)},
            input_set_generator=StaticSetGenerator(
                user_potcar_functional=self.potcar_functional,
                user_kpoints_settings=static_kpoints,
                user_incar_settings={
                    **self.global_incar,
                    "ISTART": 1,
                    "IBRION": -1,
                    "ISIF": 2,
                    "NSW": 0,
                    # Static SQS output is consumed for structure and energy;
                    # keep the Relax WAVECAR as input but do not write a new
                    # wavefunction or charge-density file.
                    "LWAVE": False,
                    "LCHARG": False,
                    "LAECHG": False,
                    **static_incar,
                },
            ),
        )

    def _make_flow(self, structure: Structure, prev_dir: Path | str | None = None) -> Flow:
        relax_job = self.relax_maker.make(structure, prev_dir=prev_dir)
        static_job = self.static_maker.make(
            relax_job.output.structure,
            prev_dir=relax_job.output.dir_name,
        )
        return Flow([relax_job, static_job], output=static_job.output)

    def get_result(self, output_job_name: str = "sqs mixing static") -> dict[str, Any] | None:
        return super().get_result(output_job_name)


__all__ = ["SqsMixingWorker"]
