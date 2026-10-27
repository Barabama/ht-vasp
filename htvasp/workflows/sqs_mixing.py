"""SQS-mixing Worker for one structure."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from jobflow import Flow
from pymatgen.core import Structure

from htvasp.model.sqs import SqsCalculationSpec
from htvasp.sqs.maker import SqsCalculationMaker
from htvasp.workflows.base import Worker


class SqsMixingWorker(Worker):
    """Execute one SQS endpoint or point ``Relax -> Static`` Flow."""

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
        calculation_spec: SqsCalculationSpec | None = None,
        calculation_id: str | None = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(
            worker_name=worker_name,
            vasp_args=vasp_args,
            potcar_functional=potcar_functional,
            global_incar=global_incar,
        )
        if calculation_spec is not None and calculation_id is not None:
            raise ValueError("Pass calculation_spec or calculation_id, not both")
        self.calculation_spec = calculation_spec
        self.calculation_id = calculation_id or (
            calculation_spec.stable_id if calculation_spec is not None else worker_name
        )
        self.flow_maker = SqsCalculationMaker(
            name="sqs calculation",
            run_vasp_kwargs=self.run_vasp_kwargs,
            potcar_functional=self.potcar_functional,
            global_incar=self.global_incar,
            relax_incar=relax_incar or {},
            static_incar=static_incar or {},
            relax_reciprocal_density=relax_reciprocal_density,
            static_reciprocal_density=static_reciprocal_density,
        )

    def _make_flow(self, structure: Structure, prev_dir: Path | str | None = None) -> Flow:
        return self.flow_maker.make(
            structure,
            calculation_spec=self.calculation_spec,
            calculation_id=self.calculation_id,
            prev_dir=prev_dir,
        )

    def get_result(self, output_job_name: str = "sqs observation") -> dict[str, Any] | None:
        return super().get_result(output_job_name)


__all__ = ["SqsMixingWorker"]
