"""DFPT-QHA Worker: the production wrapper around :class:`DfptQhaMaker`.

The Worker owns the run lifecycle only -- it builds the Maker, runs the flow
through ``run_locally_custom`` (with resume), and returns the terminal analysis
job's output.  Everything scientific lives in the Maker and its children.
"""

from __future__ import annotations

import traceback
from pathlib import Path
from typing import Any

from jobflow import Flow
from pymatgen.core import Structure

from htvasp.dfpt_qha.maker import DfptQhaMaker
from htvasp.model.dfpt_qha_recipe import Recipe
from htvasp.utils.vasp_handlers import dfpt_handlers
from htvasp.workflows.base import Worker, log


class DfptQhaWorker(Worker):
    """Execute one endpoint's DFPT-QHA flow."""

    def __init__(
        self,
        recipe: Recipe,
        worker_name: str = "dfpt-qha-worker",
        vasp_args: dict[str, Any] | None = None,
        potcar_functional: str = "PBE",
        global_incar: dict[str, Any] | None = None,
        linear_strain: tuple[float, float] = (-0.05, 0.05),
        number_of_frames: int = 6,
        supercell_matrix: Any | None = None,
        mesh: tuple[int, int, int] = (12, 12, 12),
        temperatures: list[float] | None = None,
        pressure_gpa: float = 0.0,
        **kwargs: Any,
    ) -> None:
        super().__init__(
            worker_name=worker_name,
            vasp_args=vasp_args,
            potcar_functional=potcar_functional,
            global_incar=global_incar,
        )
        self.recipe = recipe
        self.supercell_matrix = supercell_matrix
        # DFPT adds one detector on top of the base set: VASP aborts a
        # dynamical-matrix run when the perturbed cell would need a different
        # k-point set, and the only "fixes" would change the mesh or the
        # parallelisation -- both scientific inputs.  So it fails loudly instead.
        self.run_vasp_kwargs["handlers"] = (
            list(self.run_vasp_kwargs["handlers"]) + dfpt_handlers()
        )
        self.flow_maker = DfptQhaMaker(
            name=worker_name,
            recipe=recipe,
            linear_strain=linear_strain,
            number_of_frames=number_of_frames,
            run_vasp_kwargs=dict(self.run_vasp_kwargs),
            mesh=mesh,
            temperatures=temperatures,
            pressure_gpa=pressure_gpa,
        )

    def _make_flow(
        self, structure: Structure, prev_dir: Path | str | None = None
    ) -> Flow:
        return self.flow_maker.make(
            structure,
            supercell_matrix=self.supercell_matrix,
            prev_dir=prev_dir,
        )

    def get_result(
        self, output_job_name: str = "analyze DFPT-QHA"
    ) -> dict[str, Any] | None:
        """The analysis if it exists, otherwise the stability verdict behind it.

        A stopped flow has no analysis job, and returning ``None`` there would
        report a failed endpoint with no reason.  The classification job is the
        gate's own record, so it is the fallback: a soft volume and a failed
        integrity check are both scientific outcomes the caller has to see.
        """

        result = self._query_store(output_job_name)
        if result is None:
            result = self._query_stability_gate()
        self.close()
        return result

    def _query_stability_gate(self) -> dict[str, Any] | None:
        if not self.store:
            log.error("Store is not initialized. Run the flow first.")
            return None
        try:
            self.store.connect()
            jobs = self.store.query(
                criteria={"name": {"$regex": "classify DFPT volume"}},
                properties=["uuid", "index", "name"],
                sort={"index": 1},
            )
            reports: list[dict[str, Any]] = []
            for job in jobs:
                output = self.store.get_output(
                    uuid=job["uuid"], which="last", load=True
                )
                if not isinstance(output, dict):
                    continue
                reports.append(
                    {
                        "job_name": job.get("name"),
                        "index": job.get("index"),
                        "status": output.get("status"),
                        "qha_allowed": output.get("qha_allowed"),
                        "reason": output.get("reason"),
                        "report": output,
                    }
                )
        except Exception as exc:  # noqa: BLE001 - a missing store is reported
            log.error(f"Failed to read the DFPT stability gate: {exc}")
            log.error(traceback.format_exc())
            return None
        if not reports:
            log.error("No job matching 'classify DFPT volume' found")
            return None
        blocked = [r for r in reports if not r.get("qha_allowed")]
        chosen = blocked[0] if blocked else reports[-1]
        return {
            "gate": "dfpt_stability",
            "n_volumes_classified": len(reports),
            "n_volumes_blocked": len(blocked),
            "blocked_volume_indices": [r["index"] for r in blocked],
            "volumes": reports,
            "report": chosen["report"],
        }


__all__ = ["DfptQhaWorker"]
