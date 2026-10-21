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

from htvasp.dfpt_qha import names
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
        magmom: dict[str, float] | list[float] | None = None,
        incar_overrides: dict[str, Any] | None = None,
        ncore: int | None = None,
        kpar: int | None = None,
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
        # NCORE and KPAR are mutually exclusive in VASP, and the two stages use
        # different ones: the relaxation family carries NCORE and the DFPT
        # family carries KPAR.  Taking them as separate arguments keeps a
        # caller from putting both in one override dict and producing an INCAR
        # VASP rejects.
        relax_overrides: dict[str, Any] = {}
        dfpt_overrides: dict[str, Any] = {}
        if ncore is not None:
            if int(ncore) < 1:
                raise ValueError(f"ncore must be positive, got {ncore}")
            relax_overrides["NCORE"] = int(ncore)
        if kpar is not None:
            if int(kpar) < 1:
                raise ValueError(f"kpar must be positive, got {kpar}")
            dfpt_overrides["KPAR"] = int(kpar)
        self.flow_maker = DfptQhaMaker(
            name=worker_name,
            recipe=recipe,
            linear_strain=linear_strain,
            number_of_frames=number_of_frames,
            run_vasp_kwargs=dict(self.run_vasp_kwargs),
            magmom=dict(magmom) if magmom is not None else None,
            incar_overrides=dict(incar_overrides or {}),
            relax_incar_overrides=relax_overrides,
            dfpt_incar_overrides=dfpt_overrides,
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
        self, output_job_name: str = names.ANALYSIS
    ) -> dict[str, Any] | None:
        """The analysis if it exists, otherwise the stability verdict behind it.

        A stopped flow has no analysis job, and returning ``None`` there would
        report a failed endpoint with no reason.  The classification job is the
        gate's own record, so it is the fallback: a soft volume and a failed
        integrity check are both scientific outcomes the caller has to see.
        """

        result = self._query_store(output_job_name)
        if result is None:
            # Only a FAILED gate is a result.  A gate that passed means the
            # analysis should exist; if it does not, the flow died somewhere
            # else and returning the gate document would be recorded as a
            # completed endpoint with no V0, B0 or Gibbs curve.
            gate = self._query_stability_gate()
            if gate is not None and gate.get("ok") is False:
                result = gate
        if result is None:
            result = self._query_validation_failure()
        self.close()
        return result

    def _query_validation_failure(self) -> dict[str, Any] | None:
        """A volume rejected before the stability job: report which and why.

        The validator stops the flow when a DFPT document does not match the
        recipe (wrong atom count, unusable shape).  That is a different failure
        from a soft mode, and without this the endpoint would report ``None``.
        """

        if not self.store:
            return None
        try:
            self.store.connect()
            jobs = self.store.query(
                criteria={"name": {"$regex": names.VALIDATION}},
                properties=["uuid", "index", "name"],
                sort={"index": 1},
            )
            failed = []
            for job in jobs:
                output = self.store.get_output(
                    uuid=job["uuid"], which="last", load=True
                )
                if isinstance(output, dict) and output.get("ok") is False:
                    failed.append(
                        {"job_name": job.get("name"), "error": output.get("error")}
                    )
            if not failed:
                return None
            return {
                "gate": "dfpt_validation",
                "ok": False,
                "n_failed": len(failed),
                "failures": failed,
            }
        except Exception as exc:  # noqa: BLE001 - a missing store is reported
            log.error(f"Failed to read the DFPT validation gate: {exc}")
            log.error(traceback.format_exc())
            return None

    def _query_stability_gate(self) -> dict[str, Any] | None:
        """The classification job's stored output, if the flow reached it.

        A flow halted at the stability gate has no analysis job, and returning
        ``None`` there would report a failed endpoint with no reason.  The
        verdict -- status, exact-q frequencies, eigenvector subspace -- is the
        deliverable in that case.
        """

        if not self.store:
            log.error("Store is not initialized. Run the flow first.")
            return None
        try:
            self.store.connect()
            job = self.store.query_one(
                criteria={"name": {"$regex": names.STABILITY}},
                properties=["uuid", "index", "name"],
                sort={"index": -1},
            )
            if not job:
                log.error("No job matching %r found", names.STABILITY)
                return None
            return self.store.get_output(uuid=job["uuid"], which="last", load=True)
        except Exception as exc:  # noqa: BLE001 - a missing store is reported
            log.error(f"Failed to read the DFPT stability gate: {exc}")
            log.error(traceback.format_exc())
            return None


__all__ = ["DfptQhaWorker"]
