import json
import logging
import traceback
from pathlib import Path
from typing import Any, Literal

from jobflow import Flow, JobStore, Maker
from jobflow.managers.fireworks import flow_to_workflow
from maggma.stores import JSONStore, MemoryStore
from maggma.core import Store
from pymatgen.core import Structure

from htvasp.utils.run_locally import run_locally_custom

log = logging.getLogger(__name__)


class Worker:
    """
    HT-VASP - Worker Base Class

    Base class for HT-VASP workers. Supports both local execution via
    run_locally_custom and FireWorks-based distributed execution.
    """

    def __init__(
        self,
        worker_name: str,
        vasp_args: dict[str, Any],
        potcar_functional: Literal["PBE", "PBE_54", "PBE_64"] = "PBE_64",
        use_fireworks: bool = True,
        global_incar: dict[str, Any] | None = None,
        execution_mode: Literal["local", "fireworks"] = "local",
        fireworks_config_dir: Path | str | None = None,
        **kwargs,
    ):
        """
        Initialize the flow for the worker.

        Args:
            worker_name: Name of the worker
            vasp_args: VASP command and handler settings
            potcar_functional: POTCAR functional type
            global_incar: Global INCAR settings
            use_fireworks: Whether to use FireWorks for job management (legacy parameter)
            execution_mode: Execution mode - "local" or "fireworks"
            fireworks_config_dir: Directory containing FireWorks configuration files
        """
        # Store will be initialized in run_flow to allow custom paths
        self.store = None
        self.worker_name = worker_name
        self.execution_mode = execution_mode
        self.fireworks_config_dir = Path(fireworks_config_dir) if fireworks_config_dir else None
        self.fireworks_runner = None

        # Initialize FireWorks runner if in fireworks mode
        if execution_mode == "fireworks" and self.fireworks_config_dir:
            from htvasp.utils.fireworks_runner import FireWorksRunner

            lp_config = self.fireworks_config_dir / "mcmf_launchpad.yaml"
            qadapter_config = self.fireworks_config_dir / "mcmf_qadapter.yaml"

            if lp_config.exists():
                self.fireworks_runner = FireWorksRunner(
                    launchpad_config=lp_config,
                    qadapter_config=qadapter_config if qadapter_config.exists() else None,
                )
                log.info(f"FireWorks mode enabled with config dir: {self.fireworks_config_dir}")
            else:
                log.warning(f"LaunchPad config not found at {lp_config}, FireWorks mode may not work")

        # Default INCAR settings
        default_incar = {
            "ENCUT": 400,
            "ISTART": 0,
            "ICHARG": 2,
            # Electronic
            "ISMEAR": 1,
            "SIGMA": 0.2,
            "ALGO": "Fast",
            "NELM": 100,
            "NELMIN": 6,
            "NELMDL": -6,
            # Ionic
            "IBRION": 2,
            "ISIF": 2,
            "NSW": 10,
            "POTIM": 0.2,
            "EDIFF": 1e-5,
            "EDIFFG": -0.05,
            # Magnetic
            "ISPIN": 2,
            # Precision
            "ISYM": 0,
            "LREAL": "Auto",
            "PREC": "Normal",
            "SYMPREC": 1e-5,
            # Output
            "LWAVE": False,
            "LCHARG": False,
            "LORBIT": None,
            "LOPTICS": False,
            "LVTOT": False,
            "GGA": "PE",
        }
        # Override with custom settings
        if global_incar:
            default_incar.update(global_incar)
        self.global_incar = default_incar

        self.flow_maker = None
        self.flow = None  # Store the flow for later access

    def _make_flow(self, structure: Structure) -> Flow:
        raise NotImplementedError

    def close(self):
        """Close the store connection."""
        if hasattr(self, "store") and isinstance(self.store, Store):
            self.store.close()

    def run_flow(
        self,
        name: str,
        structure: Structure,
        flow_dir: Path | str,
        store_path: Path | str = "",
        ensure_success: bool = True,
        raise_immediately: bool = False,
        resume: bool = True,
        wait_for_completion: bool = False,
    ):
        """
        Run the flow for the structure.

        Args:
            name: Name of the system
            structure: Structure to run the worker on
            flow_dir: Flow directory
            store_path: Path to store the results (local mode only)
            ensure_success: Raise an error if the flow did not finish successfully
            raise_immediately: Raise an error immediately if a job fails
            resume: Resume from previous completed jobs (local mode only)
            wait_for_completion: Wait for FireWorks workflow to complete (fireworks mode only)
        """
        self.flow = self._make_flow(structure)
        flow_dir = Path(flow_dir)

        if self.execution_mode == "fireworks":
            self._run_fireworks_mode(name, flow_dir, wait_for_completion)
        else:
            self._run_local_mode(
                flow_dir,
                store_path,
                ensure_success,
                raise_immediately,
                resume,
            )

    def _run_local_mode(
        self,
        flow_dir: Path,
        store_path: Path | str,
        ensure_success: bool,
        raise_immediately: bool,
        resume: bool,
    ):
        """Execute flow using run_locally_custom (existing local behavior).

        Args:
            flow_dir: Directory for workflow execution
            store_path: Path to JSONStore file
            ensure_success: Raise error if flow fails
            raise_immediately: Raise error immediately on job failure
            resume: Resume from completed jobs
        """
        if not store_path:
            store_path = Path(flow_dir, "store.json").resolve()
        else:
            store_path = Path(store_path).resolve()

        self.store = JobStore(
            JSONStore(str(store_path), read_only=False),
            additional_stores={"data": MemoryStore()},
        )

        log.info(f"Running flow in LOCAL mode: {flow_dir}")
        try:
            run_locally_custom(
                self.flow,
                store=self.store,
                root_dir=flow_dir,
                ensure_success=ensure_success,
                raise_immediately=raise_immediately,
                resume=resume,
            )
        except Exception as e:
            log.critical(f"Flow execution failed: {e}")
            raise

    def _run_fireworks_mode(
        self,
        name: str,
        flow_dir: Path,
        wait_for_completion: bool = False,
    ):
        """Submit flow to FireWorks LaunchPad for distributed execution.

        Uses JSONStore + MemoryStore for simplicity in dual-node environments.

        Args:
            name: Workflow name identifier
            flow_dir: Directory for workflow metadata and store
            wait_for_completion: Block until workflow completes

        Raises:
            RuntimeError: If FireWorks runner is not initialized
        """
        if not self.fireworks_runner:
            raise RuntimeError(
                "FireWorks runner not initialized. "
                "Set execution_mode='fireworks' and provide fireworks_config_dir"
            )

        # Use JSONStore + MemoryStore (same as local mode)
        store_path = Path(flow_dir, "store.json").resolve()
        self.store = JobStore(
            JSONStore(str(store_path), read_only=False),
            additional_stores={"data": MemoryStore()},
        )

        workflow_name = f"{self.worker_name}-{name}"
        log.info(f"Submitting flow to FireWorks in FIREWORKS mode: {workflow_name}")

        try:
            wf_id = self.fireworks_runner.convert_and_submit(
                flow=self.flow,
                store=self.store,
                name=workflow_name,
            )

            # Save workflow metadata for later reference
            self.fireworks_runner.save_workflow_metadata(flow_dir, wf_id, name)

            log.info(f"Workflow submitted with wf_id: {wf_id}")

            if wait_for_completion:
                success = self.fireworks_runner.wait_for_completion(wf_id)
                if not success:
                    log.error(f"Workflow {wf_id} did not complete successfully")

        except Exception as e:
            log.critical(f"FireWorks submission failed: {e}")
            raise

    def get_result(self, output_job_name: str) -> dict[str, Any] | None:
        """
        Get the result from the specified job.

        Args:
            output_job_name: Name of the job to get the result from
        Returns:
            Output from the specified job or None if not found or failed
        """
        if not self.store:
            log.error("Store is not initialized. Run the flow first.")
            return None
        try:
            self.store.connect()
            job = self.store.query_one(
                criteria={"name": {"$regex": output_job_name}},
                properties=["uuid", "index", "name"],
                sort={"index": -1},
            )
            if not job:
                log.error(f"No job matching '{output_job_name}' found")
                return None
            output = self.store.get_output(uuid=job["uuid"], which="last", load=True)
            return output
        except Exception as e:
            log.error(f"Failed to get result for job '{output_job_name}': {e}")
            log.error(traceback.format_exc())
            return None
        finally:
            self.close()
