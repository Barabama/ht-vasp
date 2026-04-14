import json
import logging
import traceback
from pathlib import Path
from typing import Any, Literal

from jobflow import Flow, JobStore
from custodian.vasp.handlers import VaspErrorHandler
from maggma.stores import JSONStore, MemoryStore
from maggma.core import Store
from pymatgen.core import Structure

from htvasp.utils.local import run_locally_custom

log = logging.getLogger(__name__)


class Worker:
    """
    HT-VASP - Worker Base Class

    Base class for HT-VASP workers. Uses local execution via run_locally_custom.
    """

    def __init__(
        self,
        worker_name: str,
        vasp_args: dict[str, Any],
        potcar_functional = "PBE_64",
        global_incar: dict[str, Any] | None = None,
        **kwargs,
    ):
        """
        Initialize the flow for the worker.

        Args:
            worker_name: Name of the worker
            vasp_args: VASP command and handler settings
            potcar_functional: POTCAR functional type
            global_incar: Global INCAR settings
        """
        # Store will be initialized in run_flow to allow custom paths
        self.store = None
        self.worker_name = worker_name
        self.run_vasp_kwargs = {
            "handlers": [VaspErrorHandler()],
            "custodian_kwargs": {
                "max_errors_per_job": 3,
                "gzipped_output": False,
            },
            **vasp_args,
        }
        self.potcar_functional = potcar_functional

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
        flow_dir: Path | str = "/tmp",
        store_dir: Path | str = ".",
        ensure_success: bool = True,
        raise_immediately: bool = False,
        resume: bool = True,
    ):
        """
        Run the flow for the structure.

        Args:
            name: Name of the system
            structure: Structure to run the worker on
            flow_dir: Flow directory for running the flow (fast storage, e.g., /tmp)
            store_dir: Directory to store the results (persistent storage)
            ensure_success: Raise an error if the flow did not finish successfully
            raise_immediately: Raise an error immediately if a job fails
            resume: Resume from previous completed jobs
        """
        import shutil

        self.flow = self._make_flow(structure)

        # 生成目录路径
        flow_dir = Path(flow_dir) / f"{name}-{self.flow.uuid[:8]}"
        final_store_dir = Path(store_dir).resolve() / name

        # Resume: 从 store_dir 复制到 flow_dir
        if resume and final_store_dir.exists():
            log.info(f"Resuming from {final_store_dir}, copying to {flow_dir}")
            shutil.copytree(final_store_dir, flow_dir)

        # 初始化 store
        self.store = JobStore(
            JSONStore(str(flow_dir / "store.json"), read_only=False),
            additional_stores={"data": MemoryStore()},
        )

        log.info(f"Running flow in {flow_dir}")
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
        finally:
            # 无论成功失败，移动 flow_dir 到 store_dir
            try:
                if final_store_dir.exists():
                    shutil.rmtree(final_store_dir)
                shutil.move(str(flow_dir), str(final_store_dir))
                log.info(f"Moved results to {final_store_dir}")
            except Exception as move_error:
                log.error(f"Failed to move results to store_dir: {move_error}")

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
