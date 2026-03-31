"""
HT-VASP - OJ Workflows

Magnetic exchange calculation workflows using OstravaJ.
"""

import logging
import traceback
from pathlib import Path
from typing import Any, Literal

from pymatgen.core import Structure
from maggma.stores import JSONStore, MemoryStore
from jobflow.core.store import JobStore

from htvasp.workflows.base import Worker
from htvasp.oj.config import OJConfig
from htvasp.oj.maker import OJMaker
from htvasp.utils import run_locally_custom

log = logging.getLogger(__name__)


class OJWorker(Worker):
    """Worker for magnetic exchange calculations using OstravaJ
    
    Performs magnetic exchange interaction calculations:
    1. Generates magnetic configurations
    2. Runs VASP for each configuration
    3. Solves for exchange parameters J and Tc
    
    Example:
        >>> from pymatgen.core import Structure
        >>> from htvasp.workflows import OJWorker
        >>> from htvasp.oj import OJConfig
        >>> 
        >>> structure = Structure.from_file("POSCAR")
        >>> config = OJConfig(j_count=3, magnetic_ion_types=["Fe", "Co"])
        >>> worker = OJWorker(
        ...     worker_name="fe_co",
        ...     config=config,
        ... )
        >>> output = worker.run_flow("FeCo", structure, "./oj_work")
    """

    def __init__(
        self,
        worker_name: str,
        config: OJConfig | None = None,
        vasp_args: dict[str, Any] | None = None,
        potcar_functional: Literal["PBE", "PBE_54", "PBE_64"] = "PBE_64",
        **kwargs,
    ):
        self.store = None
        self.worker_name = worker_name
        self.config = config or OJConfig()
        self.vasp_args = vasp_args or {}
        self.potcar_functional = potcar_functional

        vasp_cmd = self.vasp_args.get("vasp_cmd", "vasp_std")
        workdir = self.vasp_args.get("workdir", "./oj_work")
        
        self._maker = OJMaker(
            name=f"{worker_name}_oj",
            config=self.config,
            vasp_cmd=vasp_cmd,
            workdir=workdir,
        )

    def run_flow(
        self,
        name: str,
        structure: Structure,
        flow_dir: Path | str,
        dir_format: str = "{name}",
        store_path: Path | str = "",
    ) -> dict[str, Any] | None:
        """Run the OJ workflow
        
        Args:
            name: Structure name
            structure: Input structure
            flow_dir: Flow directory
            dir_format: Directory format for job folders
            store_path: Path to store results
            
        Returns:
            Dictionary with J_reprs, Js, Tc_MFA, Tc_RPA, etc.
        """
        flow_dir = Path(flow_dir)
        flow_dir.mkdir(parents=True, exist_ok=True)

        if not store_path:
            store_path = Path(flow_dir, "store.json").resolve()
        else:
            store_path = Path(store_path).resolve()
        
        self.store = JobStore(
            JSONStore(store_path, read_only=False),
            additional_stores={"data": MemoryStore()},
        )

        flow = self._maker.make(structure)
        flow.name = name
        
        log.info(f"Running OJ flow for {name} in {flow_dir}")

        try:
            run_locally_custom(
                flow,
                store=self.store,
                root_dir=flow_dir,
                dir_format=dir_format,
            )

            self.store.connect()
            job = self.store.query_one(
                criteria={"name": {"$regex": "solve"}},
                properties=["uuid", "index", "name"],
                sort={"index": -1},
            )
            
            if not job:
                raise ValueError(f"No 'solve' job found in store {store_path}")

            output = self.store.get_output(uuid=job["uuid"], which="last", load=True)
            log.info(f"OJ flow for {name} completed successfully")
            return output

        except Exception as e:
            log.error(f"OJ flow for {name} failed: {e}")
            log.error(traceback.format_exc())
            return None

        finally:
            self.close()

    def create_flow(self, structure: Structure, name: str = "oj") -> tuple:
        """Create a Flow for external execution
        
        Returns:
            Tuple of (flow, maker) for use with run_locally or FireWorks
        """
        flow = self._maker.make(structure)
        flow.name = name
        return flow, self._maker
