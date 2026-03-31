from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Literal

from maggma.core import Store
from pymatgen.core import Structure


class Worker(ABC):
    """
    HT-VASP - Worker Base Class

    Base class for HT-VASP workers.
    """

    @abstractmethod
    def __init__(
        self,
        worker_name: str,
        vasp_args: dict[str, Any],
        potcar_functional: Literal["PBE", "PBE_54", "PBE_64"] = "PBE_64",
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
        pass

    @abstractmethod
    def run_flow(
        self,
        name: str,
        structure: Structure,
        flowdir: Path | str,
        dir_format: str = "{name}",
        store_path: str | None = None,
    ) -> dict[str, Any] | None:
        """
        Run the worker.

        Args:
            name: Name of the system
            structure: Structure to run the worker on
            flowdir: Flow directory
            dir_format: Directory format
            store_path: Path to store the results

        Returns:
            Results or None if failed
        """
        pass

    def close(self):
        """Close the store connection."""
        if hasattr(self, "store") and isinstance(self.store, Store):
            self.store.close()
