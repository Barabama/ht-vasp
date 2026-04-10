"""
HT-VASP - FireWorks Integration Runner

Handles conversion and submission of jobflow Flows to FireWorks Workflows.
Provides monitoring and result retrieval capabilities for FireWorks-based execution.
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path
from typing import Any

from fireworks import LaunchPad
from jobflow import Flow, JobStore
from jobflow.managers.fireworks import flow_to_workflow

log = logging.getLogger(__name__)


class FireWorksRunner:
    """Manages FireWorks workflow submission and monitoring."""

    def __init__(
        self,
        launchpad_config: Path | str,
        qadapter_config: Path | str | None = None,
    ):
        """Initialize FireWorks runner.

        Args:
            launchpad_config: Path to LaunchPad YAML config file
            qadapter_config: Path to QueueAdapter YAML config (for Slurm integration)
        """
        self.lp = LaunchPad.from_file(str(launchpad_config))
        self.qadapter_config = Path(qadapter_config) if qadapter_config else None
        log.info(f"FireWorksRunner initialized with LaunchPad: {self.lp.host}")

    def convert_and_submit(
        self,
        flow: Flow,
        store: JobStore,
        name: str = "qha-workflow",
    ) -> int:
        """Convert jobflow Flow to FireWorks Workflow and submit to LaunchPad.

        Args:
            flow: jobflow Flow object containing the computational workflow
            store: JobStore instance for result storage (should be MongoDB-backed)
            name: Descriptive name for the workflow

        Returns:
            Workflow ID (wf_id) assigned by LaunchPad
        """
        # Convert jobflow Flow to FireWorks Workflow
        wf = flow_to_workflow(flow, store)
        wf.name = name

        # Submit to LaunchPad
        wf_id = self.lp.add_wf(wf)
        log.info(f"Submitted workflow '{name}' with wf_id: {wf_id}")

        return wf_id

    def get_workflow_status(self, wf_id: int) -> dict[str, Any]:
        """Get workflow execution status from LaunchPad.

        Args:
            wf_id: Workflow ID to query

        Returns:
            Dictionary containing workflow state and metadata
        """
        wf = self.lp.get_wf_by_fw_id_lzy(wf_id)
        if not wf:
            return {"error": f"Workflow {wf_id} not found"}

        return {
            "wf_id": wf_id,
            "name": wf.name,
            "state": wf.state,
            "created_on": wf.created_on.isoformat() if wf.created_on else None,
            "updated_on": wf.updated_on.isoformat() if wf.updated_on else None,
        }

    def wait_for_completion(
        self,
        wf_id: int,
        poll_interval: int = 30,
        timeout: int | None = None,
    ) -> bool:
        """Wait for workflow to complete with periodic status polling.

        Args:
            wf_id: Workflow ID to monitor
            poll_interval: Seconds between status checks
            timeout: Maximum wait time in seconds (None for infinite wait)

        Returns:
            True if workflow completed successfully, False otherwise
        """
        start_time = time.time()
        log.info(f"Waiting for workflow {wf_id} to complete (poll interval: {poll_interval}s)")

        while True:
            status = self.get_workflow_status(wf_id)
            state = status.get("state", "UNKNOWN")

            if state == "COMPLETED":
                log.info(f"Workflow {wf_id} completed successfully")
                return True
            elif state in ["FIZZLED", "DEFUSED", "ARCHIVED"]:
                log.error(f"Workflow {wf_id} ended with state: {state}")
                return False

            # Check timeout
            if timeout and (time.time() - start_time) > timeout:
                log.warning(f"Timeout waiting for workflow {wf_id}")
                return False

            time.sleep(poll_interval)

    def save_workflow_metadata(self, flow_dir: Path, wf_id: int, name: str) -> Path:
        """Save workflow metadata to JSON file for later reference.

        Args:
            flow_dir: Directory where workflow results are stored
            wf_id: Workflow ID from LaunchPad
            name: Workflow name

        Returns:
            Path to the saved metadata file
        """
        metadata_file = flow_dir / "fireworks_metadata.json"
        metadata_file.parent.mkdir(parents=True, exist_ok=True)

        metadata = {
            "wf_id": wf_id,
            "name": name,
            "submitted_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        }

        with open(metadata_file, "w") as f:
            json.dump(metadata, f, indent=2)

        log.info(f"Saved workflow metadata to {metadata_file}")
        return metadata_file

    def load_workflow_metadata(self, flow_dir: Path) -> dict[str, Any] | None:
        """Load workflow metadata from JSON file.

        Args:
            flow_dir: Directory where workflow results are stored

        Returns:
            Metadata dictionary or None if file doesn't exist
        """
        metadata_file = flow_dir / "fireworks_metadata.json"
        if not metadata_file.exists():
            return None

        with open(metadata_file) as f:
            return json.load(f)
