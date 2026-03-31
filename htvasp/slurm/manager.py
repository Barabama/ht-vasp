"""
HT-VASP - Slurm Job Management

Python interface for submitting and managing VASP jobs on Slurm clusters.
Supports CPU/GPU partition switching.
"""

import os
import logging
import subprocess
from pathlib import Path
from typing import Any, Literal
from dataclasses import dataclass, field

log = logging.getLogger(__name__)


@dataclass
class SlurmConfig:
    """
    Slurm job configuration.

    Supports both CPU and GPU partitions with appropriate defaults.

    Example:
        >>> # CPU job
        >>> config = SlurmConfig(
        ...     job_name="relax-cpu",
        ...     partition="partCPU",
        ...     tasks_per_node=44,
        ...     memory="80G"
        ... )

        >>> # GPU job
        >>> config = SlurmConfig(
        ...     job_name="relax-gpu",
        ...     partition="partGPU",
        ...     tasks_per_node=1,
        ...     gpus_per_task=1,
        ...     memory="20G"
        ... )
    """

    # Basic settings
    job_name: str = "vasp-job"
    output_log: str = "job.log"  # Default to job.log
    error_log: str = ""  # Empty = don't use --error unless specified
    memory: str = "20G"
    nodes: int = 1
    ntasks: int = 1
    cpus_per_task: int = 1
    gpus_per_task: int = 0  # Only for GPU partition
    partition: Literal["partCPU", "partGPU"] = "partCPU"
    time_limit: str = "100:00:00"

    # VASP commands
    vasp_cmd: str = "srun vasp_std"
    vasp_gam_cmd: str = "srun vasp_gam"

    # Environment settings
    # conda_path: str = "/nfs_ssd/.conda"
    # conda_path: str = "/opt/miniconda3"

    # Extra commands
    extra_commands: list[str] = field(default_factory=list)

    def init_module(self, module_name: str) -> list[str]:
        """
        Initialize module commands.

        Returns:
            List of module load commands
        """
        # Use '.' instead of 'source' for better compatibility with /bin/sh
        return [
            f". /etc/profile.d/modules.sh",
            f"ulimit -s unlimited",
            f"module purge",
            f"module load {module_name}",
        ]

    def init_conda(self, conda_env: str) -> list[str]:
        """
        Initialize conda commands.

        Returns:
            List of conda activation commands
        """
        return [
            f". /opt/miniconda3/etc/profile.d/conda.sh",
            # f"conda info --envs",
            f"conda activate {conda_env}",
        ]


class SlurmJobManager:
    """
    Manager for submitting and monitoring Slurm jobs.

    Example:
        >>> manager = SlurmJobManager()
        >>> config = manager.get_cpu_config()
        >>> job_id = manager.submit_command("srun vasp_std")
    """

    def __init__(self):
        self.submitted_jobs: list[dict[str, Any]] = []

    def submit_command(
        self,
        command: str,
        config: SlurmConfig | None = None,
        conda_env: str  = "",
        module_name: str = "",
        workdir: Path | str = "",
    ) -> str | None:
        """
        Submit a job using sbatch --wrap (without script file).
        This method allows running commands directly without creating a script.

        Args:
            command: Shell command to run (e.g., "srun vasp_std")
            config: Slurm configuration
            workdir: Working directory (default: current directory)

        Returns:
            Job ID if successful, None otherwise
        """
        config = config or self.get_cpu_config()
        workdir = Path(workdir).resolve()
        workdir.mkdir(parents=True, exist_ok=True)
        log.info(f"Working directory: {workdir}")

        # Build the actual command to wrap (with module initialization)

        init_cmds = []
        if module_name:
            init_cmds.extend(config.init_module(module_name))
        if conda_env:
            init_cmds.extend(config.init_conda(conda_env))
        init_cmds.extend(config.extra_commands)
        full_cmd = "; ".join(init_cmds) + "; " + command

        # Build sbatch command with --wrap
        cmd_parts = [
            "sbatch",
            f"--job-name={config.job_name}",
            f"--output={config.output_log}",
            f"--mem={config.memory}",
            f"--nodes={config.nodes}",
            f"--ntasks={config.ntasks}",
            f"--cpus-per-task={config.cpus_per_task}",
            f"--partition={config.partition}",
            f"--time={config.time_limit}",
            f'--wrap="{full_cmd}"',
        ]

        # Add --error if specified
        if config.error_log:
            cmd_parts.append(f"--error={config.error_log}")

        # Add GPU-specific option
        if config.gpus_per_task > 0:
            cmd_parts.append(f"--gpus-per-task={config.gpus_per_task}")

        cmd = " ".join(cmd_parts)
        log.info(f"Submitting command: {cmd}")

        # Submit the job
        job_id = None
        try:
            result = subprocess.run(
                cmd,
                shell=True,
                cwd=workdir,
                capture_output=True,
                text=True,
                check=True,
            )
            for line in result.stdout.strip().splitlines():
                if "Submitted batch job" in line:
                    job_id = line.split()[-1]
                    break
        except subprocess.CalledProcessError as e:
            log.error(f"Error submitting job: {e}")
            log.error(f"stdout: {result.stdout}")
            return None

        if job_id is None:
            log.error(f"Failed to parse job ID from output: {result.stdout}")
            return None

        log.info(f"Job submitted successfully: {job_id}")
        self.submitted_jobs.append(
            {
                "job_id": job_id,
                "command": command,
                "config": config,
                "workdir": str(workdir),
            }
        )
        return job_id

    def get_cpu_config(
        self,
        job_name: str = "vasp-cpu",
        ntasks: int = 48,
        memory: str = "20G",
        **kwargs,
    ) -> SlurmConfig:
        """
        Get default CPU configuration.

        Args:
            job_name: Job name (default: "vasp-cpu")
            ntasks: Number of tasks (default: 48)
            memory: Memory per node (default: "20G")
            **kwargs: Additional SlurmConfig parameters

        Returns:
            SlurmConfig instance
        """
        return SlurmConfig(
            job_name=job_name,
            ntasks=ntasks,
            memory=memory,
            partition="partCPU",
            **kwargs,
        )

    def get_gpu_config(
        self,
        job_name: str = "vasp-gpu",
        ntasks: int = 1,
        gpus_per_node: int = 1,
        memory: str = "10G",
        **kwargs,
    ) -> SlurmConfig:
        """
        Get default GPU configuration.

        Args:
            job_name: Job name (default: "vasp-gpu")
            ntasks: Number of tasks (default: 1 CPU for GPU)
            gpus_per_node: Number of GPUs per node (default: 1)
            memory: Memory per node (default: "10G")
            **kwargs: Additional SlurmConfig parameters

        Returns:
            SlurmConfig instance
        """
        return SlurmConfig(
            job_name=job_name,
            ntasks=ntasks,
            gpus_per_task=gpus_per_node,
            memory=memory,
            partition="partGPU",
            **kwargs,
        )

    def get_job_status(self, job_id: str) -> dict[str, Any] | None:
        """
        Get job status from Slurm.

        Args:
            job_id: Job ID
        Returns:
            Job status dict or None
        """
        cmd = ["scontrol", "show", "job", job_id]
        status = {}
        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
            )
            for line in result.stdout.strip().splitlines():
                if "=" not in line:
                    continue
                for item in line.split():
                    if "=" not in item:
                        continue
                    key, value = item.split("=", 1)
                    status[key.strip()] = value.strip()
            return status
        except subprocess.CalledProcessError as e:
            log.error(f"Error getting job status: {e}")
            log.error(f"stdout: {result.stdout}")
            return None

    def cancel_job(self, job_id: str) -> bool:
        """
        Cancel a Slurm job.

        Args:
            job_id: Job ID
        Returns:
            True if cancelled successfully, False otherwise
        """
        cmd = ["scancel", job_id]
        try:
            subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
            )
            return True
        except subprocess.CalledProcessError as e:
            log.error(f"Error canceling job: {e}")
            return False

    def list_jobs(self, user: str = "") -> list[dict[str, Any]]:
        """
        List Slurm jobs.

        Args:
            user: Username to filter jobs (default: current user)
        Returns:
            List of job status dicts
        """
        cmd = ["squeue", "-u", user or os.getlogin()]
        jobs = []
        try:
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
            )
            lines = result.stdout.strip().splitlines()
            header = lines[0].split()
            for line in lines[1:]:
                job = dict(zip(header, line.split()))
                jobs.append(job)
            return jobs
        except subprocess.CalledProcessError as e:
            log.error(f"Error listing jobs: {e}")
            return []
