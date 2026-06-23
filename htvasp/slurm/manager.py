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
    """Slurm job configuration for CPU and GPU partitions.

    Example:
        >>> # CPU job
        >>> config = SlurmConfig(
        ...     job_name="relax-cpu",
        ...     partition="partCPU",
        ...     ntasks=32,
        ...     memory="64G",
        ... )

        >>> # GPU job
        >>> config = SlurmConfig(
        ...     job_name="relax-gpu",
        ...     partition="partGPU",
        ...     ntasks=1,
        ...     gpus_per_task=1,
        ...     memory="20G",
        ... )
    """

    job_name: str = "vasp-job"
    output_log: str = "job.log"
    error_log: str = ""  # Empty = don't use --error unless specified
    memory: str = "20G"
    nodes: int = 1
    ntasks: int = 1
    ntasks_per_node: int = 0  # don't use unless > 0
    nodelist: str = ""
    cpus_per_task: int = 1
    gpus_per_task: int = 0  # Only for GPU partition
    partition: Literal["partCPU", "partGPU"] = "partCPU"
    time_limit: str = "100:00:00"
    dependency: str = ""  # Slurm job dependency (e.g., afterok:12345)

    # VASP commands
    vasp_cmd: str = "srun vasp_std"
    vasp_gam_cmd: str = "srun vasp_gam"

    # Environment
    conda_prefix: str = "/opt/miniconda3"
    conda_env: str = ""
    module_name: str | list[str] = ""

    # Extra commands
    extra_commands: list[str] = field(default_factory=list)

    def init_module(self, module_name: str | list[str] = "") -> list[str]:
        names = module_name if isinstance(module_name, list) else [module_name]
        names = [n for n in names if n]
        cmds = [
            ". /etc/profile.d/modules.sh",
            "ulimit -s unlimited",
            "module purge",
        ]
        for name in names:
            cmds.append(f"module load {name}")
        return cmds

    def init_conda(self, conda_env: str = "") -> list[str]:
        if conda_env.startswith("./"):
            conda_env = str(Path.cwd() / conda_env[2:])
        elif conda_env.startswith("."):
            conda_env = str(Path.cwd() / conda_env[1:])

        return [
            f". {self.conda_prefix}/etc/profile.d/conda.sh",
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
        conda_env: str = "",
        module_name: str = "",
        workdir: Path | str = "",
    ) -> str | None:
        """Submit a job via sbatch --wrap.

        Args:
            command: Shell command to run (e.g., "srun vasp_std").
            config: Slurm configuration (conda_env/module_name read from here as fallback).
            conda_env: Conda environment name (overrides config.conda_env if set).
            module_name: Module name(s) to load (overrides config.module_name if set).
            workdir: Working directory.

        Returns:
            Job ID string, or None on failure.
        """
        config = config or self._default_cpu_config()
        workdir = Path(workdir).resolve()
        workdir.mkdir(parents=True, exist_ok=True)

        # Resolve env/module: explicit param > config field
        _conda_env = conda_env or config.conda_env
        _module_name = module_name or config.module_name

        init_cmds: list[str] = []
        if _module_name:
            init_cmds.extend(config.init_module(_module_name))
        if _conda_env:
            init_cmds.extend(config.init_conda(_conda_env))
        init_cmds.extend(config.extra_commands)
        full_cmd = "; ".join(init_cmds + [command])

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
        if config.error_log:
            cmd_parts.append(f"--error={config.error_log}")
        if config.ntasks_per_node:
            cmd_parts.append(f"--ntasks-per-node={config.ntasks_per_node}")
        if config.nodelist:
            cmd_parts.append(f"--nodelist={config.nodelist}")
        if config.dependency:
            cmd_parts.append(f"--dependency={config.dependency}")
        if config.gpus_per_task > 0:
            cmd_parts.append(f"--gpus-per-task={config.gpus_per_task}")

        cmd = " ".join(cmd_parts)
        log.info(f"Submitting: {cmd}")

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                cwd=workdir,
                capture_output=True,
                text=True,
                check=True,
            )
        except subprocess.CalledProcessError as e:
            log.error(f"sbatch failed: {e}")
            return None

        for line in result.stdout.strip().splitlines():
            if "Submitted batch job" in line:
                job_id = line.split()[-1]
                break
        else:
            log.error(f"Failed to parse job ID from: {result.stdout}")
            return None

        log.info(f"Job submitted: {job_id}")
        self.submitted_jobs.append(
            {
                "job_id": job_id,
                "command": command,
                "config": config,
                "workdir": str(workdir),
            }
        )
        return job_id

    def _default_cpu_config(self) -> SlurmConfig:
        return SlurmConfig(
            job_name="vasp-cpu",
            output_log="job.log",
            ntasks=48,
            memory="96G",
            partition="partCPU",
        )

    def get_cpu_config(self, **kwargs) -> SlurmConfig:
        """Get CPU config with defaults, overridden by kwargs."""
        cfg = self._default_cpu_config()
        for k, v in kwargs.items():
            if hasattr(cfg, k):
                setattr(cfg, k, v)
        return cfg

    def get_gpu_config(self, **kwargs) -> SlurmConfig:
        """Get GPU config with defaults, overridden by kwargs."""
        return SlurmConfig(
            job_name=kwargs.pop("job_name", "vasp-gpu"),
            output_log=kwargs.pop("output_log", "job.log"),
            ntasks=kwargs.pop("ntasks", 1),
            gpus_per_task=kwargs.pop("gpus_per_task", 1),
            memory=kwargs.pop("memory", "20G"),
            partition="partGPU",
            module_name=kwargs.pop("module_name", "vasp-gpu"),
            **kwargs,
        )

    def get_job_status(self, job_id: str) -> dict[str, Any] | None:
        """Query job status via scontrol."""
        try:
            result = subprocess.run(
                ["scontrol", "show", "job", job_id],
                capture_output=True,
                text=True,
                check=True,
            )
        except subprocess.CalledProcessError:
            return None

        status: dict[str, Any] = {}
        for line in result.stdout.strip().splitlines():
            for item in line.split():
                if "=" not in item:
                    continue
                key, value = item.split("=", 1)
                status[key.strip()] = value.strip()
        return status

    def is_job_completed(self, job_id: str) -> bool:
        """Check if a job has completed (finished, cancelled, or failed).

        Args:
            job_id: Job ID
        Returns:
            True if job is no longer running/pending, False otherwise
        """
        status = self.get_job_status(job_id)
        if status is None:
            # Job not found means it has completed
            return True

        job_state = status.get("JobState", "")
        # Running states: PENDING, RUNNING, CONFIGURING, COMPLETING
        # Completed states: COMPLETED, CANCELLED, FAILED, TIMEOUT, NODE_FAIL, etc.
        return job_state not in ("PENDING", "RUNNING", "CONFIGURING", "COMPLETING")

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
