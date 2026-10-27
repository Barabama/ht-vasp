"""Resolve a POTCAR directory that actually exists.

``PMG_VASP_PSP_DIR`` is read from the environment first and from a pymatgen
config file otherwise.  The two can disagree between nodes of the same
cluster: on B5070 the config points at ``/nfs_ssd/documents`` while B5071's
home directory still carries an old ``/workspace/documents``.  A job that
lands on the wrong node then dies deep inside pymatgen with
``PMG_VASP_PSP_DIR=... does not exist`` -- after the allocation has started.

The submitting process resolves the directory once, and the Slurm script
exports that value, so the job sees the same POTCAR set wherever it runs.
"""

from __future__ import annotations

import os
from pathlib import Path

from pymatgen.core import SETTINGS


def vasp_psp_dir() -> str:
    """A PMG_VASP_PSP_DIR that exists on disk, or a clear error."""

    candidates = [str(SETTINGS.get("PMG_VASP_PSP_DIR") or "")]
    for config in (Path.home() / ".config" / ".pmgrc.yaml", Path.home() / ".pmgrc.yaml"):
        if not config.is_file():
            continue
        for line in config.read_text(encoding="utf-8", errors="replace").splitlines():
            key, separator, value = line.partition(":")
            if separator and key.strip() == "PMG_VASP_PSP_DIR":
                candidates.append(value.strip().strip("'\""))

    for candidate in candidates:
        if candidate and Path(candidate).is_dir():
            return candidate
    raise RuntimeError(
        "no usable PMG_VASP_PSP_DIR: tried "
        + ", ".join(repr(c) for c in candidates if c)
        + ".  Set it in the environment or in ~/.config/.pmgrc.yaml."
    )


def psp_dir_exports() -> list[str]:
    """Shell lines that pin the directory for a Slurm job."""

    value = vasp_psp_dir()
    os.environ.setdefault("PMG_VASP_PSP_DIR", value)
    return [f"export PMG_VASP_PSP_DIR={value}"]


__all__ = ["psp_dir_exports", "vasp_psp_dir"]
