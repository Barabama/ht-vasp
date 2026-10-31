"""Choosing between several runs of one endpoint, when building the manifest.

An endpoint can accumulate more than one ``*-dfpt-qha.json``: an aborted attempt
beside a good one, a resubmission under a new run id, a re-analysis.  Which one
is the endpoint's answer has to be a decision, not an accident of sort order.

It used to be an accident.  ``gather_status.scan`` added one row per run to a
``set``, and ``build_manifest.read_status`` keeps the last line it reads, so the
alphabetically greatest applicability string won.  BCC-Ta-Ta-2 therefore
reported ``unknown`` -- from a 5-volume artifact of 2026-10-16 -- while a
complete 7-volume run from 2026-10-20 sitting beside it said
``harmonic_QHA_applicable``.  ``"unknown"`` sorts after every real verdict,
which is why that one failed; the same mechanism would have preferred
``numerical_unresolved`` over ``harmonic_QHA_applicable`` on any endpoint where
the two runs disagreed the other way.

These tests pin the choice to the artifacts' content, and pin the manifest
round-trip so a future `set` or `sorted()` cannot quietly reintroduce it.
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _artifact(root: Path, endpoint: str, run: str, *, volumes: int,
              applicability: str | None, ok: bool = True,
              mtime: float | None = None) -> Path:
    """One ``<endpoint>-dfpt-qha.json`` under the run directory the tree expects."""
    path = root / endpoint / "dfpt-qhaflow" / "runs" / run / f"{endpoint}-dfpt-qha.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "ok": ok,
        "endpoint": endpoint,
        "temperatures": list(range(0, 1001, 100)),
        "volumes": [{"index": i} for i in range(volumes)],
        "eos": {"V0_angstrom3": 70.0, "B0_GPa": 150.0},
    }
    if applicability is not None:
        payload["qha_applicability"] = applicability
    path.write_text(json.dumps(payload))
    if mtime is not None:
        os.utime(path, (mtime, mtime))
    return path


def _rows(tmp_path: Path) -> dict[str, str]:
    gather = _load("gather_status")
    rows: set[str] = set()
    gather.scan(tmp_path / "endpoints", "429pro", rows)
    return {row.split("|")[0]: row.split("|")[3] for row in rows}


def test_a_complete_run_beats_a_verdictless_partial(tmp_path):
    """The case that was actually wrong on disk: 7 volumes beat 5 + 'unknown'."""
    root = tmp_path / "endpoints"
    _artifact(root, "BCC-Ta-Ta-2", "dfpt-qha-1c5892fc7f4e68d93369",
              volumes=7, applicability="harmonic_QHA_applicable", mtime=1_700_000_000)
    _artifact(root, "BCC-Ta-Ta-2", "dfpt-qha-e868d18e1a078ea6152c",
              volumes=5, applicability=None, mtime=1_690_000_000)

    assert _rows(tmp_path)["BCC-Ta-Ta-2"] == "harmonic_QHA_applicable"


def test_a_complete_soft_run_beats_an_incomplete_harmonic_one(tmp_path):
    """Completeness is the first criterion, not the verdict.

    A 3-volume artifact cannot have fitted the same F(V,T) surface, so calling
    the endpoint "harmonic" off it would be a stronger claim than the data
    supports.  The complete run's soft-mode verdict is the endpoint's answer.
    """
    root = tmp_path / "endpoints"
    _artifact(root, "FCC-Nb-Ti-4", "dfpt-qha-aaa",
              volumes=7, applicability="finite_q_soft_mode")
    _artifact(root, "FCC-Nb-Ti-4", "dfpt-qha-zzz",
              volumes=3, applicability="harmonic_QHA_applicable")

    assert _rows(tmp_path)["FCC-Nb-Ti-4"] == "finite_q_soft_mode"


def test_equally_complete_runs_are_ranked_by_verdict(tmp_path):
    """Two complete runs: neither sort order nor file order may decide."""
    root = tmp_path / "endpoints"
    _artifact(root, "FCC-Ti-Zr-4", "dfpt-qha-zzz",
              volumes=7, applicability="numerical_unresolved", mtime=1_700_000_000)
    _artifact(root, "FCC-Ti-Zr-4", "dfpt-qha-aaa",
              volumes=7, applicability="harmonic_QHA_applicable", mtime=1_600_000_000)

    # The name that sorts first is the worse run, and it is also the older one.
    assert _rows(tmp_path)["FCC-Ti-Zr-4"] == "harmonic_QHA_applicable"


def test_a_verdict_beats_no_verdict_when_completeness_ties(tmp_path):
    """An artifact without the field predates it; it is weaker evidence."""
    root = tmp_path / "endpoints"
    _artifact(root, "FCC-Al-Ta-4", "dfpt-qha-a",
              volumes=7, applicability=None, mtime=1_700_000_000)
    _artifact(root, "FCC-Al-Ta-4", "dfpt-qha-b",
              volumes=7, applicability="gamma_soft_mode", mtime=1_600_000_000)

    assert _rows(tmp_path)["FCC-Al-Ta-4"] == "gamma_soft_mode"


def test_a_failed_run_is_not_a_candidate(tmp_path):
    """``ok=False`` is a recorded failure, not a less-good result."""
    root = tmp_path / "endpoints"
    _artifact(root, "BCC-Al-Ta-2", "dfpt-qha-bad",
              volumes=0, applicability=None, ok=False, mtime=1_700_000_000)
    _artifact(root, "BCC-Al-Ta-2", "dfpt-qha-good",
              volumes=7, applicability="finite_q_soft_mode")

    assert _rows(tmp_path)["BCC-Al-Ta-2"] == "finite_q_soft_mode"


def test_one_row_per_endpoint_not_one_per_run(tmp_path):
    """The defect was a row per run; a duplicate lets any consumer pick wrong."""
    root = tmp_path / "endpoints"
    for run, app in (("dfpt-qha-1", "finite_q_soft_mode"),
                     ("dfpt-qha-2", "harmonic_QHA_applicable"),
                     ("dfpt-qha-3", "gamma_soft_mode")):
        _artifact(root, "FCC-Zr-Ti-4", run, volumes=7, applicability=app)

    gather = _load("gather_status")
    rows: set[str] = set()
    gather.scan(root, "429pro", rows)
    assert len(rows) == 1
    assert next(iter(rows)).endswith("|harmonic_QHA_applicable")


def test_the_manifest_keeps_the_row_scan_chose(tmp_path):
    """End to end: gather_status -> build_manifest, through the real files.

    ``build_manifest`` rewrites ``status[name]`` for every line it reads, so two
    lines for one endpoint is exactly how the wrong verdict used to survive.
    """
    gather = _load("gather_status")
    build = _load("build_manifest")

    endpoint = "BCC-Ta-Ta-2"
    root = tmp_path / "endpoints"
    _artifact(root, endpoint, "dfpt-qha-1c5892fc7f4e68d93369",
              volumes=7, applicability="harmonic_QHA_applicable", mtime=1_700_000_000)
    _artifact(root, endpoint, "dfpt-qha-e868d18e1a078ea6152c",
              volumes=5, applicability=None, mtime=1_690_000_000)

    rows: set[str] = set()
    gather.scan(root, "429pro", rows)
    status_file = tmp_path / "status.txt"
    status_file.write_text("\n".join(sorted(rows)) + "\n")

    status = build.read_status([status_file])
    assert status[endpoint][2] == "harmonic_QHA_applicable"
