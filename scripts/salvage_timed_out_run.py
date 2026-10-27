#!/usr/bin/env python
"""Salvage a DFPT-QHA run that was killed by the Slurm time limit.

A flow's scratch directory is ``<flow_dir>/<name>-<flow.uuid[:8]>``, and
``jobflow`` assigns ``Flow.uuid`` with ``suid()`` -- it is **random on every
invocation**.  Resume therefore cannot find the previous scratch directory by
name; it works the other way round: ``Worker.run_flow`` copies ``store_dir``
into the new scratch directory and ``run_locally_custom`` then matches the jobs
by ``(name, index)`` inside that store.

That copy is made in a ``finally`` block, so it only happens when the flow ends
normally or raises.  A job killed at its time limit gets SIGTERM, Python exits
without unwinding, and **``store_dir`` is never created** -- the completed
volumes sit in ``/nfs_ssd/tmp`` where the next submission will not look.

This script closes that gap: it copies the newest scratch directory for an
endpoint into the ``store_dir`` the ledger expects, so a following
``--resume`` picks up at the first unfinished volume instead of starting over.

    python scripts/salvage_timed_out_run.py FCC-Nb-Zr-4 --dry-run
    python scripts/salvage_timed_out_run.py FCC-Nb-Zr-4

Run it on the cluster that ran the job: the scratch path is local to it.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRATCH = Path("/nfs_ssd/tmp")


def newest_scratch(name: str, scratch: Path = SCRATCH) -> Path | None:
    """The most recently modified ``<name>-<uuid8>`` directory holding a store."""
    candidates = [p for p in scratch.glob(f"{name}-*") if (p / "store.json").is_file()]
    if not candidates:
        return None
    return max(candidates, key=lambda p: (p / "store.json").stat().st_mtime)


def store_jobs(path: Path) -> list[dict]:
    return json.loads((path / "store.json").read_text())


def latest_run_id(ledger: Path) -> str | None:
    """The run id of the newest ledger entry, i.e. the one that just failed."""
    if not ledger.is_file():
        return None
    payload = json.loads(ledger.read_text())
    entries = payload.get("entries", [])
    return entries[-1].get("run_id") if entries else None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("names", nargs="+", help="endpoint name(s)")
    ap.add_argument("--run-id", default=None,
                    help="run id to salvage into; default: the newest ledger entry")
    ap.add_argument("--scratch", type=Path, default=SCRATCH)
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    failures = 0
    for name in args.names:
        ledger = args.root / "data" / "endpoints" / name / "dfpt-qhaflow" / "slurm_submissions.json"
        run_id = args.run_id or latest_run_id(ledger)
        if run_id is None:
            print(f"{name}: no ledger entry, so no run id to salvage into", file=sys.stderr)
            failures += 1
            continue
        source = newest_scratch(name, args.scratch)
        if source is None:
            print(f"{name}: no scratch directory with a store.json under {args.scratch}",
                  file=sys.stderr)
            failures += 1
            continue
        target = (args.root / "data" / "endpoints" / name / "dfpt-qhaflow"
                  / "runs" / run_id)
        jobs = store_jobs(source)
        completed = [j.get("name") for j in jobs]
        print(f"{name}")
        print(f"  scratch  {source}")
        print(f"  target   {target}")
        print(f"  store    {len(jobs)} jobs: {', '.join(completed[:6])}"
              + (" ..." if len(completed) > 6 else ""))
        if target.exists() and any(target.iterdir()):
            print(f"  SKIP: {target} already has content; not overwriting")
            failures += 1
            continue
        if args.dry_run:
            print("  dry-run: nothing copied")
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(source, target)
        print(f"  copied {sum(1 for _ in target.rglob('*'))} entries")
        print(f"  next: resubmit with --resume (run id {run_id})")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
