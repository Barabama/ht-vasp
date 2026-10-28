#!/usr/bin/env python
"""Mark ledger entries whose Slurm job has vanished and left no artifact.

``SlurmJobManager.reconcile_entry`` asks squeue and then sacct.  This cluster
has ``Slurm accounting storage is disabled``, so a job that has left the queue
can never be resolved: squeue no longer knows it, sacct knows nothing at all,
and the entry sticks at ``reconcile_required`` forever.  ``submit_intent``
refuses to submit anything in that state, which makes ``--resume`` unusable
for exactly the runs that most need it.

The only evidence available is therefore negative -- the job is not in squeue
and its artifact was never written -- so this script demands both before it
changes anything, and prints what it did.  Dry run by default.

    ./.conda/bin/python scripts/reconcile_vanished_jobs.py FCC-Cr-Nb-4
    ./.conda/bin/python scripts/reconcile_vanished_jobs.py --apply FCC-Cr-Nb-4

It only rewrites the bookkeeping entry: ``state`` becomes ``failed`` and
``job_id`` is cleared (``reconcile_entry`` re-queries any job id it is given,
so keeping the old one would reproduce the same dead end).  No result file,
run directory or Slurm job is touched.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TRAIL = Path("dfpt-qhaflow") / "slurm_submissions.json"
TERMINAL = {"completed", "failed"}


def in_queue(job_id: str) -> bool:
    result = subprocess.run(
        ["squeue", "-h", "-j", str(job_id), "-o", "%T"],
        capture_output=True, text=True, check=False,
    )
    return bool(result.stdout.strip())


def atomic_write(path: Path, payload: dict) -> None:
    handle, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=path.name, suffix=".tmp")
    try:
        with os.fdopen(handle, "w") as stream:
            json.dump(payload, stream, indent=1, sort_keys=True)
            stream.write("\n")
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("names", nargs="+", help="endpoint names")
    parser.add_argument("--apply", action="store_true", help="write the change")
    args = parser.parse_args()

    changed = 0
    for name in args.names:
        path = ROOT / "data" / "endpoints" / name / TRAIL
        if not path.is_file():
            print(f"{name:16s} no ledger at {path} -- nothing to reconcile")
            continue
        ledger = json.loads(path.read_text())
        touched = False
        for entry in ledger.get("entries", []):
            job_id = entry.get("job_id")
            state = entry.get("state")
            if not job_id:
                continue
            if state in TERMINAL and job_id is None:
                continue
            if in_queue(str(job_id)):
                print(f"{name:16s} {job_id}: still in squeue, leave it alone")
                continue
            artifact = entry.get("artifact_path")
            if artifact and Path(artifact).is_file():
                print(f"{name:16s} {job_id}: gone from squeue but {Path(artifact).name} "
                      "exists -- reconcile by hand, not by deleting evidence")
                continue
            print(f"{name:16s} {job_id}: not in squeue, no artifact "
                  f"(state={state}) -> failed, job_id cleared")
            entry["state"] = "failed"
            entry["job_id"] = None
            entry["updated_at"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
            entry["reconcile_note"] = (
                "job absent from squeue and accounting disabled; artifact never "
                "written, so the run is recorded as failed"
            )
            touched = True
            changed += 1
        if touched:
            if args.apply:
                atomic_write(path, ledger)
                print(f"{name:16s} ledger written")
            else:
                print(f"{name:16s} (dry run; pass --apply to write)")

    print(f"\n{changed} entr{'y' if changed == 1 else 'ies'} "
          f"{'updated' if args.apply else 'would be updated'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
