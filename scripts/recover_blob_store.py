#!/usr/bin/env python
"""Rebuild a lost ``store-data.json`` from copies the store still holds.

Why this exists: ``Worker.run_flow`` used to create the "data" additional store
as a ``MemoryStore``.  atomate2's ``@vasp_job`` passes ``data=_DATA_OBJECTS``,
so every VASP TaskDoc's ``force_constants`` and ``normalmode_eigenvecs`` were
replaced by a ``blob_uuid`` pointer whose target lived only in the writer's
process.  The pointers are in ``store.json``; the payloads are not.

A resume opens a new process, so it cannot resolve those pointers, reads the
jobs as "incomplete output", and re-runs the VASP they represent.  That is what
happened to FCC-Cr-Nb-4 on 2026-10-29: six completed DFPT volumes were
recomputed.

The payloads are not actually lost.  The per-volume collector embeds the whole
*resolved* TaskDoc in its own output, so for every ``DFPT volume i/N`` row with
a dangling pointer there is a ``collect DFPT volume i/N`` row holding the same
value, inline, at the same JSON path.  Where the VASP directory is still
archived, the value can also be re-derived from ``vasprun.xml`` -- and the two
must agree, which is what makes this a reconstruction rather than a guess.

The script only ever *creates* ``store-data.json``.  It refuses to overwrite
an existing one, and it records where each payload came from in a sidecar:

    python scripts/recover_blob_store.py <run-dir> --check
    python scripts/recover_blob_store.py <run-dir> --write
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

DOCS_NAME = "store.json"
BLOBS_NAME = "store-data.json"
SIDECAR_NAME = "store-data.recovery.json"

#: ``DFPT volume i/N`` -> ``collect DFPT volume i/N``.  The collector is the one
#: job that received the resolved TaskDoc, so its own output carries the payload
#: the volume job had offloaded.
VOLUME_PREFIX = "DFPT volume "
COLLECT_PREFIX = "collect DFPT volume "

#: Where the collector keeps the resolved TaskDoc inside its output.
TASK_DOC_KEY = "task_doc"


def _pointers(node, path=()):
    """``(path, blob_dict)`` for every ``blob_uuid`` pointer under ``node``."""
    if isinstance(node, dict):
        if isinstance(node.get("blob_uuid"), str):
            yield path, node
        for key, value in node.items():
            yield from _pointers(value, path + (key,))
    elif isinstance(node, list):
        for index, value in enumerate(node):
            yield from _pointers(value, path + (index,))


def _dig(node, path):
    for key in path:
        node = node[key]
    return node


def _donor_name(name: str) -> str | None:
    if name.startswith(VOLUME_PREFIX):
        return COLLECT_PREFIX + name[len(VOLUME_PREFIX):]
    return None


def plan(run_dir: Path) -> tuple[list[dict], list[str]]:
    """Every blob to write, with the document it came from and its digest."""
    rows = json.loads((run_dir / DOCS_NAME).read_text())
    by_name = {row.get("name"): row for row in rows}

    recovered: list[dict] = []
    problems: list[str] = []
    for name, row in by_name.items():
        for path, blob in _pointers(row.get("output")):
            donor = _donor_name(str(name))
            if donor is None or donor not in by_name:
                problems.append(
                    f"{name}: blob {blob['blob_uuid']} at {path} has no donor row "
                    f"({donor or 'no name mapping'})"
                )
                continue
            try:
                value = _dig(by_name[donor]["output"][TASK_DOC_KEY], path)
            except (KeyError, IndexError, TypeError) as exc:
                problems.append(
                    f"{name}: blob {blob['blob_uuid']} at {path} is not present in "
                    f"{donor} ({type(exc).__name__}: {exc})"
                )
                continue
            if isinstance(value, dict) and "blob_uuid" in value:
                problems.append(
                    f"{name}: blob {blob['blob_uuid']} at {path} is still a pointer "
                    f"in {donor}; nothing to copy"
                )
                continue
            recovered.append(
                {
                    "@class": blob.get("@class", ""),
                    "@module": blob.get("@module", ""),
                    "blob_uuid": blob["blob_uuid"],
                    "data": value,
                    "job_uuid": row["uuid"],
                    "job_index": row["index"],
                    "recovered_from": {
                        "job": donor,
                        "uuid": by_name[donor]["uuid"],
                        "path": [str(p) for p in (TASK_DOC_KEY,) + path],
                    },
                }
            )
    return recovered, problems


def verify_against_archive(run_dir: Path, recovered: list[dict]) -> list[str]:
    """Re-derive every payload from the archived ``vasprun.xml`` and compare.

    This is the independent check: the collector's copy and the VASP archive are
    separate records of the same run.  Disagreement means the reconstruction is
    not the payload the pointer referred to, and must not be written.
    """
    try:
        import numpy as np
        from pymatgen.io.vasp.outputs import Vasprun
    except ImportError as exc:  # pragma: no cover - the project env has these
        return [f"cannot verify against the archive: {exc}"]

    rows = json.loads((run_dir / DOCS_NAME).read_text())
    by_name = {row.get("name"): row for row in rows}
    problems: list[str] = []

    # Group the recovered payloads by the volume they belong to.
    by_volume: dict[str, dict[str, list]] = {}
    for item in recovered:
        donor = item["recovered_from"]["job"]
        if not donor.startswith(COLLECT_PREFIX):
            continue
        by_volume.setdefault(donor, {})[str(item["recovered_from"]["path"][-1])] = item["data"]

    for donor, fields in by_volume.items():
        label = donor[len(COLLECT_PREFIX):]           # e.g. "3/7"
        index = label.split("/")[0]
        candidates = sorted(run_dir.glob(f"*-DFPT_volume_{index}/*/vasprun.xml*"))
        if not candidates:
            problems.append(f"{donor}: no archived vasprun.xml to check against")
            continue
        # The collector parsed the .gz the job gzipped at the end; prefer it.
        candidates.sort(key=lambda p: (not p.name.endswith(".gz"), len(str(p))))
        path = candidates[0]
        try:
            parsed = Vasprun(path, parse_potcar_file=False, parse_dos=False,
                             parse_eigen=False)
        except Exception as exc:  # noqa: BLE001 - reported, not raised
            problems.append(f"{donor}: {path.name} could not be parsed ({exc})")
            continue
        for field, value in fields.items():
            actual = getattr(parsed, field, None)
            if actual is None:
                problems.append(f"{donor}: {path.name} carries no {field}")
                continue
            if not np.array_equal(np.asarray(actual, dtype=float),
                                  np.asarray(value, dtype=float)):
                problems.append(
                    f"{donor}: recovered {field} differs from {path.name}"
                )
    return problems


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("run_dir", type=Path, help="the run directory holding store.json")
    ap.add_argument("--write", action="store_true",
                    help="write store-data.json (refused if one already exists)")
    ap.add_argument("--check", action="store_true", default=True,
                    help="verify each payload against the archived vasprun (default)")
    ap.add_argument("--no-archive-check", dest="check", action="store_false")
    args = ap.parse_args(argv)

    run_dir = args.run_dir.resolve()
    if not (run_dir / DOCS_NAME).is_file():
        print(f"{run_dir}: no {DOCS_NAME}", file=sys.stderr)
        return 2
    blobs_path = run_dir / BLOBS_NAME

    recovered, problems = plan(run_dir)
    print(f"{run_dir.name}: {len(recovered)} blob(s) recoverable, "
          f"{len(problems)} problem(s)")
    for problem in problems:
        print(f"  ! {problem}")

    if args.check and recovered:
        archive_problems = verify_against_archive(run_dir, recovered)
        for problem in archive_problems:
            print(f"  ! {problem}")
        if archive_problems:
            problems.extend(archive_problems)
        elif recovered:
            print("  every payload matches the archived vasprun.xml")

    if problems:
        print("refusing to write: the reconstruction is incomplete or unverified",
              file=sys.stderr)
        return 1

    if blobs_path.exists():
        print(f"{blobs_path.name} already exists; nothing to do (it is never "
              f"overwritten)")
        return 0

    if not args.write:
        print("dry run: pass --write to create " + BLOBS_NAME)
        return 0

    # The pointer fields are written exactly as jobflow writes them; the two
    # recovery fields document where the payload came from and are ignored on
    # read (``get_output`` selects only blob_uuid and data).
    blobs_path.write_text(json.dumps(recovered))
    sidecar = run_dir / SIDECAR_NAME
    sidecar.write_text(json.dumps({
        "schema": "htvasp-blob-recovery/v1",
        "run": run_dir.name,
        "store_sha256": hashlib.sha256((run_dir / DOCS_NAME).read_bytes()).hexdigest(),
        "blobs": [
            {"blob_uuid": item["blob_uuid"], "job_uuid": item["job_uuid"],
             **item["recovered_from"]}
            for item in recovered
        ],
        "verified_against": "archived vasprun.xml"
        if args.check else "collector copy only",
        "reason": "the additional store was a MemoryStore; a resume could not "
                  "resolve these pointers",
    }, indent=2, sort_keys=True) + "\n")
    print(f"wrote {blobs_path.name} ({len(recovered)} blobs) and {sidecar.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
