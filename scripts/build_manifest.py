#!/usr/bin/env python
"""Rebuild ``config/dfpt_qha_batches/MANIFEST.md`` from the artifacts.

The status column comes from each endpoint's ``*-dfpt-qha.json``, not from a
booking sheet -- the previous hand-kept table had drifted far enough that it
listed finished endpoints as queued.  Run ``gather_status.py`` on each cluster
first, then this script merges the two lists.

    # on 429pro
    ./.conda/bin/python scripts/gather_status.py 429pro /tmp/status_429.txt
    # on B5070, then copy /tmp/status_507.txt back
    ./.conda/bin/python scripts/build_manifest.py \
        --status /tmp/status_429.txt --status /tmp/status_507.txt \
        --queue config/dfpt_qha_batches/batch-a-429pro.txt=429pro \
        --queue config/dfpt_qha_batches/batch-b-b5070.txt=B5070 \
        --queue config/dfpt_qha_batches/batch-c-b5070-requeued.txt=B5070 \
        --queue config/dfpt_qha_batches/batch-d-b5070-cr.txt=B5070 \
        --notes /tmp/manifest-notes.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
POOL = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_recipes.json")
STRUCTURES = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_structures.json")

COLUMNS = ("endpoint", "class", "cell", "supercell", "magnetic", "batch",
           "status", "qha_applicability")


def read_status(paths: list[Path]) -> dict[str, tuple[str, str, str]]:
    """``name|cluster|status|applicability`` lines; a later file wins."""
    status: dict[str, tuple[str, str, str]] = {}
    for path in paths:
        if not path.is_file():
            continue
        for line in path.read_text().splitlines():
            if not line.strip() or "|" not in line:
                continue
            parts = line.split("|")
            if len(parts) < 4:
                continue
            status[parts[0]] = (parts[1], parts[2], parts[3])
    return status


def read_queue(specs: list[str]) -> dict[str, str]:
    queued: dict[str, str] = {}
    for spec in specs:
        path_text, _, cluster = spec.partition("=")
        path = ROOT / path_text if not Path(path_text).is_absolute() else Path(path_text)
        if not path.is_file():
            continue
        for line in path.read_text().splitlines():
            name = line.strip()
            if name and not name.startswith("#"):
                queued.setdefault(name, cluster)
    return queued


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--status", action="append", type=Path, default=[],
                    help="output of gather_status.py; repeatable, later wins")
    ap.add_argument("--queue", action="append", default=[],
                    help="BATCHFILE=CLUSTER; repeatable")
    ap.add_argument("--notes", type=Path, default=None,
                    help="JSON {endpoint: status text} for the rows whose state "
                         "is not captured by an artifact or a batch file "
                         "(partial runs, resubmissions)")
    ap.add_argument("--pool", type=Path, default=POOL)
    ap.add_argument("--structures", type=Path, default=STRUCTURES)
    ap.add_argument("--out", type=Path,
                    default=ROOT / "config" / "dfpt_qha_batches" / "MANIFEST.md")
    args = ap.parse_args(argv)

    pool = json.loads(args.pool.read_text())["endpoints"]
    structures = json.loads(args.structures.read_text())["endpoints"]
    status = read_status(args.status)
    queued = read_queue(args.queue)
    notes = json.loads(args.notes.read_text()) if args.notes and args.notes.is_file() else {}

    rows: list[tuple] = []
    for name in sorted(pool):
        recipe = pool[name]["recipe"]
        cls = recipe["structure_class"]
        cell = len(structures[name]["structure"]["sites"])
        magnetic = recipe["magnetic_policy"]
        if name in status:
            cluster, text, applicability = status[name]
            batch = cluster
            note = text
        elif name in notes:
            note = notes[name]
            batch = note.split()[0] if note.split() else "-"
            applicability = ""
        elif name in queued:
            batch = queued[name]
            note = f"queued {batch}"
            applicability = ""
        else:
            batch, note, applicability = "-", "not submitted", ""
        rows.append((name, cls, cell, cell * 8, magnetic, batch, note,
                     applicability))

    widths = [max(len(COLUMNS[i]), max(len(str(row[i])) for row in rows))
              for i in range(len(COLUMNS))]
    lines = [
        "# DFPT-QHA endpoint manifest",
        "",
        "Derived from each endpoint's `*-dfpt-qha.json` on the cluster that produced",
        "it -- not from a booking sheet.  Regenerate with `scripts/gather_status.py`",
        "on each cluster, then `scripts/build_manifest.py`.",
        "",
        "| " + " | ".join(COLUMNS[i].ljust(widths[i]) for i in range(len(COLUMNS))) + " |",
        "|" + "|".join("-" * (widths[i] + 2) for i in range(len(COLUMNS))) + "|",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(row[i]).ljust(widths[i])
                                       for i in range(len(COLUMNS))) + " |")

    done = sum(1 for row in rows if str(row[6]).startswith("done"))
    partial = sum(1 for row in rows if "partial" in str(row[6]))
    lines += [
        "",
        f"**{done} / {len(rows)} done.** {partial} more have partial scratch state "
        "that `--resume-from` continues.",
        "",
        "## Resuming a run that hit the time limit",
        "",
        "A job killed at its Slurm limit gets SIGTERM, so `Worker.run_flow` never",
        "reaches its `finally` block and `store_dir` is never created.  The completed",
        "volumes sit in `/nfs_ssd/tmp/<name>-<uuid>` instead.  Continue from there:",
        "",
        "```bash",
        "./.conda/bin/python main_dfpt_qha.py -n <endpoint> --slurm --resume \\",
        "    --resume-from /nfs_ssd/tmp/<name>-<uuid>",
        "```",
        "",
        "`--resume-from` is threaded into the submitted command, so the job itself",
        "sees it.  Without `--resume` it is rejected rather than silently ignored.",
    ]
    args.out.write_text("\n".join(lines) + "\n")
    print(f"wrote {args.out} ({done}/{len(rows)} done, {partial} partial)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
