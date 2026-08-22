#!/usr/bin/env python3
"""Add W90 POSCAR structures to existing TB2J JSON result files.

Run from the ht-vasp repository root:

    python scripts/backfill_tb2j_structures.py

The script is idempotent by default. Use ``--force`` to replace existing
``structure`` fields and ``--dry-run`` to inspect candidate files first.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from pymatgen.core import Structure


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_ROOT = REPOSITORY_ROOT / "data" / "endmembers"


def _find_poscar(json_path: Path) -> Path:
    """Return the W90 POSCAR associated with an end-member TB2J JSON file."""
    w90_dir = json_path.parent / "tb2jflow" / "2-w90_static"
    for name in ("POSCAR", "POSCAR.gz"):
        candidate = w90_dir / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"No W90 POSCAR found below {w90_dir}")


def backfill(json_path: Path, *, dry_run: bool, force: bool) -> str:
    """Insert the serialized W90 structure into one TB2J JSON result."""
    with json_path.open(encoding="utf-8") as handle:
        payload = json.load(handle)

    if "structure" in payload and not force:
        return "skipped"

    poscar_path = _find_poscar(json_path)
    structure = Structure.from_file(poscar_path)
    expected_atoms = payload.get("num_atoms")
    if expected_atoms is not None and int(expected_atoms) != len(structure):
        raise ValueError(
            f"num_atoms={expected_atoms} disagrees with {len(structure)} sites in "
            f"{poscar_path}"
        )

    if dry_run:
        return "would-update"

    payload["structure"] = structure.as_dict()
    tmp_path = json_path.with_suffix(".json.tmp")
    with tmp_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    tmp_path.replace(json_path)
    return "updated"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root",
        type=Path,
        default=DEFAULT_DATA_ROOT,
        help=f"End-member data root (default: {DEFAULT_DATA_ROOT})",
    )
    parser.add_argument("--dry-run", action="store_true", help="Do not write JSON files")
    parser.add_argument("--force", action="store_true", help="Replace existing structure fields")
    args = parser.parse_args()

    json_paths = sorted(args.data_root.rglob("*-tb2j.json"))
    if not json_paths:
        print(f"No TB2J JSON files found below {args.data_root}", file=sys.stderr)
        return 1

    counts: dict[str, int] = {"updated": 0, "would-update": 0, "skipped": 0, "failed": 0}
    for json_path in json_paths:
        try:
            action = backfill(json_path, dry_run=args.dry_run, force=args.force)
            counts[action] += 1
            print(f"{action:12} {json_path.relative_to(args.data_root)}")
        except (FileNotFoundError, ValueError, json.JSONDecodeError) as exc:
            counts["failed"] += 1
            print(f"failed       {json_path}: {exc}", file=sys.stderr)

    print(
        "Summary: "
        + ", ".join(f"{name}={count}" for name, count in counts.items())
    )
    return 1 if counts["failed"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
