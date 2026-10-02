#!/usr/bin/env python3
"""Rebuild legacy TB2J JSON artifacts with the producer magnetic v2 contract.

The migration is metadata-only: it reads the persisted ``exchange.out`` and
the producer POSCAR, then writes a new tree.  Existing JSON files are never
modified.  ``exchange.out`` is used as the numerical source because legacy
JSON files may be incomplete or use the pre-contract schema.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from typing import Any

from pymatgen.core import Structure

from htvasp.tb2j.context import build_magnetic_context
from htvasp.tb2j.jobs import canonicalize_tb2j_result, parse_exchange_out


MIGRATION_SCHEMA_VERSION = "htvasp-tb2j-migration/v2"


def _jsonable_structure(structure: Structure, context: dict[str, Any]) -> dict[str, Any]:
    """Attach stable producer site IDs to the pymatgen structure payload."""

    payload = structure.as_dict()
    sites = payload.get("sites", [])
    context_sites = context.get("sites", [])
    if len(sites) != len(context_sites):
        raise ValueError("context site count does not match POSCAR site count")
    for site, context_site in zip(sites, context_sites):
        site["label"] = context_site["id"]
        site["id"] = context_site["id"]
    return payload


def _source_dirs(input_root: Path) -> list[Path]:
    result = []
    for directory in sorted(path for path in input_root.iterdir() if path.is_dir()):
        if (directory / "tb2jflow/3-tb2j_solve/exchange.out").is_file():
            result.append(directory)
    return result


def rebuild(
    input_root: str | Path,
    output_root: str | Path,
    *,
    names: list[str] | None = None,
) -> dict[str, Any]:
    source_root = Path(input_root).resolve()
    destination_root = Path(output_root).resolve()
    if not source_root.is_dir():
        raise ValueError(f"input root is not a directory: {source_root}")
    if destination_root.exists() and any(destination_root.iterdir()):
        raise ValueError(f"output root must be new or empty: {destination_root}")

    selected = set(names or ())
    directories = [
        directory
        for directory in _source_dirs(source_root)
        if not selected or directory.name in selected
    ]
    if selected and len(directories) != len(selected):
        missing = sorted(selected - {directory.name for directory in directories})
        raise ValueError(f"requested structures have no exchange.out: {missing}")

    records: list[dict[str, Any]] = []
    for directory in directories:
        name = directory.name
        exchange_path = directory / "tb2jflow/3-tb2j_solve/exchange.out"
        poscar_path = source_root.parent / "poscars" / f"{name}.vasp"
        destination = destination_root / name / f"{name}-tb2j.json"
        record: dict[str, Any] = {
            "name": name,
            "exchange_out": exchange_path.as_posix(),
            "poscar": poscar_path.as_posix(),
            "output": destination.as_posix(),
        }
        try:
            structure = Structure.from_file(poscar_path)
            context = build_magnetic_context(name, structure)
            parsed = parse_exchange_out(exchange_path)
            parsed["structure"] = _jsonable_structure(structure, context)
            parsed["elements"] = [str(element) for element in structure.symbol_set]
            parsed["num_atoms"] = len(structure)
            canonical = canonicalize_tb2j_result(
                parsed,
                output_path=destination,
                magnetic_context=context,
            )
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_text(
                json.dumps(canonical, indent=2, sort_keys=True, ensure_ascii=True) + "\n",
                encoding="utf-8",
            )
            observable = canonical["magnetic_observables"]
            record.update(
                {
                    "status": "pass",
                    "context_status": context.get("calphad", {}).get("status"),
                    "phase": context.get("calphad", {}).get("phase"),
                    "observable_status": observable.get("status"),
                    "transition_type": observable.get("transition_type"),
                    "TC": observable.get("TC"),
                    "TN": observable.get("TN"),
                    "BMAGN": observable.get("BMAGN"),
                    "production_ready": canonical["audit"].get("production_ready", False),
                }
            )
        except Exception as exc:  # retain per-file diagnostics and continue
            record.update(
                {
                    "status": "fail",
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                    "production_ready": False,
                }
            )
        records.append(record)

    destination_root.mkdir(parents=True, exist_ok=True)
    report = {
        "schema_version": MIGRATION_SCHEMA_VERSION,
        "operation": "legacy_tb2j_to_htvasp_magnetic_v2",
        "input_root": source_root.as_posix(),
        "output_root": destination_root.as_posix(),
        "source_policy": "exchange.out plus producer POSCAR; legacy JSON is not numerical authority",
        "counts": {
            "total": len(records),
            "passed": sum(item["status"] == "pass" for item in records),
            "failed": sum(item["status"] == "fail" for item in records),
            "production_ready": sum(item.get("production_ready", False) for item in records),
        },
        "files": records,
    }
    (destination_root / "migration-report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", default="data/endmembers")
    parser.add_argument("--output-root", default="data/endmembers-tb2j-v2")
    parser.add_argument("--name", action="append", dest="names", help="migrate one structure; repeatable")
    args = parser.parse_args()
    try:
        report = rebuild(args.input_root, args.output_root, names=args.names)
    except Exception as exc:
        parser.error(str(exc))
    print(
        f"{report['counts']['passed']}/{report['counts']['total']} migrated; "
        f"failed={report['counts']['failed']}; "
        f"production_ready={report['counts']['production_ready']}; "
        f"report={Path(report['output_root']) / 'migration-report.json'}"
    )
    return 0 if report["counts"]["failed"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
