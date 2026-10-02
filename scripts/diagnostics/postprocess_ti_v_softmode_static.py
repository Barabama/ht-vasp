#!/usr/bin/env python3
"""Audit completed soft-mode static calculations without promoting them."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from pymatgen.io.vasp.outputs import Vasprun


def archive(path: Path) -> None:
    target = path.with_name(path.name + ".gz")
    with path.open("rb") as source, gzip.open(target, "wb") as dest:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            dest.write(block)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    records = []
    errors = []
    for run in manifest["runs"]:
        path = root / run["path"]
        try:
            vasprun = Vasprun(
                str(path / "vasprun.xml"),
                parse_dos=False,
                parse_eigen=False,
                parse_projected_eigen=False,
                parse_potcar_file=False,
            )
            forces = np.asarray(vasprun.ionic_steps[-1]["forces"], dtype=float)
            record = {
                **run,
                "converged": bool(vasprun.converged),
                "energy_eV": float(vasprun.final_energy),
                "max_force_eV_A": float(np.abs(forces).max()),
                "nsteps": len(vasprun.ionic_steps),
            }
            if not vasprun.converged:
                errors.append(f"not_converged:{run['path']}")
            for filename in ("vasprun.xml", "OUTCAR", "vasp.out"):
                item = path / filename
                if not item.is_file():
                    errors.append(f"missing:{item}")
                else:
                    archive(item)
        except Exception as exc:
            errors.append(f"parse:{run['path']}:{type(exc).__name__}:{exc}")
            record = {**run, "error": f"{type(exc).__name__}: {exc}"}
        records.append(record)
    # Keep q-point cases separate.  Two distinct soft modes share the same
    # ``mid`` volume; grouping only by volume silently overwrites one case.
    summary = {}
    for case in sorted({item["case"] for item in records}):
        rows = [item for item in records if item["case"] == case and "energy_eV" in item]
        ref = next((item for item in rows if item["sign"] == "reference"), None)
        summary[case] = {
            "volume": rows[0]["volume"] if rows else None,
            "q": rows[0]["q"] if rows else None,
            "frequency_THz": rows[0]["frequency_THz"] if rows else None,
            "reference_energy_eV": ref["energy_eV"] if ref else None,
            "energies_eV": {item["sign"]: item["energy_eV"] for item in rows},
            "delta_to_reference_eV": {
                item["sign"]: item["energy_eV"] - ref["energy_eV"]
                for item in rows if ref and item["sign"] != "reference"
            },
        }
    result = {
        "schema": "ti-v-softmode-static-result/v1",
        "status": "pass" if not errors and len(records) == len(manifest["runs"]) else "failed",
        "production_consumable": False,
        "errors": errors,
        "records": records,
        "summary_by_case": summary,
        # Kept as an alias for consumers of the first diagnostic output; unlike
        # the old implementation, this mapping is case-keyed and lossless.
        "summary_by_volume": summary,
        "interpretation": "Static energy comparison only; relaxation and QHA promotion require a separate gate.",
    }
    (root / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "errors": len(errors), "summary_by_case": summary}, indent=2))
    return 0 if result["status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
