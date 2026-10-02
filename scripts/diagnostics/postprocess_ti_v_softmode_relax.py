#!/usr/bin/env python3
"""Post-process fixed-cell Ti-V soft-mode internal relaxations."""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
from pathlib import Path

import numpy as np
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Vasprun


def archive(path: Path) -> str:
    target = path.with_name(path.name + ".gz")
    if not target.exists():
        with path.open("rb") as source, gzip.open(target, "wb") as dest:
            shutil.copyfileobj(source, dest)
    return str(target.name)


def min_image_displacements(initial: Structure, final: Structure) -> np.ndarray:
    frac_delta = np.asarray(final.frac_coords) - np.asarray(initial.frac_coords)
    frac_delta -= np.rint(frac_delta)
    return np.asarray(frac_delta @ initial.lattice.matrix, dtype=float)


def lattice_metrics(initial: Structure, final: Structure) -> dict:
    initial_matrix = np.asarray(initial.lattice.matrix, dtype=float)
    final_matrix = np.asarray(final.lattice.matrix, dtype=float)
    return {
        "initial_volume_A3": float(initial.volume),
        "final_volume_A3": float(final.volume),
        "volume_delta_A3": float(final.volume - initial.volume),
        "volume_abs_delta_A3": float(abs(final.volume - initial.volume)),
        "initial_lattice_matrix_A": initial_matrix.tolist(),
        "final_lattice_matrix_A": final_matrix.tolist(),
        "lattice_max_abs_delta_A": float(np.max(np.abs(final_matrix - initial_matrix))),
        "lattice_unchanged": bool(np.allclose(final_matrix, initial_matrix, rtol=0.0, atol=1e-8)),
    }


def symmetry_label(structure: Structure) -> dict:
    try:
        analyzer = SpacegroupAnalyzer(structure, symprec=1e-3, angle_tolerance=5.0)
        return {
            "symbol": analyzer.get_space_group_symbol(),
            "number": int(analyzer.get_space_group_number()),
        }
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {exc}"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    records = []
    errors = []
    for run in manifest["runs"]:
        run_dir = root / run["path"]
        record = {**run}
        try:
            initial = Structure.from_file(run_dir / "POSCAR")
            vasprun = Vasprun(
                str(run_dir / "vasprun.xml"),
                parse_dos=False,
                parse_eigen=False,
                parse_projected_eigen=False,
                parse_potcar_file=False,
            )
            final = vasprun.final_structure
            forces = np.asarray(vasprun.ionic_steps[-1]["forces"], dtype=float)
            displacement = min_image_displacements(initial, final)
            metrics = lattice_metrics(initial, final)
            max_force = float(np.abs(forces).max())
            electronic_converged = bool(vasprun.converged)
            ionic_converged = bool(max_force <= 0.01)
            record.update(
                {
                    "electronic_converged": electronic_converged,
                    "ionic_converged": ionic_converged,
                    "converged": bool(electronic_converged and ionic_converged),
                    "final_energy_eV": float(vasprun.final_energy),
                    "delta_energy_from_static_eV": float(vasprun.final_energy - run["static_energy_eV"]),
                    "final_max_force_eV_A": max_force,
                    "n_ionic_steps": len(vasprun.ionic_steps),
                    "lattice": metrics,
                    "atomic_displacement": {
                        "max_A": float(np.linalg.norm(displacement, axis=1).max()),
                        "rms_A": float(np.sqrt(np.mean(np.sum(displacement**2, axis=1)))),
                        "per_atom_A": np.linalg.norm(displacement, axis=1).tolist(),
                    },
                    "symmetry": {
                        "initial": symmetry_label(initial),
                        "final": symmetry_label(final),
                    },
                    "final_structure_formula": final.composition.reduced_formula,
                    "final_structure_nsites": len(final),
                }
            )
            if not electronic_converged:
                errors.append(f"electronic_not_converged:{run['path']}")
            if not ionic_converged:
                errors.append(f"force_not_converged:{run['path']}:{max_force}")
            if not metrics["lattice_unchanged"]:
                errors.append(f"lattice_changed:{run['path']}")
            for filename in ("vasprun.xml", "OUTCAR", "vasp.out", "CONTCAR"):
                item = run_dir / filename
                if not item.is_file() or not item.stat().st_size:
                    errors.append(f"missing:{run['path']}/{filename}")
                else:
                    archive(item)
        except Exception as exc:
            errors.append(f"parse:{run['path']}:{type(exc).__name__}:{exc}")
            record["error"] = f"{type(exc).__name__}: {exc}"
        records.append(record)

    by_case = {}
    for record in records:
        by_case.setdefault(record["case"], {})[record["sign"]] = record
    summary = {}
    for case, signs in sorted(by_case.items()):
        minus = signs.get("minus", {})
        plus = signs.get("plus", {})
        summary[case] = {
            "q": minus.get("q", plus.get("q")),
            "frequency_THz": minus.get("frequency_THz", plus.get("frequency_THz")),
            "static_energy_eV": {
                "minus": minus.get("static_energy_eV"),
                "plus": plus.get("static_energy_eV"),
            },
            "final_energy_eV": {
                "minus": minus.get("final_energy_eV"),
                "plus": plus.get("final_energy_eV"),
            },
            "relaxation_delta_from_static_eV": {
                "minus": minus.get("delta_energy_from_static_eV"),
                "plus": plus.get("delta_energy_from_static_eV"),
            },
            "final_minus_plus_energy_difference_eV": (
                float(minus["final_energy_eV"] - plus["final_energy_eV"])
                if "final_energy_eV" in minus and "final_energy_eV" in plus
                else None
            ),
            "final_symmetry": {
                "minus": minus.get("symmetry", {}).get("final"),
                "plus": plus.get("symmetry", {}).get("final"),
            },
        }
    result = {
        "schema": "ti-v-softmode-relax-result/v1",
        "status": "pass" if not errors and len(records) == len(manifest["runs"]) else "failed",
        "production_consumable": False,
        "errors": errors,
        "records": records,
        "summary_by_case": summary,
        "interpretation": "Fixed-cell internal-coordinate diagnostic only; no production promotion.",
    }
    (root / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "errors": len(errors), "summary_by_case": summary}, indent=2))
    return 0 if result["status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
