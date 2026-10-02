#!/usr/bin/env python3
"""Post-process completed independent Ti-V static displacement jobs."""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
from phonopy import load
from pymatgen.io.vasp.outputs import Vasprun


def archive_copy(path: Path) -> str | None:
    """Create a compressed provenance copy without deleting the parsed file."""

    if not path.is_file():
        return None
    destination = path.with_name(path.name + ".gz")
    with path.open("rb") as source, gzip.open(destination, "wb", compresslevel=6) as target:
        while True:
            block = source.read(1024 * 1024)
            if not block:
                break
            target.write(block)
    return str(destination)


def process_volume(root: Path, label: str) -> dict:
    volume_dir = root / "volumes" / label
    dirs = [root / line.strip() for line in (volume_dir / "displacement_dirs.txt").read_text().splitlines() if line.strip()]
    forces = []
    for path in dirs:
        vasprun = path / "vasprun.xml"
        if not vasprun.is_file():
            raise FileNotFoundError(vasprun)
        task = Vasprun(str(vasprun), parse_dos=False, parse_eigen=False, parse_potcar_file=False)
        if not task.converged:
            raise RuntimeError(f"VASP not converged: {vasprun}")
        forces.append(np.asarray(task.ionic_steps[-1]["forces"], dtype=float))
    phonon = load(str(volume_dir / "phonopy.yaml"), produce_fc=False, is_nac=False)
    phonon.produce_force_constants(forces=forces, show_drift=True)
    phonon.save(str(volume_dir / "phonopy.yaml"), settings={"force_constants": True})
    phonon.run_mesh([20, 20, 20], is_gamma_center=True)
    phonon.run_total_dos()
    phonon.write_total_dos(str(volume_dir / "total_dos.dat"))
    # ``write_total_dos`` emits the two-column text DOS format. Do not give
    # that file a YAML suffix: downstream auditors must not mistake text DOS
    # for a phonopy YAML artifact.
    phonon.write_total_dos(str(volume_dir / "phonon_dos.dat"))
    paths = [
        [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.5, 0.25, 0.75], [0.0, 0.0, 0.0]],
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [0.5, 0.25, 0.75], [0.0, 0.0, 0.0]],
    ]
    phonon.run_band_structure(paths, labels=["Γ", "H", "N", "Γ", "Γ", "P", "N", "Γ"])
    phonon.write_yaml_band_structure(filename=str(volume_dir / "phonon_band_structure.yaml"))
    archived = []
    for path in dirs:
        run_dir = root / path
        for name in ("vasprun.xml", "OUTCAR", "vasp.out"):
            item = archive_copy(run_dir / name)
            if item:
                archived.append(item)
    return {
        "label": label,
        "displacements": len(dirs),
        "converged": len(forces),
        "archived_files": len(archived),
        "status": "pass",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    results = [process_volume(root, label) for label in ("low", "mid", "high")]
    payload = {"status": "pass", "volumes": results}
    (root / "postprocess.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    raise SystemExit(main())
