#!/usr/bin/env python3
"""Prepare an independent BCC-Nb-Ti magnetic initial-state search.

This script deliberately does not use htvasp/jobflow. It copies the archived
W90 geometry and exact POTCAR source, then creates standalone static VASP
inputs for Phase 1 magnetic initial-state tests.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shutil
from pathlib import Path
from typing import Any

from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar

STATES = {
    "fm": [3.0, 3.0],
    "afm_ti_pos_nb_neg": [3.0, -3.0],
    "afm_ti_neg_nb_pos": [-3.0, 3.0],
    "weak_fm": [0.5, 0.5],
    "weak_afm": [0.5, -0.5],
    "nm": [0.0, 0.0],
}

COMMON_INCAR = {
    "SYSTEM": "BCC-Nb-Ti magnetic initial-state search",
    "ENCUT": 520,
    "PREC": "Accurate",
    "ALGO": "Normal",
    "ISTART": 0,
    "ICHARG": 2,
    "ISPIN": 2,
    "ISMEAR": 1,
    "SIGMA": 0.2,
    "EDIFF": 1e-6,
    "NELM": 200,
    "NELMIN": 6,
    "NELMDL": -6,
    "IBRION": -1,
    "ISIF": 2,
    "NSW": 0,
    "ISYM": 0,
    "LREAL": False,
    "LASPH": True,
    "LMAXMIX": 4,
    "ENAUG": 1360,
    "LORBIT": 11,
    "LWAVE": False,
    "LCHARG": False,
    "KPAR": 1,
    "NCORE": 2,
    "NBANDS": 64,
    "GGA": "PE",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def potcar_titles(data: bytes) -> list[str]:
    text = data.decode("utf-8", errors="replace")
    return re.findall(r"TITEL\s*=\s*([^\n]+)", text)


def poscar_species(path: Path) -> list[str]:
    poscar = Poscar.from_file(path)
    return [str(x) for x in poscar.site_symbols]


def kpoints_signature(path: Path) -> dict[str, Any]:
    kp = Kpoints.from_file(path)
    return {
        "style": str(kp.style),
        "kpts": [list(map(int, row)) for row in (kp.kpts or [])],
        "num_kpts": kp.num_kpts,
    }


def write_incar(path: Path, moments: list[float]) -> None:
    values = dict(COMMON_INCAR)
    values["MAGMOM"] = moments
    Incar(values).write_file(path)


def render_wrapper(root: Path, states: list[str]) -> str:
    state_lines = "\n".join(f'  "{state}"' for state in states)
    return f'''#!/usr/bin/env bash
#SBATCH --job-name=nbti-mag
#SBATCH --partition=partCPU
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --no-requeue
#SBATCH --output={root}/logs/nbti-mag-%A_%a.out
#SBATCH --error={root}/logs/nbti-mag-%A_%a.err

set -euo pipefail
ROOT={root!s}
STATES=(
{state_lines}
)
INDEX="${{SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}}"
STATE="${{STATES[$INDEX]:?invalid array index $INDEX}}"
RUN_DIR="$ROOT/phase1/$STATE"
cd "$RUN_DIR"

for required in INCAR POSCAR KPOINTS POTCAR; do
  test -s "$required" || {{ echo "missing or empty $required" >&2; exit 2; }}
done

. /etc/profile.d/modules.sh
ulimit -s unlimited
module purge
module load vasp-cpu
command -v vasp_std
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
printf 'state=%s\\n' "$STATE" > run-environment.txt
printf 'job_id=%s\\n' "${{SLURM_JOB_ID:-unknown}}" >> run-environment.txt
printf 'array_id=%s\\n' "${{SLURM_ARRAY_TASK_ID:-unknown}}" >> run-environment.txt
printf 'host=%s\\n' "$(hostname)" >> run-environment.txt
printf 'started=%s\\n' "$(date --iso-8601=seconds)" >> run-environment.txt
srun --ntasks="$SLURM_NTASKS" --cpu-bind=cores vasp_std > vasp.out 2>&1
printf 'finished=%s\\n' "$(date --iso-8601=seconds)" >> run-environment.txt
'''


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--source-root", type=Path, default=Path("data/endmembers/BCC-Nb-Ti/tb2jflow/2-w90_static"))
    args = parser.parse_args()

    output = args.output.resolve()
    source_root = args.source_root.resolve()
    source_poscar = source_root / "POSCAR"
    source_potcar_gz = source_root / "POTCAR.gz"
    source_kpoints_gz = source_root / "KPOINTS.gz"
    for path in (source_poscar, source_potcar_gz, source_kpoints_gz):
        if not path.is_file():
            raise SystemExit(f"missing source artifact: {path}")

    if output.exists() and any(output.iterdir()):
        raise SystemExit(f"output must be new or empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    (output / "logs").mkdir()
    (output / "phase0").mkdir()
    (output / "phase1").mkdir()

    shutil.copy2(source_poscar, output / "phase0" / "POSCAR")
    with gzip.open(source_potcar_gz, "rb") as src, (output / "phase0" / "POTCAR").open("wb") as dst:
        shutil.copyfileobj(src, dst)
    with gzip.open(source_kpoints_gz, "rb") as src, (output / "phase0" / "KPOINTS").open("wb") as dst:
        shutil.copyfileobj(src, dst)

    species = poscar_species(output / "phase0" / "POSCAR")
    titles = potcar_titles((output / "phase0" / "POTCAR").read_bytes())
    if species != [re.search(r"PAW_PBE ([A-Za-z]+)", title).group(1) for title in titles]:
        raise SystemExit(f"POSCAR/POTCAR species mismatch: POSCAR={species}, POTCAR={titles}")
    if species != ["Ti", "Nb"]:
        raise SystemExit(f"unexpected source site order: {species}; expected ['Ti', 'Nb']")

    kp = kpoints_signature(output / "phase0" / "KPOINTS")
    if kp["kpts"] != [[8, 8, 8]] or "Gamma" not in kp["style"]:
        raise SystemExit(f"unexpected KPOINTS: {kp}")

    phase0_incar = dict(COMMON_INCAR)
    phase0_incar["MAGMOM"] = STATES["fm"]
    Incar(phase0_incar).write_file(output / "phase0" / "INCAR")

    phase1_records: list[dict[str, Any]] = []
    for state, moments in STATES.items():
        directory = output / "phase1" / state
        directory.mkdir()
        for name in ("POSCAR", "POTCAR", "KPOINTS"):
            shutil.copy2(output / "phase0" / name, directory / name)
        write_incar(directory / "INCAR", moments)
        record = {
            "state": state,
            "moments_muB": moments,
            "site_order": species,
            "purpose": "fixed-geometry static spin-polarized initial-state search",
            "input_files": {name: sha256(directory / name) for name in ("INCAR", "POSCAR", "KPOINTS", "POTCAR")},
        }
        (directory / "state.json").write_text(json.dumps(record, indent=2) + "\n")
        phase1_records.append(record)

    wrapper = output / "run_phase1_array.sbatch"
    wrapper.write_text(render_wrapper(output, list(STATES)))
    wrapper.chmod(0o755)
    submit = output / "submit_phase1.sh"
    submit.write_text(f'''#!/usr/bin/env bash\nset -euo pipefail\nROOT="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"\nmkdir -p "$ROOT/logs"\nexec sbatch --parsable "$ROOT/run_phase1_array.sbatch"\n''')
    submit.chmod(0o755)

    manifest = {
        "schema_version": "nbti-magnetic-search/v1",
        "created_date": "2026-09-11",
        "mode": "standalone_vasp_static_search",
        "source_root": source_root.as_posix(),
        "source_files": {
            "POSCAR": {"path": source_poscar.as_posix(), "sha256": sha256(source_poscar)},
            "POTCAR.gz": {"path": source_potcar_gz.as_posix(), "sha256": sha256(source_potcar_gz)},
            "KPOINTS.gz": {"path": source_kpoints_gz.as_posix(), "sha256": sha256(source_kpoints_gz)},
        },
        "potcar_titles": titles,
        "site_order": species,
        "phase0": {
            "geometry": "archived BCC-Nb-Ti W90 POSCAR; fixed geometry",
            "kpoints": kp,
            "incar": phase0_incar,
            "icharg_policy": "2: fresh atomic superposition; no CHGCAR dependency",
        },
        "phase1": {
            "states": phase1_records,
            "slurm": {
                "partition": "partCPU",
                "nodes": 1,
                "ntasks": 32,
                "ntasks_per_node": 32,
                "cpus_per_task": 1,
                "memory": "64G",
                "time": "02:00:00",
                "launcher": "srun --ntasks=$SLURM_NTASKS --cpu-bind=cores vasp_std",
                "module_init": ["module purge", "module load vasp-cpu"],
            },
        },
        "safety": {
            "source_tree_modified": False,
            "formal_tb2jflow_overwritten": False,
            "existing_results_overwritten": False,
        },
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (output / "README.md").write_text(f'''# BCC-Nb-Ti magnetic initial-state search\n\nGenerated 2026-09-11. This is an independent static VASP search; it does not use htvasp/jobflow and does not overwrite the archived `tb2jflow` result.\n\n- Geometry: archived `tb2jflow/2-w90_static/POSCAR`\n- POTCAR: exact decompressed `tb2jflow/2-w90_static/POTCAR.gz` (`Ti_sv`, `Nb_sv`)\n- KPOINTS: Gamma 8x8x8\n- Fresh SCF: `ISTART=0`, `ICHARG=2`, no CHGCAR dependency\n- States: {', '.join(STATES)}\n\nRun Phase 1 after audit: `./submit_phase1.sh`\n''')
    print(json.dumps({"output": output.as_posix(), "site_order": species, "potcar_titles": titles, "states": list(STATES), "wrapper": wrapper.as_posix()}, indent=2))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
