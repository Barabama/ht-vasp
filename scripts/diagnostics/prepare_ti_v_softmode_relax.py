#!/usr/bin/env python3
"""Prepare and submit fixed-cell Ti-V soft-mode internal relaxations.

This diagnostic consumes only the six +/- structures from the completed
2024-09-24 soft-mode static calculation.  It creates a new, independent
directory and never mutates production QHA or endmember data.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import textwrap
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar


PROJECT = Path("/nfs_hdd/2025/gaominliang/ht-vasp")
STATIC_ROOT = PROJECT / "data/diagnostics/BCC-Ti-V-softmode-static-20260924-v1"
STATIC_RESULT = STATIC_ROOT / "result.json"
OUTPUT = PROJECT / "data/diagnostics/BCC-Ti-V-softmode-relax-20260926-v1"
SCRIPT_ROOT = PROJECT / "scripts/diagnostics"
RECIPE = "ti-v-softmode-fixed-cell-internal-relax-v1-20260926"

TARGETS = (
    "low_q550/minus",
    "low_q550/plus",
    "mid_q550/minus",
    "mid_q550/plus",
    "mid_q555/minus",
    "mid_q555/plus",
)

INCAR_SETTINGS = {
    "SYSTEM": "BCC-Ti-V soft-mode fixed-cell internal relaxation diagnostic",
    "ENCUT": 520,
    "PREC": "Accurate",
    "ALGO": "Normal",
    "ISTART": 0,
    "ICHARG": 2,
    "ISPIN": 2,
    "ISMEAR": 1,
    "SIGMA": 0.2,
    "EDIFF": 1e-7,
    "EDIFFG": -0.01,
    "NELM": 200,
    "NELMIN": 6,
    "NELMDL": -6,
    "IBRION": 2,
    "ISIF": 2,
    "NSW": 120,
    "POTIM": 0.1,
    "ISYM": 0,
    "LREAL": False,
    "LASPH": True,
    "LMAXMIX": 4,
    "LORBIT": 11,
    "LWAVE": False,
    "LCHARG": False,
    "KPAR": 1,
    "NCORE": 2,
    "GGA": "PE",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def structure_signature(path: Path) -> dict:
    structure = Structure.from_file(path)
    return {
        "nsites": len(structure),
        "volume_A3": float(structure.volume),
        "lattice_matrix_A": np.asarray(structure.lattice.matrix, dtype=float).tolist(),
        "species": [str(site.specie) for site in structure],
    }


def static_records() -> dict[str, dict]:
    if not STATIC_RESULT.is_file():
        raise FileNotFoundError(STATIC_RESULT)
    payload = json.loads(STATIC_RESULT.read_text())
    records = {}
    for record in payload.get("records", []):
        key = f"{record['case']}/{record['sign']}"
        records[key] = record
    missing = [key for key in TARGETS if key not in records]
    if missing:
        raise RuntimeError(f"static result missing target records: {missing}")
    return records


def source_hashes() -> dict[str, str]:
    paths = [STATIC_RESULT, STATIC_ROOT / "manifest.json"]
    return {str(path): sha256(path) for path in paths if path.is_file()}


def render_array_script(root: Path) -> str:
    return textwrap.dedent(
        f"""\
        #!/usr/bin/env bash
        #SBATCH --job-name=tiv-soft-relax
        #SBATCH --partition=partCPU
        #SBATCH --nodes=1
        #SBATCH --ntasks=32
        #SBATCH --ntasks-per-node=32
        #SBATCH --cpus-per-task=1
        #SBATCH --mem=64G
        #SBATCH --time=04:00:00
        #SBATCH --no-requeue
        #SBATCH --output={root}/logs/relax-%A_%a.out
        #SBATCH --error={root}/logs/relax-%A_%a.err
        set -euo pipefail
        ROOT={root}
        mapfile -t RUNS < "$ROOT/all_run_dirs.txt"
        INDEX=${{SLURM_ARRAY_TASK_ID:?required array id}}
        RUN_DIR="$ROOT/${{RUNS[$INDEX]}}"
        cd "$RUN_DIR"
        for required in INCAR POSCAR KPOINTS POTCAR; do test -s "$required"; done
        . /etc/profile.d/modules.sh
        ulimit -s unlimited
        module purge
        module load vasp-cpu
        command -v vasp_std
        export OMP_NUM_THREADS=1
        export MKL_NUM_THREADS=1
        srun --ntasks="$SLURM_NTASKS" --cpu-bind=cores vasp_std > vasp.out 2>&1
        """
    )


def render_post_script(root: Path) -> str:
    return textwrap.dedent(
        f"""\
        #!/usr/bin/env bash
        #SBATCH --job-name=tiv-soft-relax-post
        #SBATCH --partition=partCPU
        #SBATCH --nodes=1
        #SBATCH --ntasks=1
        #SBATCH --cpus-per-task=4
        #SBATCH --mem=16G
        #SBATCH --time=01:00:00
        #SBATCH --no-requeue
        #SBATCH --output={root}/logs/post-%j.out
        #SBATCH --error={root}/logs/post-%j.err
        set -euo pipefail
        cd {PROJECT}
        exec {PROJECT}/.conda/bin/python scripts/diagnostics/postprocess_ti_v_softmode_relax.py --root {root}
        """
    )


def render_submit_script(root: Path) -> str:
    return textwrap.dedent(
        f"""\
        #!/usr/bin/env bash
        set -euo pipefail
        ROOT={root}
        cd {PROJECT}
        test -s "$ROOT/validation.json"
        grep -q '"status": "pass"' "$ROOT/validation.json"
        if test -s "$ROOT/submissions.json"; then
          echo "refusing duplicate submission: $ROOT/submissions.json" >&2
          exit 3
        fi
        N=$(wc -l < "$ROOT/all_run_dirs.txt")
        test "$N" -eq 6
        ARRAY_ID=$(sbatch --parsable --array=0-5%2 "$ROOT/run_relax_array.sbatch")
        POST_ID=$(sbatch --parsable --dependency=afterok:"$ARRAY_ID" "$ROOT/postprocess.sbatch")
        ROOT="$ROOT" ARRAY_ID="$ARRAY_ID" POST_ID="$POST_ID" N="$N" {PROJECT}/.conda/bin/python - <<'PY'
        import json, os
        from datetime import datetime, timezone
        from pathlib import Path
        root = Path(os.environ["ROOT"])
        payload = {{
            "submitted_at": datetime.now(timezone.utc).isoformat(),
            "static_array_job_id": os.environ["ARRAY_ID"],
            "postprocess_job_id": os.environ["POST_ID"],
            "array_size": int(os.environ["N"]),
            "array_spec": "0-5%2",
            "dependency": "afterok:" + os.environ["ARRAY_ID"],
            "relaxation_kind": "fixed_cell_internal_coordinate",
        }}
        (root / "submissions.json").write_text(json.dumps(payload, indent=2) + "\\n")
        manifest = json.loads((root / "manifest.json").read_text())
        manifest["submissions"] = payload
        manifest["status"] = "submitted_relax_array_with_dependent_postprocess"
        (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\\n")
        print(json.dumps(payload, indent=2))
        PY
        squeue --jobs="$ARRAY_ID,$POST_ID" --format='%.18i %.12P %.30j %.10T %.20R'
        """
    )


def prepare() -> None:
    if OUTPUT.exists():
        raise FileExistsError(
            f"refusing to overwrite existing diagnostic directory: {{OUTPUT}}"
        )
    records = static_records()
    OUTPUT.mkdir(parents=True)
    (OUTPUT / "logs").mkdir()
    run_records = []
    for key in TARGETS:
        static_record = records[key]
        source_run = STATIC_ROOT / static_record["path"]
        if not source_run.is_dir():
            raise FileNotFoundError(source_run)
        structure_path = source_run / "POSCAR"
        for filename in ("POSCAR", "KPOINTS", "POTCAR"):
            if not (source_run / filename).is_file():
                raise FileNotFoundError(source_run / filename)
        case, sign = key.split("/")
        run_dir = OUTPUT / "cases" / case / sign
        run_dir.mkdir(parents=True)
        shutil.copy2(structure_path, run_dir / "POSCAR")
        Kpoints.gamma_automatic(kpts=(4, 4, 4)).write_file(run_dir / "KPOINTS")
        shutil.copy2(source_run / "POTCAR", run_dir / "POTCAR")
        Incar(INCAR_SETTINGS).write_file(run_dir / "INCAR")
        initial = structure_signature(run_dir / "POSCAR")
        metadata = {
            "schema": "ti-v-softmode-relax-input/v1",
            "case": static_record["case"],
            "volume": static_record["volume"],
            "sign": static_record["sign"],
            "q": static_record["q"],
            "branch": static_record["branch"],
            "frequency_THz": static_record["frequency_THz"],
            "static_energy_eV": static_record["energy_eV"],
            "static_max_force_eV_A": static_record["max_force_eV_A"],
            "static_source_path": static_record["path"],
            "static_source_inputs_sha256": static_record["inputs_sha256"],
            "source_static_result_sha256": sha256(STATIC_RESULT),
            "initial_structure": initial,
            "relaxation": {
                "fixed_cell": True,
                "internal_coordinates_only": True,
                "ISIF": 2,
                "IBRION": 2,
                "NSW": 120,
                "EDIFFG_eV_A": -0.01,
                "POTIM_A": 0.1,
            },
        }
        (run_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
        input_hashes = {name: sha256(run_dir / name) for name in ("INCAR", "POSCAR", "KPOINTS", "POTCAR")}
        run_records.append(
            {
                "index": len(run_records),
                "case": static_record["case"],
                "volume": static_record["volume"],
                "sign": static_record["sign"],
                "q": static_record["q"],
                "branch": static_record["branch"],
                "frequency_THz": static_record["frequency_THz"],
                "static_energy_eV": static_record["energy_eV"],
                "static_max_force_eV_A": static_record["max_force_eV_A"],
                "static_source_path": static_record["path"],
                "path": str(run_dir.relative_to(OUTPUT)),
                "inputs_sha256": input_hashes,
                "initial_structure": initial,
            }
        )

    (OUTPUT / "all_run_dirs.txt").write_text("\n".join(item["path"] for item in run_records) + "\n")
    manifest = {
        "schema_version": "ti-v-softmode-relax/v1",
        "recipe_version": RECIPE,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "system": "BCC-Ti-V",
        "scope": "six commensurate soft-mode +/- fixed-cell internal-coordinate relaxations",
        "source_static_result": str(STATIC_RESULT),
        "source_hashes": source_hashes(),
        "software": {
            "python": sys.version,
            "pymatgen": getattr(__import__("pymatgen"), "__version__", "unknown"),
        },
        "recipe": {
            "INCAR": INCAR_SETTINGS,
            "KPOINTS": {"style": "Gamma", "mesh": [4, 4, 4]},
            "POTCAR": ["Ti_sv", "V_sv"],
            "fixed_cell": True,
            "internal_coordinates_only": True,
            "walltime": "04:00:00",
            "array": "0-5%2",
        },
        "runs": run_records,
        "production_paths_untouched": True,
        "status": "prepared_not_submitted",
    }
    (OUTPUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (OUTPUT / "README.md").write_text(
        textwrap.dedent(
            f"""\
            # BCC-Ti-V soft-mode fixed-cell internal relax diagnostic

            Recipe: `{RECIPE}`.

            This directory contains six independent fixed-cell internal-coordinate
            VASP relaxations: minus/plus structures for low_q550, mid_q550 and
            mid_q555.  The input POSCARs and provenance come from:

            `{STATIC_RESULT}`

            `ISIF=2` keeps the cell fixed while relaxing only internal coordinates.
            This is a diagnostic calculation only.  It does not modify production
            `data/endmembers/BCC-Ti-V`, `qha.py`, `main_qha.py`, or any TDB.

            Run `./submit.sh` only after `validation.json` reports `status=pass`.
            """
        )
    )
    for name, content in {
        "run_relax_array.sbatch": render_array_script(OUTPUT),
        "postprocess.sbatch": render_post_script(OUTPUT),
        "submit.sh": render_submit_script(OUTPUT),
    }.items():
        path = OUTPUT / name
        path.write_text(content)
        path.chmod(0o755)


def validate() -> dict:
    errors = []
    if not OUTPUT.is_dir():
        errors.append(f"missing_output:{OUTPUT}")
        return {"status": "fail", "errors": errors}
    try:
        manifest = json.loads((OUTPUT / "manifest.json").read_text())
    except Exception as exc:
        return {"status": "fail", "errors": [f"manifest:{exc}"]}
    if len(manifest.get("runs", [])) != 6:
        errors.append("run_count_not_6")
    expected_incar = {
        "ENCUT": 520,
        "EDIFF": 1e-7,
        "EDIFFG": -0.01,
        "IBRION": 2,
        "ISIF": 2,
        "NSW": 120,
        "POTIM": 0.1,
        "LREAL": False,
        "ISYM": 0,
        "ISPIN": 2,
    }
    for run in manifest.get("runs", []):
        run_dir = OUTPUT / run["path"]
        for name in ("INCAR", "POSCAR", "KPOINTS", "POTCAR", "metadata.json"):
            if not (run_dir / name).is_file() or not (run_dir / name).stat().st_size:
                errors.append(f"missing:{run['path']}/{name}")
        try:
            incar = Incar.from_file(run_dir / "INCAR")
            for key, expected in expected_incar.items():
                if incar.get(key) != expected:
                    errors.append(f"incar:{run['path']}:{key}={incar.get(key)!r} expected={expected!r}")
            kpoints = Kpoints.from_file(run_dir / "KPOINTS")
            if list(kpoints.kpts[0]) != [4, 4, 4]:
                errors.append(f"kpoints:{run['path']}:{kpoints.kpts}")
            structure = Structure.from_file(run_dir / "POSCAR")
            if len(structure) != run["initial_structure"]["nsites"]:
                errors.append(f"structure_nsites:{run['path']}")
            actual_hashes = {name: sha256(run_dir / name) for name in ("INCAR", "POSCAR", "KPOINTS", "POTCAR")}
            if actual_hashes != run["inputs_sha256"]:
                errors.append(f"hash_mismatch:{run['path']}")
            metadata = json.loads((run_dir / "metadata.json").read_text())
            if metadata["static_source_path"] != run["static_source_path"]:
                errors.append(f"metadata_source:{run['path']}")
        except Exception as exc:
            errors.append(f"validate:{run['path']}:{type(exc).__name__}:{exc}")
    for script in ("run_relax_array.sbatch", "postprocess.sbatch", "submit.sh"):
        if not (OUTPUT / script).is_file():
            errors.append(f"missing_script:{script}")
    status = "pass" if not errors else "fail"
    result = {
        "schema": "ti-v-softmode-relax-validation/v1",
        "status": status,
        "errors": errors,
        "run_count": len(manifest.get("runs", [])),
        "production_paths_untouched": True,
        "duplicate_submission_guard": not (OUTPUT / "submissions.json").exists(),
    }
    (OUTPUT / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--validate", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="prepare if needed, validate, and print submission command; never call sbatch")
    args = parser.parse_args()
    if not any((args.prepare, args.validate, args.dry_run)):
        parser.error("choose --prepare, --validate, or --dry-run")
    if args.prepare or args.dry_run:
        if not OUTPUT.exists():
            prepare()
        elif args.prepare and not (OUTPUT / "manifest.json").exists():
            raise RuntimeError(f"incomplete existing output directory: {OUTPUT}")
    if args.validate or args.dry_run:
        result = validate()
        print(json.dumps(result, indent=2))
        if result["status"] != "pass":
            return 2
    if args.dry_run:
        print(f"DRY-RUN: would submit {OUTPUT / 'submit.sh'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
