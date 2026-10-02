#!/usr/bin/env python3
"""Prepare and validate a minimal commensurate Ti-V soft-mode static test.

This diagnostic consumes only the completed, independent three-volume phonon
artifacts.  It creates reference and +/- eigenvector-displaced structures for
commensurate q-points representable by the existing 2x2x2 supercell.  It does
not alter the production QHA directory and does not perform relaxation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import textwrap
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from phonopy import load
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar


ROOT = Path("/nfs_hdd/2025/gaominliang/ht-vasp")
SOURCE = ROOT / "data/diagnostics/BCC-Ti-V-phonon-20260921-v3"
DEFAULT_OUTPUT = ROOT / "data/diagnostics/BCC-Ti-V-softmode-static-20260924-v1"
RECIPE = "ti-v-softmode-static-v1-20260924"

CASES = (
    {"label": "low_q550", "volume": "low", "q": [0.5, 0.5, 0.0], "branch": 0},
    {"label": "mid_q550", "volume": "mid", "q": [0.5, 0.5, 0.0], "branch": 0},
    {"label": "mid_q555", "volume": "mid", "q": [0.5, 0.5, 0.5], "branch": 0},
)
AMPLITUDE_A = 0.05

INCAR = {
    "SYSTEM": "BCC-Ti-V commensurate soft-mode static diagnostic",
    "ENCUT": 520,
    "PREC": "Accurate",
    "ALGO": "Normal",
    "ISTART": 0,
    "ICHARG": 2,
    "ISPIN": 2,
    "ISMEAR": 1,
    "SIGMA": 0.2,
    "EDIFF": 1e-7,
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
    "LORBIT": 11,
    "LWAVE": False,
    "LCHARG": False,
    "KPAR": 1,
    "NCORE": 2,
    "GGA": "PE",
    "MAGMOM": [3.0] * 16,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_hashes() -> dict[str, str]:
    paths = [
        SOURCE / "manifest.json",
        SOURCE / "postprocess.json",
        SOURCE / "static-audit.json",
        SOURCE / "volumes/low/phonopy.yaml",
        SOURCE / "volumes/mid/phonopy.yaml",
    ]
    return {str(path): sha256(path) for path in paths if path.is_file()}


def structure_from_atoms(atoms, cartesian: np.ndarray | None = None) -> Structure:
    positions = np.asarray(atoms.positions, dtype=float)
    if cartesian is not None:
        positions = cartesian
    return Structure(
        lattice=np.asarray(atoms.cell, dtype=float),
        species=list(atoms.symbols),
        coords=positions,
        coords_are_cartesian=True,
        to_unit_cell=True,
    )


def prepare(output: Path) -> dict:
    if output.exists() and any(output.iterdir()):
        raise RuntimeError(f"refusing to overwrite non-empty output: {output}")
    output.mkdir(parents=True, exist_ok=True)
    (output / "cases").mkdir()
    (output / "logs").mkdir()

    records = []
    case_records = []
    for case in CASES:
        volume = case["volume"]
        phonopy_yaml = SOURCE / "volumes" / volume / "phonopy.yaml"
        potcar = SOURCE / "volumes" / volume / "POTCAR"
        kpoints = SOURCE / "volumes" / volume / "KPOINTS"
        if not phonopy_yaml.is_file() or not potcar.is_file() or not kpoints.is_file():
            raise FileNotFoundError(f"missing phonon source for {case['label']}")

        phonon = load(str(phonopy_yaml), produce_fc=True)
        q = [float(x) for x in case["q"]]
        branch = int(case["branch"])
        phonon.run_qpoints([q], with_eigenvectors=True)
        qpoint_data = phonon.get_qpoints_dict()
        frequency = float(np.asarray(qpoint_data["frequencies"])[0, branch])

        # Phonopy's modulation amplitude is normalized by its mode convention.
        # Rescale the resulting real Cartesian modulation to an explicit 0.05 Å
        # maximum displacement so the VASP input has a directly auditable size.
        phonon.run_modulations([[2, 2, 2]], [[q, branch, 1.0, 0.0]])
        modulations, supercell = phonon.get_modulations_and_supercell()
        mode = np.real(np.asarray(modulations[0]))
        mode = np.asarray(mode, dtype=float)
        max_norm = float(np.linalg.norm(mode, axis=1).max())
        if not np.isfinite(max_norm) or max_norm <= 0.0:
            raise RuntimeError(f"zero modulation for {case['label']}")
        mode *= AMPLITUDE_A / max_norm
        reference = structure_from_atoms(supercell)

        case_dir = output / "cases" / case["label"]
        case_dir.mkdir(parents=True)
        Poscar(reference).write_file(case_dir / "reference.POSCAR")
        (case_dir / "phonopy_source.yaml").write_bytes(phonopy_yaml.read_bytes())
        (case_dir / "POTCAR").write_bytes(potcar.read_bytes())
        Kpoints.gamma_automatic(kpts=(4, 4, 4)).write_file(case_dir / "KPOINTS")
        case_record = {
            **case,
            "frequency_THz": frequency,
            "amplitude_max_A": AMPLITUDE_A,
            "mode_max_before_rescale_A": max_norm,
            "source_phonopy_yaml": str(phonopy_yaml),
            "structures": [],
        }

        for sign_name, sign in (("minus", -1.0), ("reference", 0.0), ("plus", 1.0)):
            structure = structure_from_atoms(supercell, np.asarray(supercell.positions) + sign * mode)
            run_dir = case_dir / sign_name
            run_dir.mkdir()
            Poscar(structure).write_file(run_dir / "POSCAR")
            Incar(INCAR).write_file(run_dir / "INCAR")
            (run_dir / "KPOINTS").write_bytes((case_dir / "KPOINTS").read_bytes())
            (run_dir / "POTCAR").write_bytes(potcar.read_bytes())
            metadata = {
                "case": case["label"],
                "volume": volume,
                "sign": sign_name,
                "sign_numeric": sign,
                "q": q,
                "branch": branch,
                "frequency_THz": frequency,
                "max_displacement_A": float(np.linalg.norm(sign * mode, axis=1).max()),
                "reference_volume_A3": float(reference.volume),
                "structure_nsites": len(reference),
                "source_phonopy_sha256": sha256(phonopy_yaml),
            }
            (run_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
            input_hashes = {name: sha256(run_dir / name) for name in ("INCAR", "POSCAR", "KPOINTS", "POTCAR")}
            records.append({
                "case": case["label"],
                "volume": volume,
                "sign": sign_name,
                "q": q,
                "branch": branch,
                "frequency_THz": frequency,
                "path": str(run_dir.relative_to(output)),
                "inputs_sha256": input_hashes,
            })
            case_record["structures"].append({"sign": sign_name, "path": str(run_dir.relative_to(output)), "inputs_sha256": input_hashes})
        case_records.append(case_record)

    (output / "all_run_dirs.txt").write_text("\n".join(item["path"] for item in records) + "\n")
    manifest = {
        "schema_version": "ti-v-softmode-static/v1",
        "recipe_version": RECIPE,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "system": "BCC-Ti-V",
        "scope": "commensurate soft-mode +/- static energy diagnostic",
        "source_directory": str(SOURCE),
        "source_hashes": source_hashes(),
        "software": {
            "python": sys.version,
            "phonopy": __import__("phonopy").__version__,
            "pymatgen": getattr(__import__("pymatgen"), "__version__", "unknown"),
        },
        "recipe": {
            "q_points_commensurate_with_supercell": True,
            "supercell_matrix": [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
            "amplitude_max_A": AMPLITUDE_A,
            "kpoints": {"style": "Gamma", "mesh": [4, 4, 4]},
            "INCAR": {key: value for key, value in INCAR.items() if key not in {"MAGMOM"}},
            "POTCAR": ["Ti_sv", "V_sv"],
            "static_only": True,
            "relaxation_submitted": False,
        },
        "cases": case_records,
        "runs": records,
        "production_paths_untouched": True,
        "status": "prepared_not_submitted",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    root_text = str(output)
    (output / "README.md").write_text(textwrap.dedent(f"""\
        # BCC-Ti-V commensurate soft-mode static diagnostic

        Recipe: `{RECIPE}`. This is a static-only diagnostic; no relaxation,
        QHA or TDB generation is performed. The source phonon artifacts and
        production `data/endmembers/BCC-Ti-V` directory are read-only inputs.

        Cases: low q=(0.5,0.5,0), mid q=(0.5,0.5,0), and mid q=(0.5,0.5,0.5),
        each with reference/minus/plus structures. The maximum Cartesian
        displacement is 0.05 Angstrom and the q points are commensurate with
        the existing 2x2x2 supercell.

        Run `./submit.sh` only after `validation.json` reports pass.
        """), encoding="utf-8")
    (output / "run_static_array.sbatch").write_text(textwrap.dedent(f"""\
        #!/usr/bin/env bash
        #SBATCH --job-name=tiv-soft-static
        #SBATCH --partition=partCPU
        #SBATCH --nodes=1
        #SBATCH --ntasks=32
        #SBATCH --ntasks-per-node=32
        #SBATCH --cpus-per-task=1
        #SBATCH --mem=64G
        #SBATCH --time=02:00:00
        #SBATCH --no-requeue
        #SBATCH --output={root_text}/logs/static-%A_%a.out
        #SBATCH --error={root_text}/logs/static-%A_%a.err
        set -euo pipefail
        ROOT={root_text}
        mapfile -t RUNS < "$ROOT/all_run_dirs.txt"
        RUN_DIR="$ROOT/${{RUNS[${{SLURM_ARRAY_TASK_ID:?required array id}}]}}"
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
        """), encoding="utf-8")
    (output / "postprocess.sbatch").write_text(textwrap.dedent(f"""\
        #!/usr/bin/env bash
        #SBATCH --job-name=tiv-soft-post
        #SBATCH --partition=partCPU
        #SBATCH --nodes=1
        #SBATCH --ntasks=1
        #SBATCH --cpus-per-task=4
        #SBATCH --mem=16G
        #SBATCH --time=01:00:00
        #SBATCH --no-requeue
        #SBATCH --output={root_text}/logs/post-%j.out
        #SBATCH --error={root_text}/logs/post-%j.err
        set -euo pipefail
        cd {ROOT}
        exec /nfs_hdd/2025/gaominliang/ht-vasp/.conda/bin/python scripts/diagnostics/postprocess_ti_v_softmode_static.py --root {root_text}
        """), encoding="utf-8")
    (output / "submit.sh").write_text(textwrap.dedent(f"""\
        #!/usr/bin/env bash
        set -euo pipefail
        ROOT={root_text}
        cd {ROOT}
        test -s "$ROOT/validation.json"
        grep -q '"status": "pass"' "$ROOT/validation.json"
        if test -s "$ROOT/submissions.json"; then
          echo "refusing duplicate submission: $ROOT/submissions.json" >&2
          exit 3
        fi
        N=$(wc -l < "$ROOT/all_run_dirs.txt")
        ARRAY_ID=$(sbatch --parsable --array=0-$((N-1)) "$ROOT/run_static_array.sbatch")
        POST_ID=$(sbatch --parsable --dependency=afterok:"$ARRAY_ID" "$ROOT/postprocess.sbatch")
        ROOT="$ROOT" ARRAY_ID="$ARRAY_ID" POST_ID="$POST_ID" N="$N" /nfs_hdd/2025/gaominliang/ht-vasp/.conda/bin/python - <<'PY'
        import json, os
        from datetime import datetime, timezone
        from pathlib import Path
        root = Path(os.environ["ROOT"])
        payload = {{
            "submitted_at": datetime.now(timezone.utc).isoformat(),
            "static_array_job_id": os.environ["ARRAY_ID"],
            "postprocess_job_id": os.environ["POST_ID"],
            "array_size": int(os.environ["N"]),
            "dependency": "afterok:" + os.environ["ARRAY_ID"],
            "qha_submitted": False,
            "relaxation_submitted": False,
        }}
        (root / "submissions.json").write_text(json.dumps(payload, indent=2) + "\\n")
        manifest = json.loads((root / "manifest.json").read_text())
        manifest["submissions"] = payload
        manifest["status"] = "submitted_static_array_with_dependent_postprocess"
        (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\\n")
        print(json.dumps(payload, indent=2))
        PY
        squeue --jobs="$ARRAY_ID,$POST_ID" --format='%.18i %.12P %.28j %.10T %.20R'
        """), encoding="utf-8")
    for path in (output / "run_static_array.sbatch", output / "postprocess.sbatch", output / "submit.sh"):
        path.chmod(0o755)
    return manifest


def validate(output: Path) -> dict:
    manifest = json.loads((output / "manifest.json").read_text())
    errors = []
    for run in manifest["runs"]:
        path = output / run["path"]
        for name in ("INCAR", "POSCAR", "KPOINTS", "POTCAR", "metadata.json"):
            if not (path / name).is_file() or not (path / name).stat().st_size:
                errors.append(str(path / name))
        for name, expected in run["inputs_sha256"].items():
            if sha256(path / name) != expected:
                errors.append(f"hash:{path / name}")
    result = {
        "status": "pass" if not errors else "fail",
        "errors": errors,
        "run_count": len(manifest["runs"]),
        "production_paths_untouched": True,
    }
    (output / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--validate", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    if args.prepare:
        prepare(output)
        return 0
    if args.validate:
        result = validate(output)
        print(json.dumps(result, indent=2))
        return 0 if result["status"] == "pass" else 2
    parser.error("choose --prepare or --validate")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
