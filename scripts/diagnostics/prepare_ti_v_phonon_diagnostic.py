#!/usr/bin/env python3
"""Prepare and validate an independent three-volume Ti-V phonon rerun.

This diagnostic deliberately does not import or modify ``htvasp.workflows.qha``.
It reuses the archived primitive structures and POTCARs, then uses phonopy to
make a fresh 2x2x2 finite-displacement set for static VASP force calculations.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
import textwrap
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar


RECIPE_VERSION = "ti-v-phonon-diagnostic-v1-20260921"
ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "data/diagnostics/BCC-Ti-V-phonon-20260921"
SOURCE_ROOT = ROOT / "data/endmembers/BCC-Ti-V/qhaflow"
TARGETS = {
    "low": {
        "formula_volume": 26.116757852,
        "source": "1-phonon_static_1/12_eos_deformation_2",
    },
    "mid": {
        "formula_volume": 31.118716599,
        "source": "1-phonon_static_1/12_eos_deformation_6",
    },
    "high": {
        "formula_volume": 35.262763445,
        "source": "1-phonon_static_1/12_eos_deformation_9",
    },
}
SUPERCELL = np.diag([2, 2, 2])
COMMON_INCAR = {
    "SYSTEM": "Independent BCC-Ti-V three-volume phonon diagnostic",
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
    "ENAUG": 1360,
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


def command_output(command: list[str]) -> str | None:
    try:
        return subprocess.run(command, check=False, text=True, capture_output=True).stdout.strip()
    except OSError:
        return None


def protected_hashes() -> dict[str, str]:
    paths = [
        ROOT / "data/endmembers/BCC-Ti-V",
        ROOT / "htvasp/workflows/qha.py",
        ROOT / "main_qha.py",
    ]
    result: dict[str, str] = {}
    for path in paths:
        if path.is_file():
            result[str(path)] = sha256(path)
    return result


def structure_to_phonopy(structure: Structure) -> PhonopyAtoms:
    return PhonopyAtoms(
        symbols=[str(site.specie) for site in structure.sites],
        cell=np.array(structure.lattice.matrix),
        scaled_positions=np.array(structure.frac_coords),
    )


def phonopy_to_structure(atoms: PhonopyAtoms) -> Structure:
    return Structure(
        lattice=np.array(atoms.cell),
        species=list(atoms.symbols),
        coords=np.array(atoms.scaled_positions),
        coords_are_cartesian=False,
    )


def load_source(target: dict) -> tuple[Path, Structure, bytes, dict]:
    source = SOURCE_ROOT / target["source"]
    info_path = source / "phonon_info.json"
    potcar_path = source / "POTCAR.gz"
    poscar_path = source / "POSCAR.gz"
    kpoints_path = source / "KPOINTS.gz"
    for path in (info_path, potcar_path, poscar_path, kpoints_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    info = json.loads(info_path.read_text())
    primitive = Structure.from_dict(info["original_structure"])
    with gzip.open(potcar_path, "rb") as handle:
        potcar = handle.read()
    return source, primitive, potcar, {
        "phonon_info": str(info_path),
        "archived_poscar": str(poscar_path),
        "archived_kpoints": str(kpoints_path),
        "archived_potcar": str(potcar_path),
        "archived_poscar_sha256": sha256(poscar_path),
        "archived_kpoints_sha256": sha256(kpoints_path),
        "archived_potcar_gz_sha256": sha256(potcar_path),
        "primitive_source": "phonon_info.json:original_structure",
        "archived_displaced_poscar_not_reused": True,
    }


def prepare(output: Path) -> dict:
    if (output / "manifest.json").exists() or (output / "volumes").exists():
        raise RuntimeError(f"output already contains generated inputs: {output}")
    output.mkdir(parents=True, exist_ok=True)
    (output / "logs").mkdir()
    (output / "volumes").mkdir()

    volume_records = []
    displacement_records = []
    for label, target in TARGETS.items():
        source, primitive, potcar, source_meta = load_source(target)
        actual_formula_volume = float(primitive.volume)
        if abs(actual_formula_volume - target["formula_volume"]) > 2e-6:
            raise RuntimeError(
                f"{label}: source volume {actual_formula_volume} does not match "
                f"target {target['formula_volume']}"
            )
        volume_dir = output / "volumes" / label
        displacement_dir = volume_dir / "displacements"
        displacement_dir.mkdir(parents=True)
        phonon = Phonopy(
            structure_to_phonopy(primitive),
            supercell_matrix=SUPERCELL,
            primitive_matrix=np.eye(3),
            symprec=1e-3,
            # Match the production QHA recipe's ``sym_reduce=False``.  The
            # diagnostic must not silently reduce the 12 finite-displacement
            # force calculations to the 2 symmetry representatives.
            is_symmetry=False,
        )
        phonon.generate_displacements(distance=0.01, is_plusminus="auto", is_diagonal=True)
        phonon.save(str(volume_dir / "phonopy.yaml"), settings={"force_constants": False})
        Poscar(primitive).write_file(volume_dir / "primitive.POSCAR")
        Poscar(phonopy_to_structure(phonon.supercell)).write_file(volume_dir / "supercell.POSCAR")
        (volume_dir / "POTCAR").write_bytes(potcar)
        Kpoints.gamma_automatic(kpts=(4, 4, 4)).write_file(volume_dir / "KPOINTS")
        displacement_count = len(phonon.supercells_with_displacements)
        dirs = []
        for index, displaced in enumerate(phonon.supercells_with_displacements, start=1):
            run_dir = displacement_dir / f"disp-{index:04d}"
            run_dir.mkdir()
            Poscar(phonopy_to_structure(displaced)).write_file(run_dir / "POSCAR")
            Incar(COMMON_INCAR).write_file(run_dir / "INCAR")
            shutil.copy2(volume_dir / "KPOINTS", run_dir / "KPOINTS")
            shutil.copy2(volume_dir / "POTCAR", run_dir / "POTCAR")
            dirs.append(str(run_dir.relative_to(output)))
            displacement_records.append(
                {
                    "volume": label,
                    "index": index,
                    "path": str(run_dir.relative_to(output)),
                    "inputs_sha256": {name: sha256(run_dir / name) for name in ("INCAR", "POSCAR", "KPOINTS", "POTCAR")},
                }
            )
        (volume_dir / "displacement_dirs.txt").write_text("\n".join(dirs) + "\n")
        volume_records.append(
            {
                "label": label,
                "formula_volume_A3_per_fu": actual_formula_volume,
                # Phonopy 2.48 exposes the supercell lattice through ``cell``;
                # ``cell_volume`` is not a public attribute of PhonopyAtoms.
                "supercell_volume_A3": float(abs(np.linalg.det(np.asarray(phonon.supercell.cell, dtype=float)))),
                "target_formula_volume_A3_per_fu": target["formula_volume"],
                "source_qhaflow": str(source),
                "source": source_meta,
                "supercell_matrix": [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
                "primitive_atoms": len(primitive),
                "supercell_atoms": len(phonon.supercell),
                "displacement_count": displacement_count,
                "phonopy_yaml": str((volume_dir / "phonopy.yaml").relative_to(output)),
            }
        )

    (output / "all_displacement_dirs.txt").write_text(
        "\n".join(record["path"] for record in displacement_records) + "\n"
    )
    manifest = {
        "schema_version": "ti-v-phonon-diagnostic/v1",
        "recipe_version": RECIPE_VERSION,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "system": "BCC-Ti-V",
        "scope": "independent three-volume finite-displacement static force diagnostic",
        "software": {
            "python": sys.version,
            "phonopy": __import__("phonopy").__version__,
            "pymatgen": getattr(__import__("pymatgen"), "__version__", "unknown"),
            "atomate2": getattr(__import__("atomate2"), "__version__", "unknown"),
            "vasp_module": "vasp-cpu (loaded by submit script)",
            "sbatch": command_output(["sbatch", "--version"]),
        },
        "recipe": {
            "ENCUT_eV": 520,
            "LREAL": False,
            "ISYM": 0,
            "EDIFF_eV": 1e-7,
            "ISPIN": 2,
            "POTCAR": ["Ti_sv", "V_sv"],
            "supercell_matrix": [[2, 0, 0], [0, 2, 0], [0, 0, 2]],
            "kpoints": {"style": "Gamma", "mesh": [4, 4, 4]},
            "displacement_A": 0.01,
            "static_force_workflow": "phonopy finite displacement + independent VASP static force jobs",
            "postprocess_dependency": "phonopy force constants -> band structure and total DOS",
        },
        "volumes": volume_records,
        "displacements": displacement_records,
        "protected_paths_before": protected_hashes(),
        "production_paths_untouched": True,
        "submissions": None,
        "status": "prepared_not_submitted",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=True) + "\n")
    root_text = str(output)
    (output / "README.md").write_text(textwrap.dedent(f"""\
        # Independent BCC-Ti-V three-volume phonon diagnostic

        Created by `{RECIPE_VERSION}`. This directory is independent from
        `data/endmembers/BCC-Ti-V`; it is not a production QHA artifact.

        Volumes: low, mid, high. Each volume has fresh 2x2x2 finite-displacement
        phonon inputs and static VASP force jobs with `LREAL=False`, `ISYM=0`,
        `ENCUT=520 eV`, `EDIFF=1e-7 eV`, and Gamma `4x4x4` KPOINTS.

        Workflow:

        1. `python scripts/diagnostics/prepare_ti_v_phonon_diagnostic.py --output {root_text} --validate`
        2. `./submit.sh` (submits only the static array and dependent phonon postprocess)
        3. `postprocess.sbatch` writes force constants, band structure and DOS artifacts.

        The output is diagnostic-only until imaginary-mode, convergence and provenance
        gates pass. No QHA Gibbs curve is generated by this workflow.
        """), encoding="utf-8")
    (output / "run_static_array.sbatch").write_text(textwrap.dedent(f"""\
        #!/usr/bin/env bash
        #SBATCH --job-name=tiv-ph-static
        #SBATCH --partition=partCPU
        #SBATCH --nodes=1
        #SBATCH --ntasks=32
        #SBATCH --ntasks-per-node=32
        #SBATCH --cpus-per-task=1
        #SBATCH --mem=64G
        #SBATCH --time=04:00:00
        #SBATCH --no-requeue
        #SBATCH --output={root_text}/logs/static-%A_%a.out
        #SBATCH --error={root_text}/logs/static-%A_%a.err
        set -euo pipefail
        ROOT={root_text}
        mapfile -t RUNS < "$ROOT/all_displacement_dirs.txt"
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
        printf 'job_id=%s\\narray_id=%s\\nrun_dir=%s\\nhost=%s\\nstarted=%s\\n' "$SLURM_JOB_ID" "$SLURM_ARRAY_TASK_ID" "$RUN_DIR" "$(hostname)" "$(date --iso-8601=seconds)" > run-environment.txt
        srun --ntasks="$SLURM_NTASKS" --cpu-bind=cores vasp_std > vasp.out 2>&1
        printf 'finished=%s\\n' "$(date --iso-8601=seconds)" >> run-environment.txt
        """), encoding="utf-8")
    (output / "postprocess.sbatch").write_text(textwrap.dedent(f"""\
        #!/usr/bin/env bash
        #SBATCH --job-name=tiv-ph-post
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
        exec .conda/bin/python scripts/diagnostics/postprocess_ti_v_phonon_diagnostic.py --root {root_text}
        """), encoding="utf-8")
    (output / "submit.sh").write_text(textwrap.dedent(f"""\
        #!/usr/bin/env bash
        set -euo pipefail
        ROOT={root_text}
        cd {ROOT.parent.parent.parent}
        test -s "$ROOT/validation.json"
        grep -q '"status": "pass"' "$ROOT/validation.json"
        if test -s "$ROOT/submissions.json"; then
          echo "refusing duplicate submission: $ROOT/submissions.json" >&2
          exit 3
        fi
        N=$(wc -l < "$ROOT/all_displacement_dirs.txt")
        ARRAY_ID=$(sbatch --parsable --array=0-$((N-1)) "$ROOT/run_static_array.sbatch")
        POST_ID=$(sbatch --parsable --dependency=afterok:"$ARRAY_ID" "$ROOT/postprocess.sbatch")
        ROOT="$ROOT" ARRAY_ID="$ARRAY_ID" POST_ID="$POST_ID" N="$N" .conda/bin/python - <<'PY'
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
    required = ["manifest.json", "README.md", "submit.sh", "run_static_array.sbatch", "postprocess.sbatch"]
    errors.extend(str(output / name) for name in required if not (output / name).is_file())
    for record in manifest["displacements"]:
        path = output / record["path"]
        for name in ("INCAR", "POSCAR", "KPOINTS", "POTCAR"):
            if not (path / name).is_file() or (path / name).stat().st_size == 0:
                errors.append(str(path / name))
    if protected_hashes() != manifest["protected_paths_before"]:
        errors.append("protected path hash changed")
    result = {
        "status": "pass" if not errors else "fail",
        "errors": errors,
        "displacement_count": len(manifest["displacements"]),
        "production_paths_untouched": not errors or "protected path hash changed" not in errors,
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
