import json
import time
import shutil
import logging
import argparse
import shlex
import os
import hashlib
from dataclasses import asdict, is_dataclass
from pathlib import Path
from pymatgen.core import Structure
from htvasp.model import Endmember
from htvasp.workflows import QhaWorker, StaticWorker
from htvasp.slurm import (
    SlurmJobManager,
    artifact_marker_matches,
    current_artifact_matches,
    mark_artifact_from_ledger,
    new_run_id,
    payload_fingerprint,
    stable_run_id,
)

QHA_RECIPE_VERSION = "main_qha-v3-reverse-endmember-b2"

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

VASP_ARGS = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-cpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-cpu && srun vasp_gam'",
}

# Initial magnetic moments per element (muB) - from verified submit_tb2j.sh
# 13-element initial guess table (restored; shared with main_qha.py)
MAGMOM = {
    "Al": 1.0,
    "Co": 3.0,
    "Cr": 3.0,
    "Cu": 5.0,
    "Fe": 5.0,
    "Mn": 5.0,
    "Nb": 3.0,
    "Ni": 2.0,
    "Ta": 3.0,
    "Ti": 3.0,
    "V": 3.0,
    "W": 3.0,
    "Zr": 1.0,
}
GLOBAL_INCAR = {
    "ENCUT": 520,
    # Electronic
    "ISMEAR": 1,
    "SIGMA": 0.2,
    "ALGO": "Normal",
    "NELM": 120,
    "NELMIN": 6,
    "NELMDL": -6,
    # Ionic
    "IBRION": 2,
    "ISIF": 3,
    "NSW": 80,
    "POTIM": 0.2,
    "EDIFF": 1e-6,
    "EDIFFG": -0.02,
    # Magnetism
    "ISPIN": 2,
    "MAGMOM": MAGMOM,
    "ISYM": 0,
    "LREAL": "Auto",
    "PREC": "Accurate",
    "SYMPREC": 1e-5,
    # Output
    "LWAVE": False,
    "LCHARG": False,
    "LOPTICS": False,
    "LVTOT": False,
    "GGA": "PE",
    # "AMIX": 0.1,
    # "BMIX": 1e-4,
    # "AMIX_MAG": 0.4,
    # "BMIX_MAG": 1e-4,
}

RELAX_INCAR = {"KPAR": 4, "NCORE": 2, "LORBIT": 10}
PHONON_INCAR = {"KPAR": 2, "NCORE": 4}
STATIC_INCAR = {"LORBIT": 11}

STRUCTURE_NAMES = [
    # Al-Nb-Ti-V
    "SER-Al",
    "SER-Nb",
    "SER-Ti",
    "SER-V",
    "BCC-Al-Al",
    "BCC-Al-Nb",
    "BCC-Al-Ti",
    "BCC-Al-V",
    "BCC-Nb-Nb",
    "BCC-Nb-Ti",
    "BCC-Nb-V",
    "BCC-Ti-Ti",
    "BCC-Ti-V",
    "BCC-V-V",

    # Co-Fe-Mn-Ni
    "SER-Co",
    "SER-Fe",
    "SER-Mn",
    "SER-Ni",
    "BCC-Co-Co",
    "BCC-Co-Fe",
    "BCC-Co-Mn",
    "BCC-Co-Ni",
    "BCC-Fe-Fe",
    "BCC-Fe-Mn",
    "BCC-Fe-Ni",
    "BCC-Mn-Mn",
    "BCC-Mn-Ni",
    "BCC-Ni-Ni",
    "FCC-Co-Co",
    "FCC-Co-Fe",
    "FCC-Co-Mn",
    "FCC-Co-Ni",
    "FCC-Fe-Co",
    "FCC-Fe-Fe",
    "FCC-Fe-Mn",
    "FCC-Fe-Ni",
    "FCC-Mn-Co",
    "FCC-Mn-Fe",
    "FCC-Mn-Mn",
    "FCC-Mn-Ni",
    "FCC-Ni-Co",
    "FCC-Ni-Fe",
    "FCC-Ni-Mn",
    "FCC-Ni-Ni",
]


# Keep the canonical forward targets and the derived reverse targets explicit.
# Reverse targets are deliberately not appended to STRUCTURE_NAMES: the legacy
# no-argument submission behavior must not silently expand.
FORWARD_B2_QHA_TARGETS = [
    "BCC-Al-Nb",
    "BCC-Al-Ti",
    "BCC-Al-V",
    "BCC-Nb-Ti",
    "BCC-Nb-V",
    "BCC-Ti-V",
]
REVERSE_B2_QHA_TARGETS = [
    "BCC-Nb-Al",
    "BCC-Ti-Al",
    "BCC-V-Al",
    "BCC-Ti-Nb",
    "BCC-V-Nb",
    "BCC-V-Ti",
]
B2_QHA_TARGETS = FORWARD_B2_QHA_TARGETS + REVERSE_B2_QHA_TARGETS
SUBMITTABLE_TARGETS = tuple(dict.fromkeys([*STRUCTURE_NAMES, *B2_QHA_TARGETS]))

root_dir = Path(__file__).parent.resolve()
flow_dir = Path("/nfs_ssd/tmp")
conda_env = root_dir / ".conda"
B2_REVERSE_MANIFEST = root_dir / "data" / "poscars" / "b2_reverse_manifest.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fractional_site_index(structure: Structure, coordinate: tuple[float, float, float]) -> int:
    matches = []
    for index, site in enumerate(structure):
        deltas = []
        for actual, expected in zip(site.frac_coords, coordinate):
            delta = abs((float(actual) % 1.0) - expected)
            deltas.append(min(delta, 1.0 - delta))
        if max(deltas) <= 1e-7:
            matches.append(index)
    if len(matches) != 1:
        raise ValueError(f"expected one B2 site at {coordinate}, found {matches}")
    return matches[0]


def _b2_phase_context(name: str) -> dict | None:
    """Return audited B2 metadata without inferring sublattices from names."""

    if name not in B2_QHA_TARGETS:
        return None
    poscar_path = root_dir / "data" / "poscars" / f"{name}.vasp"
    if not poscar_path.is_file():
        raise FileNotFoundError(f"missing B2 POSCAR for {name}: {poscar_path}")
    if name in REVERSE_B2_QHA_TARGETS:
        if not B2_REVERSE_MANIFEST.is_file():
            raise FileNotFoundError(f"missing reverse POSCAR manifest: {B2_REVERSE_MANIFEST}")
        manifest = json.loads(B2_REVERSE_MANIFEST.read_text(encoding="utf-8"))
        record = next((item for item in manifest.get("targets", []) if item.get("reverse_name") == name), None)
        if record is None:
            raise ValueError(f"reverse target {name} is absent from {B2_REVERSE_MANIFEST}")
        mapping = {
            key: {
                "fractional_coordinate": value["fractional_coordinate"],
                "element": value["derived_element"],
                "source_element": value["source_element"],
                "site_index": value["derived_site_index"],
            }
            for key, value in record["site_mapping"].items()
        }
        return {
            "canonical_phase": "B2_BCC",
            "ordered": True,
            "reverse_of": record["reverse_of"],
            "sublattice_mapping": mapping,
            "source_poscar": str((root_dir / record["source_path"]).resolve()),
            "source_sha256": record["source_sha256"],
            "derived_poscar": str((root_dir / record["derived_path"]).resolve()),
            "derived_sha256": record["derived_sha256"],
            "generation_rule": record["generation_rule"],
        }

    structure = Structure.from_file(poscar_path)
    if len(structure) != 2:
        raise ValueError(f"B2 POSCAR must have two sites: {poscar_path}")
    origin = _fractional_site_index(structure, (0.0, 0.0, 0.0))
    body = _fractional_site_index(structure, (0.5, 0.5, 0.5))
    return {
        "canonical_phase": "B2_BCC",
        "ordered": True,
        "sublattice_mapping": {
            "sublattice_1": {
                "fractional_coordinate": [0.0, 0.0, 0.0],
                "element": str(structure[origin].specie),
                "site_index": origin,
            },
            "sublattice_2": {
                "fractional_coordinate": [0.5, 0.5, 0.5],
                "element": str(structure[body].specie),
                "site_index": body,
            },
        },
        "source_poscar": str(poscar_path.resolve()),
        "source_sha256": _sha256(poscar_path),
        "generation_rule": "audited B2 POSCAR; no site exchange",
    }


def resolve_targets(name: str | None = None, targets: list[str] | None = None) -> list[str]:
    """Resolve a single name, an explicit list, or legacy all-target behavior."""

    if name and targets:
        raise ValueError("--name and --targets are mutually exclusive")
    requested = [name] if name else (list(targets) if targets else list(STRUCTURE_NAMES))
    if not requested:
        raise ValueError("at least one target is required")
    if len(set(requested)) != len(requested):
        raise ValueError("submission targets must be unique")
    unknown = [item for item in requested if item not in SUBMITTABLE_TARGETS]
    if unknown:
        raise ValueError(f"unknown QHA target(s): {', '.join(unknown)}")
    return requested


def _stage_paths(name: str, stage: str, run_id: str | None = None):
    base = root_dir / "data" / "endmembers" / name / f"{stage}flow"
    if run_id:
        if Path(run_id).name != run_id or not run_id.strip():
            raise ValueError(f"Invalid run_id: {run_id!r}")
        run_root = base / "runs" / run_id
        return run_root, run_root / f"{name}-{stage}.json", base / "slurm_submissions.json"
    return base, base.parent / f"{name}-{stage}.json", base / "slurm_submissions.json"


def run_static(name: str, rerun: bool = False, run_id: str | None = None):
    if rerun:
        run_id = run_id or new_run_id("qha-static")
    static_store, static_json, _ = _stage_paths(name, "static", run_id)

    if not rerun and static_json.exists() and (
        (run_id is None and not os.environ.get("HTVASP_INTENT_KEY"))
        or current_artifact_matches(_stage_paths(name, "static", run_id)[2], static_json)
    ):
        log.info(f"StaticWorker for Structure {name} already done")
        return

    log.info(f"StaticWorker for Structure {name} started")

    endmember = Endmember()
    structure = endmember.get_poscar(name, Path(root_dir) / "data" / "poscars")
    try:
        static_worker = StaticWorker(
            vasp_args=VASP_ARGS,
            potcar_functional="PBE",
            global_incar=GLOBAL_INCAR,
            relax_incar=RELAX_INCAR,
            static_incar=STATIC_INCAR,
        )
        static_worker.run_flow(
            name=name,
            structure=structure,
            flow_dir=flow_dir,
            store_dir=static_store,
            resume=not rerun,
        )
        static_output = static_worker.get_result() or {}
        static_worker.write_result(data=static_output, json_path=static_json)
        mark_artifact_from_ledger(
            _stage_paths(name, "static", run_id)[2], static_json
        )
    except Exception as e:
        log.error(f"StaticWorker for Structure {name} failed")
        log.error(e)

    log.info(f"StaticWorker for Structure {name} done")
    return


def run_qha(name: str, rerun: bool = False, run_id: str | None = None):
    if rerun:
        run_id = run_id or new_run_id("qha-qha")
    qha_store, qha_json, _ = _stage_paths(name, "qha", run_id)

    if not rerun and qha_json.exists() and (
        (run_id is None and not os.environ.get("HTVASP_INTENT_KEY"))
        or current_artifact_matches(_stage_paths(name, "qha", run_id)[2], qha_json)
    ):
        log.info(f"QhaWorker for Structure {name} already done")
        return

    log.info(f"QhaWorker for Structure {name} started")

    endmember = Endmember()
    structure = endmember.get_poscar(name, Path(root_dir) / "data" / "poscars")
    try:
        qha_worker = QhaWorker(
            vasp_args=VASP_ARGS,
            potcar_functional="PBE",
            global_incar=GLOBAL_INCAR,
            relax_incar=RELAX_INCAR,
            eos_incar=RELAX_INCAR,
            phonon_incar=PHONON_INCAR,
            temperature_range=(0, 3100, 50),
            supercell_matrix=((2, 0, 0), (0, 2, 0), (0, 0, 2)),
            reference_name=name,
            phase_context=_b2_phase_context(name),
        )
        qha_worker.run_flow(
            name=name,
            structure=structure,
            flow_dir=flow_dir,
            store_dir=qha_store,
            resume=not rerun,
        )
        qha_output = qha_worker.get_result() or {}
        qha_worker.write_result(data=qha_output, json_path=qha_json)
        mark_artifact_from_ledger(_stage_paths(name, "qha", run_id)[2], qha_json)
    except Exception as e:
        log.error(f"QhaWorker for Structure {name} failed")
        log.error(e)

    log.info(f"QhaWorker for Structure {name} done")
    return


def check_jobs(job: str):

    for name in STRUCTURE_NAMES:
        store_dir = root_dir / "data" / "endmembers" / name / f"{job}flow"
        json_path = store_dir.parent / f"{name}-{job}.json"
        if not json_path.exists():
            log.info(f"{name}: not run")
            continue
        with open(json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        # print(data.keys())
        if job == "static":
            energy = data["output"]["energy"]
            total_magnetization = data["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"]
            magnetization = data["calcs_reversed"][0]["output"]["outcar"]["magnetization"]
            log.info(f"{name}: {energy}")               # 总能量
            log.info(f"{name}: {total_magnetization}")  # 总磁矩
            # log.info(f"{magnetization}\n")              # 轨道分解 [s,p,d,tot]
        elif job == "qha":
            temp_range = (data["temperatures"][0], data["temperatures"][-1])
            len_temp = len(data["temperatures"])
            # deformations = data["deformation_energies"]
            log.info(f"{name}: {temp_range} ({len_temp} temps)")
            # log.info(f"{name}: {deformations}")




def _submission_plan(name: str, job: str, manager: SlurmJobManager, rerun: bool = False) -> dict:
    stage = job.lower()
    phase_context = _b2_phase_context(name)
    run_payload = {
        "system": name,
        "stage": stage,
        "recipe": QHA_RECIPE_VERSION,
        "phase_context": phase_context,
        "settings": {
            "vasp_args": VASP_ARGS,
            "global": GLOBAL_INCAR,
            "relax": RELAX_INCAR,
            "static": STATIC_INCAR,
            "phonon": PHONON_INCAR,
            "temperature_range": (0, 3100, 50),
            "supercell_matrix": ((2, 0, 0), (0, 2, 0), (0, 0, 2)),
        },
    }
    run_id = new_run_id(f"qha-{stage}") if rerun else stable_run_id(f"qha-{stage}", run_payload)
    run_root, result_path, ledger_path = _stage_paths(name, stage, run_id)
    command = (
        f"{shlex.quote(str(root_dir / '.conda' / 'bin' / 'python'))} "
        f"{shlex.quote(str(Path(__file__).resolve()))} -n {shlex.quote(name)} -j {shlex.quote(job)}"
    )
    command += f" --run-id {shlex.quote(run_id)}"
    if rerun:
        command += " --rerun"
    config = manager.get_cpu_config(
        job_name=f"{name}-{job}",
        output_log=f"logs/{name}-{job}.log",
        ntasks=32,
        memory="64G",
        module_name="vasp-cpu",
    )
    settings_fp = payload_fingerprint(run_payload["settings"])
    manifest_fp = payload_fingerprint({
        "system": name,
        "stage": stage,
        "recipe": QHA_RECIPE_VERSION,
        "phase_context": phase_context,
    })
    return {
        "name": name,
        "stage": stage,
        "recipe_version": QHA_RECIPE_VERSION,
        "run_id": run_id,
        "command": command,
        "config": config,
        "config_dict": asdict(config) if is_dataclass(config) else dict(vars(config)),
        "artifact_path": result_path,
        "run_root": run_root,
        "ledger_path": ledger_path,
        "settings_fingerprint": settings_fp,
        "manifest_fingerprint": manifest_fp,
        "phase_context": phase_context,
    }


def _plan_conflict(plan: dict, rerun: bool) -> str | None:
    """Reject only artifacts that this invocation could overwrite.

    An existing run directory without its final JSON is a valid resume target;
    submit_intent remains the authority for active/unknown ledger states.
    """

    if rerun or plan["stage"] != "qha":
        # Preserve the established static-entrypoint behavior.  The QHA
        # reverse-endmember safety gate below is intentionally scoped to the
        # artifacts this task is adding; static is a separate legacy stage.
        return None
    final_path = _stage_paths(plan["name"], plan["stage"])[1]
    if final_path.is_file():
        return f"existing canonical artifact: {final_path}"
    artifact_path = Path(plan["artifact_path"])
    if artifact_path.is_file() and not artifact_marker_matches(artifact_path, {
        "intent_key": "never-match",
        "system_name": plan["name"],
        "stage": plan["stage"],
        "recipe_version": plan["recipe_version"],
        "settings_fingerprint": plan["settings_fingerprint"],
        "manifest_fingerprint": plan["manifest_fingerprint"],
    }):
        return f"existing unmarked run artifact: {artifact_path}"
    return None


def _print_submission_plan(plan: dict, *, conflict: str | None = None) -> None:
    payload = {
        "target": plan["name"],
        "stage": plan["stage"],
        "recipe_version": plan["recipe_version"],
        "run_id": plan["run_id"],
        "command": plan["command"],
        "artifact_path": str(Path(plan["artifact_path"]).resolve()),
        "run_root": str(Path(plan["run_root"]).resolve()),
        "ledger_path": str(Path(plan["ledger_path"]).resolve()),
        "resources": plan["config_dict"],
        "phase_context": plan["phase_context"],
    }
    if conflict:
        payload["conflict"] = conflict
    print(json.dumps(payload, ensure_ascii=False, indent=2, default=str))


def submit_job(
    name: str,
    job: str,
    manager: SlurmJobManager,
    rerun: bool = False,
    dry_run: bool = False,
):
    plan = _submission_plan(name, job, manager, rerun=rerun)
    conflict = _plan_conflict(plan, rerun=rerun)
    if dry_run:
        _print_submission_plan(plan, conflict=conflict)
        return plan
    if conflict:
        log.error("Refusing QHA submission for %s: %s", name, conflict)
        return -1

    result_path = Path(plan["artifact_path"])
    result = manager.submit_intent(
        system_name=plan["name"],
        stage=plan["stage"],
        recipe_version=plan["recipe_version"],
        settings_fingerprint=plan["settings_fingerprint"],
        manifest_fingerprint=plan["manifest_fingerprint"],
        command=plan["command"],
        config=plan["config"],
        ledger_path=plan["ledger_path"],
        conda_env=str(conda_env),
        work_dir=root_dir,
        run_id=plan["run_id"],
        retry=rerun,
        artifact_path=result_path,
        artifact_checker=lambda entry: artifact_marker_matches(result_path, entry),
        artifact_shortcut=True,
    )
    if not result.ok:
        log.error("Refused qha job for Structure %s: %s", name, result.state)
        return -1
    log.info("QHA job for Structure %s: %s (%s)", name, result.state, result.job_id)
    return result.job_id


def submit_jobs(
    job: str,
    names: list[str] | None = None,
    rerun: bool = False,
    dry_run: bool = False,
):
    manager = SlurmJobManager()
    selected = resolve_targets(targets=names) if names is not None else list(STRUCTURE_NAMES)
    return [submit_job(name=name, job=job, manager=manager, rerun=rerun, dry_run=dry_run) for name in selected]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Qha workflow")
    parser.add_argument("--name", "-n", type=str, help="Run or submit one structure")
    parser.add_argument("--targets", nargs="+", help="Explicit target list for Slurm submission")
    parser.add_argument("--job", "-j", default="static", type=str, help="Type of job name")
    parser.add_argument("--check", "-c", action="store_true", help="Check job json")
    parser.add_argument("--slurm", "-s", action="store_true", help="Submit to Slurm")
    parser.add_argument("--dry-run", action="store_true", help="Print submission plans without calling sbatch")
    parser.add_argument("--rerun", "-r", action="store_true", help="Re-run calculations")
    parser.add_argument("--run-id", help="Operate on an explicit immutable run")
    args = parser.parse_args()

    if args.name and args.targets:
        parser.error("--name and --targets are mutually exclusive")
    if args.targets and not (args.slurm or args.dry_run):
        parser.error("--targets requires --slurm or --dry-run")
    if args.dry_run or args.slurm:
        try:
            selected = resolve_targets(name=args.name, targets=args.targets)
            submit_jobs(args.job, names=selected, rerun=args.rerun, dry_run=args.dry_run)
        except (FileNotFoundError, ValueError, json.JSONDecodeError) as exc:
            parser.error(str(exc))
    elif args.check:
        check_jobs(args.job)
    elif args.name:
        job = args.job.lower()
        if job.startswith("q"):
            run_qha(args.name, rerun=args.rerun, run_id=args.run_id)
        elif job.startswith("s"):
            run_static(args.name, rerun=args.rerun, run_id=args.run_id)
    else:
        parser.print_help()
