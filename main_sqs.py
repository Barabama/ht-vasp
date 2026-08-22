#!/usr/bin/env python3
"""Generate and run binary SQS mixing points.

The script has two deliberate stages:

1. ``--generate`` creates SQS structures and a manifest using the strict
   ``sqsgenerator`` adapter.
2. Running the script executes one independent Relax -> Static flow for each
   pure reference and SQS point, then ``--collect`` writes the portable mixing
   JSON consumed downstream.

``phase-A-B`` names are canonicalized by ``SqsSystem``; reverse declarations
are therefore not duplicated.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import shutil
from pathlib import Path
from typing import Any

from htvasp.model import (
    CalculationProvenance,
    SqsGenerator,
    SqsMixingData,
    SqsPlan,
    SqsSystem,
    TemplateRegistry,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

ROOT_DIR = Path(__file__).resolve().parent
FLOW_DIR = Path("/nfs_ssd/tmp")
DATA_DIR = ROOT_DIR / "data" / "sqs_mixing"

VASP_ARGS = {
    "vasp_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-cpu && srun vasp_std'",
    "vasp_gamma_cmd": "/bin/bash -c '. /etc/profile.d/modules.sh && module purge && module load vasp-cpu && srun vasp_gam'",
}
MAGMOM = {
    "Al": 1.0, "Co": 3.0, "Cr": 3.0, "Cu": 5.0, "Fe": 5.0, "Mn": 5.0,
    "Nb": 3.0, "Ni": 2.0, "Ta": 3.0, "Ti": 3.0, "V": 3.0, "W": 3.0, "Zr": 1.0,
}
GLOBAL_INCAR = {
    "ISMEAR": 1, "SIGMA": 0.2, "ALGO": "Normal", "NELM": 120,
    "NELMIN": 6, "NELMDL": -6,
    "ISPIN": 2, "MAGMOM": MAGMOM, "ISYM": 0,
    "LWAVE": False, "LCHARG": False, "LAECHG": False, "LORBIT": 10,
    "GGA": "PE",
}
RELAX_INCAR = {
    "ENCUT": 400,
    "PREC": "Normal",
    "EDIFF": 1e-5,
    "EDIFFG": -0.05,
    "ISIF": 3,
    "NSW": 80,
    "POTIM": 0.2,
    "LWAVE": True,
    "LCHARG": False,
    "LAECHG": False,
    "LORBIT": 10,
}
STATIC_INCAR = {
    "ENCUT": 520,
    "PREC": "Accurate",
    "EDIFF": 1e-6,
    "NSW": 0,
    "IBRION": -1,
    "LORBIT": 10,
    "LCHARG": False,
    "LAECHG": False,
    "LWAVE": False,
}

STRUCTURE_NAMES = (
    # Full production matrix is intentionally disabled for this end-to-end
    # smoke run. Restore the entries below after BCC-Al-Nb completes.
    "BCC-Al-Nb",
    # "BCC-Al-Ti", 
    # "BCC-Al-V", "BCC-Nb-Ti", "BCC-Nb-V", "BCC-Ti-V",
    # "FCC-Co-Fe", "FCC-Co-Mn", "FCC-Co-Ni", "FCC-Fe-Mn", "FCC-Fe-Ni", "FCC-Mn-Ni",
    # "HCP-Co-Ni",
)


def get_system(name: str) -> SqsSystem:
    """Parse and canonicalize a declared SQS system."""

    system = SqsSystem.from_name(name)
    if system.canonical_name not in STRUCTURE_NAMES:
        raise ValueError(
            f"Unknown SQS system {name!r}; choose one of the canonical names in STRUCTURE_NAMES"
        )
    return system


def get_template(system: SqsSystem, registry: TemplateRegistry | None = None):
    return (registry or TemplateRegistry()).get(system.phase)


def make_plan(
    system: SqsSystem,
    *,
    registry: TemplateRegistry | None = None,
    iterations: int = 10000,
    seed: int = 0,
    supercell: tuple[int, int, int] | None = None,
) -> SqsPlan:
    return SqsPlan.default(
        system,
        registry=registry,
        iterations=iterations,
        seed=seed,
        supercell=supercell,
    )


def _system_dir(system: SqsSystem, data_dir: Path = DATA_DIR) -> Path:
    return data_dir / system.canonical_name


def _manifest_path(system: SqsSystem, data_dir: Path = DATA_DIR) -> Path:
    return _system_dir(system, data_dir) / "manifest.json"


def _structure_dir(system: SqsSystem, data_dir: Path = DATA_DIR) -> Path:
    return _system_dir(system, data_dir) / "structures"


def _result_path(system: SqsSystem, data_dir: Path = DATA_DIR) -> Path:
    return _system_dir(system, data_dir) / "static-results.json"


def _mixing_path(system: SqsSystem, data_dir: Path = DATA_DIR) -> Path:
    return _system_dir(system, data_dir) / f"{system.canonical_name}-sqs-mixing.json"


def _settings_fingerprint() -> str:
    payload = json.dumps(
        {
            "global_incar": GLOBAL_INCAR,
            "relax_incar": RELAX_INCAR,
            "static_incar": STATIC_INCAR,
            "relax_reciprocal_density": 64,
            "static_reciprocal_density": 64,
            "potcar_functional": "PBE",
            "workflow": "single_relax_then_static",
        },
        sort_keys=True,
        default=str,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _calculation_provenance() -> CalculationProvenance:
    return CalculationProvenance(
        settings_fingerprint=_settings_fingerprint(),
        workflow={
            "worker": "SqsMixingWorker",
            "flow": "single_relax_then_static",
            "relax_reciprocal_density": 64,
            "static_reciprocal_density": 64,
        },
        relax_incar=RELAX_INCAR,
        static_incar=STATIC_INCAR,
        vasp={
            "command": VASP_ARGS,
            "relax_encut_eV": RELAX_INCAR["ENCUT"],
            "static_encut_eV": STATIC_INCAR["ENCUT"],
            "ispin": GLOBAL_INCAR["ISPIN"],
        },
        potcar={"functional": "PBE"},
    )


def _write_poscar(structure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    structure.to(fmt="poscar", filename=path)


def _read_plan(system: SqsSystem, *, options: argparse.Namespace | None = None) -> SqsPlan:
    data_dir = getattr(options, "data_dir", DATA_DIR)
    path = _manifest_path(system, data_dir)
    if path.is_file():
        return SqsPlan.read_manifest(path)
    if options is None:
        return make_plan(system)
    return make_plan(
        system,
        iterations=getattr(options, "iterations", 10000),
        seed=getattr(options, "seed", 0),
        supercell=tuple(options.supercell) if getattr(options, "supercell", None) else None,
    )


def generate_tick(
    name: str,
    *,
    rerun: bool = False,
    iterations: int = 10000,
    seed: int = 0,
    supercell: tuple[int, int, int] | None = None,
    data_dir: Path = DATA_DIR,
) -> SqsPlan:
    """Generate all SQS structures for one canonical binary system."""

    system = get_system(name)
    system_dir = _system_dir(system, data_dir)
    manifest = _manifest_path(system, data_dir)
    if manifest.exists() and not rerun:
        log.info("%s already generated; use --rerun to regenerate", system.canonical_name)
        return SqsPlan.read_manifest(manifest)
    if rerun and system_dir.exists():
        shutil.rmtree(system_dir)

    plan = make_plan(
        system, iterations=iterations, seed=seed, supercell=supercell
    )
    structures = _structure_dir(system, data_dir)
    for element, reference in plan.reference_structures().items():
        _write_poscar(reference, structures / f"{plan.system.reference_name(element)}.vasp")

    generator = SqsGenerator(plan.engine)
    expanded = plan.expanded_template()
    for index, point in enumerate(plan.points):
        structure, metadata = generator.generate(
            expanded,
            point,
            element_a=system.element_a,
            element_b=system.element_b,
            seed=plan.seed + index,
            iterations=plan.iterations,
        )
        _write_poscar(structure, structures / f"{point.label}.vasp")
        (structures / f"{point.label}.generation.json").write_text(
            json.dumps(metadata, indent=2, ensure_ascii=True, default=str) + "\n",
            encoding="utf-8",
        )
    plan.write_manifest(manifest)
    log.info("Generated %s (%d points, %d sites)", system.canonical_name, len(plan.points), plan.n_sites)
    return plan


def make_worker():
    """Construct the dedicated SQS worker lazily (so model-only tools stay offline)."""

    from htvasp.workflows import SqsMixingWorker

    return SqsMixingWorker(
        vasp_args=VASP_ARGS,
        potcar_functional="PBE",
        global_incar=GLOBAL_INCAR,
        relax_incar=RELAX_INCAR,
        static_incar=STATIC_INCAR,
        relax_reciprocal_density=64,
        static_reciprocal_density=64,
    )


def _run_one(worker, label: str, structure, system_dir: Path, rerun: bool) -> dict[str, Any]:
    store_dir = system_dir / "runs" / label
    if rerun and store_dir.exists():
        shutil.rmtree(store_dir)
    worker.run_flow(
        name=label,
        structure=structure,
        flow_dir=FLOW_DIR,
        store_dir=store_dir,
        resume=not rerun,
    )
    result = worker.get_result()
    if not result:
        raise RuntimeError(f"No static result returned for {label}")
    if isinstance(result, dict) and result.get("state") not in (None, "successful"):
        raise RuntimeError(f"Static task {label} returned state={result.get('state')!r}")
    # Base Worker.get_result returns the structured job output rather than the
    # surrounding task document.  Add the success envelope here so the emitted
    # contract never relies on an inferred user-provided convergence flag.
    return {"state": "successful", "output": result, "run_id": label}


def run_tick(name: str, *, rerun: bool = False, dry_run: bool = False, data_dir: Path = DATA_DIR) -> None:
    """Run references and SQS points, or only build their flows in dry-run mode."""

    system = get_system(name)
    plan = _read_plan(system, options=argparse.Namespace(data_dir=data_dir))
    structures = _structure_dir(system, data_dir)
    if not structures.is_dir():
        raise FileNotFoundError(
            f"No generated structures for {system.canonical_name}; run --generate first"
        )
    reference_structures = plan.reference_structures()
    point_structures = {}
    for point in plan.points:
        path = structures / f"{point.label}.vasp"
        if not path.is_file():
            raise FileNotFoundError(f"Missing generated SQS structure: {path}")
        from pymatgen.core import Structure

        point_structures[point.label] = Structure.from_file(path)

    worker = make_worker()
    if dry_run:
        for label, structure in {
            **{plan.system.reference_name(e): s for e, s in reference_structures.items()},
            **point_structures,
        }.items():
            flow = worker._make_flow(structure)
            log.info("%s: %s jobs for %s", system.canonical_name, len(flow.jobs), label)
        return

    results: dict[str, Any] = {"references": {}, "points": {}}
    for element, structure in reference_structures.items():
        label = plan.system.reference_name(element)
        results["references"][element] = _run_one(worker, label, structure, _system_dir(system, data_dir), rerun)
    for point, structure in point_structures.items():
        results["points"][point] = _run_one(worker, point, structure, _system_dir(system, data_dir), rerun)
    _result_path(system, data_dir).write_text(
        json.dumps(results, indent=2, ensure_ascii=True, default=str) + "\n", encoding="utf-8"
    )
    if system.phase == "BCC":
        collect_tick(name, data_dir=data_dir)
    else:
        log.info(
            "%s static results saved; formal htvasp-sqs-mixing/v1 collection is currently "
            "restricted to BCC_A2",
            system.canonical_name,
        )


def collect_tick(name: str, *, data_dir: Path = DATA_DIR) -> SqsMixingData:
    """Recompute Hmix from persisted static outputs and write the mixing JSON."""

    system = get_system(name)
    if system.phase != "BCC":
        raise ValueError(
            f"{system.canonical_name} has no formal P2 CALPHAD mapping yet; "
            "only BCC_A2 SQS JSON can be collected by this version"
        )
    plan = _read_plan(system, options=argparse.Namespace(data_dir=data_dir))
    path = _result_path(system, data_dir)
    if not path.is_file():
        raise FileNotFoundError(f"No persisted static results: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    from pymatgen.core import Structure

    structures = _structure_dir(system, data_dir)
    point_structures = {
        point.label: Structure.from_file(structures / f"{point.label}.vasp") for point in plan.points
    }
    point_provenance = {}
    for point in plan.points:
        metadata_path = structures / f"{point.label}.generation.json"
        metadata = json.loads(metadata_path.read_text(encoding="utf-8")) if metadata_path.is_file() else {}
        point_provenance[point.label] = {
            key: metadata[key]
            for key in ("engine", "engine_version", "iterations", "seed", "config", "objective", "structure_sha256")
            if key in metadata
        }
    data = SqsMixingData.from_runs(
        plan,
        payload["references"],
        payload["points"],
        point_structures=point_structures,
        reference_structures=plan.reference_structures(),
        point_provenance=point_provenance,
        dataset_id=system.canonical_name,
        calculation=_calculation_provenance(),
    )
    output = data.write_json(_mixing_path(system, data_dir))
    log.info("Wrote %s", output)
    return data


def check_tick(name: str, *, data_dir: Path = DATA_DIR) -> None:
    system = get_system(name)
    manifest = _manifest_path(system, data_dir)
    mixing = _mixing_path(system, data_dir)
    if not manifest.exists():
        log.info("%s: not generated", system.canonical_name)
        return
    plan = SqsPlan.read_manifest(manifest)
    if mixing.exists():
        data = SqsMixingData.read_json(mixing)
        log.info("%s: %s, %d points, Hmix JSON valid", system.canonical_name, data.calphad_mapping.phase_name, len(data.points))
    else:
        log.info("%s: generated (%d sites, %d points), no collected results", system.canonical_name, plan.n_sites, len(plan.points))


def run_batch(*, rerun: bool = False, dry_run: bool = False, data_dir: Path = DATA_DIR) -> None:
    for name in STRUCTURE_NAMES:
        run_tick(name, rerun=rerun, dry_run=dry_run, data_dir=data_dir)


def submit_jobs(*, rerun: bool = False, data_dir: Path = DATA_DIR) -> None:
    from htvasp.slurm import SlurmJobManager

    manager = SlurmJobManager()
    for name in STRUCTURE_NAMES:
        system = get_system(name)
        config = manager.get_cpu_config(
            job_name=f"{system.canonical_name}-sqs",
            output_log=f"logs/{system.canonical_name}-sqs.log",
            ntasks=32,
            memory="64G",
        )
        command = f"python {Path(__file__).name} --name {system.canonical_name}"
        if data_dir != DATA_DIR:
            command += f" --data-dir {data_dir}"
        if rerun:
            command += " --rerun"
        # The project environment is a path-only conda environment on this
        # cluster; there is no registered environment named ``htvasp``.
        job_id = manager.submit_command(
            command=command,
            config=config,
            conda_env=str(ROOT_DIR / ".conda"),
            work_dir=ROOT_DIR,
        )
        if job_id:
            log.info("Submitted %s as %s", system.canonical_name, job_id)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="HT-VASP binary SQS mixing workflow")
    parser.add_argument("--name", "-n", help="Canonical or reverse phase-A-B name")
    parser.add_argument("--batch", "-b", action="store_true", help="Run all canonical systems")
    parser.add_argument("--generate", action="store_true", help="Generate SQS structures")
    parser.add_argument("--collect", action="store_true", help="Collect static outputs into mixing JSON")
    parser.add_argument("--check", action="store_true", help="Validate generated/collected data")
    parser.add_argument("--slurm", "-s", action="store_true", help="Submit all systems to Slurm")
    parser.add_argument("--dry-run", action="store_true", help="Build flows without launching VASP")
    parser.add_argument("--rerun", "-r", action="store_true", help="Regenerate or rerun existing data")
    parser.add_argument("--iterations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--supercell", type=int, nargs=3, metavar=("A", "B", "C"))
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument(
        "--sync-template",
        choices=("BCC", "FCC", "HCP"),
        help="Explicitly fetch one prototype from Materials Project and update its snapshot",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    names = [args.name] if args.name else list(STRUCTURE_NAMES)
    if args.batch:
        names = list(STRUCTURE_NAMES)
    if args.sync_template:
        path = TemplateRegistry().update_from_mp(args.sync_template)
        log.info("Synchronized %s template to %s", args.sync_template, path)
        return 0
    if args.slurm:
        submit_jobs(rerun=args.rerun, data_dir=args.data_dir)
        return 0
    if args.generate:
        for name in names:
            generate_tick(
                name,
                rerun=args.rerun,
                iterations=args.iterations,
                seed=args.seed,
                supercell=tuple(args.supercell) if args.supercell else None,
                data_dir=args.data_dir,
            )
        return 0
    if args.collect:
        for name in names:
            collect_tick(name, data_dir=args.data_dir)
        return 0
    if args.check:
        for name in names:
            check_tick(name, data_dir=args.data_dir)
        return 0
    if args.batch:
        run_batch(rerun=args.rerun, dry_run=args.dry_run, data_dir=args.data_dir)
    elif args.name:
        run_tick(args.name, rerun=args.rerun, dry_run=args.dry_run, data_dir=args.data_dir)
    else:
        _parser().print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
