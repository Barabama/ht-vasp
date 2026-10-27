#!/usr/bin/env python3
"""Prepare/audit separate BCC-Ti-V legacy-DFPT and auto-kpoint pilot inputs.

This tool never submits Slurm jobs and never invokes VASP. The paired shell
script defaults to review mode; its calculation path is separately gated.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
import re
import subprocess
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

NAME = "BCC-Ti-V"
POTCAR_ROOT_DEFAULT = Path("/nfs_ssd/documents")
POTCAR_LIBRARY = "potpaw_PBE"

# VASP writes this once when it exits cleanly. It only proves the process
# finished; it says nothing about whether the SCF converged.
VASP_COMPLETION_MARKER = "General timing and accounting informations for this job"
# VASP writes this when an SCF loop is abandoned because NELM was reached
# before EDIFF. Any occurrence invalidates the run for downstream use.
VASP_UNCONVERGED_MARKER = "aborting loop EDIFF was not reached (unconverged)"
# VASP writes this for every SCF loop that did reach EDIFF.
VASP_CONVERGED_MARKER = "aborting loop because EDIFF is reached"
# DFPT (IBRION=8) prints one "Degree of freedom:   n/  M" line per perturbation
# in OUTCAR. The run only solved all perturbations when the last line is M/M.
DFPT_DEGREE_RE = re.compile(r"Degree of freedom:\s*(\d+)\s*/\s*(\d+)")
# VASP prints this on stdout (vasp.log), not in OUTCAR. Audit only.
VASP_LINEAR_RESPONSE_MARKER = "Linear response finished"

STAGES = ("initial-relax", "volume-static", "charge-seed", "dfpt")

PROFILE_MESHES = {
    "legacy-fixed-grid": {
        "initial-relax": (10, 10, 10),
        "volume-static": (10, 10, 10),
        "charge-seed": (5, 5, 5),
        "dfpt": (5, 5, 5),
    },
    "atomate2-auto-kpoints": {
        "initial-relax": None,
        "volume-static": None,
        "charge-seed": None,
        "dfpt": None,
    },
}

# These are the calculation tags from B5070's INCAR-R31, INCAR-R3f, and
# INCAR-Phonopy templates, represented as pymatgen Incar mappings. INCAR-R3f
# is intentionally static (IBRION=-1, NSW=0), not an ionic relaxation.
R3_INCAR: dict[str, Any] = {
    "SYSTEM": "BCC-Ti-V R3",
    "ENCUT": 520,
    "PREC": "Accurate",
    "EDIFF": 1e-6,
    "NELMIN": 4,
    "NELM": 300,
    # B5070 used FAST, but BCC-Ti-V starts far from its relaxed volume.
    # NORMAL avoids the observed ZBRENT bracketing failure during R3.
    "ALGO": "Normal",
    "LREAL": "Auto",
    "LWAVE": False,
    "LCHARG": False,
    "NCORE": 2,
    "KPAR": 4,
    "IBRION": 2,
    "POTIM": 0.2,
    "NSW": 100,
    "ISIF": 3,
    "EDIFFG": -0.01,
    "ISMEAR": 1,
    "SIGMA": 0.2,
    "ISTART": 0,
    "ICHARG": 2,
    "ISPIN": 2,
    "LORBIT": 11,
}

R3F_STATIC_INCAR: dict[str, Any] = {
    "SYSTEM": "BCC-Ti-V R3f static volume point",
    "ENCUT": 520,
    "PREC": "Accurate",
    "EDIFF": 1e-6,
    "NELM": 120,
    "NELMIN": 4,
    "ALGO": "Fast",
    "LREAL": False,
    "NSIM": 4,
    "NCORE": 2,
    "KPAR": 4,
    "IBRION": -1,
    "NSW": 0,
    "ISIF": 3,
    "ISMEAR": 1,
    "ISTART": 0,
    "ICHARG": 2,
    "ISPIN": 2,
    "LWAVE": False,
    "LCHARG": False,
}

# Independent SCF that produces the CHGCAR the DFPT step reads.
#
# Why this stage exists: B5070's INCAR-Phonopy uses ICHARG=1 ("read the initial
# charge density"), but the legacy driver never produced a CHGCAR for the DFPT
# step to read (R3 and R3f both run LCHARG=False). VASP therefore fell back to
# the atomic superposition, and because the template sets no MAGMOM the guess
# is all-spin-up: 200 electrons -> exactly 16 mu_B. That guess does not relax
# under ALGO=Fast/IMIX=4; job 53681 exhausted NELM=120 and job 53688 drifted
# from -132.75 eV to -125.05 eV over 244 steps with rms(c) climbing 0.004->0.29.
#
# Raising NELM does not help, so this stage supplies the charge density the
# template always assumed was there, and the DFPT INCAR is left byte-for-byte
# equal to INCAR-Phonopy. The recipe below is the one already validated on this
# exact 16-site BCC-Ti-V cell by data/diagnostics/BCC-Ti-V-phonon-20260921-v3,
# where all 36 non-magnetic static cases converged (28-58 SCF steps).
CHARGE_SEED_INCAR: dict[str, Any] = {
    "SYSTEM": "BCC-Ti-V DFPT charge seed",
    "ENCUT": 520,
    "PREC": "Accurate",
    "EDIFF": 1e-6,
    "ISMEAR": 1,
    "SIGMA": 0.2,
    "ISPIN": 2,
    # Supplying MAGMOM is the decisive difference from the DFPT template. The
    # diagnostic found the cell non-magnetic (|m| -> 0.008 mu_B), so this is a
    # starting guess, not an asserted moment. It must be a per-site list: a
    # "16*3.0" string is not understood by pymatgen's Incar.
    "MAGMOM": [3.0] * 16,
    "ALGO": "Normal",
    "ISYM": 0,
    "LREAL": False,
    "NCORE": 2,
    "KPAR": 4,
    "IBRION": -1,
    "NSW": 0,
    "ISIF": 2,
    "ISTART": 0,
    "ICHARG": 2,
    "NELM": 200,
    "NELMIN": 6,
    "NELMDL": -6,
    # ENAUG, LASPH, LMAXMIX and LORBIT are deliberately NOT set here.
    #
    # The CHGCAR this stage writes must be the ground-state density of exactly
    # the Hamiltonian DFPT will use. VASP stores the density on a grid sized by
    # ENCUT and, when set, ENAUG; LASPH changes the PAW on-site terms; LMAXMIX
    # controls how much of the one-center PAW density is written to the file.
    # Setting them here while DFPT inherits VASP's defaults is a real, measured
    # inconsistency: with them set, the seed and DFPT ground states disagreed by
    # 5.01 meV; once dropped they agree to 0.29 meV.
    #
    # It is NOT, however, the reason DFPT diverges. Jobs 53690 (seed with the
    # tags) and 53692 (seed without them, Hamiltonian matched to 0.29 meV) both
    # had their DFPT ground-state SCF abandoned at NELM while the ground state
    # plateaued at -132.7514 eV and then drifted to -125 eV. The divergence is
    # driven by ALGO=Fast (RMM-DIIS, no charge mixer - see DFPT_INCAR).
    "LWAVE": False,
    # The whole point of the stage: hand a converged charge density to DFPT.
    "LCHARG": True,
}

DFPT_INCAR: dict[str, Any] = {
    "SYSTEM": "BCC-Ti-V DFPT",
    "ENCUT": 520,
    "PREC": "A",
    "EDIFF": 1e-6,
    "EDIFFG": -0.01,
    "NELM": 120,
    "NELMIN": 10,
    # STILL B5070's ALGO=Fast, and this is the remaining open scientific
    # question. Evidence that it is the actual failure driver:
    #
    #   * With ALGO=Fast, the ground-state SCF plateaus at -132.7514 eV and then
    #     drifts monotonically to about -125 eV. Three independent runs (53681
    #     unseeded, 53690 seeded with mismatched PAW tags, 53692 seeded with
    #     tags matched to 0.29 meV) trace the same curve, so the drift does not
    #     depend on where the SCF starts.
    #   * The plateaus agree between seeded and unseeded runs (-132.751449 vs
    #     -132.751443 eV), i.e. the seed reaches the same state - it just cannot
    #     hold it.
    #   * rms(c) turns over at step ~22-24 (~0.7e-2) and grows monotonically
    #     from there; the run is at its best approximation to the ground state
    #     right where the drift starts.
    #   * With ALGO=Normal on the same 16-site cell, the same Hamiltonian and
    #     the same 5x5x5 mesh, the charge-seed stage converges in 26 steps.
    #
    # Mechanism: ALGO=Fast uses blocked Davidson only until the RMM-DIIS
    # criterion is met (OUTCAR reports IALGO=68), then switches to RMM-DIIS,
    # which is a direct-optimization scheme that updates the density from the
    # orbitals WITHOUT a charge mixer (see the VASP wiki ALGO page). On this
    # metal/gamma-centered 5x5x5 mesh that becomes an undamped fixed-point
    # iteration. Every mixer-class remedy (AMIX/BMIX/AMIN, Kerker) is bypassed
    # in that mode, and the wiki's own wording for Fast is that it may
    # "diverge or cannot find the true minimum" and should be used "to obtain a
    # first estimate". Its orbitals, "and hence forces, stress, etc., may be
    # less accurate compared to Normal" - which matters directly here, because
    # DFPT force constants come from perturbed orbitals.
    #
    # Changing it is a scientific input change to a protected template, so the
    # legacy value is kept as the default and the alternatives live in the
    # DFPT_VARIANTS table below, selectable per run.
    "ALGO": "Fast",
    "LREAL": False,
    # VASP 6.4.3 DFPT (IBRION=8) on 429pro aborts when NCORE is set.
    # Omit NCORE as a platform compatibility override; preserve it in provenance.
    "KPAR": 2,
    "IBRION": 8,
    "NSW": 1,
    "ISMEAR": 1,
    "SIGMA": 0.2,
    "ISTART": 0,
    # B5070 value, unchanged. It is only meaningful now that the shell driver
    # writes a converged CHGCAR from the charge-seed stage into this directory.
    "ICHARG": 1,
    "ISPIN": 2,
}

# PAW/one-center tags used by every converged BCC-Ti-V calculation already in the
# repository (data/diagnostics/BCC-Ti-V-phonon-20260921-v3,
# BCC-Ti-V-softmode-*, scripts/prepare_nbti_*, htvasp/workflows/tb2j.py). They
# are NOT in B5070's INCAR-Phonopy, so they are opt-in here.
#
# ENAUG is deliberately NOT in this set. VASP ignores it unless PREC is one of
# the legacy Low/Medium/High settings, and the VASP wiki marks it deprecated.
# Measured on this cell with PREC=Accurate, ENAUG=1360 vs the POTCAR default
# 526.6 eV changed nothing: identical FFT grids (NGXF=96), identical
# "Maximum index for augmentation-charges" (1186), identical CHGCAR byte count.
# Carrying it would put a tag in the manifest that looks like a physics choice
# but is not one.
PAW_TIGHT = {"LASPH": True, "LMAXMIX": 4}

# DFPT electronic-structure variants. Every variant keeps the legacy template's
# ICHARG=1 / NELM=120 / NELMIN=10 / IBRION=8 / KPAR=2 / 5x5x5 mesh; only the
# tags under test move. The charge-seed INCAR is adjusted to match, so
# check-stage-compat stays a meaningful gate in every variant.
#
# "production"  : ALGO=Normal + LASPH=True + LMAXMIX=4. The agreed production
#                 recipe. ALGO=Normal is required: with the legacy ALGO=Fast the
#                 ground-state SCF is abandoned at NELM (jobs 53681/53690/53692),
#                 while Normal converges in ~51 steps (job 53696). LASPH/LMAXMIX
#                 match every other converged Ti-V calculation in this repo.
# "algo-normal" : E1 - isolates ALGO=Fast vs Normal, no PAW tags. Converged
#                 (53696) and was used to prove ALGO is the divergence driver.
# "legacy"      : B5070 INCAR-Phonopy verbatim. Kept for reproducibility of the
#                 original failure; known to diverge (53681, 53690, 53692).
DFPT_VARIANTS: dict[str, dict[str, Any]] = {
    "production": {
        "dfpt": {"ALGO": "Normal", **PAW_TIGHT},
        "charge_seed": dict(PAW_TIGHT),
    },
    "algo-normal": {
        "dfpt": {"ALGO": "Normal"},
        "charge_seed": {},
    },
    "legacy": {
        "dfpt": {},
        "charge_seed": {},
    },
}
DEFAULT_DFPT_VARIANT = "production"

# Tags VASP provably ignores under the settings this workflow uses. They are
# excluded from the recipe and Hamiltonian comparisons, because comparing them
# would flag two physically identical runs as different.
#
# ENAUG: the VASP wiki states it "has an effect only if PREC is set to one of
# the old settings (Low, Medium or High), otherwise it is ignored", and marks
# the tag deprecated. Measured on this 16-site cell at PREC=Accurate with
# ENAUG=1360 vs the POTCAR default 526.6 eV: identical FFT grids (NGXF=96),
# identical maximum augmentation-charge index (1186), identical CHGCAR byte
# count (32273653). This matters in practice: the earlier E3 runs carry
# ENAUG=1360 while the production recipe does not, and treating that as a
# mismatch would discard verified results for no physical reason.
INERT_TAGS = frozenset({"ENAUG"})


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def versions() -> dict[str, str]:
    names = ("pymatgen", "phonopy", "atomate2", "jobflow")
    found: dict[str, str] = {}
    for name in names:
        try:
            found[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            found[name] = "unavailable"
    return found


def git_revision(root: Path) -> str | None:
    try:
        result = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def configure_potcar(root: Path) -> Path:
    # Pymatgen expects PMG_VASP_PSP_DIR to be the parent of the functional
    # folder, matching /home/mcmf429/.pmgrc.yaml on 429pro.
    if root.name == POTCAR_LIBRARY and (root / "Ti_sv" / "POTCAR").is_file():
        root = root.parent
    if not root.is_dir():
        raise FileNotFoundError(f"PMG_VASP_PSP_DIR not found: {root}")
    library = root / POTCAR_LIBRARY
    for symbol in ("Ti_sv", "V_sv"):
        if not (library / symbol / "POTCAR").is_file():
            raise FileNotFoundError(f"missing PBE POTCAR: {library / symbol / 'POTCAR'}")
    os.environ["PMG_VASP_PSP_DIR"] = str(root)
    return root


def load_model(root: Path, potcar_root: Path):
    potcar_root = configure_potcar(potcar_root)
    model_path = root / "data" / "poscars" / f"{NAME}.vasp"
    if not model_path.is_file():
        raise FileNotFoundError(f"ht-vasp model POSCAR missing: {model_path}")
    source_hash = sha256(model_path)

    # Reuse the local ht-vasp Endmember/Pymatgen construction method. The
    # required local POSCAR prevents a network lookup or replacement.
    from htvasp.model import Endmember

    structure = Endmember().get_poscar(NAME, root / "data" / "poscars")
    if sha256(model_path) != source_hash:
        raise RuntimeError(f"model POSCAR changed while loading: {model_path}")
    species = [str(site.specie) for site in structure]
    if len(structure) != 2 or species != ["Ti", "V"]:
        raise ValueError(f"unexpected BCC-Ti-V model: {len(structure)} sites {species}")
    return structure, model_path, source_hash, potcar_root


def get_generator(stage: str, incar: dict[str, Any]):
    from atomate2.vasp.sets.core import RelaxSetGenerator, StaticSetGenerator

    cls = RelaxSetGenerator if stage == "initial-relax" else StaticSetGenerator
    # The generator's own MAGMOM handling only accepts a {species: moment} dict,
    # but our templates express per-site moments (and the written INCAR is
    # replaced wholesale from `incar` below anyway). Drop it here so the
    # generator never sees a shape it cannot parse.
    generator_settings = {k: v for k, v in incar.items() if k != "MAGMOM"}
    return cls(user_potcar_functional="PBE", user_incar_settings=generator_settings)


def generator_kpoints(generator: Any, structure: Any):
    # get_input_set is intentionally called on the actual calculation cell:
    # primitive unit cell for R3/R3f, phonopy's 16-site SPOSCAR for DFPT.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Overriding the POTCAR functional.*")
        return generator.get_input_set(structure).kpoints


def mesh_product(mesh: tuple[int, int, int]) -> int:
    return math.prod(mesh)


def kpoint_payload(kpoints: Any, natoms: int) -> dict[str, Any]:
    mesh = tuple(int(value) for value in kpoints.kpts[0])
    return {
        "style": str(kpoints.style),
        "grid": list(mesh),
        "grid_points_full_mesh": mesh_product(mesh),
        "atom_count": natoms,
        "nominal_kppra_full_mesh": mesh_product(mesh) * natoms,
        "comment": kpoints.comment,
    }


def selected_mesh(profile: str, stage: str) -> tuple[int, int, int] | None:
    return PROFILE_MESHES[profile][stage]


def build_vasp_input(
    *,
    structure: Any,
    stage: str,
    profile: str,
    incar: dict[str, Any],
):
    from pymatgen.io.vasp.inputs import Incar, Kpoints

    generator = get_generator(stage, incar)
    vasp_input = generator.get_input_set(structure)
    # Keep the legacy INCAR template tags, while still relying on atomate2 /
    # pymatgen to construct the structure, POTCAR, and (auto profile) KPOINTS.
    vasp_input["INCAR"] = Incar(incar)
    mesh = selected_mesh(profile, stage)
    if mesh is not None:
        vasp_input["KPOINTS"] = Kpoints.gamma_automatic(kpts=mesh)
    return vasp_input, generator, kpoint_payload(vasp_input["KPOINTS"], len(structure))


def write_inputs(
    *,
    root: Path,
    potcar_root: Path,
    outdir: Path,
    structure: Any,
    source_structure: Path,
    stage: str,
    profile: str,
    incar: dict[str, Any],
    scale: float | None = None,
    dfpt_variant: str | None = None,
) -> dict[str, Any]:
    potcar_root = configure_potcar(potcar_root)
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"refusing to overwrite non-empty output directory: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    vasp_input, generator, kp = build_vasp_input(
        structure=structure, stage=stage, profile=profile, incar=incar
    )
    vasp_input.write_input(outdir)

    # Replace the generated POSCAR with the exact source cell for the stages
    # that run on phonopy's supercell, so the site order matches SPOSCAR (which
    # is what fixes the POTCAR ordering and the atom indexing in
    # FORCE_CONSTANTS).
    if stage in ("charge-seed", "dfpt"):
        (outdir / "POSCAR").write_bytes(source_structure.read_bytes())

    file_records = {
        path.name: {"sha256": sha256(path), "bytes": path.stat().st_size}
        for path in sorted(outdir.iterdir())
        if path.is_file()
    }
    manifest = {
        "recipe": profile,
        "recipe_family": "legacy-cxy-hmd-dfpt-qha-v1",
        "name": NAME,
        "stage": stage,
        # Recorded for every stage so a run directory is self-describing even
        # outside the run root (mirrors the run-manifest env echo).
        "dfpt_variant": dfpt_variant or DEFAULT_DFPT_VARIANT,
        "run_id": os.environ.get("RUN_ID"),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "incar_variant_overrides": (
            DFPT_VARIANTS[dfpt_variant or DEFAULT_DFPT_VARIANT]["charge_seed"]
            if stage == "charge-seed"
            else DFPT_VARIANTS[dfpt_variant or DEFAULT_DFPT_VARIANT]["dfpt"]
            if stage == "dfpt"
            else {}
        ),
        "scale": scale,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_structure": str(source_structure.resolve()),
        "source_structure_sha256": sha256(source_structure),
        "model_poscar": str((root / "data" / "poscars" / f"{NAME}.vasp").resolve()),
        "model_poscar_sha256": sha256(root / "data" / "poscars" / f"{NAME}.vasp"),
        "structure_formula": structure.composition.reduced_formula,
        "cell_volume_ang3": structure.volume,
        "site_count": len(structure),
        "site_species": [str(site.specie) for site in structure],
        "pymatgen_generator": f"{type(generator).__module__}.{type(generator).__name__}",
        "potcar_functional": "PBE",
        "PMG_VASP_PSP_DIR": str(potcar_root),
        "potcar_library": str(potcar_root / POTCAR_LIBRARY),
        "potcar_symbols": list(vasp_input["POTCAR"].symbols),
        "kpoints": kp,
        "incar": incar,
        "platform_compatibility_overrides": (
            [{"tag": "NCORE", "legacy_value": 4, "effective_value": "omitted",
              "reason": "429pro VASP 6.4.3 IBRION=8 rejects NCORE != 1"}]
            if stage == "dfpt" else []
        ),
        "software": versions(),
        "htvasp_git_revision": git_revision(root),
        "files": file_records,
    }
    manifest_path = outdir / "input-manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )
    return manifest


def stage_incar(stage: str, dfpt_variant: str | None = None) -> dict[str, Any]:
    """Return the INCAR for a stage, with an optional DFPT experiment variant.

    The variant adjusts DFPT and, where the change is Hamiltonian-affecting, the
    charge-seed stage that feeds it. ``legacy`` reproduces B5070's templates.
    """
    base = {
        "initial-relax": R3_INCAR,
        "volume-static": R3F_STATIC_INCAR,
        "charge-seed": CHARGE_SEED_INCAR,
        "dfpt": DFPT_INCAR,
    }[stage]
    variant = dfpt_variant or DEFAULT_DFPT_VARIANT
    if variant not in DFPT_VARIANTS:
        raise ValueError(
            f"unknown DFPT variant: {variant} (known: {', '.join(DFPT_VARIANTS)})"
        )
    key = "charge_seed" if stage == "charge-seed" else "dfpt"
    overrides = DFPT_VARIANTS[variant][key] if stage in ("charge-seed", "dfpt") else {}
    return {**base, **overrides}


def cmd_check_model(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    structure, source, source_hash, potcar_root = load_model(
        root, Path(args.potcar_root).resolve()
    )
    # Reuse pymatgen's PBE symbol resolver to fail before any VASP work.
    from pymatgen.io.vasp.inputs import Potcar

    potcar = Potcar(["Ti_sv", "V_sv"], functional="PBE")
    print(
        json.dumps(
            {
                "name": NAME,
                "model_path": str(source),
                "model_sha256": source_hash,
                "composition": structure.composition.reduced_formula,
                "site_species": [str(site.specie) for site in structure],
                "volume_ang3": structure.volume,
                "PMG_VASP_PSP_DIR": str(potcar_root),
                "potcar_library": str(potcar_root / POTCAR_LIBRARY),
                "potcar_symbols": potcar.symbols,
                "software": versions(),
                "htvasp_git_revision": git_revision(root),
                "slurm": {
                    "partition": "partCPU",
                    "nodes": 1,
                    "ntasks": 32,
                    "cpus_per_task": 1,
                    "memory": "64G",
                    "time_limit": "100:00:00",
                    "vasp_module": "vasp-cpu",
                },
            },
            indent=2,
            sort_keys=True,
        )
    )


def cmd_audit_kpoints(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    structure, _, _, potcar_root = load_model(root, Path(args.potcar_root).resolve())
    # Must mirror SCALES in legacy_dfpt_qha_pilot_bcc_ti_v.sh: nine points,
    # uniform in volume, V/V0 = 0.95 .. 1.13 (see the comment there).
    scales = (
        0.9830476, 0.9907480, 0.9983305, 1.0057996, 1.0131594,
        1.0204138, 1.0275664, 1.0346209, 1.0415804,
    )
    rows = []
    for scale in scales:
        cell = structure.copy()
        cell.scale_lattice(structure.volume * scale**3)
        vasp_in, _, kpr = build_vasp_input(
            structure=cell,
            stage="volume-static",
            profile="atomate2-auto-kpoints",
            incar=R3F_STATIC_INCAR,
        )
        dfpt_cell = cell.copy()
        dfpt_cell.make_supercell([2, 2, 2])
        dfpt_input, _, kpd = build_vasp_input(
            structure=dfpt_cell,
            stage="dfpt",
            profile="atomate2-auto-kpoints",
            incar=DFPT_INCAR,
        )
        legacy_dfpt = kpoint_payload(
            __import__("pymatgen.io.vasp.inputs", fromlist=["Kpoints"]).Kpoints.gamma_automatic(
                kpts=(5, 5, 5)
            ),
            len(dfpt_cell),
        )
        rows.append(
            {
                "scale": scale,
                "legacy_R3_R3f": {"grid": [10, 10, 10], "nominal_kppra_full_mesh": 2000},
                "auto_R3_R3f": kpr,
                "legacy_DFPT_supercell": legacy_dfpt,
                "auto_DFPT_actual_16_site_supercell": kpd,
                "model_sites": len(cell),
                "dfpt_sites": len(dfpt_cell),
            }
        )
    print("scale  legacy R3/R3f  legacy DFPT  auto R3/R3f (KPPRA)  auto DFPT 16-site (KPPRA)")
    for row in rows:
        auto_r3 = row["auto_R3_R3f"]
        auto_dfpt = row["auto_DFPT_actual_16_site_supercell"]
        print(
            f"{row['scale']:.2f}   10x10x10 (2000)   5x5x5 (2000)   "
            f"{'x'.join(map(str, auto_r3['grid']))} ({auto_r3['nominal_kppra_full_mesh']})   "
            f"{'x'.join(map(str, auto_dfpt['grid']))} ({auto_dfpt['nominal_kppra_full_mesh']})"
        )
    print("KPPRA = full uniform-grid point count × atom count; no symmetry reduction.")


def cmd_prepare(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    potcar_root = Path(args.potcar_root).resolve()
    profile = args.profile
    if profile not in PROFILE_MESHES:
        raise ValueError(f"unknown profile: {profile}")

    if args.stage == "initial-relax":
        structure, source, _, potcar_root = load_model(root, potcar_root)
        write_inputs(
            root=root,
            potcar_root=potcar_root,
            outdir=Path(args.outdir).resolve(),
            structure=structure,
            source_structure=source,
            stage=args.stage,
            profile=profile,
            incar=stage_incar(args.stage, args.dfpt_variant),
            dfpt_variant=args.dfpt_variant,
        )
        return

    if not args.structure:
        raise ValueError(f"--structure is required for {args.stage}")
    source = Path(args.structure).resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    from pymatgen.core import Structure

    structure = Structure.from_file(source)
    scale = None
    if args.stage == "volume-static":
        scale = float(args.scale)
        if not (0.90 <= scale <= 1.10):
            raise ValueError(f"unexpected volume scale: {scale}")
        structure.scale_lattice(structure.volume * scale**3)
    elif args.stage in ("charge-seed", "dfpt"):
        if len(structure) != 16:
            raise ValueError(
                f"{args.stage} input must be phonopy 2x2x2 16-site SPOSCAR; got {len(structure)}"
            )
        composition = structure.composition.get_el_amt_dict()
        if composition != {"Ti": 8.0, "V": 8.0}:
            raise ValueError(f"unexpected {args.stage} supercell composition: {composition}")
    write_inputs(
        root=root,
        potcar_root=potcar_root,
        outdir=Path(args.outdir).resolve(),
        structure=structure,
        source_structure=source,
        stage=args.stage,
        profile=profile,
        incar=stage_incar(args.stage, args.dfpt_variant),
        scale=scale,
        dfpt_variant=args.dfpt_variant,
    )


class VerificationError(RuntimeError):
    """Raised when a VASP run fails the fail-closed convergence check."""


def cmd_verify(args: argparse.Namespace) -> None:
    """Fail-closed scientific validation of one finished VASP run.

    A VASP job that exits cleanly (``General timing ...``) and writes a
    ``vasprun.xml`` is *not* automatically a valid result. This check keeps the
    two notions separate:

    * ``outcar_job_completion_marker`` - the process ran to completion.
    * ``outcar_no_unconverged_scf`` / ``pymatgen_electronic_convergence`` -
      the electronic structure actually converged.

    Any failed check makes the command exit non-zero so the caller (and
    ``set -e`` in the paired shell script) stops before a ``vasprun.xml`` is
    fed to ``phonopy --fc``. A machine-readable audit report is written to
    ``--report`` even when the check fails.
    """
    stage = args.stage
    vasprun_path = Path(args.vasprun).resolve()
    outcar_path = Path(args.outcar).resolve()
    report_path = Path(args.report).resolve() if args.report else None

    checks: dict[str, Any] = {}
    failures: list[str] = []
    notes: list[str] = []
    report: dict[str, Any] = {
        "stage": stage,
        "vasprun": str(vasprun_path),
        "outcar": str(outcar_path),
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "failures": failures,
        "notes": notes,
    }

    def record(name: str, ok: bool | None, *, detail: Any = None, required: bool = True) -> None:
        entry: dict[str, Any] = {"ok": ok, "required": required}
        if detail is not None:
            entry["detail"] = detail
        checks[name] = entry
        if required and ok is not True:
            failures.append(name)

    def emit() -> None:
        report["status"] = "verified" if not failures else "failed"
        text = json.dumps(report, indent=2, sort_keys=True, default=str)
        if report_path is not None:
            report_path.parent.mkdir(parents=True, exist_ok=True)
            report_path.write_text(text + "\n", encoding="utf-8")
        print(text)

    # ---- file level ------------------------------------------------------
    vasprun_ok = vasprun_path.is_file() and vasprun_path.stat().st_size > 0
    outcar_ok = outcar_path.is_file() and outcar_path.stat().st_size > 0
    record("vasprun_present", vasprun_ok, detail=str(vasprun_path))
    record("outcar_present", outcar_ok, detail=str(outcar_path))
    if not (vasprun_ok and outcar_ok):
        emit()
        raise VerificationError(
            f"missing/empty VASP output (vasprun={vasprun_ok}, outcar={outcar_ok})"
        )

    outcar_text = outcar_path.read_text(errors="replace")
    report["outcar_marker_counts"] = {
        "job_completion": outcar_text.count(VASP_COMPLETION_MARKER),
        "scf_converged": outcar_text.count(VASP_CONVERGED_MARKER),
        "scf_unconverged": outcar_text.count(VASP_UNCONVERGED_MARKER),
    }
    record(
        "outcar_job_completion_marker",
        VASP_COMPLETION_MARKER in outcar_text,
        detail="VASP process reached its final accounting block",
    )
    record(
        "outcar_no_unconverged_scf",
        VASP_UNCONVERGED_MARKER not in outcar_text,
        detail="OUTCAR must not contain 'EDIFF was not reached (unconverged)'",
    )
    record(
        "outcar_has_converged_scf",
        VASP_CONVERGED_MARKER in outcar_text,
        detail="at least one SCF loop reached EDIFF",
    )

    if stage == "dfpt":
        # DFPT solves one perturbation per "Degree of freedom:  n/  M" block.
        # The run is only usable when the last block is M/M.
        degrees = DFPT_DEGREE_RE.findall(outcar_text)
        solved = [(int(n), int(m)) for n, m in degrees]
        report["dfpt_degrees_of_freedom"] = [{"index": n, "total": m} for n, m in solved]
        completed = bool(solved) and solved[-1][0] == solved[-1][1]
        record(
            "dfpt_all_perturbations_solved",
            completed,
            detail=(
                f"last block {solved[-1][0]}/{solved[-1][1]}"
                if solved
                else "no 'Degree of freedom' block found in OUTCAR"
            ),
        )
        # ICHARG=1 is only meaningful when a converged starting charge density
        # is actually present. Without one VASP silently falls back to the
        # atomic superposition (the 53681 failure), so make that an explicit,
        # fail-closed precondition instead of a log warning.
        chgcar = outcar_path.parent / "CHGCAR"
        chgcar_ok = chgcar.is_file() and chgcar.stat().st_size > 0
        if chgcar_ok:
            chgcar_head = chgcar.read_text(errors="replace")[:8192].splitlines()
            report["chgcar_first_line"] = chgcar_head[0].strip() if chgcar_head else ""
        record(
            "chgcar_present_for_icharg1",
            chgcar_ok,
            detail=f"{chgcar} ({chgcar.stat().st_size if chgcar.is_file() else 0} bytes)",
        )
        atomic_guess = "charge density of overlapping atoms calculated" in outcar_text
        record(
            "dfpt_used_supplied_charge_density",
            not atomic_guess,
            detail=(
                "VASP reported 'charge density of overlapping atoms calculated', "
                "i.e. ICHARG=1 found no usable CHGCAR and fell back to the "
                "atomic superposition"
                if atomic_guess
                else "OUTCAR shows the supplied charge density was read"
            ),
        )
        # VASP prints the linear-response completion line on stdout only, so it
        # lands in vasp.log rather than OUTCAR. Require it whenever that file
        # was captured; otherwise record it as informational.
        vasp_log = outcar_path.parent / "vasp.log"
        if vasp_log.is_file():
            log_text = vasp_log.read_text(errors="replace")
            record(
                "vasp_log_linear_response_finished",
                VASP_LINEAR_RESPONSE_MARKER in log_text,
                detail=str(vasp_log),
            )
            record(
                "vasp_log_no_incomplete_chargedensity_warning",
                "chargedensity file is incomplete" not in log_text,
                detail=str(vasp_log),
            )
        else:
            record(
                "vasp_log_linear_response_finished",
                None,
                detail=f"no stdout capture at {vasp_log}; not evaluated",
                required=False,
            )
        notes.append(
            "DFPT (IBRION=8) concatenates the ground-state and perturbation SCF "
            "iterations into a single pymatgen ionic step, so "
            "pymatgen.converged_electronic is not a reliable gate here; it is "
            "reported for audit only and the OUTCAR SCF markers above are "
            "authoritative."
        )

    # ---- pymatgen level --------------------------------------------------
    from pymatgen.io.vasp.outputs import Vasprun

    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            parsed = Vasprun(
                str(vasprun_path),
                parse_dos=False,
                parse_eigen=False,
                parse_projected_eigen=False,
            )
        parsed_ok = True
        parse_warnings = [f"{w.category.__name__}: {w.message}" for w in caught]
    except Exception as exc:  # noqa: BLE001 - any parse failure is fatal here
        parsed, parsed_ok, parse_warnings = None, False, [f"{type(exc).__name__}: {exc}"]
    record("vasprun_parses", parsed_ok, detail=parse_warnings or "clean parse")
    if parsed is None:
        emit()
        raise VerificationError(f"vasprun.xml could not be parsed: {vasprun_path}")

    has_structure = parsed.final_structure is not None and len(parsed.final_structure) > 0
    record("final_structure_present", has_structure)
    energy = float(parsed.final_energy) if has_structure else float("nan")
    record("final_energy_finite", math.isfinite(energy), detail=energy)
    if not has_structure:
        emit()
        raise VerificationError(f"vasprun.xml lacks a final structure: {vasprun_path}")

    report["final_energy_eV"] = energy
    report["final_volume_ang3"] = float(parsed.final_structure.volume)
    report["nelm"] = parsed.parameters.get("NELM")
    report["ionic_step_count"] = len(parsed.ionic_steps)
    report["electronic_steps_in_final_ionic_step"] = (
        len(parsed.ionic_steps[-1].get("electronic_steps", [])) if parsed.ionic_steps else 0
    )

    # pymatgen's heuristic is authoritative for relax/static stages only.
    electronic_required = stage != "dfpt"
    record(
        "pymatgen_electronic_convergence",
        bool(parsed.converged_electronic),
        detail="len(final ionic-step electronic steps) < NELM",
        required=electronic_required,
    )
    record(
        "pymatgen_ionic_convergence",
        bool(parsed.converged_ionic),
        detail="ionic loop stopped before NSW",
    )

    forces = parsed.ionic_steps[-1].get("forces") if parsed.ionic_steps else None
    max_force = None
    if forces is not None and len(forces) > 0:
        max_force = max(math.sqrt(sum(float(x) ** 2 for x in row)) for row in forces)
    report["max_force_eV_per_A"] = max_force

    if stage == "charge-seed":
        # The seed stage exists only to hand a converged charge density to
        # DFPT, so a missing CHGCAR means it accomplished nothing.
        chgcar = outcar_path.parent / "CHGCAR"
        record(
            "seed_wrote_chgcar",
            chgcar.is_file() and chgcar.stat().st_size > 0,
            detail=f"{chgcar} ({chgcar.stat().st_size if chgcar.is_file() else 0} bytes)",
        )

    if args.max_force_tolerance is not None and max_force is not None:
        if max_force > args.max_force_tolerance:
            message = (
                f"max force {max_force:.6g} eV/Angstrom exceeds "
                f"{args.max_force_tolerance:.6g}; stop before DFPT and review "
                "fixed-volume ionic relaxation sensitivity"
            )
            if args.fail_on_force:
                record("max_force_within_tolerance", False, detail=message)
            else:
                record("max_force_within_tolerance", True, detail=message)
                notes.append(f"WARNING: {message}")

    emit()
    if failures:
        raise VerificationError(
            "fail-closed VASP verification failed for "
            f"{vasprun_path}: {', '.join(failures)}"
        )


def cmd_check_recipe(args: argparse.Namespace) -> None:
    """Fail when a finished stage was built with a superseded INCAR recipe.

    Reuse decisions keyed only on "outputs verify" would happily reuse a
    charge-seed CHGCAR produced by an older recipe. Job 53690 showed how bad
    that is: the seed carried ENAUG=1360/LASPH=True/LMAXMIX=4 while DFPT ran at
    VASP's defaults, which is exactly the mismatch that made DFPT's ground-state
    SCF diverge. Comparing the recorded manifest against the current template
    makes that situation visible instead of silent.
    """
    manifest_path = Path(args.manifest).resolve()
    recorded = json.loads(manifest_path.read_text(encoding="utf-8")).get("incar", {})
    current = stage_incar(args.stage, getattr(args, "dfpt_variant", None))
    changed = sorted(
        key
        for key in set(recorded) | set(current)
        if recorded.get(key) != current.get(key) and key not in INERT_TAGS
    )
    ignored = sorted(
        key
        for key in set(recorded) | set(current)
        if recorded.get(key) != current.get(key) and key in INERT_TAGS
    )
    report = {
        "manifest": str(manifest_path),
        "stage": args.stage,
        "status": "current" if not changed else "superseded",
        "changed_tags": changed,
        "ignored_inert_tags": ignored,
        "recorded_incar": recorded,
        "current_incar": current,
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    print(json.dumps(report, indent=2, sort_keys=True, default=str))
    if changed:
        raise RuntimeError(
            f"{manifest_path} was produced by a superseded {args.stage} recipe "
            f"(changed tags: {', '.join(changed)})"
        )


def parse_mesh_conf(path: Path) -> dict[str, Any]:
    """Minimal reader for the phonopy mesh.conf written by the shell script."""
    values: dict[str, Any] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, _, raw = line.partition("=")
        tokens = raw.split()
        values[key.strip().upper()] = tokens[0] if len(tokens) == 1 else tokens
    return values


# INCAR tags that determine the Hamiltonian and the charge-density
# representation, and whose VASP default depends on other tags rather than
# being a fixed number. PREC is one of these: its default is "Accurate" when
# ENCUT is > 1.3 * ENMAX and "Normal" otherwise, so PREC=A and PREC=Accurate
# resolve to the same cutoff table when both stages use the same ENCUT and
# POTCAR. Comparing the raw INCAR strings would flag that benign spelling
# difference; the effective values are compared instead, using VASP's own rule.
EFFECTIVE_DEFAULT_TAGS = frozenset({"PREC"})

# Tags with a fixed VASP default, independent of any other tag. Comparing the
# raw INCAR values is exact for these.
HAMILTONIAN_TAGS = (
    "ENCUT", "ENAUG", "PREC", "LASPH", "LMAXMIX", "LREAL", "ISPIN", "GGA",
)


def normalise_tag(tag: str, value: Any, incar: dict[str, Any]) -> Any:
    """Resolve a tag to the value VASP will actually use."""
    if value is None and tag in EFFECTIVE_DEFAULT_TAGS:
        # VASP: PREC defaults from ENCUT and the POTCAR's ENMAX.
        # ENCUT=520 eV > 1.3 * 400 eV (Ti_sv / V_sv ENMAX), so "Accurate".
        return "Accurate" if float(incar.get("ENCUT", 0)) > 520 else "Normal"
    if tag == "PREC" and isinstance(value, str):
        # "A" and "Accurate" select the same precision.
        return {"A": "Accurate", "Accurate": "Accurate"}.get(value, value)
    return value


def cmd_check_stage_compat(args: argparse.Namespace) -> None:
    """Fail-closed check that the seed and DFPT Hamiltonians agree.

    Job 53690 seeded DFPT with a density computed at ENAUG=1360 / LASPH=True /
    LMAXMIX=4 while DFPT itself ran at VASP's defaults (ENAUG=526.6, LASPH=False,
    LMAXMIX=2). Measured on this cell, the two Hamiltonians differ by 5.01 meV
    (seed -132.756450 eV vs DFPT's ground-state plateau at -132.751443 eV);
    with the tags dropped, seed and DFPT agree to 0.29 meV. That 5 meV is a
    genuine inconsistency between the density supplied and the Hamiltonian
    solved, so the check is worth keeping.

    It is not, however, what makes DFPT diverge: job 53692 ran with the matched
    Hamiltonian and its ground-state SCF still had to be abandoned at NELM. See
    DFPT_INCAR for that failure's actual driver.

    Rule: for every tag in ``HAMILTONIAN_TAGS``, the value VASP will actually
    use must agree between the two stages. A tag one stage sets while the other
    leaves to a default is rejected, because the mismatch is then invisible in
    the INCAR. Values are compared after resolving VASP's own defaults and
    spellings (see ``normalise_tag``).
    """
    seed_manifest = Path(args.seed_root).resolve() / "input-manifest.json"
    dfpt_manifest = Path(args.dfpt_root).resolve() / "input-manifest.json"
    for path in (seed_manifest, dfpt_manifest):
        if not path.is_file():
            raise FileNotFoundError(f"missing input manifest: {path}")

    seed_incar = json.loads(seed_manifest.read_text(encoding="utf-8"))["incar"]
    dfpt_incar = json.loads(dfpt_manifest.read_text(encoding="utf-8"))["incar"]

    rows: list[dict[str, Any]] = []
    failures: list[str] = []
    for tag in HAMILTONIAN_TAGS:
        seed_raw = seed_incar.get(tag, None)
        dfpt_raw = dfpt_incar.get(tag, None)
        seed_value = normalise_tag(tag, seed_raw, seed_incar)
        dfpt_value = normalise_tag(tag, dfpt_raw, dfpt_incar)
        if seed_raw is None and dfpt_raw is None:
            status = "both-unset"
        elif seed_value == dfpt_value:
            status = "match" if seed_raw == dfpt_raw else "match-default"
        else:
            status = "MISMATCH"
            failures.append(tag)
        rows.append(
            {
                "tag": tag,
                "seed": seed_value,
                "dfpt": dfpt_value,
                "seed_raw": seed_raw,
                "dfpt_raw": dfpt_raw,
                "status": status,
            }
        )

    report = {
        "seed_root": str(Path(args.seed_root).resolve()),
        "dfpt_root": str(Path(args.dfpt_root).resolve()),
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "tags": rows,
        "failures": failures,
        "status": "ok" if not failures else "failed",
        "note": (
            "Tags outside HAMILTONIAN_TAGS may differ freely: they change how the "
            "SCF is solved, not the density the seed hands over."
        ),
    }
    print(json.dumps(report, indent=2, sort_keys=True, default=str))
    if failures:
        raise RuntimeError(
            "seed/DFPT Hamiltonian mismatch on: "
            + ", ".join(failures)
            + " - the seed CHGCAR would not be the converged density of the "
            "Hamiltonian DFPT solves"
        )


def cmd_check_phonopy_inputs(args: argparse.Namespace) -> None:
    """Audit the DFPT -> phonopy structure/force-constant data flow.

    The DFPT VASP cell is phonopy's 2x2x2 ``SPOSCAR`` (16 sites), while the
    phonon calculation must be built from the 2-site ``POSCAR-unitcell`` so the
    supercell phonopy reconstructs is the same 16-site cell the forces came
    from. Feeding the 16-site ``SPOSCAR`` in as the unit cell applies ``DIM=2 2
    2`` a second time (128 sites) and compacting the 16x16 force constants then
    fails with ``IndexError: index 16 is out of bounds for axis 0 with size 16``.

    This command is read-only and never invokes VASP or phonopy.
    """
    import numpy as np
    from phonopy.interface.vasp import read_vasp
    from phonopy.structure.cells import get_supercell

    root = Path(args.dfpt_root).resolve()
    unitcell_path = root / "POSCAR-unitcell"
    sposcar_path = root / "SPOSCAR"
    fc_path = root / "FORCE_CONSTANTS"
    mesh_conf_path = root / "mesh.conf"

    checks: dict[str, Any] = {}
    failures: list[str] = []

    def record(name: str, ok: bool, detail: Any = None) -> None:
        checks[name] = {"ok": ok, "detail": detail}
        if not ok:
            failures.append(name)

    def require(path: Path) -> bool:
        ok = path.is_file() and path.stat().st_size > 0
        record(f"present::{path.name}", ok, str(path))
        return ok

    report: dict[str, Any] = {
        "dfpt_root": str(root),
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "checks": checks,
        "failures": failures,
    }

    missing = [p.name for p in (unitcell_path, sposcar_path, fc_path, mesh_conf_path) if not require(p)]
    if missing:
        report["status"] = "failed"
        print(json.dumps(report, indent=2, sort_keys=True, default=str))
        raise RuntimeError(f"missing phonopy inputs under {root}: {missing}")

    mesh_conf = parse_mesh_conf(mesh_conf_path)
    dim_raw = mesh_conf.get("DIM")
    dim = [int(x) for x in (dim_raw if isinstance(dim_raw, list) else str(dim_raw).split())]
    report["mesh_conf"] = mesh_conf
    report["dim"] = dim
    record("dim_is_3_positive_integers", len(dim) == 3 and all(d > 0 for d in dim), dim)

    unitcell = read_vasp(unitcell_path)
    sposcar = read_vasp(sposcar_path)
    report["unitcell_site_count"] = len(unitcell)
    report["sposcar_site_count"] = len(sposcar)
    record("unitcell_site_count_is_2", len(unitcell) == 2, len(unitcell))

    expected_sites = math.prod(dim) * len(unitcell)
    record(
        "sposcar_matches_dim_times_unitcell",
        len(sposcar) == expected_sites,
        f"SPOSCAR has {len(sposcar)} sites, DIM={dim} x {len(unitcell)} = {expected_sites}",
    )

    # The decisive check: does phonopy rebuild exactly the cell VASP ran on?
    rebuilt = get_supercell(unitcell, np.diag(dim))
    report["rebuilt_supercell_site_count"] = len(rebuilt)
    same_symbols = list(rebuilt.symbols) == list(sposcar.symbols)
    same_cell = np.allclose(rebuilt.cell, sposcar.cell, atol=1e-6)
    delta = np.abs(
        np.mod(rebuilt.scaled_positions - sposcar.scaled_positions + 0.5, 1.0) - 0.5
    )
    same_positions = bool(np.allclose(delta, 0.0, atol=1e-6))
    report["max_scaled_position_delta"] = float(delta.max()) if delta.size else None
    record("supercell_symbol_order_matches_sposcar", same_symbols)
    record("supercell_lattice_matches_sposcar", bool(same_cell))
    record("supercell_positions_match_sposcar", same_positions)

    header = fc_path.read_text(encoding="utf-8").splitlines()[0].split()
    fc_shape = tuple(int(x) for x in header[:2]) if len(header) >= 2 else None
    report["force_constants_shape"] = fc_shape
    record(
        "force_constants_match_supercell",
        fc_shape == (len(sposcar), len(sposcar)),
        f"FORCE_CONSTANTS {fc_shape} vs SPOSCAR {len(sposcar)} sites",
    )

    if args.vasprun:
        vasprun_path = Path(args.vasprun).resolve()
        if require(vasprun_path):
            from pymatgen.io.vasp.outputs import Vasprun

            parsed = Vasprun(
                str(vasprun_path),
                parse_dos=False,
                parse_eigen=False,
                parse_projected_eigen=False,
            )
            final = parsed.final_structure
            report["vasprun_site_count"] = len(final)
            record("vasprun_site_count_matches_sposcar", len(final) == len(sposcar), len(final))
            vp = np.array([site.frac_coords for site in final])
            record(
                "vasprun_positions_match_sposcar",
                bool(np.allclose(np.mod(vp - sposcar.scaled_positions + 0.5, 1.0) - 0.5, 0.0, atol=1e-6)),
            )
            record(
                "vasprun_species_match_sposcar",
                [str(site.specie) for site in final] == list(sposcar.symbols),
            )

    # The DFPT job reads ICHARG=1, so a live CHGCAR must sit in the DFPT
    # directory in exactly the cell VASP will run on (16-site SPOSCAR). A
    # 2-site seed-derived CHGCAR would be silently ignored (or, worse, misread).
    chgcar_path = root / "CHGCAR"
    if require(chgcar_path):
        lines = chgcar_path.read_text(errors="replace")[:8192].splitlines()
        try:
            chgcar_sites = sum(int(x) for x in lines[6].split())
        except (IndexError, ValueError):
            chgcar_sites = None
        report["chgcar_site_count"] = chgcar_sites
        report["chgcar_comment"] = lines[0] if lines else ""
        record(
            "chgcar_present_for_dfpt",
            True,
            detail=f"{chgcar_path} ({chgcar_path.stat().st_size} bytes)",
        )
        record(
            "chgcar_matches_supercell",
            chgcar_sites == len(sposcar),
            f"CHGCAR declares {chgcar_sites} sites, SPOSCAR has {len(sposcar)}",
        )

    report["status"] = "ok" if not failures else "failed"
    print(json.dumps(report, indent=2, sort_keys=True, default=str))
    if failures:
        raise RuntimeError(
            "phonopy post-processing inputs are inconsistent: " + ", ".join(failures)
        )


def cmd_record_e(args: argparse.Namespace) -> None:
    from pymatgen.io.vasp.outputs import Vasprun

    path = Path(args.vasprun).resolve()
    parsed = Vasprun(
        str(path),
        parse_dos=False,
        parse_eigen=False,
        parse_projected_eigen=False,
    )
    record = {
        "scale": float(args.scale),
        "volume_ang3": float(parsed.final_structure.volume),
        "energy_eV": float(parsed.final_energy),
        "vasprun": str(path),
        "vasprun_sha256": sha256(path),
    }
    Path(args.output).resolve().write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def cmd_write_ve(args: argparse.Namespace) -> None:
    records = [json.loads(Path(item).read_text(encoding="utf-8")) for item in args.records]
    records.sort(key=lambda item: (float(item["volume_ang3"]), str(item["vasprun"])))
    if len(records) != 9:
        raise ValueError(f"expected exactly 9 volume records; got {len(records)}")
    output = Path(args.output).resolve()
    output.write_text(
        "".join(
            f"{item['volume_ang3']:.10f} {item['energy_eV']:.12f}\n"
            for item in records
        ),
        encoding="utf-8",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("check-model")
    p.add_argument("--root", required=True)
    p.add_argument("--potcar-root", default=str(POTCAR_ROOT_DEFAULT))
    p.set_defaults(func=cmd_check_model)

    p = sub.add_parser("audit-kpoints")
    p.add_argument("--root", required=True)
    p.add_argument("--potcar-root", default=str(POTCAR_ROOT_DEFAULT))
    p.set_defaults(func=cmd_audit_kpoints)

    p = sub.add_parser("prepare")
    p.add_argument("--root", required=True)
    p.add_argument("--potcar-root", default=str(POTCAR_ROOT_DEFAULT))
    p.add_argument("--profile", choices=tuple(PROFILE_MESHES), required=True)
    p.add_argument("--stage", choices=STAGES, required=True)
    p.add_argument(
        "--dfpt-variant",
        choices=tuple(DFPT_VARIANTS),
        default=DEFAULT_DFPT_VARIANT,
        help="DFPT electronic-structure experiment (see DFPT_VARIANTS)",
    )
    p.add_argument("--outdir", required=True)
    p.add_argument("--structure")
    p.add_argument("--scale", type=float)
    p.set_defaults(func=cmd_prepare)

    p = sub.add_parser("verify")
    p.add_argument("--vasprun", required=True)
    p.add_argument("--outcar", required=True)
    p.add_argument("--stage", choices=STAGES)
    p.add_argument("--report", help="write the JSON audit report to this path")
    p.add_argument("--max-force-tolerance", type=float)
    p.add_argument("--fail-on-force", action="store_true")
    p.set_defaults(func=cmd_verify)

    p = sub.add_parser("check-phonopy-inputs")
    p.add_argument("--dfpt-root", required=True)
    p.add_argument("--vasprun")
    p.set_defaults(func=cmd_check_phonopy_inputs)

    p = sub.add_parser("check-stage-compat")
    p.add_argument("--seed-root", required=True)
    p.add_argument("--dfpt-root", required=True)
    p.set_defaults(func=cmd_check_stage_compat)

    p = sub.add_parser("check-recipe")
    p.add_argument("--manifest", required=True)
    p.add_argument("--stage", choices=STAGES, required=True)
    p.add_argument(
        "--dfpt-variant",
        choices=tuple(DFPT_VARIANTS),
        default=DEFAULT_DFPT_VARIANT,
    )
    p.set_defaults(func=cmd_check_recipe)

    p = sub.add_parser("record-e")
    p.add_argument("--vasprun", required=True)
    p.add_argument("--scale", required=True, type=float)
    p.add_argument("--output", required=True)
    p.set_defaults(func=cmd_record_e)

    p = sub.add_parser("write-ve")
    p.add_argument("--records", nargs="+", required=True)
    p.add_argument("--output", required=True)
    p.set_defaults(func=cmd_write_ve)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.command == "prepare" and args.stage != "initial-relax" and not args.structure:
        raise SystemExit("--structure is required for volume-static, charge-seed and dfpt")
    if args.command == "prepare" and args.stage == "volume-static" and args.scale is None:
        raise SystemExit("--scale is required for volume-static")
    try:
        args.func(args)
    except VerificationError as exc:
        # Distinct exit code: the files exist and the job "finished", but the
        # result is not scientifically usable.
        print(f"VERIFICATION FAILED: {exc}", file=sys.stderr)
        return 3
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())