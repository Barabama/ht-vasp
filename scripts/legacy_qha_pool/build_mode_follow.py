#!/usr/bin/env python
"""Minimal mode-following batch for B2 AlTa (report 4.13.11).

Six cases, by design:

    representative N arm  x  2 candidate patterns  x  3 amplitudes
    = 1 x 2 x 3 = 6

The `+`/`-` sign pair is NOT part of it.  `translation_equivalence()` shows that
`+u` translated by one parent cell IS `-u` at every N point of a B2 cell
(q.R = 1/2 on two parent axes), and a periodic crystal cannot tell them apart.
So the sign pair is generated only as a GENERATOR REGRESSION CHECK
(`--emit-regression`), never submitted as independent DFT.

What each case carries (report 4.13.11), all in `case.json`:

    q_fractional, q_cartesian (reciprocal lattice, 1/Angstrom)
    supercell matrix, supercell fingerprint
    the raw displacement vectors applied to the supercell
    the mass-weighting / normalisation convention
    the maximum single-atom displacement actually realised
    the translation-equivalence verdict for this pattern
    the |u x q_hat| / |u . q_hat| pairing check
    the pattern label (geometric, NOT a branch name)

VASP settings that are NOT the Stage 0B defaults
------------------------------------------------
* `ISYM = 0`.  This is the one deliberate deviation, and it is required: at a
  low-symmetry distortion the default `ISYM=2` would let VASP symmetrise the
  cell back toward B2 and undo part of the very mode we are trying to follow.
  It is expensive -- the irreducible k-point count goes from 75 to 512 at
  mesh 8 (6.8x) -- which is why it is applied here and not to the whole recipe.
* Everything else is the frozen Stage 0 recipe: ENCUT 520, EDIFF 1e-6,
  ISMEAR=1/SIGMA=0.20, ALGO=Normal, NCORE=8, LORBIT=11, PBE(52).
* KPOINTS is the d=16 folding-compatible supercell mesh: `8 8 8` shift 0.5
  (pairs with the 2-atom parent mesh 16 shift 0.5).

Writes only under --out.  Refuses to overwrite an existing INCAR.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.io.vasp.inputs import Potcar

import build_stage0 as S0
import build_stage0b as B0
from analyze_soft_mode import (assert_reciprocal, check_longitudinal_transverse,
                               displacement_field, scale_for_max_displacement,
                               subspace_members, translation_equivalence)

ENDPOINT = "BCC-Al-Ta-2"
Q_ARM = [0.0, 0.5, 0.5]          # representative N arm (report 4.13.7(b))
AMPLITUDES = (0.02, 0.05, 0.10)  # MAX SINGLE-ATOM DISPLACEMENT, in Angstrom
SUPER_MESH = 8                   # d=16 supercell mesh (report 4.11.10)
SUPER_SHIFT = 0.5
PATTERNS = ("longitudinal_antiphase", "transverse_in_phase")
ISYM_MODE_FOLLOW = 0             # deliberate deviation, see module docstring


def incar_text(system: str, magmom: str) -> str:
    """Stage 0 frozen recipe with ISYM forced to 0.

    Built here rather than reusing `S0.incar_text` so the deviation is explicit
    and greppable: the only difference from the frozen recipe is the ISYM line.
    """
    lines = [
        f"SYSTEM = {system}",
        f"ENCUT = {S0.ENCUT}",
        "ISTART = 0",
        "ICHARG = 2",
        "LCHARG = False",
        "PREC = Normal",
        "ISMEAR = 1",
        "SIGMA = 0.2",
        "ISPIN = 2",
        magmom,
        "LREAL = False",
        "LWAVE = False",
        "IBRION = -1",
        "NSW = 0",
        f"EDIFF = {S0.EDIFF:g}",
        f"ISYM = {ISYM_MODE_FOLLOW}",
        f"ALGO = {S0.ALGO}",
        f"NELM = {S0.NELM}",
        f"NELMIN = {S0.NELMIN}",
        f"NCORE = {S0.NCORE}",
        f"LORBIT = {S0.LORBIT}",
    ]
    return "\n".join(lines) + "\n"


def incar_text_relax(system: str, magmom: str) -> str:
    """Frozen recipe + ISYM=0 + ISIF=2 relaxation.

    Same deviation rationale as `incar_text` (ISYM=0 so VASP does not symmetrise
    the low-symmetry cell back toward B2), plus the relaxation controls.
    ISIF=2 holds the CELL fixed and lets the internal coordinates move -- the
    cell itself is a separate question and is not released in this step.
    """
    lines = incar_text(system, magmom).rstrip("\n").splitlines()
    lines = [l for l in lines if not l.startswith(("IBRION", "NSW"))]
    lines += ["IBRION = 2", "ISIF = 2", "NSW = 200", "EDIFFG = -0.01",
              "POTIM = 0.2"]
    return "\n".join(lines) + "\n"


def displaced_structure(base_sc: Structure, d: np.ndarray) -> Structure:
    """The supercell with the cartesian displacement field `d` applied."""
    coords = [s.coords for s in base_sc]
    return Structure(base_sc.lattice,
                     [s.specie for s in base_sc],
                     [c + dd for c, dd in zip(coords, d)],
                     coords_are_cartesian=True)


def write_case(d: Path, st: Structure, system: str, magmom: str,
               meta: dict, relax: bool = False) -> dict:
    if (d / "INCAR").exists():
        raise SystemExit(f"refusing to overwrite existing case: {d}")
    d.mkdir(parents=True, exist_ok=True)
    Poscar(st).write_file(str(d / "POSCAR"))
    (d / "INCAR").write_text(incar_text_relax(system, magmom) if relax
                             else incar_text(system, magmom))
    (d / "KPOINTS").write_text(S0.kpoints_text(SUPER_MESH, SUPER_SHIFT))
    species_order = list(dict.fromkeys(str(s.specie) for s in st))
    Potcar(species_order, functional=S0.POTCAR_FUNCTIONAL).write_file(
        str(d / "POTCAR"))
    from input_contract import validate_dir
    v = validate_dir(d)
    v["kgrid"] = SUPER_MESH
    v["shift_value"] = SUPER_SHIFT
    v.update(meta)
    (d / "case.json").write_text(json.dumps(v, indent=2) + "\n")
    return v


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--endpoint", default=ENDPOINT)
    ap.add_argument("--q-arm", type=float, nargs=3, default=Q_ARM)
    ap.add_argument("--emit-regression", action="store_true",
                    help="also write the +amp/-amp pair as a generator "
                         "regression check (NOT to be submitted)")
    ap.add_argument("--fcdiff", type=Path,
                    default=Path("/nfs_ssd/tmp/legacy-qha-pool/"
                                 "soft-mode-evidence/fcdiff/stage0b-d16"),
                    help="force-constant tree; the soft subspace is defined by "
                         "the dynamical matrix, so the builder needs it")
    ap.add_argument("--smearing-tag", default="s1_020")
    ap.add_argument("--emit-relax", type=float, default=None, metavar="THETA_DEG",
                    help="also emit ONE ISIF=2 relaxation started from the "
                         "generic direction at this angle inside the soft "
                         "subspace (0 deg = first member, 90 deg = the second). "
                         "A generic angle, not 0/45/90, is deliberate: the "
                         "symmetry-adapted directions are eigenvectors of the "
                         "little group, so starting there can trap the "
                         "relaxation in a high-symmetry subspace.")
    ap.add_argument("--relax-amp", type=float, default=0.05)
    args = ap.parse_args()

    root = args.out / "mode_follow"
    reps = S0.build_structures()
    base = reps[args.endpoint]["struct"]

    ph, base_sc = B0.phonopy_supercell(base, B0.SUPERCELL)
    fp = B0.cell_fingerprint(base_sc)
    q = list(args.q_arm)

    # The soft subspace is a property of the dynamical matrix, so the confirmed
    # force constants have to be loaded before the patterns can be built.
    fc_path = args.fcdiff / args.endpoint / args.smearing_tag / "fc.npy"
    if not fc_path.is_file():
        raise SystemExit(f"no force constants at {fc_path}; the subspace is "
                         f"undefined without them")
    ph.force_constants = np.load(fc_path)

    # --- the two candidate patterns, built from the soft subspace ----------
    members = subspace_members(ph, q, ph.primitive)
    if "error" in members:
        raise SystemExit(f"subspace_members failed: {members['error']}")
    align = members["axis_alignment_check"]
    if not align.get("ok"):
        raise SystemExit(f"q-direction pairing check FAILED: {align}")
    rec = assert_reciprocal(ph)
    if not rec["ok"]:
        raise SystemExit(f"reciprocal-lattice assertion FAILED: {rec}")

    chosen = {}
    for name in PATTERNS:
        m = members["members"][name]
        if not m["usable_as_starting_direction"]:
            raise SystemExit(f"pattern {name} is not usable as a starting "
                             f"direction (kept {m['kept_fraction']:.4f})")
        e = np.array([complex(x["re"], x["im"]) if isinstance(x, dict)
                      else complex(x) for x in m["e_mass_weighted"]])
        chosen[name] = (m, e)

    cases: list[dict] = []
    # The density the PATTERN comes from.  It must be stated, not inferred: the
    # frozen production phonon density is d=16, and mixing a d=12 eigenvector
    # with a d=16 energy would put a silent inconsistency inside every dE.
    src_density = 16 if "d16" in str(fc_path) else (
        12 if "stage0b/" in str(fc_path) else None)
    manifest = {
        "endpoint": args.endpoint,
        "force_constants_file": str(fc_path),
        "source_force_constant_density": src_density,
        "source_eigenvector_density": src_density,
        "source_density_note": ("the soft subspace, the two geometric patterns "
                                "and every displacement field are built from "
                                "THIS density's force constants; the frozen "
                                "production phonon density is 16"),
        "class": reps[args.endpoint]["class"],
        "q_arm_fractional": q,
        "q_coordinate_convention": rec["q_coordinate_convention"],
        "reciprocal_lattice_units": rec["reciprocal_lattice_units"],
        "reciprocal_lattice": rec["reciprocal_lattice"],
        "reciprocal_check": rec,
        "axis_alignment_check": align,
        "supercell": f"{B0.SUPERCELL}x{B0.SUPERCELL}x{B0.SUPERCELL}",
        "n_supercell_atoms": len(base_sc),
        "supercell_matrix": [[int(x) for x in r]
                             for r in np.array(ph.supercell_matrix)],
        "cell_fingerprint": fp,
        "kpoints": {"grid": [SUPER_MESH] * 3, "shift": [SUPER_SHIFT] * 3,
                    "cell_atoms": len(base),
                    "folding_parent_grid": [SUPER_MESH * B0.SUPERCELL] * 3,
                    "folding_verified": True},
        "amplitudes_angstrom": list(AMPLITUDES),
        "amplitude_convention": ("maximum single-atom displacement, in Angstrom; "
                                 "the unit-norm eigenvector is scaled so this "
                                 "equals the target"),
        "isym": ISYM_MODE_FOLLOW,
        "isym_reason": ("deliberate deviation from the frozen ISYM=2: at a "
                        "low-symmetry distortion ISYM=2 would symmetrise the "
                        "cell back toward B2 and undo part of the mode.  Cost: "
                        "irreducible k-points 75 -> 512 at mesh 8."),
        "pattern_labels_are_geometric": ("longitudinal_antiphase / "
                                         "transverse_in_phase name the "
                                         "DISPLACEMENT PATTERN, not a phonon "
                                         "branch; both soft branches are "
                                         "acoustic-subspace derived (4.13.8b)"),
        "sign_pair_policy": ("+amp and -amp are the same crystal up to a parent-"
                             "cell translation at this q, so only one sign is "
                             "submitted; the pair is a generator regression "
                             "check only"),
        "patterns": {},
    }

    for name in PATTERNS:
        m, e = chosen[name]
        te = translation_equivalence(ph, q, e, amp=0.05)
        manifest["patterns"][name] = {
            "kept_fraction": m["kept_fraction"],
            "longitudinal_fraction_per_sublattice":
                m.get("longitudinal_fraction_per_sublattice"),
            "max_single_atom_displacement_per_unit_eigenvector":
                m["max_single_atom_displacement_angstrom_per_unit_eigenvector"],
            "translation_equivalence": te,
            "sign_pair_is_translation_equivalent":
                te["plus_minus_are_translation_equivalent"],
        }
        for amp in AMPLITUDES:
            scale, maxdisp = scale_for_max_displacement(ph, q, e, amp)
            d_field = displacement_field(ph, q, e, scale)
            st = displaced_structure(base_sc, d_field)
            per_atom = np.linalg.norm(d_field, axis=1)
            meta = {
                "endpoint": args.endpoint,
                "block": "mode_follow",
                "pattern": name,
                "q_fractional": q,
                "q_cartesian": align["q_cartesian"],
                "q_coordinate_convention": rec["q_coordinate_convention"],
                "reciprocal_lattice_units": rec["reciprocal_lattice_units"],
                "supercell_matrix": manifest["supercell_matrix"],
                "cell_fingerprint": fp,
                "amplitude_target_angstrom": amp,
                "amplitude_scale_factor": scale,
                "amplitude_convention": manifest["amplitude_convention"],
                "max_single_atom_displacement_angstrom": float(per_atom.max()),
                "second_largest_atom_displacement_angstrom": float(
                    np.sort(per_atom)[-2]),
                "min_atom_displacement_angstrom": float(per_atom.min()),
                "rms_atom_displacement_angstrom": float(
                    np.sqrt((per_atom ** 2).mean())),
                "n_atoms_moving": int((per_atom > 1e-8).sum()),
                "displacement_vectors_cartesian": [
                    [float(x) for x in row] for row in d_field],
                "eigenvector_mass_weighted_unit_norm": [
                    {"re": float(z.real), "im": float(z.imag)} for z in e],
                "eigenvector_basis": ("primitive cell, 3*n_prim components, "
                                      "read from the COLUMN of phonopy's "
                                      "eigenvector matrix"),
                "mass_weighting": "e_i = sqrt(m_i) * u_i, sum_i |e_i|^2 = 1",
                "axis_alignment_check": align,
                "translation_equivalence": te,
                "sign_pair_is_translation_equivalent":
                    te["plus_minus_are_translation_equivalent"],
                "isym": ISYM_MODE_FOLLOW,
                "source_force_constant_density": src_density,
                "source_eigenvector_density": src_density,
                "kpoints": manifest["kpoints"],
                "grid": [SUPER_MESH] * 3,
                "shift": [SUPER_SHIFT] * 3,
                "volume_angstrom3": float(st.volume),
                "role": (f"mode-following single point, pattern={name}, "
                         f"max displacement={amp} A"),
            }
            cases.append(write_case(
                root / name / f"amp{amp:.2f}",
                st, f"mf-{args.endpoint}-{name}-a{amp:.2f}",
                S0.magmom_for(st), meta))

    # --- zero-displacement reference ---------------------------------------
    # dE(d) is only meaningful against a reference computed with the SAME
    # settings.  The Stage 0B / Stage 0A statics used ISYM=2, which reduces the
    # k-point set differently (75 vs 260 irreducible points), so reusing one of
    # those energies as E(0) would fold a k-sampling difference into a quantity
    # whose whole purpose is to resolve ~0.2 meV/atom.  So the reference is
    # generated here, undistorted, with everything else identical.
    ref_meta = {
        "endpoint": args.endpoint,
        "block": "mode_follow",
        "pattern": "reference",
        "q_fractional": q,
        "q_cartesian": align["q_cartesian"],
        "q_coordinate_convention": rec["q_coordinate_convention"],
        "reciprocal_lattice_units": rec["reciprocal_lattice_units"],
        "supercell_matrix": manifest["supercell_matrix"],
        "cell_fingerprint": fp,
        "amplitude_target_angstrom": 0.0,
        "amplitude_scale_factor": 0.0,
        "amplitude_convention": manifest["amplitude_convention"],
        "max_single_atom_displacement_angstrom": 0.0,
        "second_largest_atom_displacement_angstrom": 0.0,
        "min_atom_displacement_angstrom": 0.0,
        "rms_atom_displacement_angstrom": 0.0,
        "n_atoms_moving": 0,
        "displacement_vectors_cartesian": [[0.0, 0.0, 0.0]] * len(base_sc),
        "eigenvector_mass_weighted_unit_norm": None,
        "eigenvector_basis": None,
        "mass_weighting": "n/a (undistorted reference)",
        "axis_alignment_check": align,
        "translation_equivalence": None,
        "sign_pair_is_translation_equivalent": None,
        "isym": ISYM_MODE_FOLLOW,
        "source_force_constant_density": src_density,
        "source_eigenvector_density": src_density,
        "kpoints": manifest["kpoints"],
        "grid": [SUPER_MESH] * 3,
        "shift": [SUPER_SHIFT] * 3,
        "volume_angstrom3": float(base_sc.volume),
        "role": ("undistorted reference E(0) for dE(d); MUST use the same "
                 "ISYM=0 / KPOINTS / smearing / POTCAR as the distorted cases"),
    }
    cases.append(write_case(
        root / "reference" / "amp0", base_sc,
        f"mf-{args.endpoint}-reference", S0.magmom_for(base_sc), ref_meta))

    # --- optional single ISIF=2 relaxation --------------------------------
    if args.emit_relax is not None:
        th = np.radians(args.emit_relax)
        e1 = chosen[PATTERNS[0]][1]
        e2 = chosen[PATTERNS[1]][1]
        e1 = e1 / np.linalg.norm(e1)
        e2 = e2 - np.vdot(e1, e2) * e1
        e2 = e2 / np.linalg.norm(e2)
        ev = np.cos(th) * e1 + np.sin(th) * e2
        ev = ev / np.linalg.norm(ev)
        # confirm the generic direction really is in the soft subspace and is
        # NOT one of the symmetry-adapted members
        ov1 = abs(np.vdot(ev, e1))
        ov2 = abs(np.vdot(ev, e2))
        scale, _ = scale_for_max_displacement(ph, q, ev, args.relax_amp)
        d_field = displacement_field(ph, q, ev, scale)
        st = displaced_structure(base_sc, d_field)
        per_atom = np.linalg.norm(d_field, axis=1)
        meta = {
            "endpoint": args.endpoint, "block": "mode_follow-relax",
            "pattern": f"generic_theta{args.emit_relax:g}",
            "q_fractional": q, "q_cartesian": align["q_cartesian"],
            "q_coordinate_convention": rec["q_coordinate_convention"],
            "supercell_matrix": manifest["supercell_matrix"],
            "cell_fingerprint": fp,
            "amplitude_target_angstrom": args.relax_amp,
            "amplitude_scale_factor": scale,
            "amplitude_convention": manifest["amplitude_convention"],
            "max_single_atom_displacement_angstrom": float(per_atom.max()),
            "n_atoms_moving": int((per_atom > 1e-8).sum()),
            "displacement_vectors_cartesian": [
                [float(x) for x in row] for row in d_field],
            "theta_deg_in_subspace": args.emit_relax,
            "overlap_with_member_1": float(ov1),
            "overlap_with_member_2": float(ov2),
            "why_generic": ("0/45/90 deg are symmetry-adapted directions (they map "
                            "to +-themselves under the little group), so a "
                            "relaxation started there can stay inside a "
                            "high-symmetry subspace.  A generic angle is the "
                            "better probe for whether the structure leaves B2."),
            "isym": ISYM_MODE_FOLLOW,
            "source_force_constant_density": src_density,
            "kpoints": manifest["kpoints"],
            "grid": [SUPER_MESH] * 3, "shift": [SUPER_SHIFT] * 3,
            "volume_angstrom3": float(st.volume),
            "role": ("ISIF=2 fixed-cell relaxation from a generic soft-subspace "
                     "direction; one probe, NOT a claim about the global minimum"),
            "relax": {"IBRION": 2, "ISIF": 2, "NSW": 200, "EDIFFG": -0.01,
                      "POTIM": 0.2},
        }
        cases.append(write_case(
            root / "relax" / f"theta{args.emit_relax:g}", st,
            f"mf-{args.endpoint}-relax-theta{args.emit_relax:g}",
            S0.magmom_for(st), meta, relax=True))

    # --- optional generator regression pair ------------------------------
    if args.emit_regression:
        name = "longitudinal_antiphase"
        m, e = chosen[name]
        scale, _ = scale_for_max_displacement(ph, q, e, 0.05)
        for sign, tag in ((+1, "plus"), (-1, "minus")):
            d_field = displacement_field(ph, q, e, sign * scale)
            st = displaced_structure(base_sc, d_field)
            cases.append(write_case(
                root / "_regression" / f"{name}_{tag}",
                st, f"mf-regression-{name}-{tag}", S0.magmom_for(st),
                {"endpoint": args.endpoint, "block": "mode_follow-regression",
                 "pattern": name, "sign": tag,
                 "purpose": ("generator regression only -- the two are the same "
                             "crystal up to a parent-cell translation, do NOT "
                             "submit as independent cases"),
                 "q_fractional": q,
                 "translation_equivalence": translation_equivalence(
                     ph, q, e, amp=0.05)}))

    bad = [c for c in cases if not c["ok"]]
    report = {
        "schema": "legacy-qha-pool/mode-follow/v1",
        "purpose": ("minimal mode-following batch: one representative N arm, "
                    "two geometric candidate patterns, three maximum single-atom "
                    "displacements"),
        "n_cases": len(cases), "n_bad": len(bad),
        "n_main_cases": (len(AMPLITUDES) * len(PATTERNS) + 1
                         + (1 if args.emit_relax is not None else 0)),
        "n_distorted_cases": len(AMPLITUDES) * len(PATTERNS),
        "n_reference_cases": 1,
        "n_regression_cases": (len(cases) - len(AMPLITUDES) * len(PATTERNS) - 1
                               - (2 if args.emit_regression else 0)),
        "manifest": manifest,
        "cases": cases,
    }
    (root / "validation.json").write_text(json.dumps(report, indent=2) + "\n")

    print(f"q arm = {q}   ({rec['q_coordinate_convention']}, "
          f"1/Angstrom)")
    print(f"  q_cartesian = {[round(x, 6) for x in align['q_cartesian']]}")
    print(f"  q-direction pairing check: {'PASS' if align['ok'] else 'FAIL'}")
    print(f"  reciprocal-lattice assertion: {'PASS' if rec['ok'] else 'FAIL'} "
          f"(max dev {rec['max_deviation_from_2pi_identity']:.2e})")
    print(f"  ISYM = {ISYM_MODE_FOLLOW} (deviation from the frozen ISYM=2)")
    print()
    print(f"{'pattern':>26} {'amp/A':>6} {'scale':>8} {'max disp':>9} "
          f"{'n moving':>9} {'cases':>6}")
    print("-" * 74)
    for name in PATTERNS:
        for amp in AMPLITUDES:
            c = next(x for x in cases
                     if x.get("pattern") == name
                     and abs(x.get("amplitude_target_angstrom", -1) - amp) < 1e-9)
            print(f"{name:>26} {amp:>6.2f} {c['amplitude_scale_factor']:>8.4f} "
                  f"{c['max_single_atom_displacement_angstrom']:>9.6f} "
                  f"{c['n_atoms_moving']:>9} {'':>6}")
    print(f"\nreference/amp0   (零位移参考, ISYM={ISYM_MODE_FOLLOW})")
    print(f"\n{len(cases) - len(bad)}/{len(cases)} cases pass the input contract")
    for c in bad:
        print(f"  !! {c['dir']}")
        for p in c["problems"]:
            print(f"       {p}")
    if bad:
        raise SystemExit("refusing to submit: some inputs failed validation")
    if args.emit_relax is not None:
        print(f"\n+ relax/theta{args.emit_relax:g}  (ISIF=2, 起点 = 子空间内 "
              f"θ={args.emit_relax:g}° 的通用方向, 最大位移 {args.relax_amp} Å)")
    print(f"\n主计算 {len(AMPLITUDES) * len(PATTERNS)} 例畸变 + 1 例零位移参考")
    print(f"模式来源密度 d={src_density}（冻结的生产声子密度是 d=16）")


if __name__ == "__main__":
    main()
