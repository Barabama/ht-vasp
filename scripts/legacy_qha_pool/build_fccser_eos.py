#!/usr/bin/env python
"""Stage 0D: FCC/L1_2 and SER multi-volume static EOS (report 8.3 step 3 / 8.2 #17).

Why this exists
---------------
FCC/L1_2 (49) + SER (7) is **2/3 of the 84-endpoint pool** and had *zero* E(V)
data: Stage 0C only produced one single-volume static plus a phonon supercell
per class.  E(V) is the first brick of QHA -- without a BM3 fit there is no
QHA at all -- so this batch is the coverage step, not a refinement.

It also targets the risk that BCC/B2's smooth progress hides: report 4.2b's
proxy screen flags FCC more often than BCC (84% vs 70% of endpoints with at
least one dropped-mode marker), i.e. the *less* tested class is the *more*
suspect one.

Stages
------
    --stage relax   five relaxations (ISIF=3) to locate V0 per endpoint
    --stage spin    two zero-moment controls at the starting cell (Cr only)
    --stage eos     seven scaled statics per endpoint, built from the relaxed
                    cell of `--stage relax`

`relax` and `spin` take the cell straight from `endmember.py` (the production
generator), so this also exercises that path.  `eos` needs the relaxed cell,
so it runs after `relax`.

Recipe
------
Frozen Stage 0 recipe throughout: ENCUT 520, EDIFF 1e-6, ISMEAR=1/SIGMA=0.20,
ALGO=Normal, NELM 120, NELMIN 6, NCORE 8, LORBIT 11, ISPIN=2, PREC=Normal,
ISYM=2, PBE(52) POTCAR, and an explicit folding-compatible KPOINTS at density
d=16 (`mesh_for(V, 16)`, even + shift 0.5).

NOTE on report 5.1: its EDIFF row says 1e-5, but EVERY Stage 0 batch
(stage0, stage0b, stage0c, sigma, mode_follow) wrote 1e-6.  The freeze table
row is stale; this builder follows the code, i.e. 1e-6, and the report row is
corrected separately.

Writes only under --out.  Refuses to overwrite an existing INCAR.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp import Poscar
from pymatgen.io.vasp.inputs import Potcar

import build_stage0 as S0

DENSITY = 16
# a-axis scales for the EOS scan.  The centre is the RELAXED cell, so a tight
# range is enough; +-3% in a is +-9.3% in volume, which comfortably brackets a
# BM3 minimum.
A_SCALES = (0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03)
# relax controls: the production GLOBAL_INCAR values from main_qha.py
RELAX = {"IBRION": 2, "ISIF": 3, "NSW": 80, "POTIM": 0.2, "EDIFFG": -0.02}
# HCP has a c/a degree of freedom; at fixed volume the shape must be free for
# E(V) to mean anything.  ISIF=4 relaxes the shape at constant volume and is a
# no-op for the cubic endpoints.
HCP_ISIF = 4

TARGETS = (
    # name,         kind,  args                 note
    ("FCC-Al-Cr-4", "fcc", ("Al", "Cr"),
     "Al@1a / Cr@3c -- the direction report 4.2b flags as 'all volume points "
     "dropped'"),
    ("FCC-Cr-Al-4", "fcc", ("Cr", "Al"),
     "Cr@1a / Al@3c -- the direction 4.2b flags as 'no marker'; the pair is "
     "the L1_2 occupancy signal and is NOT a mirror of the other"),
    ("SER-Al-1", "ser", "Al", "SER, fcc, 1-atom primitive"),
    ("SER-Ta-1", "ser", "Ta", "SER, bcc, 1-atom primitive"),
    ("SER-Zr-1", "ser", "Zr", "SER, hcp, 2-atom -- the only non-cubic one"),
)
CR_TARGETS = ("FCC-Al-Cr-4", "FCC-Cr-Al-4")


def endmember_structure(kind: str, args) -> Structure:
    """Build one target through the production generator."""
    from htvasp.model.endmember import Endmember
    em = Endmember()
    if kind == "fcc":
        return em.get_fcc(*args)
    if kind == "ser":
        return em.get_ser(args)
    raise ValueError(kind)


def cell_of(st: Structure) -> dict:
    a, b, c = st.lattice.abc
    al, be, ga = st.lattice.angles
    return {"a": a, "b": b, "c": c, "alpha": al, "beta": be, "gamma": ga,
            "volume": float(st.volume), "n_atoms": len(st),
            "volume_per_atom": float(st.volume) / len(st),
            "formula": st.composition.reduced_formula}


def scaled_cell(st: Structure, s: float) -> Structure:
    """Isotropic a-axis scaling.  Keeps angles, so it works for hexagonal too."""
    return Structure(Lattice.from_parameters(
        st.lattice.a * s, st.lattice.b * s, st.lattice.c * s,
        st.lattice.alpha, st.lattice.beta, st.lattice.gamma),
        [s2.specie for s2 in st],
        [s2.frac_coords for s2 in st],
        coords_are_cartesian=False)


def incar_text(system: str, magmom: str, extra: dict) -> str:
    """Frozen Stage 0 recipe + a per-stage override block."""
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
        f"EDIFF = {S0.EDIFF:g}",
        f"ISYM = {S0.ISYM}",
        f"ALGO = {S0.ALGO}",
        f"NELM = {S0.NELM}",
        f"NELMIN = {S0.NELMIN}",
        f"NCORE = {S0.NCORE}",
        f"LORBIT = {S0.LORBIT}",
    ]
    for k, v in extra.items():
        lines.append(f"{k} = {v}")
    return "\n".join(lines) + "\n"


def write_case(d: Path, st: Structure, system: str, magmom: str,
               extra_incar: dict, meta: dict) -> dict:
    if (d / "INCAR").exists():
        raise SystemExit(f"refusing to overwrite existing case: {d}")
    d.mkdir(parents=True, exist_ok=True)
    Poscar(st).write_file(str(d / "POSCAR"))
    (d / "INCAR").write_text(incar_text(system, magmom, extra_incar))
    n = S0.mesh_for(float(st.volume), DENSITY)
    shift = 0.5 if n % 2 == 0 else 0.0
    (d / "KPOINTS").write_text(S0.kpoints_text(n, shift))
    species_order = list(dict.fromkeys(str(s.specie) for s in st))
    Potcar(species_order, functional=S0.POTCAR_FUNCTIONAL).write_file(
        str(d / "POTCAR"))
    from input_contract import validate_dir
    v = validate_dir(d)
    v["kgrid"] = n
    v["shift_value"] = shift
    v["volume_angstrom3"] = float(st.volume)
    v["n3V"] = n ** 3 * float(st.volume)
    v.update(meta)
    (d / "case.json").write_text(json.dumps(v, indent=2) + "\n")
    return v


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--stage", required=True, choices=("relax", "spin", "eos"))
    ap.add_argument("--relax-root", type=Path, default=None,
                    help="--stage eos: directory holding the relax stage")
    args = ap.parse_args()

    root = args.out / "stage0d"
    cases: list[dict] = []
    manifest: dict = {"stage": args.stage, "density": DENSITY,
                      "recipe": {"ENCUT": S0.ENCUT, "EDIFF": S0.EDIFF,
                                 "ISMEAR": 1, "SIGMA": 0.2, "ISYM": S0.ISYM,
                                 "NCORE": S0.NCORE, "PREC": "Normal",
                                 "potcar_functional": S0.POTCAR_FUNCTIONAL},
                      "a_scales": list(A_SCALES),
                      "targets": {}}

    if args.stage in ("relax", "spin"):
        for name, kind, targs, note in TARGETS:
            if args.stage == "spin" and name not in CR_TARGETS:
                continue
            st = endmember_structure(kind, targs)
            manifest["targets"][name] = {"note": note, "cell": cell_of(st),
                                         "kind": kind, "source": str(targs)}
            if args.stage == "relax":
                extra = dict(RELAX)
                if len(st) == 2 and abs(st.lattice.angles[0] - 90) > 1:
                    extra["ISIF"] = RELAX["ISIF"]   # ISIF=3 already frees c/a
                cases.append(write_case(
                    root / "relax" / name, st, f"s0d-{name}-relax",
                    S0.magmom_for(st), extra,
                    {"endpoint": name, "block": "stage0d-relax", "note": note,
                     "cell": cell_of(st), "role": "ISIF=3 relax to locate V0"}))
            else:
                # zero-moment control: if the Cr 3.0 and Cr 0.0 branches agree,
                # the FCC Cr state is non-magnetic and report 4.2b's occupancy
                # signal is not a magnetic artefact.
                cases.append(write_case(
                    root / "spin" / f"{name}_cr0", st, f"s0d-{name}-spin-cr0",
                    S0.magmom_for(st, {"zero": True}),
                    {"IBRION": -1, "NSW": 0, "ISIF": 2},
                    {"endpoint": name, "block": "stage0d-spin",
                     "magnetic_policy": "all zero (control)",
                     "note": note, "cell": cell_of(st),
                     "role": ("zero-moment control at the starting cell; "
                              "compare against the Cr=3.0 frozen table")}))

    else:  # eos
        if args.relax_root is None:
            raise SystemExit("--stage eos needs --relax-root")
        for name, kind, targs, note in TARGETS:
            c = args.relax_root / "stage0d" / "relax" / name / "CONTCAR"
            if not c.is_file():
                raise SystemExit(f"no relaxed cell at {c}; run --stage relax")
            base = Structure.from_file(str(c))
            # Does this cell have a SHAPE degree of freedom at fixed volume?
            # The first version tested `angles[0] != 90`, which is wrong in both
            # directions: the fcc primitive cell has 60 deg angles and the bcc
            # primitive cell 109.47 deg (both CUBIC crystals, shape fixed by
            # symmetry), while hcp has 90/90/120 (hexagonal, c/a genuinely
            # free).  So the test handed ISIF=4 to the two cells that did not
            # need it and ISIF=2 to the one that did.  The crystal system is
            # the right discriminator.
            from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
            need_shape = (SpacegroupAnalyzer(base, symprec=1e-3)
                          .get_crystal_system() != "cubic")
            is_hcp = need_shape
            manifest["targets"][name] = {"note": note, "relaxed_cell": cell_of(base),
                                         "is_hcp": is_hcp}
            for s in A_SCALES:
                st = scaled_cell(base, s)
                extra = ({"IBRION": 2, "ISIF": HCP_ISIF, "NSW": 80,
                          "POTIM": 0.2, "EDIFFG": -0.02} if is_hcp
                         else {"IBRION": -1, "NSW": 0, "ISIF": 2})
                cases.append(write_case(
                    root / "eos" / name / f"a{s:.2f}", st,
                    f"s0d-{name}-eos-a{s:.2f}", S0.magmom_for(st), extra,
                    {"endpoint": name, "block": "stage0d-eos", "note": note,
                     "a_scale": s,
                     "volume_scale": s ** 3,
                     "relaxed_cell": cell_of(base),
                     "is_hcp": is_hcp,
                     "role": ("constant-volume shape relax (ISIF=4) for HCP"
                              if is_hcp else "static single point")}))

    bad = [c for c in cases if not c["ok"]]
    report = {"schema": "legacy-qha-pool/stage0d/v1",
              "purpose": ("FCC/L1_2 and SER multi-volume static EOS -- the "
                          "2/3 of the pool that had no E(V) at all"),
              "n_cases": len(cases), "n_bad": len(bad),
              "manifest": manifest, "cases": cases}
    (root / args.stage).mkdir(parents=True, exist_ok=True)
    (root / args.stage / "validation.json").write_text(
        json.dumps(report, indent=2) + "\n")

    print(f"stage = {args.stage}   d={DENSITY}   EDIFF={S0.EDIFF:g}  ISYM={S0.ISYM}")
    print(f"\n{'endpoint':>16} {'n':>3} {'a/A':>9} {'V/atom':>9} {'mesh':>5} {'shift':>6}")
    print("-" * 56)
    seen = set()
    for c in cases:
        n = c.get("endpoint")
        if args.stage == "eos":
            key = (n, c["a_scale"])
        else:
            key = n
        if key in seen:
            continue
        seen.add(key)
        print(f"{n:>16} {c.get('cell', c.get('relaxed_cell', {})).get('n_atoms', 0):>3} "
              f"{c.get('cell', c.get('relaxed_cell', {})).get('a', float('nan')):>9.4f} "
              f"{c.get('cell', c.get('relaxed_cell', {})).get('volume_per_atom', float('nan')):>9.4f} "
              f"{c['kgrid']:>5} {c['shift_value']:>6}")
    print(f"\n{len(cases) - len(bad)}/{len(cases)} cases pass the input contract")
    for c in bad:
        print(f"  !! {c['dir']}")
        for p in c["problems"]:
            print(f"       {p}")
    if bad:
        raise SystemExit("refusing to submit: some inputs failed validation")


if __name__ == "__main__":
    main()
