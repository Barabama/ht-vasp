#!/usr/bin/env python
"""Stage 0A: static electronic recipe calibration on the structural classes.

Six representatives, one per class the review named:

    BCC-Ta-Ta-2    stable non-magnetic BCC/B2
    BCC-Al-Ta-2    confirmed finite-q soft mode (the control; §4.6/§4.7)
    BCC-Al-Cr-2    Cr-bearing, magnetic status to be resolved
    BCC-Hf-V-2     pool-coverage (Hf, V) -- different BCC electronic structure
    FCC-Zr-Al-4    FCC/L1_2  (== Al3Zr: Zr on 1a, Al on 3c, per Endmember.get_fcc)
    SER-Al-1       SER reference (fcc primitive, 1 atom)

Nb is deliberately absent: the Nb POTCAR decision is still open (D0-C).

KPOINTS family -- this is the part the earlier tests got wrong twice
------------------------------------------------------------------
A Monkhorst-Pack mesh with shift 0 includes Gamma only when every n_i is ODD.
For EVEN n, Gamma-centring requires shift 0.5.  So the review's phrase
"Gamma-centred, shift 0 0 0" and its density values "k=12, k=16" cannot both
hold: k=12 with shift 0 is not Gamma-centred.

There is a further constraint from the supercell.  Two meshes of a cell and its
2x2x2 supercell sample the SAME folded k-points only when

    k_super = k_cell / 2   AND   the parities line up.

Enumerating that (analyze_discriminate.grid_folds_match) gives exactly one
Gamma-centred family that also fold-matches:

    cell (even k, shift 0.5)  <->  supercell (k/2, shift 0 if odd else 0.5)

verified empirically: at V/V0 = 1, ISMEAR=1,
    2-atom k=6 s0.5 vs 16-atom k=3 s0   ->  +0.011 meV/atom
    2-atom k=8 s0.5 vs 16-atom k=4 s0.5 ->  +0.011 meV/atom
(ISMEAR=-5 does NOT agree on these pairs -- see the report; one more reason to
keep -5 out of the production path.)

So Stage 0A uses:  **even k, shift 0.5 0.5 0.5  (Gamma-centred)**, at two
densities taken from the review (d = 12 and d = 16 in 2-atom-AlTa-equivalent
units), converted to each cell by equal reciprocal spacing

    n_cell = round_to_even( 2*pi / (dk_target * a_cell) )

with dk_target = 2*pi/(a_AlTa * d).  The rounding is explicit and recorded --
unlike automatic_density, which scales and rounds silently.

Structures
----------
Cells are built here, not copied from legacy (legacy is not a trusted source,
see report §3.2).  Lattice constants come from Materials Project elemental
volumes (queried 2026-10-11) with volume-additive (Vegard) mixing for the
compounds.  These are CALIBRATION cells, not production cells: production will
relax its own.  A volume error of a few percent changes a by ~1%, i.e. it moves
the effective density far less than the d=12 -> d=16 step being measured.

Provenance of every lattice constant is written into validation.json.

Writes only under --out.  Refuses to overwrite an existing INCAR.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp import Poscar
from pymatgen.io.vasp.inputs import Potcar

from input_contract import validate_dir

# ---------------------------------------------------------------------------
# recipe
# ---------------------------------------------------------------------------
POTCAR_FUNCTIONAL = "PBE"          # PBE(52), POT_GGA_PAW_PBE -- legacy-comparable
ENCUT = 520
EDIFF = 1e-6
ISYM = 2
NCORE = 8
NELM = 120
NELMIN = 6
ALGO = "Normal"
LORBIT = 11                        # output-only: needed for per-site moments

# main_qha.py's frozen 16-element table (user decision)
MAGMOM = {"Al": 1.0, "Cr": 3.0, "Hf": 1.0, "Nb": 3.0, "Ta": 3.0,
          "Ti": 3.0, "Zr": 1.0, "V": 3.0, "Mo": 3.0, "Re": 1.0}

# reference cell for the density definition
A_ALTA = 3.2246015525754168        # legacy-relaxed B2 AlTa, Angstrom
DENSITIES = (12, 16)               # in 2-atom-AlTa-equivalent k units

SMEARINGS = {"s1": (1, 0.2), "s2": (2, 0.2), "sm5": (-5, None)}

# ---------------------------------------------------------------------------
# elemental volumes per atom (Angstrom^3), Materials Project, 2026-10-11
# ---------------------------------------------------------------------------
MP_VOLUME_PER_ATOM = {
    "Al_bcc": 16.116,   # mp-998860  (metastable BCC polymorph, Im-3m)
    "Al_fcc": 16.472,   # mp-134
    "Cr_bcc": 13.084,   # mp-90
    "Hf_bcc": 22.087,   # mp-100     (metastable BCC polymorph)
    "Ta_bcc": 18.130,   # mp-50
    "V_bcc": 13.264,    # mp-146
    "Zr_hcp": 23.500,   # mp-131
}


def cubic_a_from_vpa(vpa: float, atoms_per_cell: int) -> float:
    """Lattice constant of a cubic cell with `atoms_per_cell` at volume vpa each."""
    return (vpa * atoms_per_cell) ** (1.0 / 3.0)


# ---------------------------------------------------------------------------
# structures
# ---------------------------------------------------------------------------
def cell_bcc_2atoms(a: float, elem: str) -> Structure:
    """Conventional cubic BCC cell, 2 atoms of one element."""
    return Structure(Lattice.cubic(a), [elem, elem], [[0, 0, 0], [0.5, 0.5, 0.5]])


def cell_b2(a: float, elem1: str, elem2: str) -> Structure:
    """B2 (CsCl) conventional cubic cell: elem1 at 1a, elem2 at 1b."""
    return Structure(Lattice.cubic(a), [elem1, elem2], [[0, 0, 0], [0.5, 0.5, 0.5]])


def cell_l12(a: float, corner: str, face: str) -> Structure:
    """L1_2 conventional cubic cell: `corner` on 1a, `face` on the three 3c sites.

    Matches Endmember.get_fcc, whose template is Au_1a-Cu_3c (mp-2258), so
    FCC-Zr-Al-4 is Zr on 1a and Al on 3c, i.e. Al3Zr.
    """
    return Structure(
        Lattice.cubic(a),
        [corner, face, face, face],
        [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
    )


def cell_fcc_primitive(a_conv: float, elem: str) -> Structure:
    """1-atom primitive fcc cell (rhombohedral, alpha = 60 deg)."""
    a = a_conv / math.sqrt(2.0)
    return Structure(Lattice.from_parameters(a, a, a, 60, 60, 60),
                     [elem], [[0, 0, 0]])


def build_structures() -> dict[str, dict]:
    v = MP_VOLUME_PER_ATOM
    a_b2_alta = A_ALTA
    a_alcr = cubic_a_from_vpa((v["Al_bcc"] + v["Cr_bcc"]) / 2, 2)
    a_hfv = cubic_a_from_vpa((v["Hf_bcc"] + v["V_bcc"]) / 2, 2)
    a_zral3 = cubic_a_from_vpa((v["Zr_hcp"] + 3 * v["Al_fcc"]) / 4, 4)
    a_ta = cubic_a_from_vpa(v["Ta_bcc"], 2)
    a_al_fcc = cubic_a_from_vpa(v["Al_fcc"], 4)
    return {
        "BCC-Ta-Ta-2": {
            "struct": cell_bcc_2atoms(a_ta, "Ta"),
            "a": a_ta,
            "source": f"MP mp-50 BCC Ta, V/atom={v['Ta_bcc']} A^3 -> a=(2V)^(1/3)",
            "class": "stable non-magnetic BCC",
        },
        "BCC-Al-Ta-2": {
            "struct": cell_b2(a_b2_alta, "Al", "Ta"),
            "a": a_b2_alta,
            "source": ("legacy-relaxed B2 AlTa (continuity with report 4.6/4.7); "
                       "DFT-PBE relaxed with the same POTCAR family and ENCUT"),
            "class": "finite-q soft mode (control)",
        },
        "BCC-Al-Cr-2": {
            "struct": cell_b2(a_alcr, "Al", "Cr"),
            "a": a_alcr,
            "source": (f"Vegard: BCC Al mp-998860 {v['Al_bcc']} + BCC Cr mp-90 "
                       f"{v['Cr_bcc']} -> V/atom {(v['Al_bcc'] + v['Cr_bcc']) / 2:.3f}"),
            "class": "Cr-bearing BCC (magnetic status to resolve)",
        },
        "BCC-Hf-V-2": {
            "struct": cell_b2(a_hfv, "Hf", "V"),
            "a": a_hfv,
            "source": (f"Vegard: BCC Hf mp-100 {v['Hf_bcc']} + BCC V mp-146 "
                       f"{v['V_bcc']} -> V/atom {(v['Hf_bcc'] + v['V_bcc']) / 2:.3f}"),
            "class": "pool coverage (Hf, V)",
        },
        "FCC-Zr-Al-4": {
            "struct": cell_l12(a_zral3, "Zr", "Al"),
            "a": a_zral3,
            "source": (f"L1_2 ZrAl3 == Al3Zr (Zr 1a, Al 3c). Vegard: 1x HCP Zr mp-131 "
                       f"{v['Zr_hcp']} + 3x FCC Al mp-134 {v['Al_fcc']} -> V/atom "
                       f"{(v['Zr_hcp'] + 3 * v['Al_fcc']) / 4:.3f}"),
            "class": "FCC/L1_2",
        },
        "SER-Al-1": {
            "struct": cell_fcc_primitive(a_al_fcc, "Al"),
            "a": a_al_fcc / math.sqrt(2.0),
            "source": (f"MP mp-134 FCC Al, a_conv=(4*{v['Al_fcc']})^(1/3)={a_al_fcc:.4f}; "
                       f"1-atom primitive a=a_conv/sqrt2"),
            "class": "SER reference",
        },
    }


# ---------------------------------------------------------------------------
# k mesh: explicit, Gamma-centred (even n + shift 0.5), density-matched
# ---------------------------------------------------------------------------
def round_even(x: float) -> int:
    n = int(round(x / 2.0)) * 2
    return max(n, 2)


def mesh_for(volume_angstrom3: float, density: int) -> int:
    """Even n whose reciprocal sampling density matches `density` on the AlTa reference.

    The invariant is the number of k-points per unit reciprocal volume,

        n^3 * V_cell = const        (BZ volume is (2*pi)^3 / V_cell)

    so n = (a_AlTa * density) / V_cell^(1/3).  For a cubic cell V^(1/3) == a and
    this reduces to n = a_AlTa*density/a, but for a NON-cubic cell -- e.g. the
    1-atom primitive fcc SER cell -- using a instead of V^(1/3) under-samples it
    by 4^(1/3) = 1.59x.  That is why this takes the volume.
    """
    return round_even((A_ALTA * density) / volume_angstrom3 ** (1.0 / 3.0))


def kpoints_text(n: int, shift: float) -> str:
    return (f"Automatic kpoint scheme\n0\nMonkhorst-Pack\n{n} {n} {n}\n"
            f"{shift} {shift} {shift}\n")


# ---------------------------------------------------------------------------
# INCAR
# ---------------------------------------------------------------------------
def incar_text(system: str, ismear: int, sigma: float | None,
               magmom: str, ispin: int = 2) -> str:
    lines = [
        f"SYSTEM = {system}",
        f"ENCUT = {ENCUT}",
        "ISTART = 0",
        "ICHARG = 2",
        "LCHARG = False",
        "PREC = Normal",
        f"ISMEAR = {ismear}",
    ]
    if sigma is not None:
        lines.append(f"SIGMA = {sigma}")
    lines += [
        f"ISPIN = {ispin}",
        magmom,
        "LREAL = False",
        "LWAVE = False",
        "IBRION = -1",
        "NSW = 0",
        f"EDIFF = {EDIFF:g}",
        f"ISYM = {ISYM}",
        f"ALGO = {ALGO}",
        f"NELM = {NELM}",
        f"NELMIN = {NELMIN}",
        f"NCORE = {NCORE}",
        f"LORBIT = {LORBIT}",
    ]
    return "\n".join(lines) + "\n"


def magmom_for(struct: Structure, override: dict | None = None) -> str:
    vals = []
    for s in struct:
        el = str(s.specie)
        vals.append(0.0 if override is not None and override.get("zero") else
                    (override or {}).get(el, MAGMOM.get(el, 1.0)))
    return "MAGMOM = " + " ".join(f"{v:.1f}" for v in vals)


# ---------------------------------------------------------------------------
# writing
# ---------------------------------------------------------------------------
def write_case(d: Path, struct: Structure, n: int, shift: float, system: str,
               ismear: int, sigma: float | None, magmom: str, ispin: int,
               meta: dict) -> dict:
    if (d / "INCAR").exists():
        raise SystemExit(f"refusing to overwrite existing case: {d}")
    d.mkdir(parents=True, exist_ok=True)
    Poscar(struct).write_file(str(d / "POSCAR"))
    (d / "INCAR").write_text(incar_text(system, ismear, sigma, magmom, ispin))
    (d / "KPOINTS").write_text(kpoints_text(n, shift))
    species_order = list(dict.fromkeys(str(s.specie) for s in struct))
    Potcar(species_order, functional=POTCAR_FUNCTIONAL).write_file(str(d / "POTCAR"))
    v = validate_dir(d)
    v["kgrid"] = n
    v["shift_value"] = shift
    v.update(meta)
    (d / "case.json").write_text(json.dumps(v, indent=2) + "\n")
    return v


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)

    reps = build_structures()
    cases: list[dict] = []

    for name, spec in reps.items():
        st = spec["struct"]
        for density in DENSITIES:
            n = mesh_for(float(st.volume), density)
            dk = 2.0 * math.pi / (spec["a"] * n)
            meta = {
                "endpoint": name, "class": spec["class"],
                "a_angstrom": spec["a"], "a_source": spec["source"],
                "density_target": density, "dk_angstrom_inv": dk,
                "volume_angstrom3": float(st.volume),
            }
            for tag, (ism, sig) in SMEARINGS.items():
                cases.append(write_case(
                    out / name / f"d{density:02d}_k{n:02d}_{tag}", st, n, 0.5,
                    f"stage0-{name}-d{density}-{tag}", ism, sig,
                    magmom_for(st), 2, meta))

    # --- AlTa family probe: complete the Gamma-CENTRED odd series (shift 0) ---
    # so the even/shift-0.5 family used above can be compared against the
    # odd/shift-0 family on the SAME cell.
    st = reps["BCC-Al-Ta-2"]["struct"]
    for n in (11, 13, 15, 17):
        dk = 2.0 * math.pi / (A_ALTA * n)
        cases.append(write_case(
            out / "BCC-Al-Ta-2" / f"probe_odd_k{n:02d}_s0_s1", st, n, 0.0,
            f"stage0-AlTa-odd-k{n}", 1, 0.2, magmom_for(st), 2,
            {"endpoint": "BCC-Al-Ta-2", "class": "family probe",
             "a_angstrom": A_ALTA, "a_source": "same cell as above",
             "density_target": None, "dk_angstrom_inv": dk,
             "volume_angstrom3": float(st.volume), "family": "odd/shift0/Gamma-centred"}))

    # --- magnetic screening for the Cr-bearing representative ----------------
    st = reps["BCC-Al-Cr-2"]["struct"]
    n = mesh_for(float(st.volume), 16)
    base = {"endpoint": "BCC-Al-Cr-2", "class": "magnetic screening",
            "a_angstrom": reps["BCC-Al-Cr-2"]["a"],
            "a_source": reps["BCC-Al-Cr-2"]["source"],
            "density_target": 16, "dk_angstrom_inv": 2 * math.pi / (reps["BCC-Al-Cr-2"]["a"] * n),
            "volume_angstrom3": float(st.volume)}
    for tag, ispin, mm, desc in (
            ("ispin1", 1, "MAGMOM = 0.0 0.0", "ISPIN=1 (no spin)"),
            ("ispin2_zero", 2, "MAGMOM = 0.0 0.0", "ISPIN=2, zero initial moment"),
            ("ispin2_cr", 2, magmom_for(st), "ISPIN=2, Al 1.0 / Cr 3.0")):
        cases.append(write_case(
            out / "BCC-Al-Cr-2" / f"mag_{tag}_k{n:02d}", st, n, 0.5,
            f"stage0-AlCr-{tag}", 1, 0.2, mm, ispin,
            {**base, "magmom_policy": desc}))

    bad = [c for c in cases if not c["ok"]]
    report = {
        "schema": "legacy-qha-pool/stage0/v1",
        "purpose": "Stage 0A static electronic recipe calibration",
        "kpoints_family": "even n, shift 0.5 0.5 0.5 (Gamma-centred, fold-matchable)",
        "density_definition": ("n_cell = round_even(2*pi/(dk*a)); "
                               f"dk = 2*pi/({A_ALTA}*d), d in {DENSITIES}"),
        "reference_cell": {"name": "B2 AlTa", "a_angstrom": A_ALTA},
        "recipe": {"ENCUT": ENCUT, "EDIFF": EDIFF, "ISYM": ISYM, "NCORE": NCORE,
                   "NELM": NELM, "NELMIN": NELMIN, "ALGO": ALGO,
                   "LORBIT": LORBIT,
                   "potcar_functional": POTCAR_FUNCTIONAL,
                   "magmom_table": MAGMOM,
                   "note_LORBIT": "output-only diagnostics flag; does not change "
                                  "energies or forces"},
        "n_cases": len(cases),
        "n_bad": len(bad),
        "cases": cases,
    }
    (out / "validation.json").write_text(json.dumps(report, indent=2) + "\n")

    print(f"{'endpoint':<16} {'class':<34} {'a/A':>8} {'d':>3} {'k':>4} "
          f"{'dk/A^-1':>9} {'natoms':>6}  cases")
    print("-" * 100)
    seen: dict[tuple, int] = {}
    for c in cases:
        key = (c.get("endpoint"), c.get("density_target"))
        seen[key] = seen.get(key, 0) + 1
    for c in cases:
        if not Path(c["dir"]).name.endswith("_s1"):
            continue
        print(f"{c.get('endpoint', '?'):<16} {c.get('class', ''):<34} "
              f"{c.get('a_angstrom', 0):>8.4f} {str(c.get('density_target')):>3} "
              f"{c['kgrid']:>4} {c.get('dk_angstrom_inv', 0):>9.5f} {c['natoms']:>6}  "
              f"{seen.get((c.get('endpoint'), c.get('density_target')), 0)}")
    print(f"\n{len(cases) - len(bad)}/{len(cases)} cases pass the input contract")
    for c in bad:
        print(f"  !! {c['dir']}")
        for p in c["problems"]:
            print(f"       {p}")
    if bad:
        raise SystemExit("refusing to submit: some inputs failed validation")


if __name__ == "__main__":
    main()
