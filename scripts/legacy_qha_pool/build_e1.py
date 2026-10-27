#!/usr/bin/env python
"""Build the E1 DFPT test inputs for BCC-Al-Ta-2 (legacy cell, verbatim).

E1 asks one question: does DFPT (VASP IBRION=8) and the finite-displacement
route that atomate2/PhononMaker uses give the same phonons for the same
structure? Everything except the method is held fixed against the legacy
run, so the difference that comes out is the method difference.

The cell is the legacy one, byte-for-byte: the relaxed R3/CONTCAR of
BCC-Al-Ta-2, which is also dfpt/QHA-04/POSCAR-unitcell. Holding the geometry
fixed is the point -- if we re-relaxed, we could not tell a method difference
from a geometry difference.

Writes only under --out. Never touches the NAS or the mirror.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

# ---- the fixed part of the recipe -----------------------------------------
A0 = 3.2246015525754168          # legacy relaxed B2 lattice constant, Angstrom
SPECIES = ["Al", "Ta"]
SUPER = 2                        # 2x2x2 supercell of the 2-atom B2 cell
KGRID = [5, 5, 5]                # legacy dfpt/QHA-04/KPOINTS
KSHIFT = [0, 0, 0]
POTCAR_DIR = Path("/nfs_ssd/documents/potpaw_PBE")
# Nb has no standard POTCAR in this library (only Nb_pv / Nb_sv); noted here
# because every Nb-containing endpoint in the pool has to make this choice.
POTCAR_VARIANT = {"Nb": "Nb_pv"}
MAGMOM = {"Al": 1.0, "Ta": 3.0}  # from main_qha.py's table


def potcar_path(el: str) -> Path:
    v = POTCAR_VARIANT.get(el, el)
    p = POTCAR_DIR / v / "POTCAR"
    if not p.is_file():
        raise SystemExit(f"no POTCAR for {el} (tried {p})")
    return p


def b2_supercell(a0: float, super_cell: int) -> tuple[list[list[float]], list[list[float]], list[str]]:
    """2-atom B2 cell -> (super x super x super) supercell, species-major order."""
    a = a0 * super_cell
    lat = [[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]]
    coords: list[list[float]] = []
    symbols: list[str] = []
    n = super_cell
    for el, base in zip(SPECIES, ([0.0, 0.0, 0.0], [0.5, 0.5, 0.5])):
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    coords.append([(base[0] + i) / n, (base[1] + j) / n, (base[2] + k) / n])
                    symbols.append(el)
    return lat, coords, symbols


def write_poscar(p: Path, lat, coords, symbols, title: str) -> None:
    counts = [sum(1 for s in symbols if s == el) for el in SPECIES]
    lines = [title, "1.0"]
    lines += ["  ".join(f"{x:.16f}" for x in row) for row in lat]
    lines += ["  ".join(SPECIES), "  ".join(str(c) for c in counts), "Direct"]
    lines += ["  ".join(f"{x:.16f}" for x in c) for c in coords]
    p.write_text("\n".join(lines) + "\n")


def write_incar(p: Path, extra: dict, title: str) -> None:
    base = {
        "SYSTEM": title,
        "ENCUT": 520,          # same as legacy and main_qha
        "ISTART": 0,
        "ICHARG": 2,           # current ht-vasp policy: atomic superposition SCF
        "LCHARG": False,       # ... and do not write CHGCAR
        "PREC": "Normal",      # legacy used Prec=N
        "ISMEAR": 1,
        "SIGMA": 0.2,
        "ISPIN": 2,
        "MAGMOM": None,        # filled per structure below
        "LREAL": False,
        "LWAVE": False,
        **extra,
    }
    out = []
    for k, v in base.items():
        if v is None:
            continue
        out.append(f"{k} = {v}")
    out.append("MAGMOM = " + " ".join(f"{MAGMOM[s]:.1f}" for s in EXTRA_SYMBOLS))
    p.write_text("\n".join(out) + "\n")


EXTRA_SYMBOLS: list[str] = []


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--method", choices=["dfpt", "fdisp"], default="dfpt")
    ap.add_argument("--super", type=int, default=SUPER,
                    help="supercell repeat (default 2, the legacy value)")
    ap.add_argument("--kmult", type=float, default=None,
                    help="if set, use a supercell k-mesh of round(kmult/super) so that "
                         "super*mesh stays roughly constant (absolute k-density held)")
    args = ap.parse_args()
    out = args.out
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"{out} is not empty; refusing to overwrite")
    out.mkdir(parents=True, exist_ok=True)

    supercell = args.super
    kgrid = list(KGRID)
    if args.kmult is not None:
        kgrid = [max(2, int(round(args.kmult / supercell)))] * 3
    lat, coords, symbols = b2_supercell(A0, supercell)
    global EXTRA_SYMBOLS
    EXTRA_SYMBOLS = symbols

    write_poscar(out / "POSCAR", lat, coords, symbols,
                 f"B2 AlTa {supercell}x{supercell}x{supercell}  a={A0:.10f}")

    (out / "KPOINTS").write_text(
        "Automatic generation\n0\nMonkhorst-pack\n"
        + " ".join(str(k) for k in kgrid) + "\n"
        + " ".join(str(s) for s in KSHIFT) + "\n")

    if args.method == "dfpt":
        # DFPT. ISYM is left at VASP's default (2) on purpose: symmetry is what
        # reduces 48 perturbations to 2, and ht-vasp's global ISYM=0 would
        # disable that and make this ~24x more expensive.
        extra = {
            "IBRION": 8,
            "NSW": 1,
            "EDIFF": 1e-5,
            "ISYM": 2,
            "ALGO": "Normal",   # project recipe for 16-atom cells (see memory)
            "NELM": 120,
            "NELMIN": 10,
            # DFPT on this cluster's VASP 6.4.3 needs NPAR = number of nodes so
            # that it can change the k-point set mid-run; NCORE=4 (which is what
            # legacy used on VASP 5.4.4) makes it refuse with
            # "VASP internal routines have requested a change of the k-point set".
            # KPAR=2 with no NCORE is the recipe the BCC-Ti-V pilot verified.
            "KPAR": 2,
        }
        title = "B2-AlTa-dfpt"
    else:
        # finite displacement, matching atomate2 PhononMaker's IBRION=6/NFREE=2
        extra = {
            "IBRION": 6,
            "NFREE": 2,
            "NSW": 0,
            "ISIF": 3,
            "EDIFF": 1e-7,      # PhononMaker's static energy setting
            "ISYM": 0,
            "ALGO": "Normal",
            "NELM": 200,
            "KPAR": 2,
        }
        title = "B2-AlTa-fdisp"
    write_incar(out / "INCAR", extra, title)

    with (out / "POTCAR").open("wb") as fh:
        for el in SPECIES:
            fh.write(potcar_path(el).read_bytes())

    recipe = {
        "schema": "legacy-qha-pool/e1-recipe/v1",
        "endpoint": "BCC-Al-Ta-2",
        "method": args.method,
        "cell": {"prototype": "B2 (CsCl)", "a0_angstrom": A0,
                 "supercell": f"{supercell}x{supercell}x{supercell}",
                 "n_atoms": len(coords),
                 "a_super_angstrom": A0 * supercell},
        "cell_provenance": ("legacy relaxed cell: BCC-Al-Ta-2/R3/CONTCAR == "
                            "dfpt/QHA-04/POSCAR-unitcell, a = 3.2246015525754168 A. "
                            "Held fixed so the DFPT-vs-finite-displacement difference "
                            "is not mixed with a geometry difference."),
        "kpoints": {"scheme": "Monkhorst-Pack", "grid": kgrid, "shift": KSHIFT,
                    "provenance": ("legacy dfpt/QHA-04/KPOINTS (5x5x5 on the 2x2x2 supercell); "
                                   "for other supercells the mesh is scaled so that "
                                   "supercell*mesh stays ~10, i.e. the absolute sampling "
                                   "density of the primitive Brillouin zone is held fixed")},
        "potcar": {el: str(potcar_path(el)) for el in SPECIES},
        "potcar_variant_note": ("the pool's standard potentials are Al/Cr/Hf/Ta/Ti/Zr "
                                "(no suffix). Nb has no standard variant in "
                                "/nfs_ssd/documents/potpaw_PBE, so every Nb endpoint "
                                "must use Nb_pv (ZVAL 11) or Nb_sv (ZVAL 13); Nb_pv is "
                                "chosen here and this is a pool-level decision."),
        "magmom": {s: MAGMOM[s] for s in SPECIES},
        "deviations_from_legacy": {
            "ALGO": "Normal here vs IALGO=48 (RMM-DIIS) in legacy; per this project's "
                    "verified recipe for 16-atom cells",
            "ICHARG": "2 here vs 1 in legacy (legacy had no CHGCAR, so it fell back to "
                      "the atomic superposition anyway)",
            "parallelisation": ("KPAR=2 and no NCORE. Legacy used NCORE=4 on VASP 5.4.4; "
                                "on this cluster's VASP 6.4.3 that makes DFPT abort with "
                                "'requested a change of the k-point set'. The BCC-Ti-V "
                                "pilot hit the same thing and settled on KPAR=2."),
            "LCHARG": "False here; legacy did not keep CHGCAR either",
        },
    }
    (out / "recipe.json").write_text(json.dumps(recipe, indent=2) + "\n")
    print(f"built {out}")
    for f in sorted(out.iterdir()):
        print(f"  {f.name}  {f.stat().st_size} B")
    print(f"  POTCAR sha256 = "
          f"{hashlib.sha256((out / 'POTCAR').read_bytes()).hexdigest()[:16]}...")


if __name__ == "__main__":
    main()
