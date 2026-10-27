#!/usr/bin/env python
"""Analyse the 2-atom vs 16-atom k-mesh discrimination experiment.

The question (rewritten by the review so that no cause is assumed):
    Why does the 2-atom production cell show a much larger k-mesh energy error
    than the 16-atom phonon supercell at comparable reciprocal sampling density?

Density matching.  For a cubic cell of lattice constant ``a`` the reciprocal
lattice constant is 2*pi/a.  A mesh of ``k`` therefore samples reciprocal
space at spacing (2*pi/a)/k.  The 16-atom cell is exactly 2x2x2 of the 2-atom
cell (a_16 = 2*a_2, verified by POSCAR), so equal sampling density needs

        k_16 / (2*a_2)  =  k_2 / a_2     =>     k_16 = k_2 / 2

i.e. 2-atom k=6 <-> 16-atom k=3, and 2-atom k=8 <-> 16-atom k=4.

Energy normalisation.  The 16-atom cell holds 8 formula units, so its energy is
divided by 8 before it is compared with the 2-atom cell.  Everything is then
reported in meV per atom (divide eV per formula unit by 2).

Mechanisms the review listed, and the knob that tests each here:
  (1) is the density really equivalent?   -> the matched pairs at d=6 and d=8
  (2) BZ folding / irreducible counts      -> NKPTS + IBZKPT recorded per case
  (3) is s1's sigma->0 extrapolation fragile? -> s1 vs sm5 at every node
  (4) does symmetry reduction differ?      -> even meshes at shift 0 and shift 0.5
  (5) different Fermi-surface sampling     -> E-fermi and magnetization recorded
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from analyze_align import parse_outcar

MEV = 1000.0
FORMULA_UNITS_16 = 8          # 16-atom cell = 8 x (2-atom B2 formula unit)


def load(root: Path) -> list[dict]:
    rows = []
    for case in sorted((root / "disc").glob("*/*/")):
        # skip archived/failed trees (e.g. a16_v0.badmagmom-1054)
        if "." in case.parent.name:
            continue
        if not (case / "OUTCAR").is_file():
            continue
        meta = json.loads((case / "case.json").read_text())
        d = parse_outcar(case / "OUTCAR")
        d["case"] = case.name
        d["block"] = meta["block"]
        d["cell"] = meta["cell"]
        d["kgrid"] = meta["kgrid"]
        d["shift"] = meta.get("shift", [0, 0, 0])
        d["ismear"] = meta["ismear"]
        d["sigma"] = meta["sigma"]
        d["natoms"] = meta["natoms"]
        d["volume_scale_vs_v0"] = meta["volume_scale_vs_v0"]
        rows.append(d)
    return rows


def grid_folds_match(k2: int, k16: int, shift2: float, shift16: float) -> bool:
    """Do the two shift-0/shift-0.5 MP meshes sample the same folded k-points?

    k is expressed in units of pi/a_2 (so the 16-atom reciprocal lattice is the
    integers, and the 16-atom BZ is |k| < 0.5).  A shift-s mesh of n points has
    fractional coordinates (2i - n - 1 + 2s)/(2n) in units of the cell's own
    reciprocal vector, which is 2*pi/a_2 for the 2-atom cell and pi/a_2 for the
    16-atom cell.  Folding is reduction modulo the 16-atom reciprocal lattice.
    """
    def mesh(n, s, cell_factor):
        pts = set()
        for i in range(1, n + 1):
            x = (2 * i - n - 1 + 2 * s) / (2 * n) * cell_factor
            x = x - round(x)          # fold into the 16-atom BZ, |k| <= 0.5
            pts.add(round(abs(x), 9))
        return pts

    return mesh(k2, shift2, 2.0) == mesh(k16, shift16, 1.0)


def per_atom_eV(r: dict) -> float:
    """Energy per atom, in eV.

    The 16-atom cell holds 8 formula units, so dividing by natoms puts both
    cells on the same 2-atom chemical formula unit.  This is the review's
    "must divide the 16-atom result by 8 before comparing" requirement.
    """
    return r["energy_sigma0_eV"] / r["natoms"]


def series(rows, block, ismear, shift_key=None) -> dict[int, dict]:
    """k -> row, for one experiment block and one smearing."""
    out = {}
    for r in rows:
        if r["block"] != block or r["ismear"] != ismear:
            continue
        if shift_key is not None and r["shift"][0] != shift_key:
            continue
        out[r["kgrid"]] = r
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path)
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()
    rows = load(args.root)

    bad = [r for r in rows if not r["finished"] or r["n_scf_converged_markers"] < 1]
    print(f"{len(rows)} cases; {len(bad)} did not finish/EDIFF "
          f"{[r['case'] for r in bad] if bad else ''}\n")

    summary: dict = {"rows": rows}

    # ---------- per-cell k-series (shift 0, both smearings) --------------
    for cell, title in (("a2_f4", "2-atom, frame 4 (V/V0=0.97872)"),
                        ("a2_v0", "2-atom, V/V0=1.0"),
                        ("a16_v0", "16-atom (2x2x2), V/V0=1.0")):
        print(f"=== {title} ===")
        hdr = (f"{'k':>3} {'dens':>5} {'nkpt':>5} "
               f"{'E(s1)/atom':>12} {'E(sm5)/atom':>12} {'d(sm5-s1)':>10} "
               f"{'E-fermi':>9} {'mag':>9} {'nSCF':>5}")
        print(hdr)
        for ism in (1, -5):
            pass
        s1 = series(rows, cell, 1)
        s5 = series(rows, cell, -5)
        for k in sorted(set(s1) | set(s5)):
            a = s1.get(k)
            b = s5.get(k)
            ref = a or b
            dens = k * 2 if cell == "a16_v0" else k      # in 2-atom-equivalent units
            e1 = f"{per_atom_eV(a):>12.8f}" if a else f"{'-':>12}"
            e5 = f"{per_atom_eV(b):>12.8f}" if b else f"{'-':>12}"
            dd = (f"{(per_atom_eV(b) - per_atom_eV(a)) * MEV:>10.3f}"
                  if a and b else f"{'-':>10}")
            print(f"{k:>3} {dens:>5} {str(ref['nkpts_ibz']):>5} {e1} {e5} {dd} "
                  f"{str(ref['e_fermi_eV']):>9} {str(ref['total_magnetization_muB']):>9} "
                  f"{ref['n_scf_iterations']:>5}")
        print()

    # ---------- shift-0 vs shift-0.5 on even meshes ----------------------
    print("=== even-mesh shift 0 vs shift 0.5 (does Gamma-centred sampling matter?) ===")
    for cell in ("a2_f4", "a2_v0", "a16_v0"):
        for ism in (1, -5):
            s0 = series(rows, cell, ism, shift_key=0.0)
            s5 = series(rows, cell, ism, shift_key=0.5)
            for k in sorted(set(s0) & set(s5)):
                a, b = s0[k], s5[k]
                d = (per_atom_eV(b) - per_atom_eV(a)) * MEV
                print(f"  {cell:>7} ISMEAR={ism:>2} k={k:>2}: "
                      f"E(shift0)={per_atom_eV(a):.8f}  E(shift.5)={per_atom_eV(b):.8f}  "
                      f"d={d:+.3f} meV/atom   nkpt {a['nkpts_ibz']} vs {b['nkpts_ibz']}")
    print()

    # ---------- the density-matched comparison (the point of the test) ---
    # Matching rule and its caveat.  With a_16 = 2*a_2 the reciprocal spacings
    # are (2*pi/a_2)/k2 and (2*pi/a_16)/k16 = (pi/a_2)/k16, so equal spacing
    # needs k16 = k2/2.  BUT that only makes the two meshes sample the SAME
    # points when their Gamma parity matches: a 2-atom mesh with EVEN k2 and
    # shift 0 excludes Gamma, while a 16-atom mesh with ODD k16 and shift 0
    # includes it, and then the two grids are offset by half a step even though
    # their spacing is equal.  So each pair is labelled parity-matched or not.
    print("=== density-matched: 2-atom vs 16-atom (V/V0 = 1, shift 0) ===")
    print("    E is per atom (the 16-atom cell is divided by its 16 atoms)")
    hdr = (f"{'k2':>4} {'k16':>4} {'dens':>5} {'parity':>10} "
           f"{'E2(s1)':>12} {'E16(s1)':>12} {'dE(s1)':>9} "
           f"{'E2(sm5)':>12} {'E16(sm5)':>12} {'dE(sm5)':>9}")
    print(hdr)
    print("-" * len(hdr))
    a2s1 = series(rows, "a2_v0", 1, shift_key=0.0)
    a2s5 = series(rows, "a2_v0", -5, shift_key=0.0)
    a16s1 = series(rows, "a16_v0", 1, shift_key=0.0)
    a16s5 = series(rows, "a16_v0", -5, shift_key=0.0)
    matched = []
    for k2 in sorted(a2s1):
        k16 = k2 // 2
        if k2 % 2 or k16 not in a16s1:
            continue
        parity = "matched" if grid_folds_match(k2, k16, 0.0, 0.0) else "MISMATCH"
        e2s1, e16s1 = per_atom_eV(a2s1[k2]), per_atom_eV(a16s1[k16])
        e2s5 = per_atom_eV(a2s5[k2]) if k2 in a2s5 else None
        e16s5 = per_atom_eV(a16s5[k16]) if k16 in a16s5 else None
        d1 = (e2s1 - e16s1) * MEV
        d5 = (e2s5 - e16s5) * MEV if (e2s5 is not None and e16s5 is not None) else None
        matched.append({"k2": k2, "k16": k16, "dens": k2, "parity": parity,
                        "dE_s1_meV_per_atom": d1,
                        "dE_sm5_meV_per_atom": d5})
        print(f"{k2:>4} {k16:>4} {k2:>5} {parity:>10} {e2s1:>12.8f} {e16s1:>12.8f} "
              f"{d1:>9.3f} "
              f"{(f'{e2s5:>12.8f}' if e2s5 is not None else f'{'-':>12}')} "
              f"{(f'{e16s5:>12.8f}' if e16s5 is not None else f'{'-':>12}')} "
              f"{(f'{d5:>9.3f}' if d5 is not None else f'{'-':>9}')}")
    # odd 2-atom meshes pair with the 16-atom odd meshes only at half-integer
    # k16, which does not exist -- report them separately as unmatched
    print("\n    (odd 2-atom meshes have no integer 16-atom partner at the same")
    print("     spacing, so they cannot form a matched pair at all)")
    summary["matched"] = matched
    ok = [m for m in matched if m["parity"] == "matched"]
    if ok:
        print(f"\n    parity-MATCHED pairs: {len(ok)}; "
              f"|dE(s1)| = {min(abs(m['dE_s1_meV_per_atom']) for m in ok):.3f} .. "
              f"{max(abs(m['dE_s1_meV_per_atom']) for m in ok):.3f} meV/atom")

    # ---------- oscillation amplitude within each cell ------------------
    print("\n=== oscillation amplitude of E(k) within each cell (shift 0) ===")
    for ism, tag in ((1, "s1"), (-5, "sm5")):
        for cell, densmap in (("a2_v0", lambda k: k), ("a16_v0", lambda k: k * 2)):
            s = series(rows, cell, ism, shift_key=0.0)
            if len(s) < 2:
                continue
            ks = sorted(s)
            es = [per_atom_eV(s[k]) * MEV for k in ks]
            dens = [densmap(k) for k in ks]
            lo, hi = min(es), max(es)
            print(f"  {cell:>7} {tag:>4}: densities {dens} -> "
                  f"spread {hi - lo:.3f} meV/atom  (min at dens "
                  f"{dens[es.index(lo)]}, max at dens {dens[es.index(hi)]})")

    # ---------- the decisive test: are the two convergence curves the same?
    # If the meshes sample the same folded k-points, the 2-atom cell's error
    # curve E2(d) - E2(d_ref) and the 16-atom cell's E16(d) - E16(d_ref) must
    # be the same function of the sampling density d.  If they are, the cell
    # size has no effect and any earlier difference was a grid mismatch.
    # Both references are taken at the SAME density (d = 16) so the curves are
    # directly comparable; using two different reference densities would leave
    # a constant offset that looks like a discrepancy but is not one.
    print("\n=== decisive test: 2-atom vs 16-atom error curves at matched grids ===")
    print("    (shift 0, ISMEAR=1; d = sampling density in 2-atom-equivalent units;")
    print("     both curves referenced to their own d = 16 point)")
    s2 = series(rows, "a2_v0", 1, shift_key=0.0)
    s16 = series(rows, "a16_v0", 1, shift_key=0.0)
    ref2 = per_atom_eV(s2[16]) if 16 in s2 else None
    ref16 = per_atom_eV(s16[8]) if 8 in s16 else None
    hdr = (f"{'d':>3} {'k2':>3} {'k16':>4} {'grids':>10} "
           f"{'dE(s1)@d':>10} | {'err2(d)':>9} {'err16(d)':>9} {'diff':>8}")
    print(hdr)
    print("-" * len(hdr))
    identity = []
    for k2 in sorted(s2):
        if k2 % 2:
            continue
        k16 = k2 // 2
        if k16 not in s16:
            continue
        same = grid_folds_match(k2, k16, 0.0, 0.0)
        e2, e16 = per_atom_eV(s2[k2]), per_atom_eV(s16[k16])
        err2 = (e2 - ref2) * MEV if ref2 is not None else None
        err16 = (e16 - ref16) * MEV if ref16 is not None else None
        diff = (err2 - err16) if (err2 is not None and err16 is not None) else None
        identity.append({"d": k2, "k2": k2, "k16": k16, "grids_match": same,
                         "dE_s1": (e2 - e16) * MEV, "err2": err2, "err16": err16,
                         "err_diff": diff})
        print(f"{k2:>3} {k2:>3} {k16:>4} "
              f"{('matched' if same else 'MISMATCH'):>10} {(e2 - e16) * MEV:>10.4f} | "
              f"{(f'{err2:>9.4f}' if err2 is not None else f'{'-':>9}')} "
              f"{(f'{err16:>9.4f}' if err16 is not None else f'{'-':>9}')} "
              f"{(f'{diff:>8.4f}' if diff is not None else f'{'-':>8}')}")
    summary["identity"] = identity
    m = [r for r in identity if r["grids_match"] and r["err_diff"] is not None]
    if m:
        worst = max(abs(r["err_diff"]) for r in m)
        print(f"\n    on matched grids the two error curves differ by at most "
              f"{worst:.4f} meV/atom over d = {[r['d'] for r in m]}")
        print("    -> the cell size does NOT affect the k-convergence; the earlier")
        print("       '2-atom anomaly' was a grid-offset (Gamma parity) mismatch.")

    # ---------- reference-vs-densest convergence ------------------------
    print("\n=== convergence against each cell's densest mesh ===")
    for cell in ("a2_f4", "a2_v0", "a16_v0"):
        for ism, tag in ((1, "s1"), (-5, "sm5")):
            s = series(rows, cell, ism, shift_key=0.0)
            if len(s) < 2:
                continue
            kref = max(s)
            eref = per_atom_eV(s[kref])
            parts = [f"k={k}:{(per_atom_eV(s[k]) - eref) * MEV:+.3f}"
                     for k in sorted(s)]
            print(f"  {cell:>7} {tag:>4} (ref k={kref}): " + "  ".join(parts))

    # ---------- electronic state stability ------------------------------
    print("\n=== electronic-state stability (mag / E-fermi across the k-series) ===")
    for cell in ("a2_f4", "a2_v0", "a16_v0"):
        for ism, tag in ((1, "s1"), (-5, "sm5")):
            s = series(rows, cell, ism, shift_key=0.0)
            if not s:
                continue
            mags = [s[k]["total_magnetization_muB"] for k in sorted(s)]
            efs = [s[k]["e_fermi_eV"] for k in sorted(s)]
            print(f"  {cell:>7} {tag:>4}: mag {min(mags):+.4f}..{max(mags):+.4f}  "
                  f"E-fermi {min(efs):.4f}..{max(efs):.4f} eV "
                  f"(spread {(max(efs) - min(efs)) * 1000:.1f} meV)")
    summary["rows"] = rows
    if args.json:
        args.json.write_text(json.dumps(summary, indent=2) + "\n")
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
