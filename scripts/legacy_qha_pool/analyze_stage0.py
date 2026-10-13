#!/usr/bin/env python
"""Analyse Stage 0A: is the static electronic recipe transferable across classes?

Four questions, in the order the review asked them:

  Q1  Is d = 12 (the cheaper density) enough, or is d = 16 required?
      -> |E(d12) - E(d16)| per representative, per smearing.
  Q2  Does the smearing choice matter at the converged density?
      -> dE(sm5 - s1) and dE(s2 - s1) at d = 16.
  Q3  Is the Cr-bearing representative magnetic, non-magnetic, or unresolved?
      -> ISPIN=1 vs ISPIN=2 (zero start) vs ISPIN=2 (Cr 3.0 start).
  Q4  Which Gamma-centred family should production use -- odd n with shift 0,
      or even n with shift 0.5?  (These are the only two Gamma-centred options;
      even n with shift 0 is NOT Gamma-centred, which is the specification
      error the earlier rounds carried.)

Everything is reported per atom.  ``energy(sigma->0)`` is the metric, because
that is what pymatgen ``Vasprun.final_energy`` returns and therefore what the
production QHA consumes (report §4.6.1).
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from analyze_align import parse_outcar

MEV = 1000.0


def parse_site_moments(outcar: Path) -> list[float] | None:
    """Per-site magnetization from the LORBIT=11 table in OUTCAR.

    Layout (VASP 6.4.3)::

         magnetization (x)

        # of ion       s       p       d       tot
        ------------------------------------------
            1       -0.026  -0.020   0.000  -0.046
            2        0.029   0.047   2.579   2.655
        --------------------------------------------------
        tot          0.003   0.027   2.579   2.609

    Parsing is row-driven, not delimiter-driven: a previous version terminated
    the block on the next run of dashes, which for a single-atom cell ran past
    the table into the timing block and reported "CHARGE: cpu time ..." lines as
    site moments.  Rows must look like ``<int> <float> <float> <float> <float>``
    and the scan stops at the first line that does not.
    """
    txt = outcar.read_text(errors="replace")
    blocks = re.findall(r"magnetization \(x\).*?# of ion.*?\n-+\n", txt, re.S)
    if not blocks:
        return None
    # start scanning from the last header
    last = txt.rfind("magnetization (x)")
    if last < 0:
        return None
    tail = txt[last:].splitlines()
    vals: list[float] = []
    started = False
    for ln in tail:
        parts = ln.split()
        if len(parts) == 5 and parts[0].isdigit():
            try:
                vals.append(float(parts[-1]))
                started = True
                continue
            except ValueError:
                pass
        if started:
            break                      # table ended
    return vals or None


def load(root: Path) -> list[dict]:
    rows = []
    for case in sorted(root.rglob("case.json")):
        d = case.parent
        if d.name.startswith("."):
            continue
        if not (d / "OUTCAR").is_file():
            continue
        meta = json.loads(case.read_text())
        o = parse_outcar(d / "OUTCAR")
        o["case"] = d.name
        o["endpoint"] = meta.get("endpoint")
        o["cls"] = meta.get("class")
        o["density_target"] = meta.get("density_target")
        o["kgrid"] = meta.get("kgrid")
        o["shift_value"] = meta.get("shift_value")
        o["natoms"] = meta.get("natoms")
        # validate_dir stores INCAR tag values as strings; compare as ints
        o["ismear"] = int(meta["ismear"]) if meta.get("ismear") is not None else None
        o["sigma"] = meta.get("sigma")
        o["dir"] = str(d)
        o["site_moments"] = parse_site_moments(d / "OUTCAR")
        rows.append(o)
    return rows


def e_atom(r: dict) -> float:
    return r["energy_sigma0_eV"] / r["natoms"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path)
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()
    rows = load(args.root)

    bad = [r for r in rows if not r["finished"] or r["n_scf_converged_markers"] < 1]
    print(f"{len(rows)} cases with OUTCAR; {len(bad)} did not finish/reach EDIFF")
    if bad:
        for r in bad:
            print(f"   !! {r['endpoint']}/{r['case']}")

    summary: dict = {"n_cases": len(rows), "n_bad": len(bad)}

    # ---------------- Q1 / Q2: density and smearing, per representative ------
    print("\n=== Q1+Q2  density and smearing, per representative ===")
    print("    E is energy(sigma->0) per atom, eV; dE in meV/atom")
    hdr = (f"{'endpoint':<15} {'ISMEAR':>6} {'d=12':>13} {'d=16':>13} "
           f"{'dE(12-16)':>10} {'dE(s2-s1)@16':>13} {'dE(sm5-s1)@16':>14} "
           f"{'mag@16':>9} {'Ef@16':>8}")
    print(hdr)
    print("-" * len(hdr))
    q1: dict = {}
    for ep in sorted({r["endpoint"] for r in rows if r["density_target"]}):
        for tag, ism in (("s1", 1), ("s2", 2), ("sm5", -5)):
            a = next((r for r in rows if r["endpoint"] == ep and r["ismear"] == ism
                      and r["density_target"] == 12), None)
            b = next((r for r in rows if r["endpoint"] == ep and r["ismear"] == ism
                      and r["density_target"] == 16), None)
            if not a or not b:
                continue
            d12, d16 = e_atom(a), e_atom(b)
            q1.setdefault(ep, {})[tag] = {
                "E_d12_eV_per_atom": d12, "E_d16_eV_per_atom": d16,
                "dE_12_16_meV_per_atom": (d12 - d16) * MEV}
            base1 = next((r for r in rows if r["endpoint"] == ep and r["ismear"] == 1
                          and r["density_target"] == 16), None)
            d_s2 = d_sm5 = None
            if base1 is not None:
                if ism == 2:
                    d_s2 = (e_atom(b) - e_atom(base1)) * MEV
                if ism == -5:
                    d_sm5 = (e_atom(b) - e_atom(base1)) * MEV
            mag = b["total_magnetization_muB"]
            ef = b["e_fermi_eV"]
            print(f"{ep:<15} {ism:>6} {d12:>13.8f} {d16:>13.8f} "
                  f"{(d12 - d16) * MEV:>10.3f} "
                  f"{(f'{d_s2:>13.3f}' if d_s2 is not None else f'{'-':>13}')} "
                  f"{(f'{d_sm5:>14.3f}' if d_sm5 is not None else f'{'-':>14}')} "
                  f"{(f'{mag:>9.4f}' if mag is not None else f'{'-':>9}')} "
                  f"{(f'{ef:>8.4f}' if ef is not None else f'{'-':>8}')}")
    summary["q1_density_smearing"] = q1

    # ---------------- Q3: Cr magnetic screening -----------------------------
    print("\n=== Q3  Cr-bearing representative: magnetic status ===")
    cr = [r for r in rows if r["cls"] == "magnetic screening"]
    if cr:
        ref = next((r for r in cr if r["case"].startswith("mag_ispin1")), None)
        print(f"{'policy':<34} {'E/atom':>13} {'dE vs ISPIN1':>13} "
              f"{'tot mag':>9} {'site moments':>28}")
        for r in sorted(cr, key=lambda x: x["case"]):
            d = (e_atom(r) - e_atom(ref)) * MEV if ref else None
            sm = r["site_moments"]
            sms = (f"Al {sm[0]:+.4f} / Cr {sm[1]:+.4f}" if sm and len(sm) == 2
                   else "n/a")
            mag = r["total_magnetization_muB"]
            print(f"{r['case']:<34} {e_atom(r):>13.8f} "
                  f"{(f'{d:>13.3f}' if d is not None else f'{'-':>13}')} "
                  f"{(f'{mag:>9.4f}' if mag is not None else f'{'-':>9}')} {sms:>28}")
        if ref:
            spread = (max(e_atom(r) for r in cr) - min(e_atom(r) for r in cr)) * MEV
            mags = [abs(m) for m in (r["total_magnetization_muB"] for r in cr)
                    if m is not None]
            print(f"\n  energy spread across the three spin policies: {spread:.3f} meV/atom")
            print(f"  |total moment| range: {min(mags):.4f} .. {max(mags):.4f} muB")
            summary["q3_cr"] = {"energy_spread_meV_per_atom": spread,
                                "abs_moment_range": [min(mags), max(mags)]}

    # ---------------- Q4: the two Gamma-centred families ---------------------
    print("\n=== Q4  Gamma-centred families on the AlTa control cell ===")
    alta = [r for r in rows if r["endpoint"] == "BCC-Al-Ta-2"]
    even = sorted((r for r in alta if r["density_target"] and r["ismear"] == 1),
                  key=lambda r: r["kgrid"])
    odd = sorted((r for r in alta if r["cls"] == "family probe" and r["ismear"] == 1),
                 key=lambda r: r["kgrid"])
    print("  family A: even n, shift 0.5  (Gamma-centred, fold-matchable with a 2x supercell)")
    for r in even:
        m = r["total_magnetization_muB"]
        print(f"    n={r['kgrid']:>2} dk={2 * 3.141592653589793 / (3.2246015525754168 * r['kgrid']):.5f} "
              f"E/atom={e_atom(r):.8f}  mag={m if m is not None else 'n/a'}")
    print("  family B: odd n, shift 0  (Gamma-centred, NOT fold-matchable)")
    for r in odd:
        m = r["total_magnetization_muB"]
        print(f"    n={r['kgrid']:>2} dk={2 * 3.141592653589793 / (3.2246015525754168 * r['kgrid']):.5f} "
              f"E/atom={e_atom(r):.8f}  mag={m if m is not None else 'n/a'}")
    if even and odd:
        ref_e = even[-1]
        ref_o = odd[-1]
        print(f"\n  referenced to each family's own densest mesh "
              f"(even n={ref_e['kgrid']}, odd n={ref_o['kgrid']}):")
        print(f"    {'dk/A^-1':>9} {'even err':>10} {'odd err':>10}")
        for r in even:
            dk = 2 * 3.141592653589793 / (3.2246015525754168 * r["kgrid"])
            print(f"    {dk:>9.5f} {(e_atom(r) - e_atom(ref_e)) * MEV:>10.3f}")
        for r in odd:
            dk = 2 * 3.141592653589793 / (3.2246015525754168 * r["kgrid"])
            print(f"    {dk:>9.5f} {'':>10} {(e_atom(r) - e_atom(ref_o)) * MEV:>10.3f}")
        summary["q4_families"] = {
            "even": {r["kgrid"]: e_atom(r) for r in even},
            "odd": {r["kgrid"]: e_atom(r) for r in odd}}

    # ---------------- diagnostics ------------------------------------------
    print("\n=== diagnostics ===")
    hdr = (f"{'endpoint':<15} {'case':<24} {'ISMEAR':>6} {'k':>3} {'sh':>4} "
           f"{'nkpt':>5} {'SCF':>4} {'maxF':>8} {'mag':>9} {'TS/eV':>10}")
    print(hdr)
    for r in sorted(rows, key=lambda x: (x["endpoint"] or "", x["case"])):
        mag = r["total_magnetization_muB"]
        print(f"{r['endpoint']:<15} {r['case']:<24} {r['ismear']:>6} {r['kgrid']:>3} "
              f"{r['shift_value']:>4} {r['nkpts_ibz']:>5} {r['n_scf_iterations']:>4} "
              f"{str(r['max_force_eV_per_A']):>8} "
              f"{(f'{mag:>9.4f}' if mag is not None else f'{'-':>9}')} "
              f"{str(r['entropy_TS_eV']):>10}")

    summary["rows"] = [{k: v for k, v in r.items() if k != "site_moments"}
                       | {"site_moments": r["site_moments"]} for r in rows]
    if args.json:
        args.json.write_text(json.dumps(summary, indent=2) + "\n")
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
