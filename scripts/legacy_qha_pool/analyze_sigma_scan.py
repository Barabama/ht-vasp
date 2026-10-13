#!/usr/bin/env python
"""Analyse the SIGMA scan.

Question: which SIGMA keeps E(V) smooth AND leaves the SCF/forces sane?

Everything is per atom; ``energy(sigma->0)`` is the primary metric because that
is what pymatgen ``Vasprun.final_energy`` returns and therefore what the
production QHA consumes (report §4.6.1).  ``free energy TOTEN`` and
``energy without entropy`` are carried as controls, plus ``T*S``, forces,
stress, total and site moments.

Engineering criteria (the review's, explicitly project-internal, not physics
constants):

    range[ΔE_sigma(V)] <= 0.5 meV/atom
    AND |dV0| / V0 < 0.1%
    AND |dB0| < 1 GPa
    AND T*S controlled

ΔE_sigma(V) = E_sigma(V) - E_ref(V), with ref = ISMEAR=1/SIGMA=0.20.
The range is taken over the volume points, so it measures the VOLUME-DEPENDENT
part -- a constant offset is harmless because BM3 absorbs it (report §4.6.6).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from analyze_align import bm3, parse_outcar
from analyze_stage0 import parse_site_moments

MEV = 1000.0
REF = "s1_020"
ORDER = ("s1_005", "s1_010", "s1_020", "s2_010", "sm5")
LABEL = {"s1_005": "1 / 0.05", "s1_010": "1 / 0.10", "s1_020": "1 / 0.20",
         "s2_010": "2 / 0.10", "sm5": "-5 / ---"}
CRITERIA = {"range_dE_meV_per_atom": 0.5, "dV0_percent": 0.1, "dB0_GPa": 1.0}


def load(root: Path) -> list[dict]:
    rows = []
    for cj in sorted((root / "sigma").rglob("case.json")):
        d = cj.parent
        if any(p.startswith(".") for p in d.relative_to(root / "sigma").parts):
            continue
        if not (d / "OUTCAR").is_file():
            continue
        meta = json.loads(cj.read_text())
        o = parse_outcar(d / "OUTCAR")
        o["case"] = d.name
        o["tag"] = d.name.split("_", 1)[1]
        o["endpoint"] = meta["endpoint"]
        o["frame"] = meta["frame"]
        o["eps"] = meta["eps"]
        o["natoms"] = meta["natoms"]
        o["kgrid"] = meta["kgrid"]
        o["volume_angstrom3"] = meta["volume_angstrom3"]
        o["volume_scale_vs_v0"] = meta["volume_scale_vs_v0"]
        o["ismear"] = int(meta["ismear"])
        o["sigma"] = meta["sigma"]
        o["dir"] = str(d)
        o["site_moments"] = parse_site_moments(d / "OUTCAR")
        rows.append(o)
    return rows


def e_atom(r: dict) -> float:
    return r["energy_sigma0_eV"] / r["natoms"]


def usable(r: dict) -> bool:
    """A case can only enter an energy comparison once it has an energy."""
    return (r.get("energy_sigma0_eV") is not None and r.get("natoms")
            and r.get("finished") and r.get("n_scf_converged_markers", 0) >= 1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path)
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()
    rows = load(args.root)

    bad = [r for r in rows if not r["finished"] or r["n_scf_converged_markers"] < 1]
    print(f"{len(rows)} cases; {len(bad)} did not finish/reach EDIFF")
    for r in bad:
        print(f"   !! {r['endpoint']}/{r['case']}")
    summary: dict = {"n_cases": len(rows), "n_bad": len(bad), "criteria": CRITERIA}

    endpoints = sorted({r["endpoint"] for r in rows})

    # ---------------- 1. does SIGMA reduce the volume-dependent scatter? ----
    print("\n=== 1. volume-dependent scatter of E_sigma(V) ===")
    print("    dE(V) = E_sigma(V) - E_ref(V) with ref = ISMEAR=1/SIGMA=0.20;")
    print("    range = max-min over the 5 volume points (meV/atom)")
    hdr = (f"{'endpoint':<15} {'smearing':>10} {'range dE(V)':>12} "
           f"{'mean dE':>9} {'max|dE|':>9} {'max|F|':>8} {'T*S min':>9} {'T*S max':>9}")
    print(hdr)
    print("-" * len(hdr))
    scatter: dict = {}
    for ep in endpoints:
        ref = {r["frame"]: r for r in rows
               if r["endpoint"] == ep and r["tag"] == REF and usable(r)}
        for tag in ORDER:
            rs = sorted((r for r in rows
                         if r["endpoint"] == ep and r["tag"] == tag and usable(r)),
                        key=lambda x: x["frame"])
            if not rs or not ref:
                continue
            d = [(e_atom(r) - e_atom(ref[r["frame"]])) * MEV for r in rs
                 if r["frame"] in ref]
            ts = [r["entropy_TS_eV"] for r in rs if r["entropy_TS_eV"] is not None]
            fs = [r["max_force_eV_per_A"] for r in rs
                  if r["max_force_eV_per_A"] is not None]
            if not d:
                continue
            rec = {"range": max(d) - min(d), "mean": float(np.mean(d)),
                   "max_abs": max(abs(x) for x in d),
                   "dE_by_frame": {r["frame"]: dd for r, dd in zip(rs, d)},
                   "ts_range": [min(ts), max(ts)] if ts else None,
                   "max_force": max(fs) if fs else None}
            scatter.setdefault(ep, {})[tag] = rec
            mf = rec["max_force"]
            ts0 = rec["ts_range"][0] if rec["ts_range"] else None
            ts1 = rec["ts_range"][1] if rec["ts_range"] else None
            print(f"{ep:<15} {LABEL[tag]:>10} {rec['range']:>12.3f} "
                  f"{rec['mean']:>9.3f} {rec['max_abs']:>9.3f} "
                  f"{(f'{mf:>8.4f}' if mf is not None else f'{'-':>8}')} "
                  f"{(f'{ts0:>9.5f}' if ts0 is not None else f'{'-':>9}')} "
                  f"{(f'{ts1:>9.5f}' if ts1 is not None else f'{'-':>9}')}")
    summary["scatter"] = scatter

    # ---------------- 2. per-volume detail for the most sensitive endpoints --
    for ep in endpoints:
        print(f"\n=== 2. {ep}: E(sigma->0)/atom by volume and smearing ===")
        hdr = f"{'frame':>5} {'V/V0':>8} " + " ".join(f"{LABEL[t]:>13}" for t in ORDER)
        print(hdr)
        print("-" * len(hdr))
        ref = {r["frame"]: r for r in rows if r["endpoint"] == ep and r["tag"] == REF}
        for fr in sorted({r["frame"] for r in rows if r["endpoint"] == ep}):
            rs = {r["tag"]: r for r in rows
                  if r["endpoint"] == ep and r["frame"] == fr and usable(r)}
            if REF not in rs:
                continue
            vs = rs[REF]["volume_scale_vs_v0"]
            cells = []
            for t in ORDER:
                if t in rs:
                    cells.append(f"{e_atom(rs[t]):>13.8f}")
                else:
                    cells.append(f"{'-':>13}")
            print(f"{fr:>5} {vs:>8.5f} " + " ".join(cells))

    # ---------------- 3. EOS per smearing ----------------------------------
    print("\n=== 3. BM3 EOS per smearing (5 production volume points) ===")
    eos_all: dict = {}
    for ep in endpoints:
        eos_all[ep] = {}
        print(f"  {ep}")
        for tag in ORDER:
            rs = sorted(((r["volume_angstrom3"], e_atom(r)) for r in rows
                         if r["endpoint"] == ep and r["tag"] == tag and usable(r)),
                        key=lambda x: x[0])
            if len(rs) < 4:
                continue
            try:
                f = bm3([p[0] for p in rs], [p[1] for p in rs])
            except Exception as exc:                      # noqa: BLE001
                print(f"    {LABEL[tag]:>10}: fit failed ({exc})")
                continue
            eos_all[ep][tag] = f
            print(f"    {LABEL[tag]:>10}: V0={f['V0_angstrom3']:>9.4f} A^3  "
                  f"E0={f['E0_eV']:>13.8f} eV  B0={f['B0_GPa']:>7.2f} GPa  "
                  f"B0'={f['B0prime']:>6.3f}")
        if REF in eos_all[ep] and len(eos_all[ep]) > 1:
            r0 = eos_all[ep][REF]
            print(f"    {'vs ref':>10}:  "
                  f"{'dV0/A^3':>10} {'dV0/V0':>9} {'dB0/GPa':>9} {'dE0/meV/at':>11}")
            for tag, f in eos_all[ep].items():
                if tag == REF:
                    continue
                dv = f["V0_angstrom3"] - r0["V0_angstrom3"]
                print(f"    {LABEL[tag]:>10}:  {dv:>10.4f} "
                      f"{100 * dv / r0['V0_angstrom3']:>8.4f}% "
                      f"{f['B0_GPa'] - r0['B0_GPa']:>9.2f} "
                      f"{(f['E0_eV'] - r0['E0_eV']) * MEV:>11.3f}")
    summary["eos"] = eos_all

    # ---------------- 4. magnetic branch reproducibility -------------------
    print("\n=== 4. magnetic branch reproducibility (total and site moments) ===")
    mag_summary: dict = {}
    for ep in endpoints:
        print(f"  {ep}")
        for tag in ORDER:
            rs = sorted((r for r in rows if r["endpoint"] == ep
                         and r["tag"] == tag and r["finished"]),
                        key=lambda x: x["frame"])
            if not rs:
                continue
            tot = [r["total_magnetization_muB"] for r in rs
                   if r["total_magnetization_muB"] is not None]
            sites = [r["site_moments"] for r in rs if r["site_moments"]]
            smin = smax = None
            if sites:
                flat = list(zip(*sites))
                smin = [min(x) for x in flat]
                smax = [max(x) for x in flat]
            mag_summary.setdefault(ep, {})[tag] = {
                "total_range": [min(tot), max(tot)] if tot else None,
                "site_min": smin, "site_max": smax}
            tr = (f"{min(tot):+.4f}..{max(tot):+.4f}" if tot else "n/a")
            sr = (", ".join(f"{lo:+.3f}..{hi:+.3f}" for lo, hi in zip(smin, smax))
                  if smin else "n/a")
            print(f"    {LABEL[tag]:>10}: tot {tr:>20}   sites {sr}")
    summary["magnetic"] = mag_summary

    # ---------------- 5. verdict -------------------------------------------
    print("\n=== 5. verdict against the engineering criteria ===")
    print(f"    range dE(V) <= {CRITERIA['range_dE_meV_per_atom']} meV/atom, "
          f"|dV0|/V0 < {CRITERIA['dV0_percent']}%, |dB0| < {CRITERIA['dB0_GPa']} GPa")
    verdict: dict = {}
    for ep in endpoints:
        for tag in ORDER:
            sc = scatter.get(ep, {}).get(tag)
            if not sc:
                continue
            f = eos_all.get(ep, {}).get(tag)
            r0 = eos_all.get(ep, {}).get(REF)
            ok_range = sc["range"] <= CRITERIA["range_dE_meV_per_atom"]
            ok_v = ok_b = None
            if f and r0:
                ok_v = abs(100 * (f["V0_angstrom3"] - r0["V0_angstrom3"])
                           / r0["V0_angstrom3"]) < CRITERIA["dV0_percent"]
                ok_b = abs(f["B0_GPa"] - r0["B0_GPa"]) < CRITERIA["dB0_GPa"]
            verdict.setdefault(ep, {})[tag] = {
                "range_dE_meV_per_atom": sc["range"], "range_ok": ok_range,
                "dV0_ok": ok_v, "dB0_ok": ok_b,
                "max_force_eV_per_A": sc["max_force"]}
            flags = "".join(["R" if ok_range else "-",
                             "V" if ok_v else ("-" if ok_v is not None else "?"),
                             "B" if ok_b else ("-" if ok_b is not None else "?")])
            print(f"  {ep:<15} {LABEL[tag]:>10}: range={sc['range']:>6.3f} "
                  f"maxF={sc['max_force'] if sc['max_force'] is not None else float('nan'):>7.4f}  "
                  f"[{flags}]  (R=range,V=dV0,B=dB0 pass)")
    summary["verdict"] = verdict
    print("\n    flags: R = range dE(V) passes, V = dV0 passes, B = dB0 passes;")
    print("           '-' = fails, '?' = EOS fit unavailable")

    summary["rows"] = rows
    if args.json:
        args.json.write_text(json.dumps(summary, indent=2, default=str) + "\n")
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
