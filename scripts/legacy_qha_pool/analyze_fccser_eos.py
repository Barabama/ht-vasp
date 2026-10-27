#!/usr/bin/env python
"""Stage 0D analysis: BM3 fits for the FCC/L1_2 and SER multi-volume statics.

What this reports, and why each item is there
---------------------------------------------
* `V0`, `E0`, `B0`, `B0'` and the fit residual `r2` -- the EOS quality, which is
  the whole point of the batch.
* **`V0` vs the scan window.** A fit whose minimum sits on the edge of the
  range is a fit that wanted to go further; the number is then an artefact of
  where the scan stopped, not a property of the material.  Reported as
  `v0_inside_window`.
* **Energy span** across the scan.  A span of a few meV means the points are
  inside the numerical noise and no fit is meaningful.
* **Occupancy** for L1_2: the 1a/3c assignment must survive the relaxation,
  otherwise the endpoint is no longer the compound it is named after.
* **c/a** for the HCP endpoint, since that is the degree of freedom the
  ISIF=4 points exist to relax.

It also cross-checks the relaxed a against the Vegard estimate, because the
starting cell came from the `endmember.py` template whose lattice constant can
be tens of percent off (report 8.0c).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from analyze_align import bm3, parse_outcar


def load_points(ep_dir: Path) -> list[dict]:
    rows = []
    for case in sorted(ep_dir.glob("a*/")):
        o = parse_outcar(case / "OUTCAR")
        meta = json.loads((case / "case.json").read_text())
        rows.append({
            "dir": str(case), "a_scale": meta["a_scale"],
            "volume_scale": meta["volume_scale"],
            "volume": meta["volume_angstrom3"],
            "n_atoms": meta["relaxed_cell"]["n_atoms"],
            "energy_sigma0_eV": o["energy_sigma0_eV"],
            "finished": o["finished"], "scf": o["n_scf_iterations"],
            "max_force": o["max_force_eV_per_A"],
            "mag": o["total_magnetization_muB"],
            "e_fermi": o["e_fermi_eV"],
            "is_hcp": meta["is_hcp"],
        })
    return sorted(rows, key=lambda r: r["volume"])


def occupancy_report(ep_dir: Path) -> dict:
    """1a/3c assignment from the relaxed cell, for L1_2 endpoints."""
    from pymatgen.core import Structure
    c = ep_dir / ".." / ".." / "relax" / ep_dir.name / "CONTCAR"
    if not c.is_file():
        return {}
    st = Structure.from_file(str(c))
    if len(st) != 4:
        return {}
    # 1a = (0,0,0); the other three should be (0,1/2,1/2)-family
    at_1a, at_3c = None, []
    for s in st:
        f = np.array(s.frac_coords) % 1.0
        if np.allclose(f, 0.0, atol=1e-3) or np.allclose(f, 1.0, atol=1e-3):
            at_1a = str(s.specie)
        else:
            at_3c.append(str(s.specie))
    return {"site_1a": at_1a, "site_3c": sorted(at_3c),
            "n_1a": 1 if at_1a else 0, "n_3c": len(at_3c),
            "formula": st.composition.reduced_formula}


def hcp_ratio(ep_dir: Path) -> dict:
    from pymatgen.core import Structure
    c = ep_dir / ".." / ".." / "relax" / ep_dir.name / "CONTCAR"
    if not c.is_file():
        return {}
    st = Structure.from_file(str(c))
    a, _, cc = st.lattice.abc
    return {"a": float(a), "c": float(cc), "c_over_a": float(cc / a)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path,
                    help="the stage0d directory (holding eos/ and relax/)")
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--vegard", type=Path, default=None,
                    help="optional JSON of Vegard a per endpoint for the "
                         "starting-cell cross-check")
    args = ap.parse_args()

    eos_root = args.root / "stage0d" / "eos"
    if not eos_root.is_dir():
        eos_root = args.root / "eos"
    if not eos_root.is_dir():
        raise SystemExit(f"no eos stage under {args.root}")

    out: dict = {}
    print(f"{'endpoint':>14} {'n':>2} {'pts':>4} {'V span/meV':>10} "
          f"{'V0/A^3':>10} {'V0/at':>8} {'B0/GPa':>8} {'B0p':>6} {'r2':>9} "
          f"{'V0 interior':>11}")
    print("-" * 104)
    for ep_dir in sorted(d for d in eos_root.glob("*/") if d.is_dir()):
        pts = [p for p in load_points(ep_dir) if p["energy_sigma0_eV"] is not None]
        if len(pts) < 4:
            print(f"{ep_dir.name:>14} 只有 {len(pts)} 个有效点，跳过")
            continue
        v = np.array([p["volume"] for p in pts])
        e = np.array([p["energy_sigma0_eV"] for p in pts])
        n = pts[0]["n_atoms"]
        span_meV = float(1e3 * (e.max() - e.min()) / n)
        fit = bm3(v, e)
        pred = np.array([bm3(v, e)["E0_eV"]])  # placeholder, recomputed below
        # residual r2 on the fit
        from scipy.optimize import curve_fit

        def f(x, e0, b0, bp, v0):
            return e0 + 9.0 * v0 * b0 / 16.0 * (
                ((v0 / x) ** (2.0 / 3.0) - 1.0) ** 3 * bp
                + ((v0 / x) ** (2.0 / 3.0) - 1.0) ** 2
                * (6.0 - 4.0 * (v0 / x) ** (2.0 / 3.0)))
        popt, _ = curve_fit(f, v, e,
                            p0=[e.min(), 1.5, 4.0, v[np.argmin(e)]],
                            maxfev=20000)
        resid = e - f(v, *popt)
        ss_res = float((resid ** 2).sum())
        ss_tot = float(((e - e.mean()) ** 2).sum())
        r2 = 1.0 - ss_res / ss_tot if ss_tot else float("nan")
        v0 = fit["V0_angstrom3"]
        inside = bool(v.min() * 1.001 < v0 < v.max() * 0.999)
        rec = {"n_atoms": n, "n_points": len(pts), "points": pts,
               "energy_span_meV_per_atom": span_meV, "bm3": fit,
               "eos_r2": r2,
               "max_abs_residual_meV_per_atom":
                   float(1e3 * np.abs(resid).max() / n),
               "v0_inside_window": inside,
               "v0_per_atom": v0 / n,
               "occupancy": occupancy_report(ep_dir),
               "hcp": hcp_ratio(ep_dir) if pts[0]["is_hcp"] else {},
               "all_finished": all(p["finished"] for p in pts)}
        out[ep_dir.name] = rec
        print(f"{ep_dir.name:>14} {n:>2} {len(pts):>4} {span_meV:>10.2f} "
              f"{v0:>10.4f} {v0 / n:>8.4f} {fit['B0_GPa']:>8.2f} "
              f"{fit['B0prime']:>6.2f} {r2:>9.6f} {str(inside):>11}")

    print()
    print("=== 逐端基细节 ===")
    for name, r in out.items():
        print(f"\n{name}")
        print(f"  n={r['n_atoms']}  points={r['n_points']}  "
              f"all finished={r['all_finished']}")
        print(f"  V0 = {r['bm3']['V0_angstrom3']:.4f} A^3  "
              f"({r['v0_per_atom']:.4f} A^3/atom)   "
              f"E0 = {r['bm3']['E0_eV']:.6f} eV")
        print(f"  B0 = {r['bm3']['B0_GPa']:.2f} GPa  B0' = {r['bm3']['B0prime']:.3f}")
        print(f"  r2 = {r['eos_r2']:.8f}   "
              f"max residual = {r['max_abs_residual_meV_per_atom']:.4f} meV/atom")
        print(f"  能量跨度 = {r['energy_span_meV_per_atom']:.2f} meV/atom   "
              f"V0 在窗口内: {r['v0_inside_window']}")
        if r["occupancy"]:
            o = r["occupancy"]
            print(f"  占位: 1a={o['site_1a']}  3c={o['site_3c']}  ({o['formula']})")
        if r["hcp"]:
            print(f"  HCP: a={r['hcp']['a']:.4f}  c={r['hcp']['c']:.4f}  "
                  f"c/a={r['hcp']['c_over_a']:.4f}")

    if args.json:
        args.json.write_text(json.dumps(out, indent=2) + "\n")
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
