#!/usr/bin/env python
"""Analyse the multi-volume ISMEAR alignment test.

Answers one question:
    Is E(ISMEAR=-5) - E(ISMEAR=1), evaluated at the SAME structure, a
    volume-independent constant (which QHA can absorb), or does it vary with
    volume (in which case the pool must use a single smearing everywhere)?

The metric that matters is ``energy(sigma->0)``, because that is exactly what
the production pipeline consumes: pymatgen ``Vasprun.final_energy`` returns
``energy(sigma->0)`` (verified on the force test), and atomate2's TaskDoc
surfaces that as ``output["energy"]``, which is what QHA fits E(V) to.
``energy without entropy`` and ``free energy TOTEN`` are recorded as controls.

Thresholds (from the review):
    range[ΔE(V)] < 1 meV/atom   -> constant-shift alignment is fine
    1 .. 2 meV/atom             -> needs_review
    > 2 meV/atom                -> smearings must be unified
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

EV_PER_ATOM_TO_MEV = 1000.0


def parse_outcar(p: Path) -> dict:
    txt = p.read_text(errors="replace")
    out: dict = {}

    # ---- energies (last occurrence wins) --------------------------------
    def last_float(pattern: str) -> float | None:
        m = re.findall(pattern, txt, re.M)
        return float(m[-1]) if m else None

    out["free_energy_toten_eV"] = last_float(r"free\s+energy\s+TOTEN\s*=\s*(-?\d+\.\d+)")
    m = re.findall(r"energy\s+without entropy=\s*(-?\d+\.\d+)\s+energy\(sigma->0\)\s*=\s*(-?\d+\.\d+)", txt)
    if m:
        out["energy_without_entropy_eV"] = float(m[-1][0])
        out["energy_sigma0_eV"] = float(m[-1][1])
    else:
        out["energy_without_entropy_eV"] = None
        out["energy_sigma0_eV"] = None
    out["entropy_TS_eV"] = last_float(r"entropy\s+T\*S\s+EENTRO\s*=\s*(-?\d+\.\d+)")

    # ---- forces ---------------------------------------------------------
    m = re.findall(r"FORCES:\s+max atom, RMS\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", txt)
    if m:
        out["max_force_eV_per_A"] = float(m[-1][0])
        out["rms_force_eV_per_A"] = float(m[-1][1])
    else:
        out["max_force_eV_per_A"] = None
        out["rms_force_eV_per_A"] = None

    # ---- stress (external pressure) -------------------------------------
    out["external_pressure_kB"] = last_float(r"external pressure\s*=\s*(-?\d+\.\d+)\s+kB")

    # ---- magnetization --------------------------------------------------
    out["total_magnetization_muB"] = last_float(
        r"number of electron\s+\d+\.\d+\s+magnetization\s+(-?\d+\.\d+)")

    # ---- Fermi level ----------------------------------------------------
    out["e_fermi_eV"] = last_float(r"E-fermi\s*:\s*(-?\d+\.\d+)")

    # ---- convergence accounting ----------------------------------------
    out["n_scf_converged_markers"] = txt.count("aborting loop because EDIFF is reached")
    out["n_unconverged_markers"] = txt.count("unconverged")
    # SCF steps are counted from OSZICAR, not OUTCAR: OUTCAR's per-iteration
    # header is "Iteration  n(m)" and does not carry the energy.
    osz = p.parent / "OSZICAR"
    if osz.is_file():
        out["n_scf_iterations"] = len(re.findall(r"^(?:DAV|RMM|CG|DIIS):\s*\d+",
                                                osz.read_text(errors="replace"), re.M))
    else:
        out["n_scf_iterations"] = None
    out["has_volume_basis"] = "VOLUME and BASIS-vectors are now" in txt
    out["finished"] = "General timing and accounting informations" in txt

    # ---- k-points -------------------------------------------------------
    ibz = p.parent / "IBZKPT"
    if ibz.is_file():
        lines = ibz.read_text().splitlines()
        out["nkpts_ibz"] = int(lines[1].split()[0]) if len(lines) > 1 else None
    else:
        out["nkpts_ibz"] = None

    # ---- volume from POSCAR (authoritative) -----------------------------
    pos = p.parent / "POSCAR"
    if pos.is_file():
        pl = pos.read_text().splitlines()
        scale = float(pl[1].split()[0])
        lat = [[scale * float(x) for x in pl[i].split()] for i in (2, 3, 4)]
        a, b, c = lat
        vol = (a[0] * (b[1] * c[2] - b[2] * c[1])
               - a[1] * (b[0] * c[2] - b[2] * c[0])
               + a[2] * (b[0] * c[1] - b[1] * c[0]))
        out["volume_angstrom3"] = abs(vol)
        out["a_angstrom"] = lat[0][0]
        out["natoms"] = sum(int(x) for x in pl[6].split())
    return out


def bm3(volumes, energies):
    """Birch-Murnaghan 3rd order via scipy; returns (V0, E0, B0_GPa)."""
    import numpy as np
    from scipy.optimize import curve_fit

    v = np.asarray(volumes, float)
    e = np.asarray(energies, float)

    def f(x, e0, b0, bp, v0):
        return e0 + 9.0 * v0 * b0 / 16.0 * (
            ((v0 / x) ** (2.0 / 3.0) - 1.0) ** 3 * bp
            + ((v0 / x) ** (2.0 / 3.0) - 1.0) ** 2 * (6.0 - 4.0 * (v0 / x) ** (2.0 / 3.0)))

    b0_guess = 1.5
    p0 = [e.min(), b0_guess, 4.0, v[np.argmin(e)]]
    popt, _ = curve_fit(f, v, e, p0=p0, maxfev=20000)
    e0, b0_ev_a3, bp, v0 = popt
    b0_gpa = b0_ev_a3 * 160.21766208
    return {"V0_angstrom3": float(v0), "E0_eV": float(e0),
            "B0_GPa": float(b0_gpa), "B0prime": float(bp)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path)
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()
    root = args.root

    rows: list[dict] = []
    for case in sorted(root.glob("*/*/")):
        if not (case / "OUTCAR").is_file():
            continue
        meta = json.loads((case / "case.json").read_text()) if (case / "case.json").is_file() else {}
        d = parse_outcar(case / "OUTCAR")
        d["case"] = case.name
        d["block"] = meta.get("block", case.parent.name)
        d["frame"] = meta.get("frame")
        d["eps"] = meta.get("eps")
        d["volume_scale_vs_v0"] = meta.get("volume_scale_vs_v0")
        d["ismear"] = meta.get("ismear")
        d["sigma"] = meta.get("sigma")
        d["kgrid"] = meta.get("kgrid")
        d["dir"] = str(case)
        rows.append(d)

    # ---------- primary: ΔE(V) -------------------------------------------
    prim = [r for r in rows if r["block"] == "primary"]
    byframe: dict[int, dict[str, dict]] = {}
    for r in prim:
        byframe.setdefault(r["frame"], {})[r["case"].split("_")[-1]] = r

    n_atoms = next((r["natoms"] for r in rows if r.get("natoms")), 2)
    print(f"natoms/cell = {n_atoms};  energy metric = energy(sigma->0) "
          f"(= pymatgen final_energy, = what QHA consumes)\n")

    hdr = (f"{'frame':>5} {'V/A^3':>9} {'V/V0':>7} "
           f"{'E_sigma0(s1)':>14} {'E_sigma0(s2)':>14} {'E_sigma0(m5)':>14} "
           f"{'dE(5->1)':>10} {'dE(5->2)':>10} {'dE(2->1)':>10}")
    print(hdr)
    print("-" * len(hdr))
    deltas = {"dE_m5_s1": [], "dE_m5_s2": [], "dE_s2_s1": []}
    vols = []
    for fr in sorted(byframe):
        g = byframe[fr]
        if not all(k in g for k in ("s1", "s2", "sm5")):
            print(f"{fr:>5}  (incomplete: {sorted(g)})")
            continue
        v = g["s1"]["volume_angstrom3"]
        e1 = g["s1"]["energy_sigma0_eV"]
        e2 = g["s2"]["energy_sigma0_eV"]
        e5 = g["sm5"]["energy_sigma0_eV"]
        d51 = (e5 - e1) * EV_PER_ATOM_TO_MEV / n_atoms
        d52 = (e5 - e2) * EV_PER_ATOM_TO_MEV / n_atoms
        d21 = (e2 - e1) * EV_PER_ATOM_TO_MEV / n_atoms
        vols.append(v)
        deltas["dE_m5_s1"].append(d51)
        deltas["dE_m5_s2"].append(d52)
        deltas["dE_s2_s1"].append(d21)
        print(f"{fr:>5} {v:>9.4f} {g['s1']['volume_scale_vs_v0']:>7.5f} "
              f"{e1:>14.8f} {e2:>14.8f} {e5:>14.8f} "
              f"{d51:>10.3f} {d52:>10.3f} {d21:>10.3f}")
    print("\n  units: dE columns are meV/atom; E columns are eV/cell\n")

    summary: dict = {"natoms": n_atoms, "volumes_angstrom3": vols, "deltas_meV_per_atom": deltas}
    for key, arr in deltas.items():
        if arr:
            summary[key] = {"min": min(arr), "max": max(arr),
                            "range": max(arr) - min(arr),
                            "mean": sum(arr) / len(arr)}
    if deltas["dE_m5_s1"]:
        rng = summary["dE_m5_s1"]["range"]
        verdict = ("OK_constant_shift" if rng < 1.0
                   else "needs_review" if rng <= 2.0
                   else "must_unify_smearing")
        summary["verdict_dE_m5_s1"] = verdict
        print(f"range[dE(-5 -> s1)] = {rng:.3f} meV/atom  ->  {verdict}")
        print(f"range[dE(-5 -> s2)] = {summary['dE_m5_s2']['range']:.3f} meV/atom")
        print(f"range[dE(s2 -> s1)] = {summary['dE_s2_s1']['range']:.3f} meV/atom")
        print("  (a 'constant shift' means range << the absolute shift "
              f"|mean| = {abs(summary['dE_m5_s1']['mean']):.3f} meV/atom)")

    # ---------- EOS per smearing -----------------------------------------
    print("\nEOS fit per smearing (5 primary volumes):")
    summary["eos"] = {}
    for tag in ("s1", "s2", "sm5"):
        pts = sorted(((byframe[f][tag]["volume_angstrom3"],
                       byframe[f][tag]["energy_sigma0_eV"]) for f in byframe
                      if tag in byframe[f]), key=lambda x: x[0])
        if len(pts) < 4:
            continue
        try:
            fit = bm3([p[0] for p in pts], [p[1] for p in pts])
            summary["eos"][tag] = fit
            print(f"  {tag:>4}: V0={fit['V0_angstrom3']:.4f} A^3  "
                  f"E0={fit['E0_eV']:.6f} eV  B0={fit['B0_GPa']:.2f} GPa  "
                  f"B0'={fit['B0prime']:.3f}")
        except Exception as exc:  # noqa: BLE001
            print(f"  {tag:>4}: EOS fit failed: {exc}")
    if len(summary["eos"]) >= 2:
        v0s = [f["V0_angstrom3"] for f in summary["eos"].values()]
        b0s = [f["B0_GPa"] for f in summary["eos"].values()]
        e0s = [f["E0_eV"] for f in summary["eos"].values()]
        print(f"  spread: dV0={max(v0s) - min(v0s):.4f} A^3 "
              f"({100 * (max(v0s) - min(v0s)) / sum(v0s) * len(v0s) / 1:.3f}% of mean)  "
              f"dB0={max(b0s) - min(b0s):.2f} GPa  "
              f"dE0={1000 * (max(e0s) - min(e0s)) / n_atoms:.3f} meV/atom")
        summary["eos_spread"] = {"dV0": max(v0s) - min(v0s), "dB0": max(b0s) - min(b0s),
                                 "dE0_meV_per_atom": 1000 * (max(e0s) - min(e0s)) / n_atoms}

    # ---------- aux: k sensitivity and SIGMA sensitivity ------------------
    print("\naux: k-mesh sensitivity at the middle volume (k=11 vs k=7):")
    summary["kmesh"] = {}
    mid = byframe.get(4, {})
    for tag in ("s1", "s2", "sm5"):
        km = [r for r in rows if r["block"] == "kmesh" and r["case"].endswith(tag)]
        if tag in mid and km:
            e11, e7 = mid[tag]["energy_sigma0_eV"], km[0]["energy_sigma0_eV"]
            d = (e11 - e7) * 1000.0 / n_atoms
            summary["kmesh"][tag] = d
            print(f"  {tag:>4}: E(k11) - E(k7) = {d:+.3f} meV/atom "
                  f"({e11:.8f} vs {e7:.8f} eV)")
    print("\naux: SIGMA sensitivity, ISMEAR=1, middle volume, k=11:")
    summary["sigma"] = {}
    if "s1" in mid:
        base = mid["s1"]["energy_sigma0_eV"]
        for r in sorted((r for r in rows if r["block"] == "sigma"),
                        key=lambda x: x["sigma"]):
            d = (r["energy_sigma0_eV"] - base) * 1000.0 / n_atoms
            summary["sigma"][f"{r['sigma']}"] = d
            print(f"  SIGMA={r['sigma']:<5}: E - E(SIGMA=0.2) = {d:+.3f} meV/atom "
                  f"(T*S = {r['entropy_TS_eV']:+.6f} eV)")

    # ---------- aux: k-mesh convergence on the 2-atom cell ---------------
    # The alignment test's k=7 control exposed a mesh error far larger than the
    # smearing effect. This block puts numbers on it.
    kc = [r for r in rows if r["block"] == "kconv"]
    if kc:
        print("\naux: k-mesh convergence on the 2-atom production cell "
              "(middle volume, V=%.4f A^3):" % (kc[0]["volume_angstrom3"]))
        series: dict[str, list[tuple[int, float]]] = {}
        # include the k=7 / k=11 points already measured
        for tag, blk in (("s1", "kmesh"), ("s1", "primary"), ("s2", "kmesh"),
                         ("s2", "primary"), ("sm5", "kmesh"), ("sm5", "primary")):
            for r in rows:
                if r["block"] == blk and r["case"].endswith(tag) and r["frame"] == 4:
                    series.setdefault(tag, []).append((r["kgrid"], r["energy_sigma0_eV"]))
        for r in kc:
            tag = r["case"].split("_")[-1]
            series.setdefault(tag, []).append((r["kgrid"], r["energy_sigma0_eV"]))
        summary["kconv"] = {}
        nkmap = {(r["kgrid"], r["case"].split("_")[-1]): r["nkpts_ibz"]
                 for r in rows if r["frame"] == 4}
        for tag in ("s1", "s2", "sm5"):
            pts = sorted(set(series.get(tag, [])))
            if not pts:
                continue
            ref = max(p[0] for p in pts)
            eref = dict(pts)[ref]
            print(f"  ISMEAR={tag}: reference k={ref}")
            for k, e in pts:
                nk = nkmap.get((k, tag), "?")
                print(f"     k={k:>2} (nkpt={nk}): E={e:.8f}  "
                      f"E - E(k{ref}) = {(e - eref) * 1000.0 / n_atoms:+.3f} meV/atom")
            summary["kconv"][tag] = {int(k): float(e) for k, e in pts}
        # does dE(-5 -> s1) shrink as the mesh converges?
        if "s1" in summary["kconv"] and "sm5" in summary["kconv"]:
            common = sorted(set(summary["kconv"]["s1"]) & set(summary["kconv"]["sm5"]))
            print("  dE(-5 -> s1) vs mesh:")
            summary["dE_m5_s1_vs_k"] = {}
            for k in common:
                d = (summary["kconv"]["sm5"][k] - summary["kconv"]["s1"][k]) * 1000.0 / n_atoms
                summary["dE_m5_s1_vs_k"][k] = d
                print(f"     k={k:>2}: {d:+.3f} meV/atom")

    # ---------- diagnostics table ----------------------------------------
    print("\ndiagnostics (convergence, forces, stress, magnetization):")
    dh = (f"{'case':>20} {'ISMEAR':>7} {'k':>3} {'nkpt':>5} {'SCF':>4} "
          f"{'conv':>4} {'maxF':>9} {'P/kB':>8} {'mag':>9} {'TS/eV':>10} {'fin':>4}")
    print(dh)
    for r in sorted(rows, key=lambda x: (x["block"], x["frame"] or 0, x["case"])):
        print(f"{r['case']:>20} {r['ismear']:>7} {r['kgrid']:>3} {r['nkpts_ibz']:>5} "
              f"{r['n_scf_iterations']:>4} {r['n_scf_converged_markers']:>4} "
              f"{str(r['max_force_eV_per_A']):>9} "
              f"{str(r['external_pressure_kB']):>8} "
              f"{str(r['total_magnetization_muB']):>9} "
              f"{str(r['entropy_TS_eV']):>10} {str(r['finished'])[:1]:>4}")
    bad = [r for r in rows if not r["finished"] or r["n_scf_converged_markers"] < 1]
    if bad:
        print(f"\n!! {len(bad)} case(s) did not finish or did not reach EDIFF: "
              f"{[r['case'] for r in bad]}")
    else:
        print(f"\nall {len(rows)} cases finished and reached EDIFF (unconverged markers: "
              f"{sum(r['n_unconverged_markers'] for r in rows)})")

    # ---------- aux: the same 5 volumes at k=15 --------------------------
    k15 = [r for r in rows if r["block"] == "k15"]
    if k15:
        print("\naux: the same 5 volumes at k=15 (s1 and sm5 only) -- "
              "re-takes the verdict at a denser mesh:")
        byf15: dict[int, dict[str, dict]] = {}
        for r in k15:
            byf15.setdefault(r["frame"], {})[r["case"].split("_")[-1]] = r
        hdr15 = (f"{'frame':>5} {'V/A^3':>9} {'V/V0':>7} "
                 f"{'dE(5->1)@k11':>13} {'dE(5->1)@k15':>13} {'diff':>9} "
                 f"{'E_s1(k15)-E_s1(k11)':>21}")
        print(hdr15)
        print("-" * len(hdr15))
        d15, d11, dk = [], [], []
        for fr in sorted(byf15):
            g15 = byf15[fr]
            g11 = byframe.get(fr, {})
            if not all(k in g15 for k in ("s1", "sm5")) or "s1" not in g11:
                continue
            a = (g15["sm5"]["energy_sigma0_eV"] - g15["s1"]["energy_sigma0_eV"]) * 1000.0 / n_atoms
            b = (g11["sm5"]["energy_sigma0_eV"] - g11["s1"]["energy_sigma0_eV"]) * 1000.0 / n_atoms
            c = (g15["s1"]["energy_sigma0_eV"] - g11["s1"]["energy_sigma0_eV"]) * 1000.0 / n_atoms
            d15.append(a)
            d11.append(b)
            dk.append(c)
            print(f"{fr:>5} {g15['s1']['volume_angstrom3']:>9.4f} "
                  f"{g15['s1']['volume_scale_vs_v0']:>7.5f} {b:>13.3f} {a:>13.3f} "
                  f"{a - b:>9.3f} {c:>21.3f}")
        if d15:
            summary["k15"] = {"dE_m5_s1": d15, "dE_m5_s1_k11": d11, "kerr_s1": dk,
                              "range_dE_m5_s1_k15": max(d15) - min(d15),
                              "range_dE_m5_s1_k11": max(d11) - min(d11)}
            print(f"  range[dE(-5->s1)]: k=11 {max(d11) - min(d11):.3f} -> "
                  f"k=15 {max(d15) - min(d15):.3f} meV/atom")
            print(f"  |E_s1(k15) - E_s1(k11)| = {min(abs(x) for x in dk):.3f} .. "
                  f"{max(abs(x) for x in dk):.3f} meV/atom (k-mesh error on the same volumes)")
        # EOS at k=15 vs k=11, ISMEAR=1
        print("  EOS fit, ISMEAR=1:")
        eos2 = {}
        for tag, blk in (("k11", "primary"), ("k15", "k15")):
            pts = sorted(((r["volume_angstrom3"], r["energy_sigma0_eV"]) for r in rows
                          if r["block"] == blk and r["case"].endswith("s1")), key=lambda x: x[0])
            if len(pts) >= 4:
                try:
                    f = bm3([p[0] for p in pts], [p[1] for p in pts])
                    eos2[tag] = f
                    print(f"    {tag:>4}: V0={f['V0_angstrom3']:.4f} A^3  "
                          f"E0={f['E0_eV']:.6f} eV  B0={f['B0_GPa']:.2f} GPa  "
                          f"B0'={f['B0prime']:.3f}")
                except Exception as exc:  # noqa: BLE001
                    print(f"    {tag:>4}: fit failed: {exc}")
        if len(eos2) == 2:
            dv = abs(eos2["k15"]["V0_angstrom3"] - eos2["k11"]["V0_angstrom3"])
            db = abs(eos2["k15"]["B0_GPa"] - eos2["k11"]["B0_GPa"])
            de = abs(eos2["k15"]["E0_eV"] - eos2["k11"]["E0_eV"]) * 1000.0 / n_atoms
            print(f"    k11 -> k15: dV0={dv:.4f} A^3 ({100 * dv / eos2['k15']['V0_angstrom3']:.3f}%)  "
                  f"dB0={db:.2f} GPa  dE0={de:.3f} meV/atom")
            summary["eos_k11_vs_k15"] = {"dV0": dv, "dB0": db, "dE0_meV_per_atom": de}

    summary["rows"] = rows

    # ---------- provenance: every case must share POTCAR / ENCUT / EDIFF ----
    val = root / "validation.json"
    kcval = root / "kconv" / "validation.json"
    k15val = root / "k15" / "validation.json"
    prov_cases = []
    for f in (val, kcval, k15val):
        if f.is_file():
            blob = json.loads(f.read_text())
            prov_cases += blob["cases"] if isinstance(blob, dict) else blob
    if prov_cases:
        pot_hashes = {c["potcar_sha256"] for c in prov_cases}
        print(f"\nprovenance: {len(prov_cases)} generated cases; "
              f"unique POTCAR sha256 = {sorted(pot_hashes)}")
        print("  (a single POTCAR hash across all cases is what makes the "
              "ISMEAR and k comparisons legitimate)")
        summary["provenance"] = {
            "n_cases_generated": len(prov_cases),
            "unique_potcar_sha256": sorted(pot_hashes),
            "unique_poscar_per_case": len({c["poscar_sha256"] for c in prov_cases}),
            "n_bad": sum(1 for c in prov_cases if not c["ok"]),
        }
        if len(pot_hashes) != 1:
            print("  !! MORE THAN ONE POTCAR HASH -- the comparison is not clean")

    if args.json:
        args.json.write_text(json.dumps(summary, indent=2) + "\n")
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
