#!/usr/bin/env python
"""Compare Stage 0B phonons at d=12 (supercell mesh 6) and d=16 (mesh 8).

Stage 0B originally ran at mesh 6 -- density 12 -- while the E(V) path uses
density 16.  build_stage0b_d16.py supplies the density-16 partner.  This script
puts the two side by side so the k-density dependence of the phonon quantities
is a measured number rather than an assumption.

Both batches are the SAME cells, displacements, smearings and amplitude; only
the supercell mesh differs.  So any difference here is purely the k-mesh.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from analyze_stage0b import build_fc, phonons, zpe, read_forces
from analyze_align import parse_outcar
from analyze_stage0 import parse_site_moments

SMEARINGS = ("s1_010", "s1_020", "s2_010", "s2_020")
LABEL = {"s1_010": "1/0.10", "s1_020": "1/0.20", "s2_010": "2/0.10", "s2_020": "2/0.20"}
ENDPOINTS = ("BCC-Ta-Ta-2", "BCC-Al-Ta-2", "BCC-Al-Cr-2")
# the soft mode is the quantity most likely to move with mesh
SOFT_Q = {"BCC-Al-Ta-2": "N", "BCC-Ta-Ta-2": "Gamma", "BCC-Al-Cr-2": "Gamma"}


def collect(root: Path, tag_dir: str) -> dict:
    out: dict = {}
    base = root / tag_dir
    for ep in ENDPOINTS:
        for tag in SMEARINGS:
            dirs = sorted(d for d in (base / ep).glob(f"a0.01_d*_{tag}")
                          if (d / "OUTCAR").is_file())
            if not dirs:
                continue
            info = build_fc(ep, dirs, Path(f"/tmp/fcdiff/{tag_dir}/{ep}/{tag}"))
            if "error" in info:
                out[(ep, tag)] = {"error": info["error"]}
                continue
            fc = np.load(info["fc_file"])
            r = phonons(ep, fc)
            r["zpe_kJ_per_mol"] = zpe(fc, ep)
            r.update({"max_abs_Phi": info["max_abs_Phi"],
                      "asr_rel": info["asr_rel"],
                      "permutation": info["permutation_residual"],
                      "n_disp": info["n_displacements"]})
            out[(ep, tag)] = r
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("root", type=Path)
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()

    d12 = collect(args.root, "stage0b")
    d16 = collect(args.root, "stage0b-d16")

    print("=== 力常数与声子：d=12（超胞 mesh 6）vs d=16（mesh 8）===\n")
    summary: dict = {}
    for ep in ENDPOINTS:
        print(f"--- {ep}")
        hdr = (f"    {'smearing':>8} {'max|Phi| d12':>13} {'d16':>13} {'diff':>8} "
               f"{'ASR_rel d16':>12} {'perm d16':>10} "
               f"{'softmode d12':>13} {'d16':>9} {'diff':>8} {'neg d12/d16':>12}")
        print(hdr)
        for tag in SMEARINGS:
            a = d12.get((ep, tag))
            b = d16.get((ep, tag))
            if not a or not b or "error" in a or "error" in b:
                print(f"    {LABEL[tag]:>8}   (missing: "
                      f"{a.get('error') if a and 'error' in a else ''} "
                      f"{b.get('error') if b and 'error' in b else ''})")
                continue
            q = SOFT_Q[ep]
            sa = min(a[q])
            sb = min(b[q])
            dp = 100 * (b["max_abs_Phi"] / a["max_abs_Phi"] - 1)
            print(f"    {LABEL[tag]:>8} {a['max_abs_Phi']:>13.5f} "
                  f"{b['max_abs_Phi']:>13.5f} {dp:>7.2f}% "
                  f"{b['asr_rel']:>12.2e} {b['permutation']:>10.2e} "
                  f"{sa:>13.4f} {sb:>9.4f} {sb - sa:>8.4f} "
                  f"{a['mesh_n_negative']:>5}/{b['mesh_n_negative']:<5}")
            summary.setdefault(ep, {})[tag] = {
                "max_abs_Phi_d12": a["max_abs_Phi"], "max_abs_Phi_d16": b["max_abs_Phi"],
                "dPhi_percent": dp,
                "soft_q": q, "soft_d12_THz": sa, "soft_d16_THz": sb,
                "soft_delta_THz": sb - sa,
                "mesh_min_d12": a["mesh_min_THz"], "mesh_min_d16": b["mesh_min_THz"],
                "n_neg_d12": a["mesh_n_negative"], "n_neg_d16": b["mesh_n_negative"],
                "zpe_d12": a["zpe_kJ_per_mol"], "zpe_d16": b["zpe_kJ_per_mol"],
                "gamma_acoustic_d12": min(a["Gamma"]),
                "gamma_acoustic_d16": min(b["Gamma"]),
            }
        print()

    print("=== 分类是否随密度改变 ===")
    for ep in ENDPOINTS:
        for tag in SMEARINGS:
            r = summary.get(ep, {}).get(tag)
            if not r:
                continue
            c12 = "stable" if r["n_neg_d12"] == 0 else "soft"
            c16 = "stable" if r["n_neg_d16"] == 0 else "soft"
            flag = "OK" if c12 == c16 else "**CHANGED**"
            print(f"  {ep:<15} {LABEL[tag]:>8}: d12 {c12:<6} d16 {c16:<6} {flag}")

    # E1 是同一个胞、同一个位移，只有 k 网格不同（5x5x5，等价 d=10）
    E1_N = -3.6327      # report 4.1, e1-fdisp-al-ta, verified bit-for-bit in 4.11.7
    print("\n=== AlTa 的 N 点软模随 k 密度（同一胞、同一 0.01 Å 位移）===")
    print(f"    {'来源':<28} {'网格':>10} {'等价密度':>9} {'N(1/2,1/2,0) / THz':>20}")
    print(f"    {'E1（独立交付结果）':<28} {'5x5x5':>10} {'d=10':>9} {E1_N:>20.4f}")
    for tag in SMEARINGS:
        a = d12.get((ep := "BCC-Al-Ta-2", tag))
        if a and "error" not in a:
            print(f"    {'Stage 0B（ISMEAR=' + LABEL[tag] + '）':<28} {'6x6x6':>10} "
                  f"{'d=12':>9} {min(a['N']):>20.4f}")
    for tag in SMEARINGS:
        b = d16.get(("BCC-Al-Ta-2", tag))
        if b and "error" not in b:
            print(f"    {'Stage 0B d16（ISMEAR=' + LABEL[tag] + '）':<28} {'8x8x8':>10} "
                  f"{'d=16':>9} {min(b['N']):>20.4f}")
    vals = [E1_N] + [min(d12[("BCC-Al-Ta-2", t)]["N"]) for t in SMEARINGS
                     if ("BCC-Al-Ta-2", t) in d12 and "error" not in d12[("BCC-Al-Ta-2", t)]] \
                  + [min(d16[("BCC-Al-Ta-2", t)]["N"]) for t in SMEARINGS
                     if ("BCC-Al-Ta-2", t) in d16 and "error" not in d16[("BCC-Al-Ta-2", t)]]
    print(f"\n    → 全部 9 个取值落在 {min(vals):.4f} … {max(vals):.4f} THz，"
          f"跨度 {max(vals) - min(vals):.4f} THz ({100 * (max(vals) - min(vals)) / abs(np.mean(vals)):.1f}%)")
    print("    → 三者非单调：d=10 给 %.4f，d=12 给 %.4f，d=16 给 %.4f"
          % (E1_N, min(d12[("BCC-Al-Ta-2", "s1_020")]["N"]),
             min(d16[("BCC-Al-Ta-2", "s1_020")]["N"])))
    print("      **d=12 是离群的那个**；声子层的 k 收敛仍未做。")

    print("\n=== ZPE 的密度依赖 ===")
    for ep in ENDPOINTS:
        for tag in SMEARINGS:
            r = summary.get(ep, {}).get(tag)
            if not r:
                continue
            d = 100 * (r["zpe_d16"] / r["zpe_d12"] - 1)
            print(f"  {ep:<15} {LABEL[tag]:>8}: {r['zpe_d12']:.5f} -> "
                  f"{r['zpe_d16']:.5f} kJ/mol  ({d:+.2f}%)")

    if args.json:
        args.json.write_text(json.dumps(summary, indent=2, default=str) + "\n")
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
