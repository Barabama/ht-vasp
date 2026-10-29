#!/usr/bin/env python3
"""Closed offline loop: wannier90 replay (both spins) -> TB2J -> J / TC.

Takes an archived ``2-w90_static`` directory, re-runs ``wannier90.x`` on the
archived ``.mmn/.amn/.eig/.chk`` with a chosen parameter set, feeds the
resulting ``_hr.dat`` + ``_centres.xyz`` to TB2J's ``wann2J``, and reports the
exchange alongside the archived values.  No VASP, no Slurm.

The archived ``_hr.dat`` is never touched: each variant gets its own scratch dir.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPLAY = HERE / "w90replay.py"
CONDA = Path("/nfs_hdd/2025/gaominliang/ht-vasp/.conda/bin/python")


def parse_exchange_out(path: Path) -> dict:
    txt = path.read_text(errors="replace")
    pairs = []
    for line in txt.splitlines():
        m = re.match(r"\s*([A-Z][a-z]?\d+)\s+([A-Z][a-z]?\d+)\s+\(\s*(-?\d+),\s*(-?\d+),\s*(-?\d+)\)\s+"
                     r"([-+0-9.eE]+)\s+\(\s*([-+0-9.eE]+),\s*([-+0-9.eE]+),\s*([-+0-9.eE]+)\)\s+([-+0-9.eE]+)", line)
        if m:
            pairs.append({"i": m.group(1), "j": m.group(2),
                          "R": [int(m.group(3)), int(m.group(4)), int(m.group(5))],
                          "J_iso_meV": float(m.group(6)), "dist": float(m.group(10))})
    moments = []
    for line in txt.splitlines():
        m = re.match(r"\s*([A-Z][a-z]?\d+)\s+([-+0-9.eE]+)\s*$", line)
        if m and "atom" not in line.lower():
            moments.append(float(m.group(2)))
    return {"n_pairs": len(pairs), "pairs": pairs, "moments": moments}


def summarize(res: dict) -> dict:
    pairs = res["pairs"]
    if not pairs:
        return {**res, "max_abs_J_meV": None, "mean_abs_J_meV": None}
    js = [abs(p["J_iso_meV"]) for p in pairs]
    top = max(pairs, key=lambda p: abs(p["J_iso_meV"]))
    return {**res, "max_abs_J_meV": max(js), "mean_abs_J_meV": sum(js) / len(js),
            "strongest": {"i": top["i"], "j": top["j"], "R": top["R"],
                          "J_iso_meV": top["J_iso_meV"], "dist": top["dist"]}}


def mfa_tc(summ: dict, n_sites_per_sublattice: dict) -> dict:
    """Crude single-sublattice MFA: k_B T_C = (2/3) * sum_j |J_ij| over one site."""
    if not summ.get("max_abs_J_meV"):
        return {}
    n_sites = max(n_sites_per_sublattice.values() or [1])
    lam = sum(abs(p["J_iso_meV"]) for p in summ["pairs"]) / max(n_sites, 1)
    k_B = 0.08617333262
    return {"lambda_max_meV": lam, "tc_mfa_K": (2.0 / 3.0) * lam / k_B}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("source")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--efermi", type=float, required=True)
    ap.add_argument("--elements", required=True)
    ap.add_argument("--kmesh", default="9 9 9")
    ap.add_argument("--set", action="append", default=[])
    ap.add_argument("--skip-replay", action="store_true")
    ap.add_argument("--baseline-json", default=str(HERE / "baseline_BCC-Nb-Nb.json"))
    args = ap.parse_args()

    src = Path(args.source).resolve()
    work = Path(args.workdir).resolve()
    work.mkdir(parents=True, exist_ok=True)
    w90_stage = work / "w90"

    if not args.skip_replay:
        if w90_stage.exists():
            shutil.rmtree(w90_stage)
        w90_stage.mkdir(parents=True)
        poscar = src / "POSCAR"
        if not poscar.exists():
            import gzip
            (w90_stage / "POSCAR").write_text(gzip.open(src / "POSCAR.gz", "rt").read())
        else:
            shutil.copy(poscar, w90_stage / "POSCAR")
        for spin in (1, 2):
            cmd = [str(CONDA), "-I", str(REPLAY), str(src), "--spin", str(spin),
                   "--workdir", str(w90_stage / f"s{spin}"), "--keep"]
            for item in args.set:
                cmd += ["--set", item]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0 and not (w90_stage / f"s{spin}" / f"wannier90.{spin}_hr.dat").exists():
                print(f"spin {spin} replay failed:\n{proc.stdout[-1500:]}\n{proc.stderr[-800:]}")
                return 1

    tb2j_dir = work / "tb2j"
    if tb2j_dir.exists():
        shutil.rmtree(tb2j_dir)
    tb2j_dir.mkdir()
    shutil.copy(w90_stage / "POSCAR", tb2j_dir / "POSCAR")
    for spin in (1, 2):
        for name in (f"wannier90.{spin}_hr.dat", f"wannier90.{spin}_centres.xyz"):
            p = w90_stage / f"s{spin}" / name
            if p.exists():
                (tb2j_dir / name).symlink_to(p)

    code = (
        "import sys; sys.path.insert(0, %r); import nc_stub\n"
        "sys.argv = ['wann2J','--path','.','--posfile','POSCAR','--efermi',%r,"
        "'--kmesh',%r,%r,%r,'--elements',*%r,'--prefix_up','wannier90.1','--prefix_down','wannier90.2',"
        "'--output_path','out']\n"
        "from TB2J.scripts.wann2J import run_wann2J\nrun_wann2J()\n"
    ) % (str(HERE), str(args.efermi), *args.kmesh.split(), tuple(args.elements.split(",")))
    env = dict(os.environ, PYTHONPATH=str(HERE))
    proc = subprocess.run([str(CONDA), "-I", "-c", code], cwd=tb2j_dir,
                          capture_output=True, text=True, env=env, timeout=7200)
    (tb2j_dir / "tb2j.log").write_text(proc.stdout + proc.stderr)
    ex = tb2j_dir / "out" / "exchange.out"
    if not ex.exists():
        print("TB2J produced no exchange.out:\n" + (proc.stdout + proc.stderr)[-2000:])
        return 1

    summary = summarize(parse_exchange_out(ex))
    out = {"source": str(src), "overrides": args.set,
           "converged": {},
           "n_pairs": summary["n_pairs"],
           "max_abs_J_meV": summary["max_abs_J_meV"],
           "mean_abs_J_meV": summary["mean_abs_J_meV"],
           "strongest": summary.get("strongest"),
           "moments": summary["moments"]}
    for spin in (1, 2):
        wp = w90_stage / f"s{spin}" / f"wannier90.{spin}.wout"
        if wp.exists():
            t = wp.read_text(errors="replace")
            dis = [l.split() for l in t.splitlines()
                   if "<-- DIS" in l and l.strip() and l.strip()[0].isdigit()]
            out["converged"][f"spin{spin}"] = {
                "satisfied": "criteria satisfied" in t and "not satisfied" not in t,
                "iters": int(dis[-1][0]) if dis else None,
                "last_delta": float(dis[-1][3]) if dis else None,
            }
            m = re.search(r"Final Omega_I\s+([-\d.]+)", t)
            if m:
                out["converged"][f"spin{spin}"]["omega_i"] = float(m.group(1))
    out.update(mfa_tc(summary, {"0": 1}))
    print(json.dumps(out, indent=1, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
