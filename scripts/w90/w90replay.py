#!/usr/bin/env python3
"""Replay wannier90.x offline on an archived VASP Wannier90 run.

Reads mmn/amn/eig/chk from an archived ``2-w90_static`` directory, rebuilds the
win file (cell, kpoints and mp_grid come from the .chk), applies parameter
overrides, and runs wannier90.x standalone.  No VASP, no Slurm.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import struct
import subprocess
import sys
from pathlib import Path

W90 = "/nfs_ssd/softwares/wannier90/bin/wannier90.x"
MKL = "/opt/intel/oneapi/2025.3/lib"


def _records(data: bytes):
    off, out = 0, []
    while off + 4 <= len(data):
        n = struct.unpack("<i", data[off:off + 4])[0]
        if n < 0 or off + 8 + n > len(data):
            break
        out.append((off + 4, n))
        off += 8 + n
    return out


def parse_chk(path: Path) -> dict:
    b = path.read_bytes()
    recs = _records(b)
    if len(recs) < 10:
        raise ValueError(f"{path}: only {len(recs)} records, not a wannier90 .chk")
    o_num_bands, _ = recs[1]
    num_bands = struct.unpack("<i", b[o_num_bands:o_num_bands + 4])[0]
    o_real, _ = recs[4]
    # The .chk stores real_lattice in Fortran column-major order
    # (``((real_lattice(i,j), i=1,3), j=1,3)``), so the 9 doubles are
    # (a1_x, a2_x, a3_x, a1_y, ...).  Reading them row-major silently
    # transposes the cell: for a cubic cell that is invisible, for HCP/FCC it
    # makes wannier90 rebuild a different k-mesh neighbourhood and abort with
    # "mmn has not the right number of nearest neighbours".
    real_col = struct.unpack("<9d", b[o_real:o_real + 72])
    real = [real_col[r + 3 * c] for r in range(3) for c in range(3)]
    o_mp, _ = recs[7]
    mp_grid = struct.unpack("<3i", b[o_mp:o_mp + 12])
    o_kpt, n_kpt = recs[8]
    num_kpts = n_kpt // 24
    kpt = struct.unpack("<%dd" % (3 * num_kpts), b[o_kpt:o_kpt + 24 * num_kpts])
    return {
        "num_bands": num_bands,
        "num_kpts": num_kpts,
        "mp_grid": mp_grid,
        "real_lattice": [list(real[i * 3:(i + 1) * 3]) for i in range(3)],
        "kpoints": [list(kpt[i * 3:(i + 1) * 3]) for i in range(num_kpts)],
    }


def read_wout(path: Path) -> dict:
    if not path.exists():
        return {"converged": None, "error": "no .wout produced"}
    txt = path.read_text(errors="replace")
    if "Disentanglement convergence criteria satisfied" in txt:
        converged = True
    elif "Disentanglement convergence criteria not satisfied" in txt:
        converged = False
    else:
        converged = None
    dis = [ln.split() for ln in txt.splitlines()
           if "<-- DIS" in ln and ln.strip() and ln.strip()[0].isdigit()]
    out: dict = {"converged": converged, "n_dis_rows": len(dis)}
    if dis:
        out["last_iter"] = int(dis[-1][0])
        out["last_delta"] = float(dis[-1][3])
        out["omega_i"] = [float(dis[-1][1]), float(dis[-1][2])]
    m = re.search(r"Final Omega_I\s+([-\d.]+)", txt)
    if m:
        out["final_omega_i"] = float(m.group(1))
    mo = re.search(r"Outer:\s*([-\d.]+)\s*to\s*([-\d.]+)", txt)
    mi = re.search(r"Inner:\s*([-\d.]+)\s*to\s*([-\d.]+)", txt)
    if mo:
        out["outer_window"] = [float(mo.group(1)), float(mo.group(2))]
    if mi:
        out["inner_window"] = [float(mi.group(1)), float(mi.group(2))]
    return out


def parse_poscar(path: Path) -> tuple[list[str], list[list[float]]]:
    """Return (species_per_site, fractional coordinates) from a POSCAR/CONTCAR."""

    lines = path.read_text().splitlines()
    scale = float(lines[1].split()[0])
    lattice = [[float(x) * scale for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    symbols = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    start = 7
    if lines[start].strip().lower().startswith("s"):
        start += 1
    if lines[start].strip()[0].lower() == "d":
        start += 1
    coords = []
    total = sum(counts)
    for line in lines[start:start + total]:
        coords.append([float(x) for x in line.split()[:3]])
    species = [s for s, n in zip(symbols, counts) for _ in range(n)]
    # POSCAR cartesian coordinates are in Angstrom * scale; win wants fractions
    if "c" in lines[start - 1].strip().lower()[:1] or lines[start - 1].strip().lower().startswith("k"):
        import numpy as np
        A = np.array(lattice)
        coords = [list(np.linalg.solve(A.T, np.array(c))) for c in coords]
    return species, coords


def _check_cell_matches_poscar(real_lattice: list[list[float]], poscar: Path, seed: str) -> None:
    """Abort if the .chk cell disagrees with the POSCAR the run used."""

    if not poscar.exists():
        return
    lines = poscar.read_text().splitlines()
    scale = float(lines[1].split()[0])
    pos = [[float(x) * scale for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    worst = max(abs(a - b) for ra, rb in zip(real_lattice, pos) for a, b in zip(ra, rb))
    if worst > 1e-4:
        raise SystemExit(
            f"{seed}.chk cell disagrees with POSCAR by {worst:.2e} Angstrom; "
            "refusing to replay on a possibly-corrupted checkpoint"
        )


def build_win(chk: dict, num_wann: int, projections: str, overrides: dict,
              species: list[str] | None = None,
              fracs: list[list[float]] | None = None) -> str:
    lines = [f"num_wann = {num_wann}", f"num_bands = {chk['num_bands']}",
             f"mp_grid = {chk['mp_grid'][0]} {chk['mp_grid'][1]} {chk['mp_grid'][2]}"]
    for k, v in overrides.items():
        lines.append(f"{k} = {v}")
    lines += ["begin unit_cell_cart", "Ang"]
    lines += ["  %.12f %.12f %.12f" % tuple(r) for r in chk["real_lattice"]]
    lines += ["end unit_cell_cart"]
    if species and fracs:
        lines += ["begin atoms_frac"]
        lines += ["  %s %.12f %.12f %.12f" % (el, *f) for el, f in zip(species, fracs)]
        lines += ["end atoms_frac"]
    lines += ["begin kpoints"]
    lines += ["  %.12f %.12f %.12f" % tuple(k) for k in chk["kpoints"]]
    lines += ["end kpoints", "begin projections", projections, "end projections"]
    for key, value in (("write_hr", ".true."), ("write_xyz", ".true.")):
        if key not in overrides:
            lines.append(f"{key} = {value}")
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("source", help="archived 2-w90_static directory")
    ap.add_argument("--spin", type=int, default=1)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--num-wann", type=int, default=None)
    ap.add_argument("--projections", default=None)
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="win keyword override, repeatable")
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args()

    src = Path(args.source).resolve()
    seed = f"wannier90.{args.spin}"
    work = Path(args.workdir)
    if work.exists() and not args.keep:
        shutil.rmtree(work)
    work.mkdir(parents=True, exist_ok=True)

    # Only .mmn/.amn/.eig are inputs, and wannier90 treats them as read-only, so
    # symlinking them is safe.
    #
    # The .chk is deliberately NOT provided.  wannier90 opens it with
    # STATUS='REPLACE' and writes its own state back at the end of every run, so
    # a symlink writes straight through into the archive -- that is how five
    # archived .chk files were overwritten on 2026-10-28.  It is also
    # unnecessary: the disentanglement starts from the projections in .amn, and
    # replaying without a .chk reproduces the archived Omega_I exactly
    # (SER-Ni 4.88224762, BCC-Nb-Nb 12.29390697, BCC-Al-Al 12.23976667).
    # The cell and k-mesh that used to be read from .chk now come from the
    # POSCAR and a regenerated Gamma-centred mp_grid.
    for ext in ("mmn", "amn", "eig"):
        p = src / f"{seed}.{ext}"
        if not p.exists():
            raise SystemExit(f"missing {p}")
        link = work / f"{seed}.{ext}"
        if not link.exists():
            os.symlink(p, link)

    chk = parse_chk(src / f"{seed}.chk")
    # The .chk supplies the k-point ORDER (wannier90 matches .mmn row by row),
    # so it has to be read.  Cross-check its cell against the POSCAR: the .chk
    # is a file wannier90 rewrites, and a stale or corrupted one would silently
    # produce a different k-mesh neighbourhood.
    _check_cell_matches_poscar(chk["real_lattice"], src / "POSCAR", seed)

    num_wann = args.num_wann
    if num_wann is None:
        old = (src / f"{seed}.win").read_text().splitlines()
        num_wann = next(int(l.split("=")[1]) for l in old if l.startswith("num_wann"))
    projections = args.projections
    if projections is None:
        old = (src / f"{seed}.win").read_text().splitlines()
        i = next(n for n, l in enumerate(old) if l.strip() == "begin projections")
        j = next(n for n, l in enumerate(old) if l.strip() == "end projections")
        projections = "\n".join(old[i + 1:j])

    # Carry over every dis_* / guiding / iter setting the archived run used,
    # then let --set override them.  Without this the replay silently falls
    # back to wannier90's defaults (dis_win = full eigenvalue range), which is
    # NOT the archived configuration and produces a different DIS trajectory.
    overrides = {}
    for line in (src / f"{seed}.win").read_text().splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "!", "begin", "end")):
            continue
        k, sep, v = line.partition("=")
        if not sep:
            continue
        k = k.strip()
        if k in {"num_wann", "num_bands", "mp_grid", "write_hr", "write_xyz"}:
            continue
        overrides[k] = v.strip()
    overrides.setdefault("write_hr", ".true.")
    overrides.setdefault("write_xyz", ".true.")
    for item in args.set:
        k, _, v = item.partition("=")
        overrides[k.strip()] = v.strip()

    poscar = src / "POSCAR"
    if not poscar.exists():
        poscar = src / "POSCAR.gz"
        if poscar.exists():
            import gzip
            (work / "POSCAR").write_text(gzip.open(poscar, "rt").read())
            poscar = work / "POSCAR"
    species, fracs = parse_poscar(poscar) if poscar.exists() else ([], [])
    win = build_win(chk, num_wann, projections, overrides, species, fracs)
    (work / f"{seed}.win").write_text(win)

    env = dict(os.environ)
    env["LD_LIBRARY_PATH"] = MKL + (":" + env["LD_LIBRARY_PATH"] if env.get("LD_LIBRARY_PATH") else "")
    proc = subprocess.run([W90, seed], cwd=work, capture_output=True, text=True, env=env, timeout=7200)
    (work / "stdout.log").write_text(proc.stdout + proc.stderr)

    result = {
        "source": str(src), "spin": args.spin, "num_wann": num_wann,
        "overrides": overrides, "returncode": proc.returncode,
        "hr_written": (work / f"{seed}_hr.dat").exists(),
        **read_wout(work / f"{seed}.wout"),
    }
    if proc.returncode != 0 and result.get("converged") is None:
        tail = (proc.stdout + proc.stderr).strip().splitlines()[-4:]
        result["error"] = "\n".join(tail)
    print(json.dumps(result, indent=1, ensure_ascii=False))
    return 0 if result.get("converged") is not None else 1


if __name__ == "__main__":
    sys.exit(main())
