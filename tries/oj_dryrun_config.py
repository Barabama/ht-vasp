#!/usr/bin/env python3
"""OJ dry-run config report: predict magnetic-config count & solvability WITHOUT VASP.

This is the actionable pure-dry-run flow for OstravaJ. It runs ONLY the
`generate` stage of OstravaJ (geometry/algebra only — no VASP, no total energies),
and reports the exact numbers that decide whether a real OJ+VASP run is worthwhile:

    N_k       = number of requested exchange interactions (J pairs)
    n_flips   = number of independent magnetic configs OJ would generate
               (== number of flip* dirs == number of VASP static jobs needed)
    full-rank = whether the coefficient matrix reached rank N_k+1 ("Success")

Usage:
    python oj_dryrun_config.py --name FCC-Co-Fe --extend 2 2 2 --jcount 4
    python oj_dryrun_config.py --name BCC-Co-Fe  --extend 4 4 2 --dist 6.0
    python oj_dryrun_config.py --name BCC-Fe-Fe  --extend 4 4 2 --jcount 6

Evidence it is VASP-free:
  * ostravaj/jmixer/cli.py `runGenerate` hardcodes s.energies.append(0) and only
    uses matrix RANK (which does not depend on the RHS energies) for
    isConfigIndependent / full-rank checks.
  * ostravaj/jmixer/generator.py searches 2^N configs using only spin products
    and J-pair geometry; config.py reads only POSCAR+INCAR+OJ.conf (no POTCAR).
"""

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

from pymatgen.core import Lattice, Structure

OJ_DIR = Path("/nfs_hdd/2025/gaominliang/ht-vasp/ostravaj")
PYTHON = "/nfs_hdd/2025/gaominliang/ht-vasp/.conda/bin/python"
CLI = OJ_DIR / "jmixer" / "cli.py"

BASE_SPIN = {"Co": 2.0, "Fe": 2.5, "Mn": 2.5, "Ni": 1.0}

# convenience builders; in production use data/poscars/*.vasp via Endmember
BUILDERS = {
    "BCC-Fe-Fe": lambda: Structure(Lattice.cubic(2.866), ["Fe", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]]),
    "BCC-Co-Fe": lambda: Structure(Lattice.cubic(2.860), ["Co", "Fe"], [[0, 0, 0], [0.5, 0.5, 0.5]]),
    "FCC-Co-Fe": lambda: Structure(Lattice.cubic(3.735), ["Co", "Fe", "Fe", "Fe"],
                                   [[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]]),
    "FCC-Co-Ni": lambda: Structure(Lattice.cubic(3.520), ["Co", "Ni", "Ni", "Ni"],
                                   [[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]]),
}


def magnetic_types(name: str) -> list[str]:
    seen, out = set(), []
    for tok in name.replace("-", " ").split():
        if tok in BASE_SPIN and tok not in seen:
            seen.add(tok)
            out.append(tok)
    return out


def write_input(struct, mag_types, extend, jcount, dist, outdir: Path) -> Path:
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)
    struct.to(fmt="poscar", filename=outdir / "POSCAR")
    (outdir / "INCAR").write_text("ISPIN = 2\n")  # minimal; POTCAR not needed for generate
    jspec = f"J_count {jcount}" if jcount is not None else f"dist_cutoff {dist}"
    spins = " ".join(str(BASE_SPIN.get(t, 2.0)) for t in mag_types)
    (outdir / "OJ.conf").write_text("\n".join([
        jspec,
        f"magnetic_ion_types {' '.join(mag_types)}",
        "noncollinear 0",
        f"base_spin {spins}",
        f"extend_poscar {extend[0]} {extend[1]} {extend[2]}",
    ]) + "\n")
    return outdir


def run_generate(in_dir: Path) -> str:
    env = dict(__import__("os").environ)
    env["PYTHONPATH"] = str(OJ_DIR) + (":" + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    p = subprocess.run([PYTHON, str(CLI), "generate", "-d", "-i", str(in_dir), "-r", str(in_dir)],
                       capture_output=True, text=True, timeout=3600, env=env)
    return p.stdout + p.stderr


def parse(text: str) -> dict:
    info = {}
    m = re.search(r"Found Js:\s*(\[.*\])", text)
    if m:
        info["Js"] = m.group(1)
        info["N_k"] = len(re.findall(r"\('", m.group(1)))
    m = re.search(r"Found (\d+) suitable magnetic configurations", text)
    info["n_suitable"] = int(m.group(1)) if m else None
    m = re.search(r"matrix of J coefficients:\n(\[\[.*?\]\])\n", text, re.DOTALL)
    if m:
        info["n_flips"] = len(re.findall(r"(?<=\[)\s*-?\d", m.group(1)))
    info["full_rank"] = "Success." in text and "Ran out of configs." not in text
    return info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--extend", nargs=3, type=int, required=True)
    ap.add_argument("--jcount", type=int, default=None)
    ap.add_argument("--dist", type=float, default=6.0)
    args = ap.parse_args()

    if args.name not in BUILDERS:
        sys.exit(f"unknown structure {args.name}; builders: {sorted(BUILDERS)}")
    struct = BUILDERS[args.name]()
    mag = magnetic_types(args.name)
    extend = tuple(args.extend)
    work = Path("/tmp") / f"oj_dryrun_{args.name}_{args.jcount or args.dist}"
    in_dir = write_input(struct, mag, extend, args.jcount, args.dist, work)
    text = run_generate(in_dir)
    info = parse(text)
    n_atoms = len(struct) * extend[0] * extend[1] * extend[2]
    odd = [d for d, f in zip("abc", extend) if f % 2 == 1]

    print(f"structure    : {args.name} ({mag})")
    print(f"supercell    : extend={extend} -> {n_atoms} atoms")
    print(f"spec         : {'J_count=' + str(args.jcount) if args.jcount else 'dist_cutoff=' + str(args.dist)}")
    print(f"N_k          : {info.get('N_k')}  Js={info.get('Js')}")
    print(f"n_suitable   : {info.get('n_suitable')}")
    print(f"n_flips      : {info.get('n_flips')}  (== number of VASP static jobs)")
    print(f"full_rank    : {info.get('full_rank')}  (needs {info.get('N_k', 0) + 1} configs)")
    advice = []
    if odd:
        advice.append(f"ODD supercell dimension(s) {odd}: observed collapse of usable configs; use even dimensions")
    if info.get("full_rank") is False and info.get("N_k") is not None:
        advice.append("ran out of configs: reduce J_count (fewer shells) or try another even supercell")
    if n_atoms > 80:
        advice.append(f"{n_atoms} atoms: the config search itself can take many minutes; run dry-run on compute node")
    for a in advice:
        print(f"  -> {a}")
    if not info.get("full_rank") and not info.get("N_k"):
        print(text[-3000:])


if __name__ == "__main__":
    main()
