#!/usr/bin/env python
"""Read the converged magnetic state out of the archived DFPT-QHA runs.

The DFPT-QHA analysis JSON carries no magnetic information at all -- no moment,
no BMAGN, no magnetic verdict.  But every run archived its
``1-EOS_equilibrium_relaxation/OUTCAR.gz``, and that file's last
``magnetization (x)`` block *is* the converged per-site moment of the cell the
EOS was built from.  So the magnetic state of the endpoints that already exist
can be recovered without a single new VASP run.

Two things make this worth a script rather than a grep:

* **The moment is not in the JSON, so nothing downstream can see it.**  A TDB
  needs ``BMAGN`` (and a ``TC`` sign) per phase.  This writes
  ``data/endpoints/<name>/magnetic-state.json`` beside the analysis artifact,
  with the source OUTCAR's SHA-256 so the number is traceable, and does not
  touch the immutable run directory.
* **The moment has to be interpreted, not just printed.**  A net moment of zero
  can mean "nonmagnetic" or "antiferromagnetic with the sublattices cancelling",
  and those two are opposite answers for CALPHAD: the first has no magnetic
  Gibbs term, the second has one with a negative TC.  ``mean|m|`` separates
  them; the net moment does not.

The EOS relaxation is the 2-atom (BCC/B2) or 4-atom (FCC/L1_2) cell, so the
per-site list is short and the site order is the recipe's.

    ./.conda/bin/python scripts/extract_magnetic_state.py --pool
    ./.conda/bin/python scripts/extract_magnetic_state.py -n BCC-Cr-Cr-2 --dry-run

For an endpoint whose run directory exists only on the other cluster, pass
``--remote-root`` to read it from a mounted copy, or copy the single
``OUTCAR.gz`` (about 100 kB) into ``data/remote/<cluster>/`` first.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import logging
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
POOL = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_recipes.json")
MAGNETIZATION_TOL = 0.05   # mu_B; below this the cell is treated as nonmagnetic

log = logging.getLogger("magnetic-state")

#: ``magnetization (x)`` followed by the per-ion table.  ``tot`` is the fourth
#: column; s/p/d are the orbital projections and are kept only for the record.
_BLOCK = re.compile(
    r"magnetization \(x\)\s*\n\s*\n"
    r"# of ion\s+s\s+p\s+d\s+tot\s*\n"
    r"-+\n"
    r"((?:\s*\d+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s*\n)+)"
)

def _display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path)

def read_text(path: Path) -> str:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", errors="replace") as handle:
            return handle.read()
    return path.read_text(errors="replace")

def per_site_moments(outcar: Path) -> list[float] | None:
    """The ``tot`` column of the last magnetization block, or None."""
    blocks = _BLOCK.findall(read_text(outcar))
    if not blocks:
        return None
    rows = [line.split() for line in blocks[-1].strip().splitlines()]
    return [float(row[4]) for row in rows]

def initial_moments(outcar: Path) -> str | None:
    incar = outcar.with_name("INCAR.gz")
    if not incar.is_file():
        incar = outcar.with_name("INCAR")
    if not incar.is_file():
        return None
    match = re.search(r"^\s*MAGMOM\s*=\s*(.+)$", read_text(incar), re.M)
    return match.group(1).strip() if match else None

def classify(moments: list[float]) -> dict:
    """Net and absolute averages, plus the magnetic verdict they imply.

    ``nonmagnetic`` needs both averages to vanish: a zero net moment with a
    large mean|m| is antiferromagnetic, and IHJ represents that as a positive
    BMAGN with a negative TC.
    """
    net = sum(moments)
    mean_signed = net / len(moments)
    mean_abs = sum(abs(m) for m in moments) / len(moments)
    if mean_abs < MAGNETIZATION_TOL:
        kind = "nonmagnetic"
    elif abs(mean_signed) < MAGNETIZATION_TOL:
        kind = "antiferromagnetic"
    else:
        kind = "ferromagnetic"
    return {
        "n_sites": len(moments),
        "moments_muB": moments,
        "net_moment_muB": net,
        "mean_signed_muB": mean_signed,
        "mean_abs_muB": mean_abs,
        "kind": kind,
        # IHJ writes BMAGN as the magnitude that enters the magnetic Gibbs term,
        # so an antiferromagnet still carries one; the sign lives on TC.
        "bmagn_ihj_muB": mean_abs,
        "tc_sign": {"nonmagnetic": 0, "ferromagnetic": 1,
                    "antiferromagnetic": -1}[kind],
    }

def eos_outcar(run_dirs: list[Path]) -> Path | None:
    """The newest ``1-EOS_equilibrium_relaxation/OUTCAR.gz`` among the runs."""
    candidates = [d / "OUTCAR.gz" for d in run_dirs]
    candidates = [p for p in candidates if p.is_file()]
    return max(candidates, key=lambda p: p.stat().st_mtime) if candidates else None

def run_dirs(name: str, roots: list[Path]) -> list[Path]:
    found: list[Path] = []
    for root in roots:
        base = root / name / "dfpt-qhaflow" / "runs"
        if base.is_dir():
            found += [p / "1-EOS_equilibrium_relaxation" for p in base.iterdir()
                      if (p / "1-EOS_equilibrium_relaxation").is_dir()]
    return found

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-n", "--name", action="append", default=[],
                    help="endpoint name; repeatable")
    ap.add_argument("--pool", action="store_true",
                    help="every endpoint whose recipe magnetic_policy is not nonmagnetic")
    ap.add_argument("--all", action="store_true",
                    help="every endpoint with a run directory, magnetic or not")
    ap.add_argument("--remote-root", type=Path, action="append", default=[],
                    help="extra endpoints/ tree to read (e.g. a B5070 mirror); "
                         "the local tree is always read first")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    recipes = json.loads(POOL.read_text())["endpoints"]
    names = list(args.name)
    if args.pool:
        names += [n for n, v in recipes.items()
                  if v["recipe"]["magnetic_policy"] != "nonmagnetic"]
    if args.all:
        names += sorted(recipes)
    names = sorted(dict.fromkeys(names))
    if not names:
        ap.error("nothing to do: pass -n, --pool or --all")

    roots = [ROOT / "data" / "endpoints"] + [Path(r).resolve() for r in args.remote_root]
    written = 0
    missing: list[str] = []
    for name in names:
        outcar = eos_outcar(run_dirs(name, roots))
        if outcar is None:
            missing.append(name)
            continue
        moments = per_site_moments(outcar)
        if moments is None:
            missing.append(f"{name} (OUTCAR has no magnetization block)")
            continue
        record = {
            "schema": "htvasp-magnetic-state/v1",
            "endpoint": name,
            # relative to the repository when possible, so the record survives
            # being read on the other cluster or after the tree is moved
            "source_outcar": _display_path(outcar),
            "source_sha256": hashlib.sha256(outcar.read_bytes()).hexdigest(),
            "initial_moments": initial_moments(outcar),
            "magnetic_policy": recipes.get(name, {}).get("recipe", {})
                                     .get("magnetic_policy"),
            **classify(moments),
        }
        policy = record["magnetic_policy"]
        if policy == "nonmagnetic" and record["kind"] != "nonmagnetic":
            record["warning"] = (
                f"the recipe says nonmagnetic but the converged cell carries "
                f"{record['mean_abs_muB']:.3f} mu_B/atom -- either the recipe is "
                "wrong or the run left the intended branch"
            )
        if policy and policy != "nonmagnetic" and record["kind"] == "nonmagnetic":
            record["warning"] = (
                "the recipe expects a magnetic branch but the run converged "
                "nonmagnetic.  This does NOT establish that the run is wrong: the "
                "SCF may have found a genuinely lower nonmagnetic state at this "
                "volume.  What it does establish is that the endpoint's magnetic "
                "state is unresolved, because the only evidence for the magnetic "
                "branch came from a different cell or volume -- so a comparison "
                "at THIS relaxed volume is needed before the curve is used."
            )
        target = ROOT / "data" / "endpoints" / name / "magnetic-state.json"
        log.info("%-14s %-17s net=%+6.3f mean|m|=%.3f %s",
                 name, record["kind"], record["net_moment_muB"],
                 record["mean_abs_muB"], "-> " + record["warning"] if "warning" in record else "")
        if not args.dry_run:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(json.dumps(record, indent=2) + "\n")
            written += 1

    print(f"wrote {written} magnetic-state.json" if not args.dry_run
          else f"dry-run: {len(names) - len(missing)} readable")
    if missing:
        print(f"{len(missing)} without a readable EOS OUTCAR:")
        for item in missing:
            print(f"  {item}")
    return 0

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    raise SystemExit(main())
