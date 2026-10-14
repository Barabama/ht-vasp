#!/usr/bin/env python
"""Run the VASP input gate over a tree of prepared case directories.

The gate itself lives in `htvasp.dfpt_qha.gate`; this is the CLI the report
documents.  Exit status is 1 if any directory fails, so a submission script can
chain on it.

    gate_cli.py DIR [DIR ...]                 # static checks
    gate_cli.py --root DIR                    # every case under DIR
    gate_cli.py --root DIR --json out.json    # machine-readable verdicts
    gate_cli.py --root DIR --dryrun --scratch /nfs_ssd/tmp/gate

`--dryrun` needs a compute allocation: each case is copied to a scratch
directory and run with NELM=1, NSW=0 so VASP parses the inputs but does not
relax anything.  Cases the static layer already rejected are not dry-run --
there is nothing to learn from asking VASP to parse inputs known to be
malformed.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from htvasp.dfpt_qha.gate import iter_cases, run_gate, validate_dir


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="*", type=Path)
    ap.add_argument("--root", type=Path, default=None,
                    help="validate every directory under here that has an INCAR")
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--dryrun", action="store_true",
                    help="also let VASP parse each case (needs a compute allocation)")
    ap.add_argument("--scratch", type=Path, default=None,
                    help="scratch root for --dryrun copies")
    ap.add_argument("--ranks", type=int, default=4)
    args = ap.parse_args()

    targets: list[Path] = list(args.dirs)
    if args.root:
        targets += iter_cases(args.root)
    if not targets:
        ap.error("give directories or --root")
    if args.dryrun and args.scratch is None:
        ap.error("--dryrun needs --scratch")

    if args.dryrun:
        args.scratch.mkdir(parents=True, exist_ok=True)
        results = run_gate(targets, scratch=args.scratch, np_ranks=args.ranks)
    else:
        results = [validate_dir(d) for d in targets]

    if not args.quiet:
        for result in results:
            extra = ""
            if result.mesh:
                extra = (f" mesh={list(result.mesh)} shift={list(result.shift)} "
                         f"gamma={result.gamma_centered} ISMEAR={result.ismear}")
            print(f"{'OK  ' if result.ok else 'FAIL'} {result.dir}{extra}")
            for problem in result.problems:
                print(f"       !! {problem}")
            for warning in result.warnings:
                print(f"       ~  {warning}")
            if result.dryrun is not None and result.dryrun.nions is not None:
                dry = result.dryrun
                print(f"       VASP read: NIONS={dry.nions} NKPTS={dry.nkpts} "
                      f"ISPIN={dry.ispin}")

    n_bad = sum(1 for result in results if not result.ok)
    if not args.quiet:
        print(f"\n{len(results) - n_bad}/{len(results)} directories pass the contract")
    if args.json:
        args.json.write_text(
            json.dumps([r.model_dump(mode="json") for r in results], indent=2) + "\n"
        )
        print(f"wrote {args.json}")
    if n_bad:
        sys.exit(1)


if __name__ == "__main__":
    main()
