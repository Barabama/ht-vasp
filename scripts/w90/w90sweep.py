#!/usr/bin/env python3
"""Apply one disentanglement configuration to every archived W90 endpoint.

Re-runs ``wannier90.x`` offline on the archived ``mmn/amn/eig/chk`` (via
``w90replay.py``) for both spins, so a parameter set can be scored across the
whole pool without any VASP time.

Configurations
--------------
``BC``  Fixed and endpoint-independent: ``dis_num_iter = 2000`` plus
        ``dis_conv_tol = 1e-6``.  This is the only combination that converged
        on BCC-Nb-Nb in the 2026-10-28 probe.
``A2``  Per-endpoint windows, because an E_F-relative window has to be
        computed per endpoint::

            dis_win_min  = E_F - 25
            dis_win_max  = E_F + 45
            dis_froz_min = E_F - 10
            dis_froz_max = min_k E(band=num_wann) - 1e-5   (the auto value)

        The bare scheme A (``dis_win = E_F +/- 15``) is not offered: it cannot
        hold ``num_wann`` bands on any of the 44 endpoints and wannier90 aborts
        with "Energy window contains fewer states than number of target WFs".

Usage::

    python w90sweep.py BC --workdir /nfs_ssd/tmp/w90sweep --workers 16
    python w90sweep.py A2 --workdir /nfs_ssd/tmp/w90sweep --workers 16

Results are written incrementally to ``<workdir>/sweep-<config>.jsonl`` and a
channel that already has a result is skipped, so an interrupted sweep resumes.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPLAY = HERE / "w90replay.py"
CONDA = Path("/nfs_hdd/2025/gaominliang/ht-vasp/.conda/bin/python")
ROOT = Path("/nfs_hdd/2025/gaominliang/ht-vasp")
ENDPOINTS = ROOT / "data" / "endmembers"


def read_efermi(static: Path) -> float | None:
    """Last ``E-fermi`` in OUTCAR(.gz) -- the converged SCF value."""
    for name, opener in (("OUTCAR.gz", gzip.open), ("OUTCAR", open)):
        p = static / name
        if not p.exists():
            continue
        last = None
        with opener(p, "rt", errors="replace") as fh:
            for line in fh:
                m = re.search(r"E-fermi\s*:\s*([-\d.]+)", line)
                if m:
                    last = float(m.group(1))
        return last
    return None


def read_num_wann(static: Path) -> int | None:
    p = static / "wannier90.1.wout"
    if not p.exists():
        return None
    m = re.search(r"Number of Wannier Functions\s*:\s*(\d+)", p.read_text(errors="replace"))
    return int(m.group(1)) if m else None


def read_auto_froz_max(static: Path, num_wann: int) -> float | None:
    """min over k of E(band = num_wann), minus 1e-5 -- VASP's ESTIMATE_AUTO_WINDOW."""
    p = static / "wannier90.1.eig"
    if not p.exists():
        return None
    bands: dict[int, list[float]] = collections.defaultdict(list)
    for line in p.read_text(errors="replace").splitlines():
        parts = line.split()
        if len(parts) != 3:
            continue
        bands[int(parts[0])].append(float(parts[2]))
    if num_wann not in bands:
        return None
    return min(bands[num_wann]) - 1e-5


def config_params(config: str, static: Path) -> list[str] | None:
    if config == "BC":
        return ["dis_num_iter=2000", "dis_conv_tol=1e-6"]
    if config == "A2":
        efermi = read_efermi(static)
        num_wann = read_num_wann(static)
        auto = read_auto_froz_max(static, num_wann) if num_wann else None
        if efermi is None or auto is None:
            return None
        return [
            f"dis_win_min={efermi - 25:.4f}",
            f"dis_win_max={efermi + 45:.4f}",
            f"dis_froz_min={efermi - 10:.4f}",
            f"dis_froz_max={auto:.6f}",
            "dis_num_iter=1000",
        ]
    raise SystemExit(f"unknown config {config!r}")


def run_channel(config: str, endpoint: str, spin: int, workdir: Path) -> dict:
    static = ENDPOINTS / endpoint / "tb2jflow" / "2-w90_static"
    wd = workdir / f"{config}-{endpoint}-s{spin}"
    out = {"config": config, "endpoint": endpoint, "spin": spin, "ok": False}
    params = config_params(config, static)
    if params is None:
        out["error"] = "missing E_F or num_wann/eig"
        return out
    out["params"] = params
    cmd = [str(CONDA), "-I", str(REPLAY), str(static), "--spin", str(spin), "--workdir", str(wd)]
    for p in params:
        cmd += ["--set", p]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
    except subprocess.TimeoutExpired:
        out["error"] = "timeout"
        return out
    text = proc.stdout
    brace = text.find("{")
    if brace < 0:
        out["error"] = (proc.stderr or text)[-200:]
        return out
    try:
        d = json.loads(text[brace:])
    except json.JSONDecodeError:
        out["error"] = "unparsable replay output"
        return out
    out.update(
        converged=d.get("converged"),
        iters=d.get("last_iter"),
        last_delta=d.get("last_delta"),
        omega_i=d.get("final_omega_i"),
        hr=d.get("hr_written"),
    )
    out["ok"] = d.get("converged") is not None
    if d.get("error"):
        out["error"] = d["error"]
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("config", choices=["BC", "A2"])
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--endpoints", nargs="*", default=None)
    args = ap.parse_args()

    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    ledger = workdir / f"sweep-{args.config}.jsonl"

    endpoints = args.endpoints or sorted(
        p.parent.parent.name for p in ENDPOINTS.glob("*/tb2jflow/2-w90_static")
    )
    done = set()
    if ledger.exists():
        for line in ledger.read_text().splitlines():
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if rec.get("ok"):
                done.add((rec["endpoint"], rec["spin"]))
    todo = [(e, s) for e in endpoints for s in (1, 2) if (e, s) not in done]
    print(f"{args.config}: {len(endpoints)} endpoints, {len(todo)} channels to run "
          f"({len(done)} already done)", flush=True)

    lock_fh = open(ledger, "a", buffering=1)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(run_channel, args.config, e, s, workdir) for e, s in todo]
        for i, fut in enumerate(futures, 1):
            rec = fut.result()
            lock_fh.write(json.dumps(rec, ensure_ascii=False) + "\n")
            mark = "OK " if rec.get("converged") else ("-- " if rec.get("converged") is False else "ERR")
            print(f"[{i}/{len(todo)}] {mark} {rec['endpoint']:14s} s{rec['spin']} "
                  f"iters={rec.get('iters')} omega={rec.get('omega_i')} "
                  f"{rec.get('error', '')}", flush=True)
    lock_fh.close()

    rows = [json.loads(line) for line in ledger.read_text().splitlines() if line.strip()]
    ok = [r for r in rows if r.get("converged") is True]
    print(f"\n{args.config}: {len(ok)}/{len(rows)} channels converged")
    per_endpoint: dict[str, list[bool]] = collections.defaultdict(list)
    for r in rows:
        per_endpoint[r["endpoint"]].append(r.get("converged") is True)
    both = [e for e, v in per_endpoint.items() if len(v) == 2 and all(v)]
    print(f"endpoints with BOTH spins converged: {len(both)}/{len(per_endpoint)}")
    print(", ".join(sorted(both)) if both else "(none)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
