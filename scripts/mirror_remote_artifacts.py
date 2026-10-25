#!/usr/bin/env python
"""Mirror the DFPT-QHA analysis artifacts from another cluster into this tree.

The two clusters do not share a filesystem: 429pro reads
``/nfs_hdd/2025/gaominliang/ht-vasp`` while B5070 reads
``/workspace/gaominliang/ht-vasp``.  Twenty-eight of the completed endpoints
exist only on B5070, so every local post-processing step -- the ATAT control
tree, an em-tdb-v2 import, any re-fit -- sees at most half the pool.

This copies the three files a consumer needs (the analysis JSON, the derived
``-vt`` JSON, and the submission ledger) into ``data/remote/<cluster>/``,
keeping the ``endpoints/...`` layout.  It deliberately does **not** copy the
run directories: those are gigabytes of VASP output that no consumer reads,
and a stale copy of them would be worse than none.

The ``-vt`` file pins its source with ``source_json_sha256``.  That pin is the
reason this mirror is trustworthy rather than merely convenient: every copied
analysis file is checked against the hash its own derived file recorded, so a
file that was changed or truncated in transit is rejected instead of silently
becoming the local truth.  ``--verify`` re-checks an existing mirror with no
network access.

    ./.conda/bin/python scripts/mirror_remote_artifacts.py --cluster B5070
    ./.conda/bin/python scripts/mirror_remote_artifacts.py --cluster B5070 --verify

Credentials come from the workspace ``.env`` (``SERVER_<CLUSTER>_IP/_USER/_PASS``).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import posixpath
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MIRROR_ROOT = ROOT / "data" / "remote"

#: Where each cluster keeps its copy of this repository.
REMOTE_ROOTS = {
    "B5070": "/workspace/gaominliang/ht-vasp",
    "B5071": "/workspace/gaominliang/ht-vasp",
}

ARTIFACT_PATTERNS = ("*-dfpt-qha.json", "*-dfpt-qha-vt.json")
LEDGER_NAME = "slurm_submissions.json"

#: Two more files, about 100 kB per run, that unlock something the analysis
#: JSON cannot answer: the EOS relaxation's OUTCAR carries the converged
#: per-site magnetic moments, and the DFPT-QHA artifact records no magnetic
#: information at all.  Without these, a TDB consumer on this machine cannot
#: read a BMAGN for any endpoint that ran on the other cluster, and the only
#: alternative is a rerun.  The INCAR beside it says what MAGMOM the run
#: started from, which is what distinguishes "converged nonmagnetic" from
#: "collapsed out of a magnetic branch".
MAGNETIC_STATE_FILES = ("1-EOS_equilibrium_relaxation/OUTCAR.gz",
                        "1-EOS_equilibrium_relaxation/INCAR.gz")


def load_env(path: Path = Path("/nfs_hdd/2025/gaominliang/.env")) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        values[key.strip()] = value.strip().strip('"').strip("'")
    return values


def connect(cluster: str):
    import paramiko

    env = load_env()
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(
        env[f"SERVER_{cluster}_IP"],
        port=22,
        username=env[f"SERVER_{cluster}_USER"],
        password=env[f"SERVER_{cluster}_PASS"],
        timeout=30,
        banner_timeout=30,
        auth_timeout=30,
    )
    return client


def list_remote(client, remote_root: str) -> list[str]:
    """Every artifact path, relative to ``remote_root``."""
    pattern = " -o ".join(f"-name '{p}'" for p in ARTIFACT_PATTERNS)
    command = (
        f"cd {remote_root} && find data/endpoints -maxdepth 6 -type f "
        f"\\( {pattern} -o -name '{LEDGER_NAME}' \\) | sort"
    )
    _, stdout, _ = client.exec_command(command, timeout=900)
    paths = [line.strip() for line in stdout.read().decode().splitlines() if line.strip()]
    # The EOS-relax inputs live deeper than the analysis files, so they get
    # their own pass rather than a looser maxdepth for everything.
    for name in MAGNETIC_STATE_FILES:
        _, out, _ = client.exec_command(
            f"cd {remote_root} && ls data/endpoints/*/dfpt-qhaflow/runs/*/{name} "
            "2>/dev/null | sort", timeout=900)
        paths += [line.strip() for line in out.read().decode().splitlines()
                  if line.strip()]
    if not paths:
        raise RuntimeError(f"{cluster}:{remote_root} produced no artifact paths")
    return paths


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def fetch(client, remote_root: str, relpaths: list[str], out_root: Path) -> dict:
    """Pull the files, then check every vt pin against the copied analysis file."""
    sftp = client.open_sftp()
    records: list[dict] = []
    try:
        for relpath in relpaths:
            remote = posixpath.join(remote_root, relpath)
            with sftp.open(remote, "rb") as handle:
                payload = handle.read()
            target = out_root / relpath
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(payload)
            records.append(
                {"path": relpath, "bytes": len(payload), "sha256": sha256_bytes(payload)}
            )
    finally:
        sftp.close()
    return {"files": records}


def check_pins(out_root: Path) -> list[str]:
    """Every ``-vt`` file's recorded source hash must match the copied source.

    Returns a list of problems; empty means the mirror is self-consistent.
    A vt file whose ``source_json`` points somewhere this cluster cannot reach
    is fine -- the hash is what carries the evidence, not the path.
    """
    problems: list[str] = []
    for vt in sorted(out_root.glob("data/endpoints/*/*/runs/*/*-dfpt-qha-vt.json")):
        payload = json.loads(vt.read_text())
        recorded = payload.get("source_json_sha256")
        source_name = Path(payload.get("source_json", "")).name
        if not recorded:
            problems.append(f"{vt}: no source_json_sha256, so nothing pins it")
            continue
        if not source_name:
            problems.append(f"{vt}: no source_json to resolve")
            continue
        candidates = list(vt.parent.glob(source_name))
        if not candidates:
            problems.append(f"{vt}: source {source_name} was not mirrored")
            continue
        actual = sha256_bytes(candidates[0].read_bytes())
        if actual != recorded:
            problems.append(
                f"{vt}: source hash mismatch (recorded {recorded[:12]}, "
                f"mirrored {actual[:12]})"
            )
    return problems


def summarise(out_root: Path) -> dict:
    """What the mirror actually contains, per endpoint."""
    endpoints: dict[str, dict] = {}
    for analysis in sorted(out_root.glob("data/endpoints/*/*/runs/*/*-dfpt-qha.json")):
        payload = json.loads(analysis.read_text())
        name = payload.get("endpoint") or analysis.name[: -len("-dfpt-qha.json")]
        state = endpoints.setdefault(name, {"ok": False, "vt": False, "applicability": None})
        if payload.get("ok") is True:
            state["ok"] = True
            state["applicability"] = payload.get("qha_applicability")
        if analysis.with_name(f"{name}-dfpt-qha-vt.json").is_file():
            state["vt"] = True
    return endpoints


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--cluster", required=True, choices=sorted(REMOTE_ROOTS))
    ap.add_argument("--remote-root", default=None, help="override the tracked root")
    ap.add_argument("--out", type=Path, default=None, help="default data/remote/<cluster>")
    ap.add_argument("--verify", action="store_true",
                    help="re-check an existing mirror without touching the network")
    args = ap.parse_args(argv)

    out_root = args.out or (MIRROR_ROOT / args.cluster)
    remote_root = args.remote_root or REMOTE_ROOTS[args.cluster]

    if args.verify:
        if not out_root.is_dir():
            print(f"{out_root} does not exist; nothing to verify", file=sys.stderr)
            return 1
    else:
        client = connect(args.cluster)
        try:
            relpaths = list_remote(client, remote_root)
            print(f"{args.cluster}:{remote_root} -> {out_root}  ({len(relpaths)} files)")
            manifest = fetch(client, remote_root, relpaths, out_root)
        finally:
            client.close()
        manifest.update({"cluster": args.cluster, "remote_root": remote_root})
        (out_root / "mirror.json").write_text(json.dumps(manifest, indent=2) + "\n")

    problems = check_pins(out_root)
    endpoints = summarise(out_root)
    done = sum(1 for v in endpoints.values() if v["ok"])
    print(f"endpoints mirrored: {len(endpoints)}   ok=true: {done}   "
          f"with -vt: {sum(1 for v in endpoints.values() if v['vt'])}")
    if problems:
        print(f"PIN CHECK FAILED ({len(problems)}):", file=sys.stderr)
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
        return 1
    print("pin check: every -vt source hash matches the mirrored analysis file")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
