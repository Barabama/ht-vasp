"""Write one line per completed endpoint: name|cluster|status|applicability.

Scans the local ``data/endpoints`` tree and, by default, every cluster mirror
under ``data/remote/<cluster>/data/endpoints``.  The mirror is how a single
invocation on 429pro can describe the whole pool: ``build_manifest.py`` used to
need a round trip per cluster, which meant the manifest was only ever as fresh
as the last manual copy.

    # on 429pro: local artifacts, plus whatever B5070 mirror is present
    ./.conda/bin/python scripts/gather_status.py 429pro /tmp/status_429.txt
    # a mirror only, tagging its rows with the mirror's cluster
    ./.conda/bin/python scripts/gather_status.py --mirrors-only /tmp/status_507.txt

An endpoint present in both places is reported once per cluster; the manifest
takes the later line, and ``--local-first`` puts the local row last so the
cluster that is still writing to its tree wins.
"""
import glob
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def scan(root: Path, cluster: str, rows: set[str]) -> int:
    found = 0
    for path in glob.glob(str(root / "*" / "dfpt-qhaflow" / "runs" / "*"
                             / "*-dfpt-qha.json")):
        payload = json.load(open(path))
        if payload.get("ok") is not True:
            continue
        name = payload.get("endpoint") or Path(path).name.split("-dfpt-qha")[0]
        rows.add("%s|%s|done %s|%s" % (name, cluster, cluster,
                                       payload.get("qha_applicability") or "unknown"))
        found += 1
    return found


def mirror_roots() -> list[tuple[Path, str]]:
    roots: list[tuple[Path, str]] = []
    for mirror in sorted((ROOT / "data" / "remote").glob("*")):
        tree = mirror / "data" / "endpoints"
        if tree.is_dir():
            roots.append((tree, mirror.name))
    return roots


def main(argv: list[str]) -> int:
    args = [a for a in argv if not a.startswith("-")]
    flags = {a for a in argv if a.startswith("-")}
    mirrors_only = "--mirrors-only" in flags
    local_first = "--local-first" in flags

    rows: set[str] = set()
    mirrors = mirror_roots()
    local = (ROOT / "data" / "endpoints", args[0] if args else "local")

    order = (mirrors + [local]) if local_first else ([local] + mirrors)
    if mirrors_only:
        order = mirrors
    if not order:
        print("nothing to scan: no local tree and no data/remote/<cluster>",
              file=sys.stderr)
        return 1

    for root, cluster in order:
        print(f"{cluster:8s} {root}  {scan(root, cluster, rows)} artifacts",
              file=sys.stderr)
    out = Path(args[1]) if len(args) > 1 else Path("status.txt")
    out.write_text("\n".join(sorted(rows)) + "\n")
    print(len(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
