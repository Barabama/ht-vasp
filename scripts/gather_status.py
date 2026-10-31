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

**One line per endpoint per cluster, not one per run.**  An endpoint that was
run twice (an aborted attempt beside a good one, a resubmission, a re-analysis)
has several ``*-dfpt-qha.json`` files, and which of them is the endpoint's
answer is a decision -- not something to leave to sort order.  This used to add
a row per run to a ``set`` and let ``build_manifest`` keep the last line, i.e.
the alphabetically greatest applicability string won: BCC-Ta-Ta-2 reported
``unknown`` from a 5-volume artifact of 2026-10-16 while a complete 7-volume run
from 2026-10-20 beside it said ``harmonic_QHA_applicable``.
"""
import glob
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


#: Worst-first, matching ``htvasp.dfpt_qha.flow_jobs.APPLICABILITY_ORDER``.
#: Duplicated rather than imported: this script is run from a checkout on either
#: cluster and must not drag the workflow package (and therefore atomate2) in.
APPLICABILITY_RANK = {
    "harmonic_QHA_applicable": 0,
    "incommensurate_soft_mode_suspected": 1,
    "finite_q_soft_mode": 2,
    "gamma_soft_mode": 3,
    "numerical_unresolved": 4,
}
#: An artifact with no verdict at all is the weakest evidence there is: it
#: predates the ``qha_applicability`` field.  Ranked after every real verdict.
UNKNOWN_RANK = len(APPLICABILITY_RANK)


def best_first(payload: dict, path: str) -> tuple:
    """Sort key for choosing between several runs of one endpoint.

    Completeness comes first: a 7-volume run is the endpoint's answer even if it
    found a soft mode, because a 5-volume artifact cannot have fitted as much of
    the F(V,T) surface -- and reporting "harmonic" off a 3-volume run would be
    worse than reporting "soft mode" off a complete one.  A verdict that exists
    then beats one that does not (an artifact without ``qha_applicability``
    predates that field), the quality of the verdict breaks the remaining ties,
    and modification time is last so the key is total: on a genuine tie the
    newer artifact wins rather than an arbitrary one.
    """
    applicability = payload.get("qha_applicability")
    return (
        -len(payload.get("volumes") or []),
        0 if applicability else 1,
        APPLICABILITY_RANK.get(applicability, UNKNOWN_RANK),
        -Path(path).stat().st_mtime,
    )


def scan(root: Path, cluster: str, rows: set[str]) -> int:
    """Add one row per endpoint, built from that endpoint's best run."""
    candidates: dict[str, list[tuple[tuple, dict]]] = {}
    found = 0
    for path in glob.glob(str(root / "*" / "dfpt-qhaflow" / "runs" / "*"
                             / "*-dfpt-qha.json")):
        payload = json.load(open(path))
        if payload.get("ok") is not True:
            continue
        name = payload.get("endpoint") or Path(path).name.split("-dfpt-qha")[0]
        candidates.setdefault(name, []).append((best_first(payload, path), payload))
        found += 1

    for name, runs in candidates.items():
        best = min(runs, key=lambda item: item[0])[1]
        rows.add("%s|%s|done %s|%s" % (name, cluster, cluster,
                                       best.get("qha_applicability") or "unknown"))
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
