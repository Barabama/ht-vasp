#!/usr/bin/env python
"""Write pool_structures.json for the 84-endpoint DFPT-QHA pool.

The inventory itself lives in `htvasp.dfpt_qha.pool`; this is the CLI that
runs it and prints the summary a human reads before submitting anything.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from htvasp.dfpt_qha.pool import build_pool

DEFAULT_CACHE = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0c/endmember_cache")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path,
                    help="directory to write pool_structures.json into")
    ap.add_argument("--cache", type=Path, default=DEFAULT_CACHE,
                    help="SER structure cache from endmember.get_ser")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    manifest = build_pool(args.cache)
    report = manifest.dump()
    (args.out / "pool_structures.json").write_text(json.dumps(report, indent=2) + "\n")

    endpoints = report["endpoints"]
    print(f"端基总数 = {len(endpoints)}   {report['by_class']}")
    print()
    print(f"{'class':>10} {'n':>4} {'cell atoms':>11} {'Vegard?':>8} {'proxy?':>7}")
    print("-" * 46)
    for endpoint_class in ("BCC_B2", "FCC_L12", "SER"):
        subset = [e for e in endpoints.values() if e["class"] == endpoint_class]
        counts = sorted({e["n_atoms"] for e in subset})
        vegard = all("Vegard" in e["lattice"] for e in subset)
        proxies = sum(1 for e in subset if e["is_proxy"])
        print(f"{endpoint_class:>10} {len(subset):>4} {str(counts):>11} "
              f"{str(vegard):>8} {proxies:>7}")

    print()
    print("=== 同一物理状态出现多次的 ===")
    for state, names in sorted(report["duplicate_states"].items()):
        first, second = (endpoints[name]["a"] for name in names[:2])
        print(f"  {state:>10}: {names[0]} (a={first:.5f})  vs  "
              f"{names[1]} (a={second:.5f})  Δa = {100 * (second / first - 1):+.5f}%")

    print()
    proxies = [e for e in endpoints.values() if e["is_proxy"]]
    print(f"=== 用了基态体积代理的端基: {len(proxies)} 个 ===")
    by_element: dict[str, int] = {}
    for endpoint in proxies:
        for element, source in (endpoint["vpa_source"] or {}).items():
            if "PROXY" in source:
                by_element[element] = by_element.get(element, 0) + 1
    for element, count in sorted(by_element.items()):
        print(f"  {element}: {count} 个端基用了代理体积")

    print()
    print(f"wrote {args.out / 'pool_structures.json'}")


if __name__ == "__main__":
    main()
