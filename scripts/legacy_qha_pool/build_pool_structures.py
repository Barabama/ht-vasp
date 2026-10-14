#!/usr/bin/env python
"""The 84-endpoint pool: structure inventory + the recipe decisions behind it.

Report 8.3 step 5 ("SER primitive vs conventional") turned out to be a bigger
question than its name: the pool's binary endpoints and its SER entries were
being built from DIFFERENT cells, and for 4 elements they are the same crystal.

The decisions this module encodes (report 8.0c / 8.0d):

1. REPRESENTATION = conventional for everything.
   B2 and L1_2 are simple-cubic Bravais lattices with a 2- or 4-atom basis, so
   their cells ARE primitive already -- `use_symmetrized_structure` is a no-op
   for all 77 binary endpoints.  The only cells it would inflate are the SER
   entries whose element is fcc or bcc (Al, Nb, Ta), and those are inflated to
   the conventional cell on purpose.

2. LATTICE CONSTANT = Vegard from elemental volumes, NOT the `endmember.py`
   template.  The template reuses one lattice constant (AlNi 2.8597 for BCC,
   Cu3Au 3.7354 for FCC) for every composition, which is up to -47% in volume
   (report 4.12.5 correction).  Vegard is within a few percent.

3. The 4 "duplicate" pairs (Al fcc, Cr/Nb/Ta bcc) collapse.  Because Vegard for
   an A=A endpoint IS the element's own volume, `BCC-Ta-Ta-2` and `SER-Ta-1`
   (converted to conventional) come out the SAME CELL -- measured agreement
   <0.002% in a.  So they are one calculation with two CALPHAD roles, not two
   independent data points (report 8.0c(f)).

The occupancy/chemistry still comes from `endmember.py` -- only the lattice
constant is overridden.  That keeps the tested generator in the loop while
removing the part of it that was wrong.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

POOL = ("Al", "Cr", "Hf", "Nb", "Ta", "Ti", "Zr")

# Ground-state phase and volume per atom, from the MP entries `endmember.py`
# already names (report 4.12.6 / 8.0c).  Hf is mp-103 per the user.
GROUND = {
    "Al": ("FCC", 16.4717, "mp-134"),
    "Cr": ("BCC", 13.0845, "mp-90"),
    "Hf": ("HCP", 22.2640, "mp-103"),
    "Nb": ("BCC", 18.2581, "mp-75"),
    "Ta": ("BCC", 18.1300, "mp-50"),
    "Ti": ("HCP", 17.3194, "mp-46"),
    "Zr": ("HCP", 23.4997, "mp-131"),
}

# Element volumes in a phase that is NOT its ground state.  Only the entries
# that have been looked up appear here; everything else falls back to the
# ground-state volume and is flagged `proxy` so a reader can see it.
POLYMORPH = {
    ("Al", "BCC"): (16.1160, "mp-998860"),
    ("Hf", "BCC"): (22.0870, "mp-100"),
    # fcc polymorphs of the bcc/hcp metals have no looked-up entry yet
}


def vpa_of(el: str, phase: str) -> tuple[float, str, bool]:
    """(volume per atom, source string, is_proxy)."""
    if (el, phase) in POLYMORPH:
        v, src = POLYMORPH[(el, phase)]
        return v, f"polymorph {src}", False
    gs_phase, v, src = GROUND[el]
    if gs_phase == phase:
        return v, f"ground state {src}", False
    return v, f"PROXY: ground state {src} is {gs_phase}, used for {phase}", True


def cubic_a(vpa: float, n_atoms: int) -> float:
    return (vpa * n_atoms) ** (1.0 / 3.0)


def bcc_2atoms(a: float, e1: str, e2: str) -> Structure:
    """B2 / A2 conventional cell: e1 at 1a (0,0,0), e2 at 1b (1/2,1/2,1/2)."""
    return Structure(Lattice.cubic(a), [e1, e2],
                     [[0, 0, 0], [0.5, 0.5, 0.5]], coords_are_cartesian=False)


def fcc_4atoms(a: float, e1: str, e2: str) -> Structure:
    """L1_2 / A1 conventional cell: e1 at 1a (0,0,0), e2 at 3c (0,1/2,1/2)..."""
    return Structure(Lattice.cubic(a), [e1, e2, e2, e2],
                     [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
                     coords_are_cartesian=False)


# `endmember.elem_mpid` has 13 keys: 7 elements outside this pool and no Hf, so
# `get_ser("Hf")` raises KeyError (report 4.12.2).  The mapping is passed in
# HERE rather than patched into htvasp/ -- the production default still lacks
# it, and this generator records that it had to supply the value.
ELEM_MPID_OVERRIDE = {"Hf": "mp-103"}


def ser_structure(el: str, cache: Path) -> tuple[Structure, dict]:
    """The element's MP ground state, fetched through `endmember.py`."""
    from htvasp.model.endmember import Endmember
    f = cache / f"ser_{el}.json"
    if f.is_file():
        return Structure.from_dict(json.loads(f.read_text())), {"source": "cache"}
    base = Endmember()
    table = dict(base.elem_mpid)
    note = None
    for k, v in ELEM_MPID_OVERRIDE.items():
        if k not in table:
            note = f"supplied {k}={v} at the call site (htvasp default lacks it)"
        table[k] = v
    em = Endmember(elem_mpid=table)
    s = em.get_ser(el)
    cache.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(s.as_dict()))
    return s, {"source": "endmember.get_ser", "elem_mpid_note": note}


def conventional(st: Structure) -> Structure:
    return SpacegroupAnalyzer(st, symprec=1e-3).get_conventional_standard_structure()


def site_role(st: Structure) -> dict:
    """Which element sits on 1a and which on 3c, from the fractional coords."""
    if len(st) not in (2, 4):
        return {}
    at_1a, at_3c = None, []
    for s in st:
        f = np.array(s.frac_coords) % 1.0
        if np.allclose(f, 0.0, atol=1e-4) or np.allclose(f, 1.0, atol=1e-4):
            at_1a = str(s.specie)
        else:
            at_3c.append(str(s.specie))
    return {"site_1a": at_1a, "site_3c": sorted(at_3c) if at_3c else []}


def enum_pool() -> list[dict]:
    out = []
    for a, b in itertools.combinations_with_replacement(POOL, 2):
        out.append({"name": f"BCC-{a}-{b}-2", "class": "BCC_B2",
                    "elements": [a, b], "phase": "BCC"})
    for a in POOL:
        for b in POOL:
            out.append({"name": f"FCC-{a}-{b}-4", "class": "FCC_L12",
                        "elements": [a, b], "phase": "FCC",
                        "site_note": f"{a}@1a / {b}@3c"})
    for a in POOL:
        out.append({"name": f"SER-{a}-1", "class": "SER",
                    "elements": [a], "phase": "SER"})
    return out


def physical_state_id(rec: dict) -> str | None:
    """Name of the physical crystal this endpoint IS, for the duplicate pairs.

    Only an A=A binary endpoint can coincide with an SER entry, and only when
    the element's ground state IS that phase.  The four that do, per report
    8.0c, are Al (fcc) and Cr/Nb/Ta (bcc).
    """
    el = rec["elements"]
    if rec["class"] == "SER":
        return f"{GROUND[rec['elements'][0]][0].lower()}-{rec['elements'][0]}"
    if rec["class"] in ("BCC_B2", "FCC_L12") and el[0] == el[1]:
        phase = rec["phase"]
        if GROUND[el[0]][0] == phase:
            return f"{phase.lower()}-{el[0]}"
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--cache", type=Path,
                    default=Path("/nfs_ssd/tmp/legacy-qha-pool/stage0c/endmember_cache"))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    recs = enum_pool()
    assert len(recs) == 84, len(recs)
    by_class = {}
    for r in recs:
        by_class[r["class"]] = by_class.get(r["class"], 0) + 1
    assert by_class == {"BCC_B2": 28, "FCC_L12": 49, "SER": 7}, by_class

    inventory: dict = {}
    states: dict[str, list[str]] = {}
    for r in recs:
        name, cls, els, phase = r["name"], r["class"], r["elements"], r["phase"]
        if cls == "SER":
            st, src = ser_structure(els[0], args.cache)
            conv = conventional(st)
            rec = {"structure": conv.as_dict(), "n_atoms": len(conv),
                   "lattice": "conventional of the MP ground state",
                   "a": float(conv.lattice.a),
                   "volume_per_atom": float(conv.volume) / len(conv),
                   "source": src["source"],
                   "elem_mpid_note": src.get("elem_mpid_note"),
                   "vpa_source": f"MP ground state {GROUND[els[0]][2]}",
                   "is_proxy": False,
                   "detour": (f"MP gives {len(st)} atoms; conventional is "
                              f"{len(conv)}" if len(st) != len(conv) else None)}
        else:
            a_el, b_el = els
            va, sa, pa = vpa_of(a_el, phase)
            vb, sb, pb = vpa_of(b_el, phase)
            if cls == "BCC_B2":
                vpa = (va + vb) / 2.0
                n = 2
                st = bcc_2atoms(cubic_a(vpa, n), a_el, b_el)
            else:
                vpa = (va + 3.0 * vb) / 4.0
                n = 4
                st = fcc_4atoms(cubic_a(vpa, n), a_el, b_el)
            rec = {"structure": st.as_dict(), "n_atoms": n,
                   "lattice": "Vegard (NOT the endmember.py template)",
                   "a": float(st.lattice.a),
                   "volume_per_atom": vpa,
                   "source": "Vegard from elemental volumes",
                   "vpa_source": {a_el: sa, b_el: sb},
                   "is_proxy": bool(pa or pb),
                   "site_role": site_role(st)}
        rec.update({k: r[k] for k in ("name", "class", "elements", "phase")
                    if k in r})
        rec["physical_state_id"] = physical_state_id(r)
        inventory[name] = rec
        if rec["physical_state_id"]:
            states.setdefault(rec["physical_state_id"], []).append(name)

    dupes = {k: v for k, v in states.items() if len(v) > 1}
    report = {"schema": "legacy-qha-pool/pool-structures/v1",
              "n_endpoints": len(inventory), "by_class": by_class,
              "decisions": {
                  "representation": "conventional for all 84",
                  "lattice_constant": "Vegard from elemental volumes",
                  "template_lattice_used": False,
                  "duplicate_pairs_share_one_curve": True},
              "physical_states": states,
              "duplicate_states": dupes,
              "endpoints": inventory}
    (args.out / "pool_structures.json").write_text(
        json.dumps(report, indent=2) + "\n")

    print(f"端基总数 = {len(inventory)}   {by_class}")
    print()
    print(f"{'class':>10} {'n':>4} {'cell atoms':>11} {'Vegard?':>8} {'proxy?':>7}")
    print("-" * 46)
    for cls in ("BCC_B2", "FCC_L12", "SER"):
        sub = [v for v in inventory.values() if v["class"] == cls]
        ns = sorted({v["n_atoms"] for v in sub})
        veg = all("Vegard" in v["lattice"] for v in sub)
        prox = sum(1 for v in sub if v["is_proxy"])
        print(f"{cls:>10} {len(sub):>4} {str(ns):>11} {str(veg):>8} {prox:>7}")
    print()
    print("=== 同一物理状态出现多次的 ===")
    for st, names in sorted(dupes.items()):
        a1 = inventory[names[0]]["a"]
        a2 = inventory[names[1]]["a"]
        print(f"  {st:>10}: {names[0]} (a={a1:.5f})  vs  {names[1]} (a={a2:.5f})  "
              f"Δa = {100*(a2/a1-1):+.5f}%")
    print()
    prox = [k for k, v in inventory.items() if v["is_proxy"]]
    print(f"=== 用了基态体积代理的端基: {len(prox)} 个 ===")
    by_el = {}
    for k in prox:
        for el, s in inventory[k]["vpa_source"].items():
            if "PROXY" in s:
                by_el.setdefault(el, 0)
                by_el[el] += 1
    for el, n in sorted(by_el.items()):
        print(f"  {el}: {n} 个端基用了代理体积")
    print()
    print(f"wrote {args.out / 'pool_structures.json'}")


if __name__ == "__main__":
    main()
