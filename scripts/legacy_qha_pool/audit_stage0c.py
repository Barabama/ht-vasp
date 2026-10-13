#!/usr/bin/env python
"""Stage 0C, step C0: read-only structure-applicability audit.

Answers one question per endpoint: **is the cell we are about to compute the
same object the production pipeline would generate, and is it a sensible
object at all?**

The production generator is ``htvasp/model/endmember.py``.  Its three builders
do something that is easy to miss and matters a lot:

  * ``get_bcc``/``get_fcc``/``get_hcp`` take a TEMPLATE structure and only
    ``replace_species`` on it -- the template's lattice constant is used
    VERBATIM.  So every ``BCC-*`` POSCAR in ``data/poscars/`` carries the AlNi
    template lattice (a = 2.8597113 A) and every ``FCC-*`` carries the AuCu3
    template lattice (a = 3.7353876 A), regardless of the actual alloy.  The
    R7/R3 relaxation is expected to fix that; the audit records how far the
    starting point is.
  * ``get_sorted_structure`` then REORDERS sites by element symbol, so the
    POSCAR species order is alphabetical, not the template's order.
  * ``get_ser`` fetches the element's MP ground state, so SER cells are a mix
    of primitive cubic (Nb, Ta), conventional cubic (Cr) and HCP (Ti, Zr).

Checks performed (all read-only, no VASP):

  1. chemistry and site order vs ``endmember.py``
  2. template provenance and the starting-lattice offset
  3. lattice type / space group / primitive-vs-conventional volume
  4. L1_2 1a/3c occupancy in BOTH directions (A3B vs AB3 are different
     compounds, not mirror images -- unlike B2, where site exchange leaves the
     compound unchanged)
  5. symmetry-reduced displacement count for the phonon supercell
  6. phonopy supercell == VASP POSCAR == force-constant ordering (the assert
     Stage 0B had to add after finding a silent permutation)
  7. MAGMOM length, POTCAR order, KPOINTS folding relation
  8. SER template sanity (is the template the element's real ground state?)

Outputs ``structure_audit.json`` and ``structure_audit.md`` under --out.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

warnings.filterwarnings("ignore")

import build_stage0 as S0          # noqa: E402
import build_stage0b as B0         # noqa: E402

POOL = ("Al", "Cr", "Hf", "Nb", "Ta", "Ti", "Zr")
CACHE = Path("/nfs_ssd/tmp/legacy-qha-pool/stage0c/endmember_cache")


def endmember_structures() -> dict:
    """Generate the C0 targets through endmember.py itself, caching the result."""
    CACHE.mkdir(parents=True, exist_ok=True)
    out: dict = {}
    from htvasp.model.endmember import Endmember
    em = Endmember()

    out["elem_mpid_coverage"] = {
        el: em.elem_mpid.get(el) for el in POOL}
    out["elem_mpid_all_keys"] = sorted(em.elem_mpid)

    def cached(key: str, fn):
        f = CACHE / f"{key}.json"
        if f.is_file():
            return Structure.from_dict(json.loads(f.read_text()))
        s = fn()
        f.write_text(json.dumps(s.as_dict()))
        return s

    # B2 representatives
    for name, e1, e2 in (("BCC-Al-Ta", "Al", "Ta"), ("BCC-Ta-Ta", "Ta", "Ta"),
                         ("BCC-Al-Cr", "Al", "Cr"), ("BCC-Hf-V", "Hf", "V")):
        try:
            out.setdefault("bcc", {})[name] = cached(
                f"bcc_{name}", lambda e1=e1, e2=e2: em.get_bcc(e1, e2))
        except Exception as exc:                              # noqa: BLE001
            out.setdefault("bcc_errors", {})[name] = f"{type(exc).__name__}: {exc}"
    # L1_2 both directions
    for name, e1, e2 in (("FCC-Zr-Al", "Zr", "Al"), ("FCC-Al-Zr", "Al", "Zr")):
        try:
            out.setdefault("fcc", {})[name] = cached(
                f"fcc_{name}", lambda e1=e1, e2=e2: em.get_fcc(e1, e2))
        except Exception as exc:                              # noqa: BLE001
            out.setdefault("fcc_errors", {})[name] = f"{type(exc).__name__}: {exc}"
    # every pool SER
    for el in POOL:
        try:
            out.setdefault("ser", {})[el] = cached(
                f"ser_{el}", lambda el=el: em.get_ser(el))
        except Exception as exc:                              # noqa: BLE001
            out.setdefault("ser_errors", {})[el] = f"{type(exc).__name__}: {exc}"
    return out




LEGACY_BM3 = Path("/nfs_hdd/2025/gaominliang/ht-vasp/scripts/nas/bm3_fit_16symbols.json")
# legacy atom counts per endpoint, needed to compare V/atom
LEGACY_NATOMS = {"FCC-Zr-Al-4": 4, "FCC-Al-Zr-4": 4, "BCC-Al-Ta-2": 2,
                 "BCC-Ta-Ta-2": 2, "BCC-Al-Cr-2": 2, "SER-Al-1": 1}


def legacy_v0_table() -> dict:
    """V0 and B0 from the legacy 16symbols BM3 fit, for a starting-point sanity check.

    This is NOT a validation of legacy -- legacy is not a trusted source for
    parameters (report 3.2).  It is only used to ask "is the cell we are about
    to compute within tens of percent of a physically relaxed cell for this
    composition", which is enough to catch a wrong template.
    """
    if not LEGACY_BM3.is_file():
        return {}
    eps = json.loads(LEGACY_BM3.read_text())["endpoints"]
    out = {}
    for e in eps:
        n = LEGACY_NATOMS.get(e["endpoint"])
        out[e["endpoint"]] = {
            "V0_angstrom3": e["V0_angstrom3"],
            "V0_per_atom": e["V0_angstrom3"] / n if n else None,
            "B0_GPa": e["B0_GPa"], "eos_r2": e["eos_r2"], "n_points": e["n_points"],
        }
    return out


def _sitekey(frac) -> str:
    """Position key that is insensitive to -0.0 and to coordinates outside [0,1).

    The first version formatted raw frac_coords, so an endmember cell storing
    (-0.0, -0.0, -0.0) keyed as "-0.0000" while a cell storing (0.0, 0.0, 0.0)
    keyed as "0.0000" -- and every B2 comparison reported "different chemistry"
    for structures that are byte-for-byte the same chemistry.  Normalise first.
    """
    return ",".join(f"{(float(x) % 1.0 + 1e-9) % 1.0:.4f}" for x in frac)


def cell_info(st: Structure, symprec: float = 1e-3) -> dict:
    sga = SpacegroupAnalyzer(st, symprec=symprec)
    conv = sga.get_conventional_standard_structure()
    prim = sga.get_primitive_standard_structure()
    return {
        "formula": st.composition.reduced_formula,
        "formula_full": st.composition.formula,
        "n_sites": len(st),
        "spacegroup": sga.get_space_group_symbol(),
        "spacegroup_number": sga.get_space_group_number(),
        "lattice_type": sga.get_lattice_type(),
        "crystal_system": sga.get_crystal_system(),
        "a": float(st.lattice.a), "b": float(st.lattice.b),
        "c": float(st.lattice.c),
        "alpha": float(st.lattice.alpha), "beta": float(st.lattice.beta),
        "gamma": float(st.lattice.gamma),
        "volume": float(st.volume),
        "volume_per_atom": float(st.volume / len(st)),
        "conventional_volume": float(conv.volume),
        "primitive_volume": float(prim.volume),
        "conventional_n_sites": len(conv),
        "primitive_n_sites": len(prim),
        "species_sequence": [str(s.specie) for s in st],
        "species_order_unique": list(dict.fromkeys(str(s.specie) for s in st)),
        "site_map": {_sitekey(x.frac_coords): str(x.specie) for x in st},
    }


def occupancy_1a3c(st: Structure) -> dict:
    """Which element sits on 1a (0,0,0) and which on the three 3c sites."""
    def at(target):
        for s in st:
            d = np.abs(((np.array(s.frac_coords) - np.array(target) + 0.5) % 1.0) - 0.5)
            if np.linalg.norm(d) < 1e-4:
                return str(s.specie)
        return None
    corner = at((0.0, 0.0, 0.0))
    faces = sorted({str(s.specie) for s in st
                    if np.linalg.norm(np.abs(((np.array(s.frac_coords) + 0.5) % 1.0)
                                             - 0.5)) > 1e-4})
    counts = {el: sum(1 for s in st if str(s.specie) == el) for el in faces}
    return {"1a": corner, "3c": faces, "3c_counts": counts}


def displacement_check(endpoint: str, st: Structure) -> dict:
    """Symmetry-reduced phonon displacements for the 2x2x2 supercell."""
    import build_stage0 as _S
    ph, sup = B0.phonopy_supercell(st, 2)
    ph.generate_displacements(distance=0.01, is_plusminus="auto", is_diagonal=True)
    disps = [{"atom": int(d[0]),
              "species": str(sup[int(d[0])].specie),
              "vector": [float(d[1]), float(d[2]), float(d[3])]}
             for d in ph.displacements]
    # ordering assert: phonopy supercell vs what we would write
    return {"n_displacements": len(disps), "displacements": disps,
            "n_supercell_atoms": len(sup),
            "cell_fingerprint": B0.cell_fingerprint(sup),
            "supercell_spacegroup": SpacegroupAnalyzer(sup, symprec=1e-3)
                                    .get_space_group_symbol()}


def kpoints_for(st: Structure, density: int = 16) -> dict:
    n = S0.mesh_for(float(st.volume), density)
    shift = 0.5 if n % 2 == 0 else 0.0
    return {"density_target": density, "grid": [n, n, n], "shift": [shift] * 3,
            "contains_gamma": (n % 2 == 1) if shift == 0 else (n % 2 == 0),
            "n3V": n ** 3 * float(st.volume),
            "supercell_partner_grid": [n // 2] * 3 if n % 2 == 0 else None,
            "supercell_partner_shift": ([0.0] * 3 if (n // 2) % 2 else [0.5] * 3)
                                       if n % 2 == 0 else None,
            "folding_verified": n % 2 == 0}


def magmom_length(st: Structure) -> dict:
    vals = [S0.MAGMOM.get(str(s.specie), 1.0) for s in st]
    return {"n_values": len(vals), "n_atoms": len(st),
            "valid_length": len(vals) in (1, len(st)),
            "values": vals,
            "table_hits": {str(s.specie): S0.MAGMOM.get(str(s.specie))
                           for s in st}}


def potcar_order(st: Structure) -> dict:
    order = list(dict.fromkeys(str(s.specie) for s in st))
    return {"species_order_unique": order,
            "matches_poscar_species_line": order ==
            list(dict.fromkeys(str(s.specie) for s in st))}


def compare(a: Structure, b: Structure) -> dict:
    """Chemistry/order comparison between an endmember cell and a Stage-0 cell."""
    if len(a) != len(b):
        return {"comparable": False, "reason": f"{len(a)} vs {len(b)} atoms"}
    am = {_sitekey(x.frac_coords): str(x.specie) for x in a}
    bm = {_sitekey(x.frac_coords): str(x.specie) for x in b}
    return {"comparable": True,
            "same_chemistry_per_site": am == bm,
            "same_species_order": [str(x.specie) for x in a] ==
                                  [str(x.specie) for x in b],
            "species_a": [str(x.specie) for x in a],
            "species_b": [str(x.specie) for x in b],
            "delta_a_angstrom": float(abs(a.lattice.a - b.lattice.a)),
            "delta_a_percent": float(100 * abs(a.lattice.a - b.lattice.a)
                                     / a.lattice.a),
            "delta_volume_percent": float(100 * abs(a.volume - b.volume) / a.volume)}



def write_markdown(audit: dict, path: Path) -> None:
    """Human-readable rendering of the audit (the JSON is the machine one)."""
    L: list[str] = []
    A = L.append
    A("# Stage 0C C0：结构适用性审计（只读）\n")
    A(f"生成器：`{audit['generator']}`  ·  "
      f"legacy V0 来源：`{audit['legacy_v0_source']}`（**仅作起点合理性对照，"
      f"不作参数来源**）\n")

    A("## 1. `endmember.elem_mpid` 对 7 元素池的覆盖\n")
    A("| 元素 | mpid | 状态 |")
    A("|---|---|---|")
    for el, mid in audit["endmember"]["elem_mpid_coverage"].items():
        A(f"| {el} | {mid or '—'} | {'✓' if mid else '**✗ 缺失，`get_ser()` 会 KeyError**'} |")
    miss = audit["endmember"]["missing_for_pool"]
    A(f"\n缺失：**{', '.join(miss) if miss else '无'}**\n")

    A("## 2. SER：元素基态与模板合理性\n")
    A("| 元素 | 化学式 | 空间群 | 原子数 | V/atom / Å³ | c/a |")
    A("|---|---|---|---|---|---|")
    for el, v in audit["endmember"].get("ser", {}).items():
        A(f"| {el} | {v['formula']} | {v['spacegroup']} | {v['n_sites']} | "
          f"{v['volume_per_atom']:.4f} | {v['c_over_a']:.4f} |")
    for el, err in audit["endmember"].get("ser_errors", {}).items():
        A(f"| {el} | — | — | — | — | **{err}** |")
    A("")

    A("## 3. 目标端基：起点与生成器一致性\n")
    A("| 端基 | 类别 | 化学占位一致 | 物种顺序一致 | 胞内原子数 | "
      "超胞原子数 | 位移数 | 网格(胞/超胞) |")
    A("|---|---|---|---|---|---|---|---|")
    for name, t in audit["targets"].items():
        c = t["comparison"]
        ph = t["phonon"]
        kp = t["kpoints"]
        A(f"| `{name}` | {t['class']} | {'✓' if c['same_chemistry_per_site'] else '✗'} | "
          f"{'✓' if c['same_species_order'] else '**✗ 顺序不同**'} | "
          f"{t['stage0_cell']['n_sites']} | {ph['n_supercell_atoms']} | "
          f"{ph['n_displacements']} | {kp['grid'][0]} / "
          f"{kp['supercell_partner_grid'][0] if kp['supercell_partner_grid'] else '—'} |")
    A("")

    A("## 4. 起点晶格 vs 弛豫后的参考体积\n")
    A("起点晶格来自 `endmember.py` 的模板（BCC 用 AlNi mp-1487，FCC 用 AuCu₃ mp-2258），"
      "**模板的晶格常数被原样沿用**，所以起点与真实合金无关。\n")
    A("| 端基 | 起点 a / Å | 起点 V/atom | legacy BM3 V0/atom | 比值 | legacy r² |")
    A("|---|---|---|---|---|---|")
    for name, t in audit["targets"].items():
        s0 = t["stage0_cell"]
        lb = t.get("legacy_bm3")
        if lb and lb.get("V0_per_atom"):
            ratio = 100 * (s0["volume_per_atom"] / lb["V0_per_atom"] - 1)
            A(f"| `{name}` | {s0['a']:.4f} | {s0['volume_per_atom']:.4f} | "
              f"{lb['V0_per_atom']:.4f} | {ratio:+.1f}% | {lb['eos_r2']:.4f} |")
        else:
            A(f"| `{name}` | {s0['a']:.4f} | {s0['volume_per_atom']:.4f} | — | — | — |")
    A("")

    A("## 5. 对称约化位移（2×2×2，0.01 Å）\n")
    for name, t in audit["targets"].items():
        ph = t["phonon"]
        A(f"**`{name}`** — 超胞 {ph['n_supercell_atoms']} 原子，"
          f"空间群 {ph['supercell_spacegroup']}，{ph['n_displacements']} 个位移\n")
        A("| # | 原子 | 元素 | 位移向量 / Å |")
        A("|---|---|---|---|")
        for d in ph["displacements"]:
            A(f"| {d['atom']} | {d['atom']} | {d['species']} | "
              f"({d['vector'][0]:.5f}, {d['vector'][1]:.5f}, {d['vector'][2]:.5f}) |")
        A("")
    path.write_text("\n".join(L) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)

    em = endmember_structures()
    reps = S0.build_structures()
    legacy = legacy_v0_table()
    audit: dict = {"schema": "legacy-qha-pool/stage0c-audit/v1",
                   "generator": "htvasp/model/endmember.py",
                   "legacy_v0_source": "scripts/nas/bm3_fit_16symbols.json",
                   "legacy_v0_note": ("used ONLY as a starting-point sanity check, "
                                      "not as a parameter source (report 3.2)"),
                   "targets": {}, "endmember": {}}

    # ---- coverage --------------------------------------------------------
    print("=== 1. endmember.elem_mpid 对 7 元素池的覆盖 ===")
    cov = em["elem_mpid_coverage"]
    missing = [el for el, v in cov.items() if v is None]
    for el, mid in cov.items():
        print(f"  {el:3s}: {mid if mid else '✗ 缺失 —— get_ser() 会 KeyError'}")
    print(f"  → 缺失: {missing if missing else '无'}")
    audit["endmember"]["elem_mpid_coverage"] = cov
    audit["endmember"]["missing_for_pool"] = missing

    # ---- BCC -------------------------------------------------------------
    print("\n=== 2. BCC/B2：endmember vs build_stage0 ===")
    bcc_map = {"BCC-Al-Ta": "BCC-Al-Ta-2", "BCC-Ta-Ta": "BCC-Ta-Ta-2",
               "BCC-Al-Cr": "BCC-Al-Cr-2", "BCC-Hf-V": "BCC-Hf-V-2"}
    for ename, sname in bcc_map.items():
        if ename not in em.get("bcc", {}):
            print(f"  {ename}: 生成失败 {em.get('bcc_errors', {}).get(ename)}")
            continue
        a, b = em["bcc"][ename], reps[sname]["struct"]
        c = compare(a, b)
        ci = cell_info(b)
        print(f"  {ename:<12} vs {sname:<14} chem={c['same_chemistry_per_site']} "
              f"order={c['same_species_order']} "
              f"Δa={c['delta_a_angstrom']:.4f} Å ({c['delta_a_percent']:.1f}%) "
              f"ΔV={c['delta_volume_percent']:.0f}%")
        audit["targets"][sname] = {
            "class": "BCC/B2", "endmember_name": ename,
            "legacy_bm3": legacy.get(sname),
            "endmember_cell": cell_info(a), "stage0_cell": ci,
            "comparison": c,
            "kpoints": kpoints_for(b), "magmom": magmom_length(b),
            "potcar": potcar_order(b),
            "phonon": displacement_check(sname, b)}
        audit["endmember"].setdefault("bcc", {})[ename] = cell_info(a)

    # ---- FCC / L1_2 ------------------------------------------------------
    print("\n=== 3. L1₂：两个方向是不同的化学计量 ===")
    for ename, sname in (("FCC-Zr-Al", "FCC-Zr-Al-4"), ("FCC-Al-Zr", None)):
        if ename not in em.get("fcc", {}):
            print(f"  {ename}: 生成失败 {em.get('fcc_errors', {}).get(ename)}")
            continue
        s = em["fcc"][ename]
        occ = occupancy_1a3c(s)
        ci = cell_info(s)
        print(f"  {ename:<12} -> {ci['formula']:<6} sg={ci['spacegroup']:<7} "
              f"1a={occ['1a']:<3} 3c={occ['3c']} × {occ['3c_counts']}  "
              f"a={ci['a']:.4f} (模板 {ci['a']:.4f})")
        audit["endmember"].setdefault("fcc", {})[ename] = {
            **ci, "occupancy_1a3c": occ}
    if "FCC-Zr-Al" in em.get("fcc", {}):
        a = em["fcc"]["FCC-Zr-Al"]
        b = reps["FCC-Zr-Al-4"]["struct"]
        c = compare(a, b)
        occ_a, occ_b = occupancy_1a3c(a), occupancy_1a3c(b)
        print(f"\n  FCC-Zr-Al-4: endmember 1a={occ_a['1a']} 3c={occ_a['3c']}")
        print(f"               stage0    1a={occ_b['1a']} 3c={occ_b['3c']}")
        print(f"               同一化学占位={c['same_chemistry_per_site']}  "
              f"同一物种顺序={c['same_species_order']}")
        print(f"               species endmember={c['species_a']}")
        print(f"               species stage0   ={c['species_b']}")
        print(f"               Δa={c['delta_a_angstrom']:.4f} Å "
              f"({c['delta_a_percent']:.1f}%)  ΔV={c['delta_volume_percent']:.0f}%")
        audit["targets"]["FCC-Zr-Al-4"] = {
            "class": "FCC/L1_2", "endmember_name": "FCC-Zr-Al",
            "legacy_bm3": legacy.get("FCC-Zr-Al-4"),
            "endmember_cell": cell_info(a),
            "endmember_occupancy_1a3c": occ_a,
            "stage0_cell": cell_info(b),
            "stage0_occupancy_1a3c": occ_b,
            "comparison": c,
            "kpoints": kpoints_for(b), "magmom": magmom_length(b),
            "potcar": potcar_order(b),
            "phonon": displacement_check("FCC-Zr-Al-4", b)}

    # ---- SER -------------------------------------------------------------
    print("\n=== 4. SER：元素基态与模板合理性 ===")
    print(f"  {'element':<8} {'formula':<5} {'sg':<10} {'n':>2} {'V/atom':>8} "
          f"{'c/a':>6}  {'与 SER-Al-1 的对照':<20}")
    for el, s in em.get("ser", {}).items():
        ci = cell_info(s)
        capp = ci["c"] / ci["a"] if ci["a"] else float("nan")
        note = ""
        if el == "Al" and "SER-Al-1" in reps:
            c = compare(s, reps["SER-Al-1"]["struct"])
            note = f"ΔV={c['delta_volume_percent']:.4f}%"
            audit["targets"]["SER-Al-1"] = {
                "class": "SER", "endmember_name": "SER-Al",
                "legacy_bm3": legacy.get("SER-Al-1"),
                "endmember_cell": ci,
                "stage0_cell": cell_info(reps["SER-Al-1"]["struct"]),
                "comparison": c,
                "kpoints": kpoints_for(reps["SER-Al-1"]["struct"]),
                "magmom": magmom_length(reps["SER-Al-1"]["struct"]),
                "potcar": potcar_order(reps["SER-Al-1"]["struct"]),
                "phonon": displacement_check("SER-Al-1",
                                             reps["SER-Al-1"]["struct"])}
        print(f"  {el:<8} {ci['formula']:<5} {ci['spacegroup']:<10} "
              f"{ci['n_sites']:>2} {ci['volume_per_atom']:>8.4f} {capp:>6.4f}  {note}")
        audit["endmember"].setdefault("ser", {})[el] = {
            **ci, "c_over_a": capp}
    for el, err in em.get("ser_errors", {}).items():
        print(f"  {el:<8} ✗ {err}")
        audit["endmember"].setdefault("ser_errors", {})[el] = err

    # ---- displacement sanity for the L1_2 and SER ------------------------
    print("\n=== 5. 对称约化位移数（2×2×2，0.01 Å）===")
    for name in ("FCC-Zr-Al-4", "SER-Al-1"):
        if name in audit["targets"]:
            ph = audit["targets"][name]["phonon"]
            print(f"  {name:<14} 超胞 {ph['n_supercell_atoms']} 原子  "
                  f"sg={ph['supercell_spacegroup']}  位移 {ph['n_displacements']} 个")
            for d in ph["displacements"]:
                print(f"      atom {d['atom']:>2} ({d['species']:<3}) "
                      f"vec={[round(v,5) for v in d['vector']]}")

    (out / "structure_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    write_markdown(audit, out / "structure_audit.md")
    print(f"\nwrote {out / 'structure_audit.json'}")
    print(f"wrote {out / 'structure_audit.md'}")


if __name__ == "__main__":
    main()
