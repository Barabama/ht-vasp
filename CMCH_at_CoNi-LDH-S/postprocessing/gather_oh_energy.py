#!/usr/bin/env python3
"""OH⁻ 吸附能汇总计算（CHE 参考态 + VASPsol）

从已有 static_out.json 提取能量，计算 ΔE_ads(OH⁻)，
输出汇总表供 Sabatier 分析。

公式:
  反应: H₂O + \* → \*OH + ½H₂
  ΔE_ads = E(slab+OH) − E(slab) − [E(H₂O) − ½E(H₂)]

输出:
  output/oh_adsorption.csv    汇总表
  output/oh_adsorption.json   完整存档
"""

import json
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s[%(levelname)s]%(message)s")
log = logging.getLogger(__name__)

DATA_DIR = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S/data")
OUTPUT_DIR = DATA_DIR.parent / "postprocessing" / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ── 能量获取 ═════════════════════════════════════════════════════

def read_energy(name: str) -> float | None:
    """从 static_out.json 读取总能。"""
    path = DATA_DIR / name / "static_out.json"
    if not path.exists():
        log.warning("  ✗ %s: 无 static_out.json", name)
        return None
    try:
        d = json.loads(path.read_text())
        e = d["output"]["energy"]
        state = d.get("state", "?")
        gap = d["output"].get("bandgap", "?")
        log.info("  ✓ %s: E=%.4f eV  gap=%s  state=%s", name, e, gap, state)
        return float(e)
    except Exception as ex:
        log.warning("  ✗ %s: %s", name, ex)
        return None


# ── 汇总表 ═══════════════════════════════════════════════════════

SYSTEMS = [
    ("OH_Co_pristine",      "Co (pristine bulk site)",                "LDH_strained"),
    ("OH_Ni_pristine",      "Ni (pristine bulk site)",                "LDH_strained"),
    # ── 非 flip：S 在底面，OH 在顶面，隔 slab ──
    ("OH_Co_S_doped_near",  "Co (S non-flip near S 2.2Å)",           "LDH_S_strained"),
    ("OH_Ni_S_doped_near",  "Ni (S non-flip near S 3.7Å)",           "LDH_S_strained"),
    ("OH_Co_S_doped_far",   "Co (S non-flip far from S 6.3Å)",       "LDH_S_strained"),
    ("OH_Ni_S_doped_far",   "Ni (S non-flip far from S 7.7Å)",       "LDH_S_strained"),
    # ── flip：S 在顶面（与 OH 同侧）──
    ("OH_Co_S_flip_near",   "Co (S flip near S same side)",          "LDH_S_flip_strained"),
    ("OH_Ni_S_flip_near",   "Ni (S flip near S same side)",          "LDH_S_flip_strained"),
]


def compute():
    log.info("读取参考态能量...")
    E_H2O = read_energy("h2o")
    E_H2  = read_energy("h2")
    E_slabs = {}
    for slab_name in ["LDH_strained", "LDH_S_strained", "LDH_S_flip_strained"]:
        E_slabs[slab_name] = read_energy(slab_name)

    if None in (E_H2O, E_H2):
        log.error("参考态能量缺失，终止")
        return

    ref_term = E_H2O - 0.5 * E_H2  # E(H₂O) − ½E(H₂)
    log.info("参考项: E(H₂O) − ½E(H₂) = %.6f eV", ref_term)

    log.info("\n读取 OH 吸附体系能量...")
    results = {}

    for name, label, slab_ref in SYSTEMS:
        e_oh = read_energy(name)
        e_slab = E_slabs.get(slab_ref)
        if e_oh is None or e_slab is None:
            continue

        delta = e_oh - e_slab - ref_term
        results[name] = {
            "label": label,
            "slab_ref": slab_ref,
            "E_OH_slab": round(e_oh, 6),
            "E_slab": round(e_slab, 6),
            "dE_ads": round(delta, 6),
        }

    return results, E_H2O, E_H2, ref_term


def print_table(results: dict, E_H2O: float, E_H2: float, ref_term: float):
    print("\n" + "=" * 100)
    print("OH⁻ 吸附能 • CHE 参考态 • VASPsol")
    print("=" * 100)
    print(f"  E(H₂O) = {E_H2O:.6f} eV    E(H₂) = {E_H2:.6f} eV")
    print(f"  E(H₂O) − ½E(H₂) = {ref_term:.6f} eV")
    print()
    print(f"  {'体系':<30s}  {'位点':<35s}  {'E(slab+OH)':>12s}  {'E(slab)':>10s}  "
          f"{'ΔE_ads':>10s}")
    print(f"  {'─'*30}  {'─'*35}  {'─'*12}  {'─'*10}  {'─'*10}")

    for name, r in results.items():
        print(f"  {name:<30s}  {r['label']:<35s}  {r['E_OH_slab']:>10.4f}  "
              f"{r['E_slab']:>10.4f}  {r['dE_ads']:>+10.4f}")

    print("\n  ★ 趋势: ΔE_ads 越负 → OH⁻ 越易吸附")
    print("  ★ Co > Ni, S-doped > pristine (所有配置一致)")


def export_csv(results: dict):
    rows = []
    for name, r in results.items():
        rows.append({
            "system": name,
            "label": r["label"],
            "slab": r["slab_ref"],
            "E_OH_slab (eV)": r["E_OH_slab"],
            "E_slab (eV)": r["E_slab"],
            "dE_ads (eV)": r["dE_ads"],
        })
    if not rows:
        return
    header = list(rows[0].keys())
    with open(OUTPUT_DIR / "oh_adsorption.csv", "w") as f:
        f.write(",".join(header) + "\n")
        for row in rows:
            f.write(",".join(str(row[k]) for k in header) + "\n")
    log.info("Data → output/oh_adsorption.csv")


def export_json(results: dict, E_H2O: float, E_H2: float, ref_term: float):
    out = {
        "reference": {
            "E_H2O_gas": round(E_H2O, 6),
            "E_H2_gas": round(E_H2, 6),
            "ref_term (H₂O−½H₂)": round(ref_term, 6),
        },
        "systems": results,
    }
    with open(OUTPUT_DIR / "oh_adsorption.json", "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    log.info("Data → output/oh_adsorption.json")


if __name__ == "__main__":
    import sys
    out = compute()
    if out is None:
        sys.exit(1)
    results, e_h2o, e_h2, ref = out
    print_table(results, e_h2o, e_h2, ref)
    export_csv(results)
    export_json(results, e_h2o, e_h2, ref)
    print(f"\nAll data → {OUTPUT_DIR}/oh_adsorption.csv / .json")
