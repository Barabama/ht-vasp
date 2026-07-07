#!/usr/bin/env python3
"""在 LDH slab 每个金属顶位生成 OH 吸附 POSCAR（全部遍历）

遍历 LDH_strained 和 LDH_S_strained 的每个 Co/Ni 顶位，
在金属正上方放置 OH 吸附质，输出到 data/poscars/OH_*/。

用法:
  python build_oh_adsorption.py          # 生成全部 24 个 POSCAR
  python build_oh_adsorption.py --check  # 检查已生成的 POSCAR
"""

import argparse
import json
from pathlib import Path

import numpy as np
from pymatgen.core import Structure, Lattice, Molecule

# ── 路径 ──────────────────────────────────────────────────────────────
ROOT = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S")
DATA = ROOT / "data"
OUTDIR = DATA / "poscars"
OUTDIR.mkdir(exist_ok=True)

# ── 吸附质参数 ────────────────────────────────────────────────────────
OH_BOND_LEN = 0.96          # O-H 键长 (Å)
INITIAL_M_O_DIST = 2.0      # 金属-O 初始距离 (Å)
OH_TILT_ANGLE = 15.0        # OH 键轴与表面法线夹角 (°)

# ── Level 2 参数 ──────────────────────────────────────────────────────
TOTAL_C = 28.0
VAC_BOTTOM = 2.0
DEFAULT_TTF_FRAC = 0.20


# ═══════════════════════════════════════════════════════════════════════
# 1. 加载 slab
# ═══════════════════════════════════════════════════════════════════════

def load_slab(name: str) -> Structure:
    """从 static_out.json 加载弛豫后的 slab 结构。"""
    path = DATA / name / "static_out.json"
    with open(path) as f:
        d = json.load(f)
    return Structure.from_dict(d["output"]["structure"])


def get_metal_sites(slab: Structure) -> list[dict]:
    """提取所有金属位点 (Co/Ni) 的信息。"""
    sites = []
    for i, site in enumerate(slab.sites):
        if site.specie.symbol in ("Co", "Ni"):
            sites.append({
                "index": i,
                "species": site.specie.symbol,
                "frac_coords": site.frac_coords.tolist(),
                "cart_coords": slab.cart_coords[i].tolist(),
            })
    return sites


def get_s_position(slab: Structure) -> np.ndarray | None:
    """找到 S 原子的笛卡尔坐标（如有）。"""
    for i, site in enumerate(slab.sites):
        if site.specie.symbol == "S":
            return slab.cart_coords[i]
    return None


# ═══════════════════════════════════════════════════════════════════════
# 2. 构建 OH 吸附模型
# ═══════════════════════════════════════════════════════════════════════

def build_oh_on_site(slab: Structure, metal_idx: int,
                     ttf_frac: float = DEFAULT_TTF_FRAC) -> Structure:
    """在指定金属位点上方放置 OH 吸附质。

    策略:
      1. 取金属的 (x, y) 坐标
      2. O 放在 top_surface_H 的 z + 1.5 Å 处（金属正上方）
      3. H 沿表面法线（+z）方向偏移 OH_BOND_LEN，略倾斜避免 H-H 重叠
      4. 让 VASP 弛豫找到最终平衡位置

    Args:
        slab: 已弛豫的 LDH slab
        metal_idx: 金属位点在 slab 中的索引
        ttf_frac: TTF 分数坐标阈值

    Returns:
        含吸附 OH 的新结构
    """
    metal_site = slab.sites[metal_idx]
    metal_cart = slab.cart_coords[metal_idx]
    mx, my = metal_cart[0], metal_cart[1]

    # top surface H 的最大 z（所有 H 中 z 最大的即为顶面）
    h_indices = [i for i, s in enumerate(slab.sites) if s.specie.symbol == "H"]
    top_z = max(slab.cart_coords[i][2] for i in h_indices)

    # O 放在 top_H 上方 1.5 Å，正对金属 xy 位置
    oz = top_z + 1.5
    ox = mx
    oy = my

    # H 沿 z 方向偏移，轻微 x 倾斜以避免与表面 H 重叠
    tilt = np.radians(OH_TILT_ANGLE)
    hx = ox + OH_BOND_LEN * np.sin(tilt)
    hy = oy
    hz = oz + OH_BOND_LEN * np.cos(tilt)

    # 合并结构
    new_species = list(slab.species) + ["O", "H"]
    new_coords = np.vstack([
        slab.cart_coords,
        [[ox, oy, oz]],
        [[hx, hy, hz]],
    ])

    result = Structure(
        slab.lattice, new_species, new_coords,
        coords_are_cartesian=True,
    )

    # TTF 约束（底层固定 + 吸附质自由）
    sd = []
    for i in range(len(result)):
        if result.frac_coords[i, 2] < ttf_frac:
            sd.append([True, True, False])   # TTF
        else:
            sd.append([True, True, True])    # TTT
    result.add_site_property("selective_dynamics", sd)

    return result


# ═══════════════════════════════════════════════════════════════════════
# 3. 主流程
# ═══════════════════════════════════════════════════════════════════════

def run_all():
    """遍历所有金属位点，生成吸附 POSCAR。"""
    # ── 加载 slab ──
    slab_pristine = load_slab("LDH_strained")
    slab_s_doped = load_slab("LDH_S_strained")
    slab_s_flip = load_slab("LDH_S_flip_strained")  # S在顶面

    configs = [
        ("pristine", slab_pristine, None),
        ("S_doped",   slab_s_doped,  get_s_position(slab_s_doped)),
        ("S_flip",    slab_s_flip,   get_s_position(slab_s_flip)),   # ← 新增
    ]

    summary = []

    for slab_label, slab, s_xyz in configs:
        metal_sites = get_metal_sites(slab)
        print(f"\n{'='*70}")
        print(f"  {slab_label}: {len(metal_sites)} 个金属位点")
        print(f"{'='*70}")

        for site in metal_sites:
            idx = site["index"]
            species = site["species"]
            cx, cy, cz = site["cart_coords"]

            # 与 S 的距离（仅 S_doped）
            dist_s = None
            if s_xyz is not None:
                dist_s = np.linalg.norm(np.array(site["cart_coords"]) - s_xyz)

            # 生成 POSCAR
            oh_struct = build_oh_on_site(slab, idx)

            # 命名: OH_Co3_pristine.vasp  or  OH_Co3_S_doped.vasp
            # 按同元素编号
            same_species_indices = [
                s["index"] for s in metal_sites if s["species"] == species
            ]
            local_num = same_species_indices.index(idx) + 1
            fname = f"OH_{species}{local_num}_{slab_label}.vasp"
            outpath = OUTDIR / fname

            oh_struct.to(str(outpath), fmt="poscar")

            entry = {
                "file": fname,
                "species": species,
                "metal_idx": idx,
                "cart_xyz": [round(cx, 3), round(cy, 3), round(cz, 3)],
                "dist_to_S": round(dist_s, 3) if dist_s is not None else None,
            }
            summary.append(entry)

            dist_str = f"  d(S)={dist_s:.2f} Å" if dist_s is not None else ""
            print(f"  ✅ {fname:<35s}  "
                  f"({species}{local_num} @ z={cz:.2f}){dist_str}")

    # ── 汇总表 ──
    print(f"\n{'='*70}")
    print(f"  汇总: {len(summary)} 个 POSCAR 已生成")
    print(f"{'='*70}")

    for slab_label in ["pristine", "S_doped", "S_flip"]:
        entries = [e for e in summary if slab_label in e["file"]]
        print(f"\n  ── {slab_label} ──")
        print(f"  {'文件':<35s}  {'元素':>4s}  {'z(Å)':>6s}  {'d(S)(Å)':>8s}")
        print(f"  {'─'*35}  {'─'*4}  {'─'*6}  {'─'*8}")
        for e in entries:
            d_s = f"{e['dist_to_S']:.2f}" if e["dist_to_S"] is not None else "  —"
            print(f"  {e['file']:<35s}  {e['species']:>4s}  "
                  f"{e['cart_xyz'][2]:>6.2f}  {d_s:>8s}")

    # S_doped 中按距离排序的推荐
    s_doped_entries = [e for e in summary if "S_doped" in e["file"] and e["dist_to_S"] is not None]
    if s_doped_entries:
        s_doped_entries.sort(key=lambda e: e["dist_to_S"])
        print(f"\n  S_doped 位点按与 S 的距离排序:")
        for i, e in enumerate(s_doped_entries):
            marker = "← 最近" if i == 0 else ("← 最远" if i == len(s_doped_entries)-1 else "")
            print(f"    {e['file']:<35s}  d(S)={e['dist_to_S']:.2f} Å  {marker}")

    print(f"\n输出目录: {OUTDIR}")


def check():
    """检查已生成的 OH 吸附 POSCAR。"""
    oh_files = sorted(OUTDIR.glob("OH_*.vasp"))
    if not oh_files:
        print("❌ 未找到 OH_* 吸附 POSCAR")
        return

    print(f"找到 {len(oh_files)} 个 OH 吸附 POSCAR:\n")
    print(f"  {'文件':<35s}  {'原子数':>6s}  {'O-metal(Å)':>10s}")
    print(f"  {'─'*35}  {'─'*6}  {'─'*10}")

    for f in oh_files:
        s = Structure.from_file(str(f))
        # 找吸附 O（最后添加的 O，z 最大的 O）
        o_indices = [i for i, site in enumerate(s.sites) if site.specie.symbol == "O"]
        max_o_idx = max(o_indices, key=lambda i: s.cart_coords[i][2])
        max_o_cart = s.cart_coords[max_o_idx]
        # 找最近的金属
        metal_indices = [i for i, site in enumerate(s.sites) if site.specie.symbol in ("Co", "Ni")]
        min_dist = min(
            np.linalg.norm(max_o_cart - s.cart_coords[m])
            for m in metal_indices
        )
        print(f"  {f.name:<35s}  {len(s):>6d}  {min_dist:>10.2f}")


# ═══════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="OH 吸附 POSCAR 遍历生成")
    parser.add_argument("--check", action="store_true", help="检查已生成的 POSCAR")
    args = parser.parse_args()

    if args.check:
        check()
    else:
        run_all()
