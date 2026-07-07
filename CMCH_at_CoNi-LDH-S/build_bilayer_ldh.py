#!/usr/bin/env python3
"""构建 LDH bilayer slab（双层锚点验证）

从已弛豫的 LDH monolayer slab 出发，用体相层间距堆叠双层，
施加 Level 2 应变（与 build_heterostructure.py 一致的公共晶格）。

输出:
  data/poscars/LDH2_strained.vasp   — 本征双层 LDH

依赖:
  data/LDH_strained/static_out.json   — 已弛豫的 monolayer slab
  data/CoNiOH2/dos_out.json           — 体相 LDH（提供层间距）

用法:
  python build_bilayer_ldh.py                # 生成 POSCAR
  python build_bilayer_ldh.py --check        # 验证已生成的 POSCAR
  python build_bilayer_ldh.py --ttf-frac 0.24  # 手动指定 TTF 阈值
"""

import argparse
import json
from pathlib import Path

import numpy as np
from pymatgen.core import Structure, Lattice

# ── 路径 ──────────────────────────────────────────────────────────────
ROOT = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S")
DATA = ROOT / "data"
OUTDIR = DATA / "poscars"
OUTDIR.mkdir(exist_ok=True)

# ── Level 2 公共晶格参数（与 build_heterostructure.py 一致）────────
AVG_A = 10.3255    # Å
AVG_B = 9.4602     # Å
TOTAL_C = 28.0     # Å（真空层 + slab）
VAC_BOTTOM = 2.0   # Å（底部预留真空）

# ── 默认 TTF 阈值 ─────────────────────────────────────────────────
#  底层 O-M-O 最高 z = VAC_BOTTOM + monolayer_thickness
#  = 2.0 + 4.08 = 6.08 Å → frac = 6.08/28 ≈ 0.217
#  顶层底部 H 起始 z = VAC_BOTTOM + monolayer_thickness + interlayer_gap
#  = 2.0 + 4.08 + 0.435 = 6.52 Å → frac = 6.52/28 ≈ 0.233
#  取 0.22 恰好覆盖底层全部 60 原子，不触及顶层
DEFAULT_TTF_FRAC = 0.22


# ═══════════════════════════════════════════════════════════════════════
# 1. 加载已有结构
# ═══════════════════════════════════════════════════════════════════════

def load_monolayer() -> Structure:
    """加载已弛豫的 LDH monolayer slab。"""
    path = DATA / "LDH_strained" / "static_out.json"
    with open(path) as f:
        d = json.load(f)
    s = Structure.from_dict(d["output"]["structure"])
    print(f"  Monolayer: atoms={len(s)}  a={s.lattice.a:.4f}  b={s.lattice.b:.4f}")
    return s


def get_interlayer_gap() -> float:
    """从体相 LDH 推算层间间距。

    体相 CoNiOH₂ (120 atoms, c=9.03 Å) 包含 2 个完整 O-M-O 层，
    层间距 = c/2 − monolayer_slab_thickness。
    """
    bulk_path = DATA / "CoNiOH2" / "dos_out.json"
    with open(bulk_path) as f:
        d = json.load(f)
    bulk = Structure.from_dict(d["output"]["structure"])
    layer_spacing = bulk.lattice.c / 2  # ≈ 4.515 Å

    mono = load_monolayer()
    mono_thickness = mono.cart_coords[:, 2].max() - mono.cart_coords[:, 2].min()  # ≈ 4.08 Å

    gap = layer_spacing - mono_thickness
    print(f"  Bulk c={bulk.lattice.c:.3f} Å, 层间距=c/2={layer_spacing:.3f} Å")
    print(f"  Monolayer 厚度={mono_thickness:.3f} Å")
    print(f"  层间间隙={gap:.3f} Å")
    return gap


# ═══════════════════════════════════════════════════════════════════════
# 2. 构建双层 slab
# ═══════════════════════════════════════════════════════════════════════

def build_bilayer(mono: Structure, interlayer_gap: float) -> Structure:
    """堆叠两层 monolayer slab 为双层。

    层序（z 增加方向）:
      Layer 1 (bottom): monolayer copy, 底部原子将被 TTF 约束
      Layer 2 (top):    monolayer copy, 全自由

    Args:
        mono: 已弛豫的 monolayer slab（60 atoms）
        interlayer_gap: 两层 O-M-O 之间的最小间距 (Å)

    Returns:
        双层 slab，c = TOTAL_C，底部起始位于 VAC_BOTTOM。
    """
    cart_z = mono.cart_coords[:, 2]
    z_min = cart_z.min()
    z_max = cart_z.max()
    thickness = z_max - z_min   # 单层厚度 ≈ 4.08 Å

    # --- Layer 1: 底部归零 ---
    coords1 = mono.cart_coords.copy()
    coords1[:, 2] -= z_min      # layer 1 bottom at z = 0

    # --- Layer 2: 整体上移 thickness + gap ---
    shift = thickness + interlayer_gap
    coords2 = coords1.copy()
    coords2[:, 2] += shift      # layer 2 bottom at z = shift

    # --- 合并 ---
    all_species = list(mono.species) * 2
    all_coords = np.vstack([coords1, coords2])

    # 平移到 VAC_BOTTOM
    all_coords[:, 2] += VAC_BOTTOM

    # 新晶格: a,b 暂用 monolayer 的（后面再应变），c = TOTAL_C
    lat = Lattice.from_parameters(
        mono.lattice.a, mono.lattice.b, TOTAL_C,
        90.0, 90.0, 90.0,
    )
    return Structure(lat, all_species, all_coords, coords_are_cartesian=True)


# ═══════════════════════════════════════════════════════════════════════
# 3. 应变与约束（复用 build_heterostructure.py 的逻辑）
# ═══════════════════════════════════════════════════════════════════════

def apply_strain(slab: Structure, target_a: float, target_b: float) -> Structure:
    """面内应变到公共晶格。"""
    new_lat = Lattice.from_parameters(
        target_a, target_b, slab.lattice.c,
        slab.lattice.alpha, slab.lattice.beta, slab.lattice.gamma,
    )
    return Structure(new_lat, slab.species, slab.frac_coords,
                     site_properties=slab.site_properties,
                     coords_are_cartesian=False)


def apply_ttf(slab: Structure, ttf_frac: float) -> Structure:
    """对 frac_z < ttf_frac 的原子施加 TTF（xy 自由，z 固定）约束。"""
    sd = []
    for i in range(len(slab)):
        if slab.frac_coords[i, 2] < ttf_frac:
            sd.append([True, True, False])   # TTF
        else:
            sd.append([True, True, True])    # TTT
    slab.add_site_property("selective_dynamics", sd)
    return slab


# ═══════════════════════════════════════════════════════════════════════
# 4. 验证
# ═══════════════════════════════════════════════════════════════════════

def check_poscar(struct: Structure, label: str = "") -> None:
    """打印验证摘要。"""
    z = struct.cart_coords[:, 2]
    sd = struct.site_properties.get("selective_dynamics", [])
    ttf = sum(1 for d in sd if list(d) == [True, True, False])
    ttt = len(sd) - ttf

    # 按元素统计
    elem_counts = {}
    for site in struct.sites:
        sym = site.specie.symbol
        elem_counts[sym] = elem_counts.get(sym, 0) + 1

    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")
    print(f"  原子数:   {len(struct)}")
    print(f"  元素组成: {elem_counts}")
    print(f"  晶格:     a={struct.lattice.a:.4f}  b={struct.lattice.b:.4f}  c={struct.lattice.c:.4f}")
    print(f"  z 范围:   [{z.min():.2f}, {z.max():.2f}] Å  厚度={z.max()-z.min():.2f} Å")
    print(f"  TTF: {ttf}  TTT: {ttt}")

    # 分层分析：找到层间间隙
    sorted_z = np.sort(z)
    diffs = np.diff(sorted_z)
    max_diff_idx = np.argmax(diffs)
    gap_z = sorted_z[max_diff_idx]
    print(f"  层间间隙位置: z≈{gap_z:.2f} Å  (Δ={diffs[max_diff_idx]:.3f} Å)")
    print(f"  底层: z=[{z.min():.2f}, {gap_z:.2f}]  顶层: z=[{sorted_z[max_diff_idx+1]:.2f}, {z.max():.2f}]")


# ═══════════════════════════════════════════════════════════════════════
# 5. 主函数
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="构建 LDH bilayer slab")
    parser.add_argument("--ttf-frac", type=float, default=DEFAULT_TTF_FRAC,
                        help=f"TTF 分数坐标阈值 (默认 {DEFAULT_TTF_FRAC})")
    parser.add_argument("--check", action="store_true",
                        help="验证已生成的 POSCAR")
    args = parser.parse_args()

    outpath = OUTDIR / "LDH2_strained.vasp"

    # ── 验证模式 ──
    if args.check:
        if not outpath.exists():
            print(f"❌ {outpath} 不存在，请先生成。")
            return
        s = Structure.from_file(str(outpath))
        check_poscar(s, "LDH2_strained (验证)")
        return

    # ── 加载 ──
    print("=== 加载 monolayer + 体相层间距 ===")
    mono = load_monolayer()
    interlayer_gap = get_interlayer_gap()

    # ── 堆叠 ──
    print("\n=== 构建 bilayer ===")
    bilayer = build_bilayer(mono, interlayer_gap)
    print(f"  堆叠完成: atoms={len(bilayer)}")

    # ── 应变 ──
    print(f"\n=== 施加 Level 2 应变 (a={AVG_A:.4f}, b={AVG_B:.4f}) ===")
    ea = (AVG_A - mono.lattice.a) / mono.lattice.a * 100
    eb = (AVG_B - mono.lattice.b) / mono.lattice.b * 100
    print(f"  应变: εa={ea:+.2f}%  εb={eb:+.2f}%")
    bilayer = apply_strain(bilayer, AVG_A, AVG_B)

    # ── TTF 约束 ──
    print(f"\n=== TTF 约束 (frac_z < {args.ttf_frac}) ===")
    bilayer = apply_ttf(bilayer, args.ttf_frac)

    # ── 保存 ──
    bilayer.to(str(outpath), fmt="poscar")
    print(f"\n✅ 已保存: {outpath}")

    # ── 验证 ──
    check_poscar(bilayer, "LDH2_strained (最终)")

    # ── 与 monolayer 对比 ──
    print(f"\n{'='*60}")
    print(f"  与 monolayer 对比")
    print(f"{'='*60}")
    print(f"  {'':<20s}  {'Monolayer':>12s}  {'Bilayer':>12s}")
    print(f"  {'原子数':<20s}  {len(mono):>12d}  {len(bilayer):>12d}")
    print(f"  {'a (Å)':<20s}  {mono.lattice.a:>12.4f}  {bilayer.lattice.a:>12.4f}")
    print(f"  {'b (Å)':<20s}  {mono.lattice.b:>12.4f}  {bilayer.lattice.b:>12.4f}")
    print(f"  {'c (Å)':<20s}  {'28.0':>12s}  {'28.0':>12s}")
    mono_thick = mono.cart_coords[:,2].max() - mono.cart_coords[:,2].min()
    bi_thick = bilayer.cart_coords[:,2].max() - bilayer.cart_coords[:,2].min()
    print(f"  {'Slab 厚度 (Å)':<20s}  {mono_thick:>12.2f}  {bi_thick:>12.2f}")

    # 预计机时
    n_relax = len(bilayer) - sum(1 for d in bilayer.site_properties["selective_dynamics"]
                                 if list(d) == [True, True, False])
    print(f"\n  预计计算量:")
    print(f"    弛豫自由原子: {n_relax}/{len(bilayer)}")
    print(f"    预计 GPU 时: relax ~4-5h, static+DOS ~2-3h")


if __name__ == "__main__":
    main()
