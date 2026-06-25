#!/usr/bin/env python3
"""构建 gap=1.0 异质结 + 匹配晶格 slab (Level 2)

输出 (data/poscars/):
  3 个异质结: hetero_intrinsic, hetero_S_doped, hetero_S_exposed
  4 个应变slab: CMCH_strained, LDH_strained, LDH_S_strained, LDH_S_flip_strained

所有 7 个体系使用同一公共晶格 (avg_a, avg_b)，保证 NGXF/NGYF 一致。

用法:
  python build_heterostructure.py                    # 全部生成
  python build_heterostructure.py --gap 1.0          # 指定gap
  python build_heterostructure.py --ttf-frac 0.18    # TTF阈值(分数坐标)
  python build_heterostructure.py --check            # 验证POSCAR
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.core.interface import Interface

# ── 路径 ──────────────────────────────────────────────────────────────────
ROOT = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S")
DATA = ROOT / "data"
OUTDIR = DATA / "poscars"
OUTDIR.mkdir(exist_ok=True)

# ── 默认参数 ──────────────────────────────────────────────────────────────
DEFAULT_GAP = 1.0
DEFAULT_TTF_FRAC = 0.20   # 分数坐标 z < 0.20 的 CMCH 原子 → TTF
VAC_BOTTOM = 2.0          # 底部真空 (Å)
VAC_TOP = 12.0            # 顶部真空 (Å)
TOTAL_C = 28.0            # 总 c 轴 (Å)

# ── slab 加载配置 ──────────────────────────────────────────────────────────
SLAB_CONFIGS = {
    "CMCH": {
        "file": "CoMnH2CO5-slab",
        "strained_out": "CMCH_strained.vasp",
        "substrate": True,
    },
    "LDH": {
        "file": "CoNiOH2-slab",
        "strained_out": "LDH_strained.vasp",
        "substrate": False,
    },
    "LDH_S": {
        "file": "CoNiOH2S-noH-slab",
        "strained_out": "LDH_S_strained.vasp",
        "substrate": False,
    },
    "LDH_S_flip": {
        "file": "CoNiOH2S-noH-slab-flip",
        "strained_out": "LDH_S_flip_strained.vasp",
        "substrate": False,
        "swap_ab": True,
    },
}

HETERO_CONFIGS = [
    {
        "name": "intrinsic",
        "film": "LDH",
        "out": "hetero_intrinsic.vasp",
        "label": "本征异质结",
    },
    {
        "name": "S_doped",
        "film": "LDH_S",
        "out": "hetero_S_doped.vasp",
        "label": "S掺杂 (S在真空侧)",
    },
    {
        "name": "S_exposed",
        "film": "LDH_S_flip",
        "out": "hetero_S_exposed.vasp",
        "label": "S暴露 (S在界面侧)",
    },
]


# ═══════════════════════════════════════════════════════════════════════
# 1. 加载 slab
# ═══════════════════════════════════════════════════════════════════════

def load_slabs() -> dict[str, Structure]:
    """从 static_out.json 读取所有优化后 slab。"""
    slabs = {}
    for key, cfg in SLAB_CONFIGS.items():
        path = DATA / cfg["file"] / "static_out.json"
        with open(path) as f:
            d = json.load(f)
        s = Structure.from_dict(d["output"]["structure"])

        if cfg.get("swap_ab"):
            s = _swap_ab(s)

        slabs[key] = s
        print(f"  加载 {cfg['file']:<30s}  a={s.lattice.a:.4f}  b={s.lattice.b:.4f}  atoms={len(s)}")
    return slabs


def _swap_ab(slab: Structure) -> Structure:
    """交换 a/b 轴。"""
    old = slab.lattice
    new_lat = Lattice.from_parameters(
        old.b, old.a, old.c,
        old.alpha, old.beta, old.gamma,
    )
    frac = slab.frac_coords.copy()
    frac[:, [0, 1]] = frac[:, [1, 0]]
    return Structure(new_lat, slab.species, frac,
                     site_properties=slab.site_properties,
                     coords_are_cartesian=False)


# ═══════════════════════════════════════════════════════════════════════
# 2. 公共晶格
# ═══════════════════════════════════════════════════════════════════════

def compute_common_lattice(slabs: dict) -> tuple[float, float]:
    """用 CMCH + LDH 计算公共晶格。"""
    cmch = slabs["CMCH"]
    ldh = slabs["LDH"]
    avg_a = (cmch.lattice.a + ldh.lattice.a) / 2
    avg_b = (cmch.lattice.b + ldh.lattice.b) / 2
    return avg_a, avg_b


def strain_to(slab: Structure, target_a: float, target_b: float) -> Structure:
    """应变 slab 面内 a/b 到目标值 (保持分数坐标)。"""
    old = slab.lattice
    new_lat = Lattice.from_parameters(
        target_a, target_b, old.c,
        old.alpha, old.beta, old.gamma,
    )
    return Structure(new_lat, slab.species, slab.frac_coords,
                     site_properties=slab.site_properties,
                     coords_are_cartesian=False)


# ═══════════════════════════════════════════════════════════════════════
# 3. 构建异质结
# ═══════════════════════════════════════════════════════════════════════

def build_hetero(substrate: Structure, film: Structure,
                 gap: float, avg_a: float, avg_b: float,
                 ttf_frac: float) -> Structure:
    """构建一个异质结。

    Args:
        substrate: CMCH 基底
        film: LDH 薄膜
        gap: 初始界面间距 (Å)
        avg_a, avg_b: 公共晶格参数
        ttf_frac: 分数坐标 TTF 阈值

    Returns:
        含 selective_dynamics 的结构。
    """
    cmch_m = strain_to(substrate, avg_a, avg_b)
    film_m = strain_to(film, avg_a, avg_b)

    iface = Interface.from_slabs(
        substrate_slab=cmch_m,
        film_slab=film_m,
        gap=gap,
        vacuum_over_film=VAC_BOTTOM + VAC_TOP,
        center_slab=True,
    )

    zs = iface.cart_coords[:, 2]
    labels = np.array(iface.site_properties["interface_label"])
    shift = VAC_BOTTOM - zs[labels == "substrate"].min()

    coords_new = iface.cart_coords.copy()
    coords_new[:, 2] += shift
    new_lat = Lattice.from_parameters(
        iface.lattice.a, iface.lattice.b, TOTAL_C,
        iface.lattice.alpha, iface.lattice.beta, iface.lattice.gamma,
    )
    result = Structure(new_lat, iface.species, coords_new,
                       site_properties=iface.site_properties,
                       coords_are_cartesian=True)

    # 选择性动力学: TTF
    zs_final = result.cart_coords[:, 2]
    labels_final = np.array(result.site_properties["interface_label"])
    sub_mask = labels_final == "substrate"

    sd = []
    for i in range(len(result)):
        if sub_mask[i] and result.frac_coords[i, 2] < ttf_frac:
            sd.append([True, True, False])    # TTF
        else:
            sd.append([True, True, True])     # TTT

    result.add_site_property("selective_dynamics", sd)
    return result


# ═══════════════════════════════════════════════════════════════════════
# 4. 构建应变 slab
# ═══════════════════════════════════════════════════════════════════════

def build_strained_slab(slab: Structure, avg_a: float, avg_b: float,
                        ttf_frac: float) -> Structure:
    """将 slab 应变到公共晶格，加真空，设 TTF 约束。"""
    s = strain_to(slab, avg_a, avg_b)
    z_min = s.cart_coords[:, 2].min()

    new_lat = Lattice.from_parameters(
        s.lattice.a, s.lattice.b, TOTAL_C,
        s.lattice.alpha, s.lattice.beta, s.lattice.gamma,
    )
    coords_new = s.cart_coords.copy()
    coords_new[:, 2] += VAC_BOTTOM - z_min

    result = Structure(new_lat, s.species, coords_new,
                       coords_are_cartesian=True)

    sd = []
    for i in range(len(result)):
        if result.frac_coords[i, 2] < ttf_frac:
            sd.append([True, True, False])    # TTF
        else:
            sd.append([True, True, True])     # TTT

    result.add_site_property("selective_dynamics", sd)
    return result


# ═══════════════════════════════════════════════════════════════════════
# 5. 验证
# ═══════════════════════════════════════════════════════════════════════

def check_poscars() -> bool:
    """验证所有 POSCAR 的晶格和 TTF 约束。"""
    files = [
        "hetero_intrinsic.vasp", "hetero_S_doped.vasp", "hetero_S_exposed.vasp",
        "CMCH_strained.vasp", "LDH_strained.vasp",
        "LDH_S_strained.vasp", "LDH_S_flip_strained.vasp",
    ]

    print(f"\n{'='*65}")
    print(f"  POSCAR 验证")
    print(f"{'='*65}")
    print(f"  {'文件':<32s}  {'a':>8s}  {'b':>8s}  {'atoms':>6s}  {'TTF':>6s}  {'TTT':>6s}")
    print(f"  {'─'*32}  {'─'*8}  {'─'*8}  {'─'*6}  {'─'*6}  {'─'*6}")

    all_ok = True
    ref_a = ref_b = None

    for fname in files:
        path = OUTDIR / fname
        if not path.exists():
            print(f"  ❌ {fname}: 文件不存在")
            all_ok = False
            continue

        s = Structure.from_file(path)
        a, b = s.lattice.a, s.lattice.b

        if ref_a is None:
            ref_a, ref_b = a, b
        elif abs(a - ref_a) > 0.001 or abs(b - ref_b) > 0.001:
            print(f"  ⚠️ {fname}: 晶格不一致! a={a:.4f} b={b:.4f} vs ref a={ref_a:.4f} b={ref_b:.4f}")
            all_ok = False

        if "selective_dynamics" in s.site_properties:
            sd = s.site_properties["selective_dynamics"]
            ttf = sum(1 for d in sd if list(d) == [True, True, False])
            ttt = sum(1 for d in sd if list(d) == [True, True, True])
            other = len(sd) - ttf - ttt
            if other:
                print(f"  ❌ {fname}: 有非 TTF/TTT 约束 (other={other})")
                all_ok = False
        else:
            ttf, ttt = 0, len(s)

        print(f"  {'✅' if path.exists() else '❌'} {fname:<30s}  {a:>8.4f}  {b:>8.4f}  {len(s):>6d}  {ttf:>6d}  {ttt:>6d}")

    if all_ok and ref_a:
        print(f"\n  ✅ 全部通过 — 所有体系晶格一致 (a={ref_a:.4f}, b={ref_b:.4f})")
    else:
        print(f"\n  ⚠️ 有错误需要修正")

    return all_ok


# ═══════════════════════════════════════════════════════════════════════
# 6. 主函数
# ═══════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="构建 gap=1.0 异质结 + 匹配晶格 slab")
    parser.add_argument("--gap", type=float, default=DEFAULT_GAP,
                        help=f"初始界面间距 (默认 {DEFAULT_GAP})")
    parser.add_argument("--ttf-frac", type=float, default=DEFAULT_TTF_FRAC,
                        help=f"TTF 分数坐标阈值 (默认 {DEFAULT_TTF_FRAC})")
    parser.add_argument("--check", action="store_true",
                        help="验证已生成的 POSCAR")
    args = parser.parse_args()

    if args.check:
        check_poscars()
        return

    # ── 加载 slab ──
    print("=== 加载优化后 slab ===")
    slabs = load_slabs()

    # ── 公共晶格 ──
    avg_a, avg_b = compute_common_lattice(slabs)
    print(f"\n公共晶格: avg_a={avg_a:.4f}  avg_b={avg_b:.4f}")
    for key, s in slabs.items():
        ea = (avg_a - s.lattice.a) / s.lattice.a * 100
        eb = (avg_b - s.lattice.b) / s.lattice.b * 100
        print(f"  {key:<20s}  εa={ea:+.2f}%  εb={eb:+.2f}%")

    # ── 构建 3 个异质结 ──
    print(f"\n=== 构建异质结 (gap={args.gap}) ===")
    for cfg in HETERO_CONFIGS:
        substrate = slabs["CMCH"]
        film = slabs[cfg["film"]]
        struct = build_hetero(substrate, film, args.gap, avg_a, avg_b, args.ttf_frac)

        sub_zs = struct.cart_coords[:, 2][
            np.array(struct.site_properties["interface_label"]) == "substrate"
        ]
        cmch_top = sub_zs.max()
        film_zs = struct.cart_coords[:, 2][
            np.array(struct.site_properties["interface_label"]) == "film"
        ]
        true_gap = film_zs.min() - cmch_top

        outpath = OUTDIR / cfg["out"]
        struct.to(str(outpath), fmt="poscar")
        print(f"  ✅ {cfg['label']:<20s} → {outpath.name}")
        print(f"     界面间隙={true_gap:.2f} Å  atoms={len(struct)}")

    # ── 构建 4 个应变 slab ──
    print(f"\n=== 构建应变 slab ===")
    for key, cfg in SLAB_CONFIGS.items():
        s = slabs[key]
        struct = build_strained_slab(s, avg_a, avg_b, args.ttf_frac)
        z_min = struct.cart_coords[:, 2].min()
        z_max = struct.cart_coords[:, 2].max()

        outpath = OUTDIR / cfg["strained_out"]
        struct.to(str(outpath), fmt="poscar")
        print(f"  ✅ {key:<20s} → {outpath.name}")
        print(f"     z=[{z_min:.2f}, {z_max:.2f}]  height={z_max-z_min:.2f}  atoms={len(struct)}")

    # ── 最终验证 ──
    print(f"\n{'='*55}")
    check_poscars()
    print(f"\n完成！")
    print(f"输出目录: {OUTDIR}")


if __name__ == "__main__":
    main()
