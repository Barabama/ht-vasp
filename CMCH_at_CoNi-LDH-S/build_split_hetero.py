#!/usr/bin/env python3
"""将异质结 CONTCAR 按 z=0.3 (分数坐标) 拆分为 LDH 和 CMCH 两部分

上半部分 (frac_z > 0.3): LDH
下半部分 (frac_z < 0.3): CMCH

输出到 data/poscars/

用法:
  python build_split_hetero.py                    # 拆分全部 3 个异质结
  python build_split_hetero.py --intrinsic        # 仅拆分 intrinsic
  python build_split_hetero.py --s-doped          # 仅拆分 S_doped
  python build_split_hetero.py --s-exposed        # 仅拆分 S_exposed
  python build_split_hetero.py --check            # 验证已生成的 POSCAR
"""

import argparse
import gzip
from pathlib import Path

import numpy as np
from pymatgen.core import Structure

# ── 路径 ──────────────────────────────────────────────────────────────────
ROOT = Path("/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S")
DATA = ROOT / "data"
OUTDIR = DATA / "poscars"
OUTDIR.mkdir(exist_ok=True)

# ── 输入配置 ──────────────────────────────────────────────────────────────
HETERO_CONFIGS = [
    {
        "name": "intrinsic",
        "input": "hetero_intrinsic/3-static/CONTCAR.gz",
        "ldh_out": "LDH_intrinsic.vasp",
        "cmch_out": "CMCH_intrinsic.vasp",
        "label": "本征异质结",
    },
    {
        "name": "S_doped",
        "input": "hetero_s_doped/3-static/CONTCAR.gz",
        "ldh_out": "LDH_S_doped.vasp",
        "cmch_out": "CMCH_S_doped.vasp",
        "label": "S掺杂异质结",
    },
    {
        "name": "S_exposed",
        "input": "hetero_s_exposed/3-static/CONTCAR.gz",
        "ldh_out": "LDH_S_exposed.vasp",
        "cmch_out": "CMCH_S_exposed.vasp",
        "label": "S暴露异质结",
    },
]

# ── 拆分参数 ──────────────────────────────────────────────────────────────
SPLIT_FRAC_Z = 0.30   # 分数坐标 z = 0.3 处拆分
DEFAULT_TTF_FRAC = 0.20  # TTF 阈值


def load_structure(gz_path: Path) -> Structure:
    """从 .gz 文件读取 Structure。"""
    with gzip.open(gz_path, "rt") as f:
        content = f.read()
    return Structure.from_str(content, fmt="poscar")


def split_structure(struct: Structure, split_frac_z: float) -> tuple[Structure, Structure]:
    """按分数坐标 z 拆分结构，保持原始坐标和晶格不变。

    Args:
        struct: 输入结构
        split_frac_z: 拆分分数坐标 (0-1)

    Returns:
        (upper_struct, lower_struct): 上半部分(LDH), 下半部分(CMCH)
    """
    frac_z = struct.frac_coords[:, 2]
    upper_mask = frac_z > split_frac_z
    lower_mask = frac_z < split_frac_z

    # 上半部分 (LDH)
    upper_species = [struct.species[i] for i in range(len(struct)) if upper_mask[i]]
    upper_coords = struct.cart_coords[upper_mask]
    upper_sd = struct.site_properties.get("selective_dynamics", [None]*len(struct))
    upper_sd = [upper_sd[i] for i in range(len(struct)) if upper_mask[i]] if upper_sd[0] is not None else None

    # 下半部分 (CMCH)
    lower_species = [struct.species[i] for i in range(len(struct)) if lower_mask[i]]
    lower_coords = struct.cart_coords[lower_mask]
    lower_sd = struct.site_properties.get("selective_dynamics", [None]*len(struct))
    lower_sd = [lower_sd[i] for i in range(len(struct)) if lower_mask[i]] if lower_sd[0] is not None else None

    # 保持原始晶格不变
    lattice = struct.lattice

    upper_struct = Structure(lattice, upper_species, upper_coords, coords_are_cartesian=True)
    lower_struct = Structure(lattice, lower_species, lower_coords, coords_are_cartesian=True)

    # --- 重新设置 TTF 约束 (frac_z < DEFAULT_TTF_FRAC 为 TTF) ---
    for s, sd_list in [(upper_struct, upper_sd), (lower_struct, lower_sd)]:
        if sd_list is not None:
            new_sd = []
            for i in range(len(s)):
                if s.frac_coords[i, 2] < DEFAULT_TTF_FRAC:
                    new_sd.append([True, True, False])  # TTF
                else:
                    new_sd.append([True, True, True])   # TTT
            s.add_site_property("selective_dynamics", new_sd)
        else:
            # 如果原结构没有 selective_dynamics，按新的 frac_z 设置
            sd = []
            for i in range(len(s)):
                if s.frac_coords[i, 2] < DEFAULT_TTF_FRAC:
                    sd.append([True, True, False])
                else:
                    sd.append([True, True, True])
            s.add_site_property("selective_dynamics", sd)

    return upper_struct, lower_struct


def check_poscars() -> bool:
    """验证已生成的 POSCAR。"""
    files = [
        ("LDH_intrinsic.vasp", "CMCH_intrinsic.vasp"),
        ("LDH_S_doped.vasp", "CMCH_S_doped.vasp"),
        ("LDH_S_exposed.vasp", "CMCH_S_exposed.vasp"),
    ]

    print(f"\n{'='*70}")
    print(f"  POSCAR 验证")
    print(f"{'='*70}")
    print(f"  {'文件':<30s}  {'a':>8s}  {'b':>8s}  {'c':>8s}  {'atoms':>6s}  {'TTF':>6s}  {'TTT':>6s}")
    print(f"  {'─'*30}  {'─'*8}  {'─'*8}  {'─'*8}  {'─'*6}  {'─'*6}  {'─'*6}")

    all_ok = True
    ref_a = ref_b = None

    for ldh_f, cmch_f in files:
        for fname in [ldh_f, cmch_f]:
            path = OUTDIR / fname
            if not path.exists():
                print(f"  ❌ {fname:<30s}: 文件不存在")
                all_ok = False
                continue

            s = Structure.from_file(path)
            a, b, c = s.lattice.a, s.lattice.b, s.lattice.c

            if ref_a is None:
                ref_a, ref_b = a, b
            elif abs(a - ref_a) > 0.001 or abs(b - ref_b) > 0.001:
                print(f"  ⚠️ {fname:<30s}: 晶格不一致! a={a:.4f} b={b:.4f} vs ref a={ref_a:.4f} b={ref_b:.4f}")
                all_ok = False

            if "selective_dynamics" in s.site_properties:
                sd = s.site_properties["selective_dynamics"]
                ttf = sum(1 for d in sd if list(d) == [True, True, False])
                ttt = sum(1 for d in sd if list(d) == [True, True, True])
                other = len(sd) - ttf - ttt
                if other:
                    print(f"  ❌ {fname:<30s}: 有非 TTF/TTT 约束 (other={other})")
                    all_ok = False
            else:
                ttf, ttt = 0, len(s)

            z_range = f"[{s.cart_coords[:,2].min():.2f}, {s.cart_coords[:,2].max():.2f}]"
            print(f"  {'✅' if path.exists() else '❌'} {fname:<30s}  {a:>8.4f}  {b:>8.4f}  {c:>8.4f}  {len(s):>6d}  {ttf:>6d}  {ttt:>6d}  z={z_range}")

    if all_ok and ref_a:
        print(f"\n  ✅ 全部通过 — 所有体系晶格一致 (a={ref_a:.4f}, b={ref_b:.4f})")
    else:
        print(f"\n  ⚠️ 有错误需要修正")

    return all_ok


def run_split(config: dict) -> tuple[Structure, Structure]:
    """执行单个异质结的拆分。"""
    input_path = DATA / config["input"]
    print(f"\n=== 拆分 {config['label']} ===")
    print(f"  输入: {input_path}")

    struct = load_structure(input_path)
    print(f"  原始结构: {struct.formula}  atoms={len(struct)}")
    print(f"  晶格: a={struct.lattice.a:.4f} b={struct.lattice.b:.4f} c={struct.lattice.c:.4f}")
    print(f"  frac_z 范围: [{struct.frac_coords[:,2].min():.4f}, {struct.frac_coords[:,2].max():.4f}]")

    ldh_struct, cmch_struct = split_structure(struct, SPLIT_FRAC_Z)

    # 保存
    ldh_path = OUTDIR / config["ldh_out"]
    cmch_path = OUTDIR / config["cmch_out"]

    ldh_struct.to(str(ldh_path), fmt="poscar")
    cmch_struct.to(str(cmch_path), fmt="poscar")

    print(f"  ✅ LDH:  {config['ldh_out']}  atoms={len(ldh_struct)}  z=[{ldh_struct.cart_coords[:,2].min():.2f}, {ldh_struct.cart_coords[:,2].max():.2f}]")
    print(f"  ✅ CMCH: {config['cmch_out']}  atoms={len(cmch_struct)}  z=[{cmch_struct.cart_coords[:,2].min():.2f}, {cmch_struct.cart_coords[:,2].max():.2f}]")

    return ldh_struct, cmch_struct


def main():
    parser = argparse.ArgumentParser(description="拆分异质结为 LDH 和 CMCH")
    parser.add_argument("--intrinsic", action="store_true", help="仅拆分 intrinsic")
    parser.add_argument("--s-doped", action="store_true", help="仅拆分 S_doped")
    parser.add_argument("--s-exposed", action="store_true", help="仅拆分 S_exposed")
    parser.add_argument("--check", action="store_true", help="验证已生成的 POSCAR")
    args = parser.parse_args()

    if args.check:
        check_poscars()
        return

    # 确定要处理的配置
    if args.intrinsic or args.s_doped or args.s_exposed:
        targets = []
        if args.intrinsic:
            targets.append("intrinsic")
        if args.s_doped:
            targets.append("S_doped")
        if args.s_exposed:
            targets.append("S_exposed")
        configs = [c for c in HETERO_CONFIGS if c["name"] in targets]
    else:
        configs = HETERO_CONFIGS

    print(f"{'='*70}")
    print(f"  异质结拆分 (frac_z = {SPLIT_FRAC_Z})")
    print(f"  上半部分 → LDH, 下半部分 → CMCH")
    print(f"  输出目录: {OUTDIR}")
    print(f"{'='*70}")

    for cfg in configs:
        run_split(cfg)

    print(f"\n{'='*70}")
    print(f"  完成！")
    print(f"{'='*70}")

    # 最终验证
    check_poscars()


if __name__ == "__main__":
    main()