#!/usr/bin/env python3
"""构建 CMCH(010) @ CoNi-LDH-S(001) 异质结

CMCH = 基底(substrate, 下层), LDH = 薄膜(film, 上层)

Interface.from_slabs 会 flip film (film_coords[:,2] *= -1.0):
  film原底面 → 异质结顶面(真空侧)
  film原顶面 → 异质结底面(界面侧)

因此:
  slab-001 (S在原底面 z≈2.7): flip后 S在顶面 → S远离CMCH → 掺杂异质结 ✅
  slab-00-1 (S在原顶面 z≈5.4): flip后 S在底面 → S靠近CMCH → 暴露掺杂异质结 ✅
"""

import json
from pathlib import Path

import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.core.interface import Interface

DATA = Path("/workspace/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S/data")
OUTDIR = DATA / "heterostructures"
OUTDIR.mkdir(exist_ok=True)


def load_optimized(name):
    with open(DATA / name / "static_out.json") as f:
        d = json.load(f)
    return Structure.from_dict(d["output"]["structure"])


def swap_ab(slab):
    """交换 slab 的 a/b 轴"""
    old = slab.lattice
    new_lat = Lattice.from_parameters(
        old.b, old.a, old.c,
        old.alpha, old.beta, old.gamma,
    )
    frac = slab.frac_coords.copy()
    frac[:, [0, 1]] = frac[:, [1, 0]]
    return Structure(
        new_lat, slab.species, frac,
        site_properties=slab.site_properties,
        coords_are_cartesian=False,
    )


def strain_to(slab, target_a, target_b):
    """应变 slab 面内 a/b 到目标值 (保持分数坐标)"""
    old = slab.lattice
    new_lat = Lattice.from_parameters(
        target_a, target_b, old.c,
        old.alpha, old.beta, old.gamma,
    )
    return Structure(
        new_lat, slab.species, slab.frac_coords,
        site_properties=slab.site_properties,
        coords_are_cartesian=False,
    )


# -----------------------------------------------------------
# 1. 读取优化后 slab
# -----------------------------------------------------------
cmch = load_optimized("CoMnH2CO5-slab")           # 基底 CMCH(010)
ldh = load_optimized("CoNiOH2-slab")               # 薄膜 LDH(001) 本征
ldh_s = load_optimized("CoNiOH2S-noH-slab")        # 薄膜 LDH-S(001) S在原底面
ldh_s_flip = load_optimized("CoNiOH2S-noH-slab-flip")  # 薄膜 LDH-S(00-1) S在原顶面

# 修复: slab-00-1 的 a/b 被 MS 交换了，换回来对齐 CMCH
ldh_s_flip = swap_ab(ldh_s_flip)

print("=== 原始 slab (修正后) ===")
for label, s in [
    ("CMCH 基底", cmch), ("LDH 本征", ldh),
    ("LDH-S(001)", ldh_s), ("LDH-S(00-1) a/b已换", ldh_s_flip),
]:
    print(f"  {label}: a={s.lattice.a:.3f} b={s.lattice.b:.3f} γ={s.lattice.gamma:.2f}°")

# -----------------------------------------------------------
# 2. 平均晶格匹配
# -----------------------------------------------------------
avg_a = (cmch.lattice.a + ldh.lattice.a) / 2
avg_b = (cmch.lattice.b + ldh.lattice.b) / 2

print(f"\n=== 晶格匹配 ===")
print(f"  a: {cmch.lattice.a:.3f} / {ldh.lattice.a:.3f} → avg={avg_a:.3f}  (各~6%应变)")
print(f"  b: {cmch.lattice.b:.3f} / {ldh.lattice.b:.3f} → avg={avg_b:.3f}  (各~0.1%应变)")

cmch_m = strain_to(cmch, avg_a, avg_b)
ldh_m = strain_to(ldh, avg_a, avg_b)
ldh_s_m = strain_to(ldh_s, avg_a, avg_b)
ldh_s_flip_m = strain_to(ldh_s_flip, avg_a, avg_b)

# -----------------------------------------------------------
# 3. 构建 3 个异质结
# -----------------------------------------------------------
GAP = 2.0
VAC_BOTTOM = 2.0    # 基底底部真空
VAC_TOP = 6.0       # 薄膜顶部真空
TOTAL_C = 20.0

# 计算实际需要的 vacuum_over_film (Interface.from_slabs 的 c = material + gap + vacuum_over_film)
# material = sub_thick + film_thick
sub_thick = cmch_m.cart_coords[:, 2].max() - cmch_m.cart_coords[:, 2].min()
film_thick = ldh_m.cart_coords[:, 2].max() - ldh_m.cart_coords[:, 2].min()
extra_vac = VAC_BOTTOM + VAC_TOP  # 额外真空 (Interface.from_slabs center后上下均分)

configs = {
    "intrinsic": {
        "film": ldh_m,
        "label": "本征异质结",
        "out": "hetero_intrinsic.vasp",
    },
    "s_doped": {
        "film": ldh_s_m,
        "label": "掺杂异质结 (S在真空侧, 远离CMCH)",
        "out": "hetero_S_doped.vasp",
    },
    "s_exposed": {
        "film": ldh_s_flip_m,
        "label": "暴露掺杂异质结 (S在界面侧, 靠近CMCH)",
        "out": "hetero_S_exposed.vasp",
    },
}

for key, cfg in configs.items():
    iface = Interface.from_slabs(
        substrate_slab=cmch_m,
        film_slab=cfg["film"],
        gap=GAP,
        vacuum_over_film=extra_vac,
        center_slab=True,
    )

    # 调整: shift 使基底底部 = VAC_BOTTOM + 设总 c = TOTAL_C
    zs = iface.cart_coords[:, 2]
    labels = np.array(iface.site_properties["interface_label"])

    sub_zs = zs[labels == "substrate"]
    film_zs = zs[labels == "film"]

    shift = VAC_BOTTOM - sub_zs.min()
    print(f"    z min: sub = {sub_zs.min()}, film = {film_zs.min()}, shift = {shift}")

    coords_new = iface.cart_coords.copy()
    coords_new[:, 2] += shift

    new_lat = Lattice.from_parameters(
        iface.lattice.a, iface.lattice.b, TOTAL_C,
        iface.lattice.alpha, iface.lattice.beta, iface.lattice.gamma,
    )
    iface_final = Structure(
        new_lat, iface.species, coords_new,
        site_properties=iface.site_properties,
        coords_are_cartesian=True,
    )

    # 选择性动力学: CMCH 底部 ~half 固定(F F F), 其余可移动(T T T)
    zs_final = iface_final.cart_coords[:, 2]
    labels_final = np.array(iface_final.site_properties["interface_label"])
    sub_mask = labels_final == "substrate"
    sub_zs = zs_final[sub_mask]

    # CMCH z 中点: 下半固定, 上半可动
    sub_z_mid = (sub_zs.min() + sub_zs.max()) / 2
    n_fixed = int(np.sum(sub_zs < sub_z_mid))
    n_free = len(iface_final) - n_fixed

    sd = []
    for i in range(len(iface_final)):
        if sub_mask[i] and zs_final[i] < sub_z_mid:
            sd.append([False, False, False])  # F F F = 固定
        else:
            sd.append([True, True, True])     # T T T = 可移动

    iface_final.add_site_property("selective_dynamics", sd)

    outpath = OUTDIR / cfg["out"]
    iface_final.to(outpath, fmt="poscar")

    # 报告
    film_zs = zs_final[~sub_mask]
    true_gap = film_zs.min() - sub_zs.max()
    print(f"\n  ✓ {cfg['label']}")
    print(f"    POSCAR: {outpath}")
    print(f"    atoms={len(iface_final)}, formula={iface_final.composition.formula}")
    print(f"    CMCH z: {sub_zs.min():.2f}-{sub_zs.max():.2f} (thick={sub_zs.max()-sub_zs.min():.2f})")
    print(f"    LDH  z: {film_zs.min():.2f}-{film_zs.max():.2f} (thick={film_zs.max()-film_zs.min():.2f})")
    print(f"    界面间隙(O/H终端): {true_gap:.2f} Å")
    print(f"    底真空: {sub_zs.min():.2f}, 顶真空: {TOTAL_C - film_zs.max():.2f}")
    print(f"    选择性动力学: {n_fixed} 固定(F F F) + {n_free} 可动(T T T)")
    print(f"    固定范围: CMCH z < {sub_z_mid:.2f} Å (底部下半)")

    # S 验证
    s_idx = [i for i, s in enumerate(iface_final) if s.species_string == "S"]
    if s_idx:
        s_z = zs_final[s_idx[0]]
        film_mid = (film_zs.min() + film_zs.max()) / 2
        d_to_interface = s_z - sub_zs.max()
        if s_z < film_mid:
            print(f"    S z={s_z:.2f} (距界面{d_to_interface:.2f}Å) → **界面侧** (靠近CMCH) ✓")
        else:
            print(f"    S z={s_z:.2f} → **真空侧** (远离CMCH) ✓")

print(f"\n=== 完成 ===")
