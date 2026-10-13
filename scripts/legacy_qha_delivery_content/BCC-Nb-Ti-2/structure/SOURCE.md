# BCC-Nb-Ti-2 / structure — 本目录里没有任何东西来自 NAS

本端基的镜像里**没有**结构文件（无 `POSCAR` / `CONTCAR` / `POSCAR-unitcell`）。
所以这里放的不是 legacy 数据，而是**从 ht-vasp 导入的工作假设**，
并且每一件都标明了来源。

## from-ht-vasp/POSCAR-unitcell

- 生成自：`ht-vasp/data/endmembers/BCC-Nb-Ti/BCC-Nb-Ti-static.json`
  （ht-vasp 的 atomate2 `staticflow` 弛豫结果）
- 生成方式：`scripts/build_legacy_qha_delivery.py --external-structure-json`，
  用 `pymatgen.core.Structure.from_dict(...).to(fmt="poscar")` 导出
- 内容：2 原子立方胞，Nb (1,1,1) ≡ (0,0,0) 与 Ti (½,½,½)
  （pymatgen 导出时把 Nb 坐标写成 0.9999999895，等价于 0）
- a = 3.266876449678747 Å，V = 34.865689687333266 Å³
- **它没有任何 legacy provenance。** 它落在 legacy 体积区间 33.16–37.23 Å³ 内、
  比 legacy 的离散极小体积 34.91 Å³ 低 0.127 %，这只能作为工作假设。
- 用途：如果要复现 legacy 曲线，这是目前唯一有依据的起点；
  但"复现"需要新跑 VASP，见 `../FINDINGS.md` §7。

## from-ht-vasp/BCC-Nb-Ti.vasp

- 来源：`ht-vasp/data/poscars/BCC-Nb-Ti.vasp`
- **不要把它当成这个端基的胞。** 它的晶格常数 2.85971129 Å 是 ht-vasp
  `data/poscars/` 下**所有 BCC-* 模板共用的占位值**：
  `BCC-Ti-V.vasp`、`BCC-Al-Nb.vasp`、`BCC-Nb-Nb.vasp`、`BCC-Ti-Ti.vasp`
  的第一行晶格常数与它**逐位相同**。它给出 V = 23.39 Å³（2 原子），
  即 11.7 Å³/atom，对 bcc Ti-Nb 来说过密，与 legacy 的 16.6–18.6 Å³/atom 不符。
- 另外它的标题行写 `Ti1 Nb1` 而物种行写 `Nb Ti`，两者不一致（模板残留）。
- 它放进来只是为了记录"曾经考虑过这个来源"以及为什么排除它。

## 一句话

- **legacy 提供的结构：无。**
- **ht-vasp 提供的：一个弛豫胞（工作假设）+ 一个被排除的模板。**
- **本端基自己的锚点：不可恢复。**
