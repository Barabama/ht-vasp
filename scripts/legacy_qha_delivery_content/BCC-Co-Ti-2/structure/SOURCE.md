# BCC-Co-Ti-2 / structure — 全部来自 legacy 镜像，逐字节

本端基的镜像**自带结构**：根目录有 `POSCAR` 与 `POSCAR-org`，`R3/` 有弛豫的
`CONTCAR`，每个体积点目录里有 `POSCAR-unitcell`。本目录是它们的副本。

> 镜像**根目录没有** `POSCAR-unitcell`。锚点胞因此以 `relaxed-CONTCAR`
> （= 根目录 `POSCAR`）为准，逐体积的胞在 `volumes/` 下。

## relaxed-CONTCAR 与 POSCAR

- 来源：`legacy-as-received/R3/CONTCAR` 与 `legacy-as-received/POSCAR`
- 两者逐位相同：a = 2.9712807178297407 Å，V = 26.231978868335762 Å³
- 2 原子立方胞，Ti (0,0,0) + Co (½,½,½)
- 这是 `ISIF=3, IBRION=2, NSW=50, EDIFF=1e-5` 弛豫的收敛结果
  （`R3/OUTCAR` 中有一次 `reached required accuracy`）
- 标题行是 `Cl Cs` —— 模板残留，但恰好命名了正确的 CsCl 原型；
  `R3/CONTCAR` 首行写 `BCC-A2`，**与内容不符**：2 原子有序胞不是无序 A2，是 B2/CsCl 型
- 物种顺序由 POTCAR 决定（Ti, Co），**不是**由标题行决定：
  `R3/R31` 打印过 `WARNING: type information on POSCAR and POTCAR are incompatible`，
  VASP 用 POTCAR 覆盖了 POSCAR 的物种名
- 弛豫过程不干净：重启过至少两次，`R3/R31` 报 `chargedensity file is incomplete`，
  `R3/R32` 报 `charge-density read from file: FeAl-bcc`（读入别的体系的电荷密度），
  `R3/R48` 报 `dimensions on CHGCAR file are different`。
  几何本身是收敛的，路径不是。

## POSCAR-org

- 来源：`legacy-as-received/POSCAR-org`
- 弛豫的**起点**，a = 2.850 Å（V = 23.15 Å³），比收敛值小 4.1 %
- 保留它是因为它说明 R3 弛豫的搜索范围，不是因为它有物理意义

## volumes/POSCAR-unitcell-QHA-01 … -10

- 来源：`legacy-as-received/dfpt/QHA-NN/POSCAR-unitcell`
- 这 10 个是**每个体积点 DFPT 实际使用的胞**，也就是对应 `FORCE_CONSTANTS` 的胞
- 与 `curves/v-e.dat` 的体积列一致（v-e.dat 记录的是 a³ 四舍五入到两位小数）
- 体积网格规律已还原到 9 位小数：
  `a_i = a0·(V/V0)^(1/3)`，`V/V0 = linspace(0.95, 1.10, 10)`，
  `a0 = 2.9712807178297407 Å`（即弛豫胞）

| 体积点 | V / Å³ (v-e.dat) | a / Å | V/V0 |
|---|---|---|---|
| QHA-01 | 24.92 | 2.9209102968535019 | 0.950000000 |
| QHA-02 | 25.36 | 2.9378927150550966 | 0.966666667 |
| QHA-03 | 25.79 | 2.9546810378526338 | 0.983333333 |
| QHA-04 | 26.23 | 2.9712807178297407 | 1.000000000 |
| QHA-05 | 26.67 | 2.9876969665442803 | 1.016666667 |
| QHA-06 | 27.11 | 3.0039347689395477 | 1.033333333 |
| QHA-07 | 27.54 | 3.0199988966763711 | 1.050000000 |
| QHA-08 | 27.98 | 3.0358939204826902 | 1.066666667 |
| QHA-09 | 28.42 | 3.0516242216071681 | 1.083333333 |
| QHA-10 | 28.86 | 3.0671940024545625 | 1.100000000 |

## 超胞

DFPT 在 **2×2×2 超胞**上做（16 原子，a_super = 2·a_unit）：
`dfpt/QHA-NN/POSCAR` 是超胞，`POSCAR-unitcell` 是原胞，`phonopy.yaml` 里
`supercell_matrix` 是 `diag(2,2,2)`。超胞 POSCAR 没有复制到本交付，
但它在 `raw/RAW_INDEX.json` 里有指针。
