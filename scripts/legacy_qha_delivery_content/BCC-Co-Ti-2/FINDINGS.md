# BCC-Co-Ti-2 — 结论与证据

**一句话：数据齐全、可完整重建，但重建出来的曲线带确定性缺陷，不能当作可直接消费的生产观测。**

本目录是从 NAS `E:\nas-shared\Endmembers\16symbols-cxy-hmd\BCC-Co-Ti-2` 逐字节镜像
（487 文件 / 32,521,023 B，内容摘要 `f284934c…`）派生出来的交付层。镜像本身在
`../../legacy-dfpt-qha-pilot/BCC-Co-Ti/legacy-as-received/`，未被本交付改动一个字节。

---

## 1. 先说能用的部分

| 项 | 结论 | 证据 |
|---|---|---|
| E(V) 曲线 | 可重建，且**逐行可验证** | `curves/v-e.dat` 的每一行 = 该体积点 `R3f/OUTCAR` 的 `energy(sigma->0)`，10/10 点吻合到 <1e-3 meV |
| 声子层 | 可重建，且**逐点可验证** | 用 `dfpt/QHA-*/FORCE_CONSTANTS` 重算 thermal properties，10/10 点的 `zero_point_energy` 与随附 yaml 吻合到 <1e-6 kJ/mol；`num_modes`/`num_integrated_modes` 全部一致 |
| 体积网格规律 | 已还原到 9 位小数 | `a_i = a0·(V/V0)^(1/3)`，`V/V0 = linspace(0.95, 1.10, 10)`，`a0 = 2.9712807178297407 Å`；`v-e.dat` 记录的是 `a_i³` 四舍五入到两位小数 |
| 结构 | 已确定 | 2 原子立方胞，Ti (0,0,0) + Co (½,½,½)，即 CsCl(B2) 原型；`a0` 即 `R3/CONTCAR`（ISIF=3 弛豫收敛） |
| 配方 | 已确定 | `recipe/` 下：`IBRION=8`、`NSW=1`、`ENCUT=520`、`EDIFF=1e-5`、`ISMEAR=1`、`ISPIN=2`、`IALGO=48`；DFPT 超胞 `2×2×2`；DFPT k 网格 `MP 5×5×5`；R3 弛豫 k 网格 `MP 7×7×7`；热力学网格 `MP 24×24×24` |
| 方法 | 已确定 | DFPT，VASP `IBRION=8`，16 原子超胞，对称性约化到 2 个微扰（原子 1 与 9 的 x 方向），两个都解完 |
| 生产者 | 已确定 | VASP `5.4.4.18Apr17-6-g9f103f2a35`；phonopy `1.11.0`（`dfpt/QHA-*/phonopy.yaml`） |

也就是说：**这个端基是可以真正复算的**——所有输入都在镜像里，重算结果与随附 yaml 逐点吻合。
这与 BCC-Nb-Ti-2 完全相反（后者没有任何 VASP 层）。

## 2. 再说不能用的部分

`audit/audit.json`（`legacy-qha-audit/v1`，工具 `ht-vasp/scripts/legacy_qha_audit.py`）
状态为 **fail，37 条 fail-closed 失败**。按缺陷类别归并：

### 2.1 全部 10 个体积点：静态与 DFPT 都退回了原子叠加初猜

- 输入写的是 `ICHARG=1`（读 CHGCAR），但镜像里**没有任何 CHGCAR 文件**。
- `R3f/OUTCAR` 逐点打印 `charge density of overlapping atoms calculated`，
  `log-DFPT` 逐点打印 `WARNING: chargedensity file is incomplete`，
  DFPT 的 `OUTCAR` 打印 `initial charge density was supplied:` +
  `charge density of overlapping atoms calculated` + `keeping initial charge density in first step`。
- 后果：产生 `v-e.dat` 的 R3f 静态（`IBRION=-1, NSW=0, EDIFF=1e-4`）**从未使用过 R3 收敛的电荷密度**。
  10/10 点都是如此，因此这不是个别点的事故，是整条链的配方缺陷。
- 上游 `R3` 弛豫本身也不干净：`R3/R31` 报 `chargedensity file is incomplete`，
  `R3/R32` 报 `charge-density read from file: FeAl-bcc`（读进了别的体系的电荷密度），
  `R3/R48` 报 `dimensions on CHGCAR file are different`。

### 2.2 QHA-08（V = 27.98 Å³）：DFPT 基态 SCF 发散，声子层报废

`dfpt/QHA-08/log-DFPT` 的基态 SCF：RMM 1–12 正常收敛到 −142.7999 eV，
RMM 13 起电荷密度更新后开始震荡，RMM 56–60 彻底崩掉，最后一步
`RMM: 60  E = +0.435946909617E+04 eV, mag = 16.9657`，撞 NELM=60 退出。

四条互相独立的失败信号：

| 信号 | QHA-08 实测 | 健康体积实测 | 阈值 |
|---|---|---|---|
| DFPT 基态能量 vs 同体积 R3f 静态 | **+280 265 meV/atom** | 4.6 – 28.7 meV/atom | 100 meV/atom |
| `FORCE_CONSTANTS` 最大 \|Φ\| | **229.05 eV/Å²** | 6.08 – 10.15 eV/Å² | 50 eV/Å² |
| 声学求和规则相对残差 | **36.2 %** | 0.8 – 2.4 % | 10 % |
| 24³ 网格上虚频加权占比 | **50.000 %**（最低 −18.5362 THz） | 0 % | > 0 |

注意 `E(V)` 的这一行**本身没问题**：`v-e.dat` 第 8 行 −15.595737 与它自己的 `R3f/OUTCAR`
吻合到 5e-5 meV，而该 R3f 静态正常收敛（16 步，`aborting loop because EDIFF is reached`）。
坏的只是这个体积的声子层。但因为 QHA 要同时用 E(V) 和 F_ph(V,T)，
QHA-08 的声子层一坏，整条 Gibbs 曲线在高 T 就失去物理意义。

### 2.3 phonopy 会静默丢弃虚频，记录的 ZPE 只覆盖部分模式

QHA-08 的 `thermal_properties-from-dfpt-08.yaml`：`num_modes: 82944`，
`num_integrated_modes: 41472`。差的 41472 就是被丢掉的虚频加权数。
phonopy 不报错、不警告，`cutoff_frequency` 也不出现在 1.11.0 的输出里。
于是那一行的 `zero_point_energy: 18.42 kJ/mol` 是**残缺模式集**的 ZPE，
与其它体积的 ~6 kJ/mol 不可比。

### 2.4 QHA-01 / 03 / 04：DFPT 基态 SCF 撞 NELM=60

`dfpt/QHA-01` 的基态 SCF 在第 12 步（固定电荷密度阶段）曾达到 dE = 2.1e-6 eV，
但第 13 步电荷密度一更新就跳到 +65 eV，之后在 −124.662 ± 0.002 eV 附近游走 48 步，
到 NELM=60 仍未满足 `EDIFF=1e-5`（末步 dE = 1.4e-4 eV，最后 10 步能量跨度 1.5 meV）。
QHA-03 / 04 同类（跨度 0.5 / 0.6 meV）。
这三点与 QHA-08 在**数量级上完全不同**（meV vs 1e6 meV），但同属"基态未达 EDIFF"。

> **这一条对 A 线本身很重要：**
> VASP 5.4.4 在电子步耗尽 NELM 时**不打印任何 `unconverged` 标记**，
> 所以 BCC-Ti-V 试点 `cmd_verify` 那套"数标记"的 fail-closed 门
> （`converged_scf_marker_count > 0`、`unconverged_scf_marker_count == 0`）
> **看不见 QHA-08 这种失败**——它的 3 个 `aborting loop because EDIFF is reached`
> 全部来自线性响应阶段，基态发散的那一段一个标记都没有。
> 本审计新增的"DFPT 基态能量 vs 同体积 R3f 静态"和"基态 SCF 是否撞 NELM"两道门才抓得住。

### 2.5 报告出来的 Gibbs 曲线不光滑

`curves/gibbs-temperature.dat`（0–2000 K，步长 10 K，201 点）：

- 1920 → 1930 K 处跳变 **0.515431788536624 eV/cell**
- 另有 14 处相邻跳变 > 0.2 eV/cell（合计 15 处）

配合 §2.2 看：高 T 段 Gibbs 逐点选不同的最优体积，而 QHA-08 的 F_ph 是坏的，
曲线就在体积之间跳。**这条曲线只能当诊断对象，不能当交付曲线。**

### 2.6 同一几何下 R3 与 R3f 相差 51.7 meV/cell

| 计算 | 设置 | energy(sigma→0) |
|---|---|---|
| `R3/OUTCAR`（最终离子步，a = 2.9712807178297407 Å） | ISMEAR=0, EDIFF=1e-5, ICHARG=1 | −15.59210606 eV |
| `dfpt/QHA-04/R3f/OUTCAR`（同几何） | ISMEAR=1, SIGMA=0.1, EDIFF=1e-4, ICHARG=1 | −15.64379742 eV |
| 差 | | **51.69 meV/cell** |

`v-e.dat` 用的是后者。这不是相对量的小修正，是同一几何两套设置
（含 smearing 类型不同）的系统差；在把 `v-e.dat` 的绝对能量接到任何参考态之前必须先处理它。

### 2.7 输入模板残留（不影响数字，但影响"按名字推断语义"）

| 文件 | 内容 | 说明 |
|---|---|---|
| `dfpt/QHA-*/INCAR`、`recipe/INCAR-Phonopy` | `SYSTEM=laves-C14` | 模板残留 |
| `R3/INCAR`、`recipe/INCAR-R31` | `SYSTEM=FeAl-bcc` | 模板残留 |
| `dfpt/QHA-*/mesh.conf`、`recipe/mesh.conf` | `ATOM_NAME= Nb` | 模板残留 |
| `dfpt/QHA-*/POSCAR-unitcell`、`R3/CONTCAR` | 标题行 `Cl Cs` | 恰好描述了正确的 CsCl 原型，但它是模板名 |
| `R3/CONTCAR` | 首行 `BCC-A2` | **与内容不符**：2 原子有序胞不是无序 A2 |
| `POTCAR` | 含 Ti、Co、**W** 三个块 | W 未被读取（POSCAR 只声明 2 个物种） |

`R3/R31` 还打印了 `WARNING: type information on POSCAR and POTCAR are incompatible`，
VASP 用 POTCAR 覆盖了 POSCAR 的物种名。真实物种顺序由 POTCAR 决定（Ti, Co），
与坐标一致，所以数字是对的。

---

## 3. 给 B 线的适用性建议

- `applicability`：**`degraded_defective`**，默认 **`needs_review`**，不得直接进 alignment / fit / TDB。
- 若要使用，最低限度需要：
  1. 排除 QHA-08（或整条声子层重算），否则 Gibbs 曲线在 T ≳ 1900 K 不可信；
  2. 承认 `v-e.dat` 的 E 来自"无 CHGCAR 的原子叠加初猜 + EDIFF=1e-4"的静态，
     且与 `R3/OUTCAR` 在同几何下相差 **51.69 meV/cell**（见 §2.6）；
  3. `reference_role` 与 `energy_reference_state` 在本端基内**无法解析**，
     必须由 profile 或用户显式给定，不得由目录名或格式推定。
- 全部 10 个点的 `E(V)` 与 `thermal_properties` 的**内部一致性**是可靠的
  （§1 两行）；不可靠的是它们的**物理质量**。二者不要混为一谈。

## 4. 复算所需的东西都在

`raw/RAW_INDEX.json` 列出 240 个原始 VASP 层文件（12.88 MiB）的路径与 sha256，
包括每个体积点的 `INCAR` / `POSCAR` / `POSCAR-unitcell` / `POTCAR` / `OUTCAR` /
`vasprun.xml` / `log-DFPT` / `FORCE_CONSTANTS` / `mesh.yaml` / `disp.yaml` / `DYNMAT`，
以及 `R3/` 的整套。本交付**故意不复制**这些字节（只给指针），以免出现两份可能漂移的副本；
但配方层（`INCAR*` / `KPOINTS*` / `mesh.conf`）已复制进 `recipe/`，
因为它们是"这条路线的 producer 身份如何证明"的直接证据。
