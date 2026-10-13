# 传统 DFPT-QHA 目录布局：A 线提案（**未冻结**）

> **本文件不是契约。**
>
> `em-tdb-v2` 计划 v2 §B.5.7 明确：具体 DFPT 数据结构**尚未产生**，
> 目录布局、列顺序、能量单位、温度网格、参考态、原子数归一化方式
> **一律不得提前定死**；在 §B.5.7 待确认清单 1–9 全部由用户确认之前，
> 任何 adapter 都不得被标记为已实现，也不得用旧 `em-tdb` 的行为冒充已冻结契约。
>
> 所以本文件做的是三件事，而不是定义标准：
>
> 1. 把 A 线**实际交付的目录**摊开，逐项说明每个选择；
> 2. 对 §B.5.7 的 9 个待确认项，给出**A 线手上有的证据**（有就给，没有就说没有）；
> 3. 明确指出哪些项**只能由用户决定**，A 线不替用户定。
>
> 同时遵守 §B.5.8：**不得按目录名推断物理语义**。本交付里每一个
> 物理断言都带自己的证据；`endpoint_label` 只作为 NAS 来源标识记录。

---

## 1. 交付树实际长什么样

```text
data/legacy-dfpt-qha-delivery/
├── README.md                      交接说明（给 B 线）
├── DELIVERY.json                  机器可读总索引
├── LAYOUT_PROPOSAL.md             本文件
└── <endpoint_label>/              = NAS 目录名，仅作来源标识
    ├── DATASET.json               逐字段 claim + evidence + status
    ├── FINDINGS.md                人能读的结论与缺陷
    ├── INVENTORY.json             机器生成：本目录所有文件的路径/大小/sha256
    ├── audit/audit.json           fail-closed 审计报告原文
    ├── curves/                    *.dat 逐字节副本（扁平化，同名取根目录那份）
    ├── phonons/                   thermal_properties-*.yaml 逐字节副本
    ├── recipe/                    INCAR / KPOINTS / mesh.conf 逐字节副本 + POTCAR 身份摘录
    ├── raw/RAW_INDEX.json         原始 VASP 层的**指针**（路径 + sha256），不复制字节
    └── structure/                 胞；逐体积胞放在 structure/volumes/
```

设计取向：

- **镜像是不可变事实源。** `legacy-as-received/` 一个字节都不改；本交付里的每个
  副本都在生成时与镜像自己的 `source-manifest.json` 对过 sha256
  （`INVENTORY.json.files[].matches_mirror_manifest`）。
- **原始 VASP 层只给指针。** 复制会产生第二份可能漂移的副本，
  且 Co-Ti 那份有 12.88 MiB、Nb-Ti 那份根本没有。
- **配方层要复制。** `INCAR*` / `KPOINTS*` / `mesh.conf` 很小，而且是
  §B.5.7 第 9 项（producer 身份）的直接证据，所以进 `recipe/`。
  `POTCAR` 不复制（大且有许可），只在 `recipe/POTCAR_IDENTITY.txt` 里记
  `TITEL` / `POMASS` / `ZVAL` / `VRHFIN` 行，原文留在镜像里。
- **审计报告随交付走。** B 线不必重跑 phonopy 就能拿到 per-volume 的适用性状态。
- **人能读的与机器可读的分开。** `FINDINGS.md` 讲判断，`DATASET.json` 讲字段。

## 2. 逐项回答 §B.5.7 待确认清单

下表的"状态"用词：

- `已有证据` — A 线有文件级证据，可复现；
- `部分` — 有部分证据，其余不可恢复；
- `用户决策` — A 线不替用户定。

### 1. 目录命名与层级（是否沿用 `{PHASE}-{elem}-{elem}[-{atom_count}]`？）

| | |
|---|---|
| 本交付做法 | 顶层目录用 **NAS 目录名原文**（`BCC-Nb-Ti-2` / `BCC-Co-Ti-2`），并在 `DATASET.json.endpoint_label_is_a_nas_directory_name = true` 标注它是来源标识而非物理断言。内部层级固定为 `curves/ phonons/ recipe/ structure/ audit/ raw/`。 |
| 证据 | 两个 NAS 路径：`…\16symbols-cxy-hmd\BCC-Nb-Ti-2`、`…\BCC-Co-Ti-2`（见各自 `source-manifest.json`）。 |
| 状态 | **用户决策**。A 线刻意不让目录名承载语义，因为 §B.5.8 明确禁止由目录名推断体系；而两个端基的实际情况正好证明了这条禁令的必要性——`BCC-Nb-Ti-2` 的文件里**一个元素名都没有**。 |
| 附带发现 | legacy 自己的**文件标题行**也不可信：Co-Ti 的 `POSCAR-unitcell` 标题写 `Cl Cs`、`R3/CONTCAR` 首行写 `BCC-A2`（与内容不符）、`INCAR` 写 `SYSTEM=FeAl-bcc` 和 `laves-C14`、`mesh.conf` 写 `ATOM_NAME= Nb`；`POTCAR` 里还多一个从未被读取的 W 块。**按名字推断语义在这批数据上会直接出错。** |

### 2. 能量文件格式与列语义

| | |
|---|---|
| 本交付做法 | 不改格式，逐字节保留 `v-e.dat`（2 列空白分隔）与 `gibbs-temperature.dat`（2 列空白分隔）。列语义记在 `DATASET.json.claims.units`。 |
| 证据 | Co-Ti：`v-e.dat` 第 i 行 = 第 i 个体积点 `R3f/OUTCAR` 的 `energy(sigma->0)`，10/10 点吻合到 <1e-3 meV（`audit.json.volumes[*].v_e_row_matches_r3f_static_meV`）。**注意是 `energy(sigma->0)`，不是 `energy without entropy`**，两者在同一行上相差约 0.5 meV。Nb-Ti：无 VASP 层，无法交叉验证，故其单位在 `DATASET.json` 里标为 `partially_resolved`（降级）。 |
| 状态 | **已有证据**（Co-Ti）/ **部分**（Nb-Ti） |
| 遗留 | 其余 `.dat`（`helmholtz-volume.dat`、`Cv-volume.dat`、`entropy-volume.dat` 等）的列语义 A 线**没有独立验证**，只是照原样交付。 |

### 3. 能量归一化基准（per supercell / per atom / per reduced formula unit）

| | |
|---|---|
| 本交付做法 | 记录为**每 2 原子胞**（= 每 TiCo 或 TiNb 化学式单元）。不换算成别的基准。 |
| 证据 | Co-Ti：`phonopy.yaml` 与 `POSCAR-unitcell` 都是 2 原子胞；`thermal_properties-*.yaml` 写 `natom: 2`；`v-e.dat` 的体积 24.92–28.86 Å³ 只有按 2 原子胞解释才落在 bcc 金属的合理区间（12.5–14.4 Å³/atom）。Nb-Ti：只有 `natom: 2` 与体积区间（16.6–18.6 Å³/atom），同样一致。 |
| 状态 | **已有证据**（2 原子胞这一点）/ **用户决策**（下游要不要改成 per-atom 或 per-mole） |
| 注意 | §B.5.6 指出 `composition_reduced` 是 v2 的硬门槛（`qha_composition_missing` 直接 error）。Co-Ti 的组成有文件级证据；**Nb-Ti 的组成在本端基内无法证明**（无元素名），必须由 profile 或用户给出。 |

### 4. 温度网格与有效温区

| | |
|---|---|
| 本交付做法 | 原样交付，并在 `DATASET.json.claims.temperature_grid` 里逐文件记录 `[first, last, step, n_points]`。**不裁剪**到旧 `em-tdb` 硬编码的 100–2900 K。 |
| 证据 | Co-Ti：`gibbs-temperature.dat` 0–2000 K / 步 10 / 201 点；`thermal_properties-*.yaml` 0–3000 K / 步 10 / 301 点。Nb-Ti：`gibbs-temperature.dat` 0–2960 K / 步 10 / 297 点；`thermal_properties-*.yaml` 0–3000 K / 步 10 / 301 点。 |
| 状态 | **已有证据**（网格本身）/ **用户决策**（有效温区上限） |
| 附带发现 | **两条网格的上限互不相同，而且两个端基之间也不同**（Gibbs 2000 vs 2960 K）。`thermal_properties` 到 3000 K 而 Gibbs 只到 2000/2960 K，说明 Gibbs 扫描另设了一个上限。这两个上限都不是物理熔点判据：Co-Ti 体系熔点远低于 2000 K，Nb-Ti 体系也低于 2960 K。**有效温区必须由用户给，不能由文件上限推。** |

### 5. 参考态（是否同样以 unary50 GHSER 为 SER 参考）

| | |
|---|---|
| 本交付做法 | `DATASET.json.claims.reference_role` 与 `energy_reference_state` 在两个端基里**都写 `null` / `unresolved`**。 |
| 证据 | 无。两个端基里**没有任何** SER / GHSER / 元素参考的声明。Co-Ti 的 `v-e.dat` 是绝对 DFT 总能量（约 −15.6 eV/2 原子），Nb-Ti 的约 −17.7 eV/2 原子，都没有相对量。 |
| 状态 | **用户决策**（且 A 线认为这是最需要用户拍板的一项） |
| 为什么不能猜 | 若按旧 `em-tdb` 的 unary50 GHSER 处理，需要先解决 Co-Ti 的 51.69 meV/cell 系统差（同一几何下 `R3` 与 `R3f` 两套设置的能量差，见 `BCC-Co-Ti-2/structure/SOURCE.md` 末节与 `BCC-Co-Ti-2/FINDINGS.md` §2.6）。在参考态定下来之前，这个差是悬空的。 |

### 6. 与 atomate2 路线共存时的优先级与冲突消解

| | |
|---|---|
| 本交付做法 | A 线**不替用户定优先级**，但把冲突摆出来。 |
| 证据 | ht-vasp 里存在**同名端基的独立现代 QHA 数据**：`data/endmembers/BCC-Nb-Ti/qhaflow/` 与 `BCC-Nb-Ti-qha.json`。其 Gibbs(T=0) = −17.83702657957205 eV/f.u.，而 legacy NAS 端基为 −17.679267504318961 eV/f.u.，**相差 0.158 eV/f.u.** |
| 状态 | **用户决策**（对应计划里的 R13 风险） |
| A 线建议 | 两者**不是同一次计算**，差 0.158 eV/f.u. 远超后处理差异所能解释。不要把 legacy 端基当作 atomate2 路线的"降级版"来填补缺失字段，也不要反过来。若两者都要用，必须显式记录来源并分别进入 provenance 链。 |

### 7. 声子虚频 / 结构失败如何表达为 applicability 状态

| | |
|---|---|
| 本交付做法 | `audit.json` 用三态 `status`：`pass` / `fail` / `unverifiable`。`fail` 只在**有失败**时给；`unverifiable` 在**没有 VASP 层可查**时给——没有东西可查**不等于**通过。逐体积的量化状态放在 `volumes[*]` 与 `warnings[]`。 |
| 证据 | Co-Ti：`fail`，37 条失败，含 QHA-08 的 DFPT 发散（4 条独立信号：基态能量偏离 280 265 meV/atom、`FORCE_CONSTANTS` 最大 \|Φ\| 229.05 eV/Å²、声学求和规则残差 36.2 %、24³ 网格上 50.000 % 虚频）。Nb-Ti：`unverifiable`，0 条失败、1 条警告（6/8 个点的 `num_integrated_modes < num_modes`，即 phonopy 静默丢弃了虚频）。 |
| 状态 | **已有证据**（A 线能给到 per-volume 量化状态） |
| 关键教训（A 线自己的） | VASP 5.4.4 在电子步耗尽 `NELM` 时**不打印任何 `unconverged` 标记**。BCC-Ti-V 试点 `cmd_verify` 那套"数标记"的门（`converged_scf_marker_count > 0` 且 `unconverged_scf_marker_count == 0`）**看不见 QHA-08**——它的 3 个 `aborting loop because EDIFF is reached` 全来自线性响应阶段，基态发散那一段一个标记都没有。新审计改用两道物理门："DFPT 基态能量 vs 同体积 R3f 静态"与"基态 SCF 是否撞 NELM"。**B 线的 applicability 判定不要只依赖标记计数。** |

### 8. 与 §B.5.6 的曲线候选链是否共用同一个曲线 reader

| | |
|---|---|
| 本交付做法 | 交付层刻意把曲线文件**扁平化**到 `<endpoint>/curves/<basename>`，同名文件取镜像根目录那份，并在 `INVENTORY.json.curve_alternate_copies` 里记录被跳过的那份及其 sha256 与"是否逐字节相同"。 |
| 证据 | Co-Ti 里 `v-e.dat` 同时存在于镜像根目录与 `QHA-TiCo/`，两者 sha256 相同（`d73ce544…`）；`QHA-TiCo/thermal_properties-0N.yaml` 与 `dfpt/QHA-NN/thermal_properties-NN.yaml` 也逐字节相同。 |
| 状态 | **用户决策**（§B.5.7 第 8 项建议"共用 reader、分开 adapter"，A 线同意但无权定） |
| A 线建议 | 共用 reader 是合理的（两列的空白分隔数字，格式一致）；但 adapter 必须分开，因为**元数据来源不同**：legacy 端基的组成/方法/参考态要靠 profile 或用户，atomate2 路线的这些字段来自 `qha.json` 自身。 |

### 9. 该路线的 producer / workflow 身份如何证明

| | |
|---|---|
| 本交付做法 | 在 `DATASET.json.claims.producer` 与 `claims.calculation_method` 里逐项给证据，并且**只在有证据时给值**。 |
| 证据（Co-Ti） | 方法 = DFPT：`recipe/INCAR-Phonopy` 与 `dfpt/QHA-*/INCAR` 里 `IBRION=8, NSW=1`；`OUTCAR` 里 `Linear response parameters`、`Linear response progress:`、`Degree of freedom: 1/2` 与 `2/2`；`log-DFPT` 里 `Linear response reoptimize wavefunctions to high precision`；`FORCE_CONSTANTS` 是完整 16×16×3×3。生产者 = VASP `5.4.4.18Apr17-6-g9f103f2a35` + phonopy `1.11.0`（来自 `phonopy.yaml` 的 `version: 1.11.0`）。 |
| 证据（Nb-Ti） | **无。** 方法、生产者、版本全部 `unresolved`。yaml 字段集合与 Co-Ti 的 phonopy 1.11.0 输出一致，但**格式相似不等于身份相同**，A 线不据此填值。 |
| 状态 | **已有证据**（Co-Ti）/ **未解析**（Nb-Ti） |
| 与 §B.5.6 第 9 条的关系 | 该条要求"方法保持 `unknown`，除非有独立 provenance 证据"。Co-Ti **有**这种证据（INCAR / OUTCAR / log 三处独立），所以可以写 `DFPT`；Nb-Ti **没有**，所以写 `unknown`。这就是"不得由格式或目录名推定"的实际执行。 |

## 3. A 线明确没有做的事

- **没有**替用户冻结任何一项。本文件所有"做法"都是**这批数据的实际形态**，不是标准。
- **没有**按目录名或文件标题行推断物理语义。凡无证据的字段一律 `unresolved`。
- **没有**提交任何 VASP / Slurm 作业，也**没有**改动任何科学输入。
  BCC-Nb-Ti-2 的"复现"需要新跑计算，其前置条件列在
  `BCC-Nb-Ti-2/FINDINGS.md` §7，等用户决定。
- **没有**把两个端基的缺陷"修好"再交付。缺陷照原样交付并量化，因为
  "这份 legacy 数据到底能不能用"正是 A 线要向 B 线回答的问题。

## 4. 建议的下一步（供用户/B 线参考，非 A 线自行执行）

1. 用户就 §B.5.7 清单第 1、3、5、6 项给出决定——这四项 A 线无法自行推进。
2. B 线可先按 §B.5.6 的"曲线候选发现 → 证据配对 → 元数据决议 → fail-closed"
   框架把这两个端基接进来，`DATASET.json` 的 `claims[*].status` 可直接映射成
   逐字段的 `metadata_source` 与 `ApplicabilityStatus`。
3. Co-Ti 若要变成可用数据，最小改动是重跑声子层（至少排除 QHA-08）
   并给出 CHGCAR；这属于新的科学决策，需用户确认。
