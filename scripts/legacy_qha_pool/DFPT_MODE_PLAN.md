> ⚠️ **本文件是过程记录，部分内容已过时。**
> 当前结论以 **[docs/legacy-qha-rebuild-report.md](../../docs/legacy-qha-rebuild-report.md)** 为准。
> 已知过时点：`PARAM_FREEZE.md` 里的 `KPAR=2` **是错的**（应为 `NCORE=8`，见 D0-D）；
> `REVIEW_RESPONSE_AND_D0.md` §3.3 的成本数字全部虚高 7.6 倍；
> `E1_DFPT_vs_FD.md` 的成本对比需按 `NCORE=8` 重算。

# ht-vasp 的 DFPT 模式：可行性与设计计划（v2）

**日期：** 2026-10-10 · **状态：** 设计草案，未实施 · **依据：** E1（job 55208/55209）、D0 进行中（job 55210）
**v2 修订：** 采纳外部评审的 6 项修改（见 §7 修订记录）

---

## 0. 三条最重要的结论

1. **DFPT 与有限位移在本体系等价**（E1 实测差 ≤2%），两种方法可以互换。
2. **真正的主风险不是"DFPT 还是有限位移"，而是"哪些端基根本不适用简谐 QHA"。**
   池内 14 个可直接审计的端基里，**6 个有真实有限 q 软模**（−1.4 到 −4.4 THz），
   且位置都在 N 点附近。
3. **ht-vasp 目前没有 DFPT 能力**（atomate2 0.1.4 零支持），实现它要写新代码。

---

## 1. DFPT 比有限位移快吗？—— 必须实测，不能推算

### 1.1 E1 的实测（16 原子胞，k=5×5×5，48 核）

| 路线 | 作业数 | wall | 核时 | 相对 |
|---|---|---|---|---|
| **DFPT**（`IBRION=8`，2 个微扰） | 1 | 97.5 s | **1.30 core·h** | 1.0× |
| 有限位移 `sym_reduce=True`（2 个位移） | 2 | 2×28.5 s | **0.76 core·h** | **0.58×** |
| 有限位移 `sym_reduce=False`（当前 ht-vasp，12 个位移） | 12 | 12×28.5 s | **4.56 core·h** | 3.5× |

### 1.2 成本模型更正（v1 这里是错的）

> **v1 写的**："DFPT 的成本主要由原胞的电子结构决定，超胞只决定 q 采样密度，不进入成本公式。"
>
> **这句话是错的，已删除。** 实际情况：本工作流的 DFPT 是在**超胞上**跑的
> （`IBRION=8` 作用于 VASP 实际收到的那张 POSCAR，E1 里是 16 原子），
> 所以电子结构、能带数、微扰数、线性响应迭代**都随超胞增大而变贵**。
> 超胞进入成本公式，而且很可能是主导项。

**正确的说法：**

- 超胞增大时，**有限位移的位移数**可能保持不变（E1 实测 2³/3³/4³ 都是 2），
  但**每个位移的静态成本**随原子数上升；
- DFPT 的**微扰数**也可能保持较少，但**每个 DFPT 任务本身**同样随胞增大而变贵；
- 因此**谁更快只能由相同超胞、相同 k 密度、相同资源下的实测 wall/core-hour 决定**。

### 1.3 D0-A 实测（已完成 2×2×2 与 3×3×3）

**结论：DFPT 在 3×3×3 上仍然没有优势，而且两者的差距在缩小。**

| | 2×2×2（16 原子） | 3×3×3（54 原子） | 增长倍数 |
|---|---|---|---|
| **DFPT** wall | 97.5 s | **2285 s**（38.1 min） | **23.4×** |
| **FD**（2 位移） wall | 57.0 s | **1741 s**（2×871 s） | **30.5×** |
| 原子数 | 16 | 54 | 3.4× |
| DFPT 微扰数 | 2 | 2（不变） | — |
| FD 位移数 | 2 | 2（不变） | — |
| 比值 DFPT/FD | 1.71× | **1.31×** | — |

**两法都随超胞超线性增长**（3.4 倍原子 → 23–30 倍时间），
DFPT 的指数约 2.6，FD 约 2.8。

> **v1 里我预测"3×3×3 时有限位移可能反超 DFPT 10 倍"——这个预测是错的。**
> 实测只有 1.3 倍，而且 FD 的标度指数反而更大，
> 所以在更大超胞上 FD 的相对优势不会扩大。

**幂律外推（只有两点，粗糙，仅作规划用）**：4×4×4（128 原子）
DFPT ≈ 6.0 h、FD ≈ 5.5 h。**两者在 4×4×4 上仍然接近。**

### 1.4 一个真实的机制差别

DFPT 的基态 SCF **只做一次**，所有微扰在 VASP 内部一次解完；
有限位移的**每个位移都要从头做一次 SCF**。
这个优势被"单次线性响应更贵 + DFPT 任务本身随胞增长更快"抵消了。
**在本体系、本资源下，两者成本始终同量级，DFPT 从未占优。**

---

## 2. ht-vasp 现在能做 DFPT 吗？不能

1. **atomate2 0.1.4 完全没有 DFPT 支持。** 我查了 `PhononMaker` 的全部参数与源码：
   `DFPT` / `IBRION` / `linear response` / `LEPSILON` / `LPHON` **零命中**。
   它是纯有限位移（`IBRION=6, NFREE=2`），位移由 `generate_phonon_displacements()`
   生成，里面写死 `Phonon(..., is_symmetry=sym_reduce)`。**没有开关能切过去。**
2. `htvasp/workflows/qha.py` 走的就是 `PhononMaker`。
3. **但 pilot 有一个已验证的 DFPT-QHA 驱动**
   （`scripts/legacy_dfpt_qha_pilot_bcc_ti_v.sh` 706 行 + `.py` 1357 行），
   做的是完整 DFPT 流程，不使用 jobflow 的 flow 结构。E1 已验证它在集群上可用。

**结论：实现 DFPT = 写新代码。**

---

## 3. 池内软模：现在有直接证据了

### 3.1 v1 的表述过强，已收紧

> **v1 写的**："池内 84 个端基里 66 个（79%）有软模。"
>
> **改为**："池内 84 个 legacy 端基里，66 个的 `thermal_properties` 出现模式未被积分的
> 现象（`num_modes > num_integrated_modes`）；在 cutoff=0 下这提示存在非正频率模式，
> 但**位置、幅度、本征矢尚未逐一确认**。"

理由：`num_modes - num_integrated_modes` 是**热性质汇总层的代理指标**，
它不告诉你负频在 Γ 还是有限 q、是声学数值噪声还是真实失稳、幅度多少。

### 3.2 但对能直接审计的 14 个端基，代理指标 14/14 正确

池内 84 个端基里只有 **14 个有 VASP 层**（即有 `mesh.yaml` 可以做直接 q 点审计）。
我把这 14 个的 `dfpt/QHA-04/mesh.yaml`（平衡体积）逐个解析：

| 端基 | mesh 最低频 / THz | 含负频 q 点数 | 判定 |
|---|---|---|---|
| `BCC-Al-Al-2` | −3.3983 | 190 | **finite_q_soft_mode** |
| `BCC-Al-Ta-2` | −3.6234 | 284 | **finite_q_soft_mode** |
| `BCC-Al-Ti-2` | −3.6543 | 128 | **finite_q_soft_mode** |
| `BCC-Ta-Ti-2` | −1.9407 | 155 | **finite_q_soft_mode** |
| `BCC-Ti-Ti-2` | −4.4343 | 317 | **finite_q_soft_mode** |
| `SER-Ta-2` | −1.3859 | 141 | finite_q_soft_mode（**疑为模板伪影，见下**） |
| `BCC-Al-Cr-2` | +0.4616 | 0 | stable |
| `BCC-Cr-Cr-2` | +0.6770 | 0 | stable |
| `BCC-Cr-Ta-2` | +0.3794 | 0 | stable |
| `BCC-Cr-Ti-2` | +0.2913 | 0 | stable |
| `BCC-Ta-Ta-2` | +0.2674 | 0 | stable |
| `SER-Al-1` | +0.1646 | 0 | stable |
| `SER-Cr-1` | +0.3842 | 0 | stable |
| `SER-Ti-8` | +0.1234 | 0 | stable |

**代理指标与直接审计完全一致**：6 个 `10/10 体积点丢模式` 的全部有有限 q 软模；
8 个 `0/10` 的全部干净（最低频 +0.12 到 +0.68 THz，全部为正）。

**所以：代理指标在可验证的子集里 14/14 有效；对另外 70 个端基它仍然是代理，
因为那 70 个根本没有 VASP 层、没有 `mesh.yaml`。**

### 3.3 D0-A 的科学副产品：软模是一片区域，不是单个 q 点

用 3×3×3 的力常数（54×54）算，并与 2×2×2 对比：

| q 点 | 2×2×2 DFPT | 3×3×3 DFPT | 3×3×3 FD |
|---|---|---|---|
| (0,0,0) | −0.0002 | −0.0045 | +0.0153 |
| **(1/3,1/3,0)** | 不可表示 | **−3.4629** | **−3.5034** |
| **(1/2,1/2,0)** | **−3.5588** | **−3.5891** | **−3.5917** |
| (1/2,1/2,1/2) | +2.4376 | +2.4814 | +2.6259 |
| Γ 光学 | 8.5204 | 8.5571 | 8.5748 |

两条结论：

1. **软模沿 [110] 延展成一片**：`(1/3,1/3,0)` 有 −3.46 THz，`(1/2,1/2,0)` 有 −3.59 THz。
   2×2×2 **看不见** `(1/3,1/3,0)`（不可表示），所以 **2×2×2 是欠采样的**——
   换更大超胞不是"消除软模"，而是**发现更多**。
2. **两法在 3×3×3 上仍然等价**：软模差 1.2%（`(1/3,1/3,0)`）和 0.07%（`(1/2,1/2,0)`），
   Γ 光学差 0.2%。最大差异 0.42 THz 出现在 `(1/3,1/3,1/3)` 的 8.78 THz 模上（4.8%）。

### 3.4 两条新发现

**(a) 软模位置是同一个 q 点。** 5 个 BCC 端基的最低频全部在
`q = (0.4792, 0.4792, 0.0208)`，距 N 点 (½,½,0) 只差 0.036（24³ 网格的步长是 0.0417）。
**这是 BCC 铝化物/β 相在 N 点的系统性失稳，不是随机噪声。**

**(b) `SER-Ta-2` 的软模疑为模板伪影。** 它的胞是

```
a = 3.2821 Å,  c = 4.0169 Å,  c/a = 1.224
```

**理想 HCP 的 c/a = 1.633，它偏低 25.1%；而 Ta 的基态是 BCC，不是 HCP。**
也就是说 legacy 拿一个严重压缩的 HCP 模板去算纯 Ta，失稳是必然的。
这与审计报告里"SER 不是标准参考态"（`SER-Mo-2` 用 CsCl、`SER-Ta-2` 用 HCP）是同一件事。

### 3.5 这对规划的含义

不是"新库有 79% 不能做 QHA"，而是：

> **legacy 的 `thermal_properties` 不能用作新库的稳定性先验，
> E(V) 拟合质量也不能替代声子适用性判断。**

正确路线是三层筛查：

```
第一层：legacy 只作为风险先验（代理指标，14/14 有效但只覆盖 14 个）
第二层：新管线在参考体积直接做 q 点声子适用性审计
第三层：只有 stable / 批准的 gamma_acoustic_noise 才进入 10 体积 QHA
```

---

## 4. 软模分级：阈值是 provisional 的

```
stable                     -> 进生产 QHA
gamma_acoustic_noise       -> 进生产 QHA（仅在 Γ 附近、幅度小于阈值）
finite_q_soft_mode         -> diagnostic_only，默认不进生产
phonon_incomplete          -> 不进下游
```

**阈值（provisional，非普适）：**

```json
{
  "gamma_acoustic_tolerance_THz": 0.1,
  "tolerance_status": "provisional",
  "basis": "E1 + D0, B2 AlTa, VASP 6.4.3, ENCUT 520, MP 5x5x5, ISMEAR 1/0.2, PBE 标准 POTCAR",
  "not_globally_frozen": true
}
```

依据：E1 实测 Γ 声学噪声 0.0001–0.0231 THz，真实有限 q 软模 1.4–4.4 THz，**相差两个数量级**。
但阈值会随元素质量、力常数尺度、k 点、ENCUT、EDIFF、超胞、ASR 修正、磁性状态、POTCAR 变化，
所以**只作为初始审计阈值，不冻结为普适门禁**。每个批次必须记录该阈值与其依据。

**分级必须在写 `thermal_properties` 之前做**，不能只靠 `QhaMaker(ignore_imaginary_modes=False)`
——那个开关在 QHA 后处理阶段，管不到上游 phonopy 已经丢弃负频模式。

---

## 5. 设计计划

### 5.1 架构：并行 Worker，不改造现有 Worker

**v2 修订**：D1 不直接改 `QhaWorker`，先新增独立的 `DFPTQhaWorker`。

```
QhaWorker                       ← 现有有限位移路径，一行不动
DFPTQhaWorker                   ← 新增
  ├─ prepare_structure
  ├─ generate_sposcar
  ├─ run_dfpt
  ├─ check_dfpt_completion      ← pilot 的 cmd_verify --stage dfpt
  ├─ build_force_constants      ← phonopy --fc vasprun.xml
  ├─ audit_phonon_stability     ← 软模分级，在 thermal 之前
  ├─ generate_thermal_properties
  └─ analyze_qha
```

理由：atomate2 的 `PhononMaker` 输出、`QhaMaker` 输入、store 生命周期都是围绕有限位移设计的。
一开始就塞进同一个 Worker 会导致 store 语义混淆、resume 无法区分方法、
失败 attempt 被覆盖、以及"上游已产出 thermal_properties 才发现软模"。

等 DFPT 端到端通过 E1 三方复核后，再加开关：

```python
QhaWorker(..., phonon_method="dfpt" | "finite_displacement")   # 默认保持现状
```

### 5.2 必须从 pilot 移植、不要重写的三段

1. **SPOSCAR ↔ VASP 胞 ↔ FORCE_CONSTANTS 的映射核对**
   （`cmd_check_phonopy_inputs`，`.py:1083`）：验证 `POSCAR-unitcell` + `DIM`
   重建的超胞 == VASP 实际跑的 SPOSCAR == FORCE_CONSTANTS 的 N×N。
   **不做这一步，力常数会静默映射到错的胞上。**
   已知坑：DFPT 目录里的 `POSCAR` 是 SPOSCAR，对它再套 `DIM=2 2 2` 会建出 128 原子胞。
2. **DFPT 完成性检查**（`cmd_verify --stage dfpt`，`.py:680`）：作业完成标记、
   `aborting loop because EDIFF is reached`、无 `unconverged` 标记、
   `dfpt_all_perturbations_solved`（最后一个 `Degree of freedom: n/M` 必须是 M/M）、
   `Linear response finished`（在 **stdout**，不在 OUTCAR）。
3. **不可变 attempt 归档**（`archive_dfpt_attempt`）：失败/被取代的尝试 `mv` 走，
   不删除、不覆盖。来自 53701 事故。

### 5.3 门禁：区分"诊断指标"与"失败条件"

**v2 修订**：v1 把几个门写成了普适 fail 条件，改为分级。

| 检查 | 类型 | 判据 |
|---|---|---|
| 基态 SCF 撞 NELM | **fail** | 迭代数 ≥ NELM；同时报末段能量跨度以区分"发散"与"只是没到 EDIFF" |
| 基态能量 vs 同体积静态 | **强异常**（recipe-specific） | > 100 meV/atom；`gate_basis = "Ti-V pilot / E1 recipe"` |
| `max\|Φ\|` | **诊断指标** | 不作为单独 fail 条件；记录原值 |
| ASR 残差 | **fail**（同时记绝对值与相对值） | 相对 > 10% **且** 绝对值 > 0.01 eV/Å² |
| permutation 残差 | 诊断指标 | 记录 |
| 软模分级 | **fail（分类型）** | 见 §4 |

### 5.4 配方

| 项 | 值 | 条件 |
|---|---|---|
| `IBRION` / `NSW` | 8 / 1 | |
| `KPAR` | **2，不写 `NCORE`** | E1 实测：`NCORE=4` 在 VASP 6.4.3 上 abort |
| `ISYM` | **2**，但**只对经过结构与磁性对称性审计的端基** | 见 §5.5 |
| `ALGO` | `Normal` | 项目已验证的 16 原子胞配方 |
| `NELM` / `NELMIN` | 120 / 10 | pilot |
| `ENCUT` / `ISMEAR` / `SIGMA` | 520 / 1 / 0.2 | legacy 与 `main_qha` 一致（**收敛未验证**） |
| `EDIFF` | 1e-5 | legacy DFPT |
| `ICHARG` / `LCHARG` | 2 / False | 当前 ht-vasp 策略（明写） |
| 超胞 / q 网格 | **待 D0** | |

### 5.5 `ISYM` 与 `sym_reduce` 是条件化策略，不是全局开关

| 结构类型 | `ISYM` | `sym_reduce` |
|---|---|---|
| 精确有序 BCC/B2、FCC/L1₂，且磁矩不破坏对称性 | 2 | True（需先做 symmetry audit） |
| 磁性排列可能破坏空间对称性 | 需审计 | 按联合对称性定 |
| 低对称 / 内部坐标自由 | 0 或按审计 | False |
| SQS / 随机结构 | 0（通常 P1） | False，除非验证了非平凡对称性 |

E1 证明 `ISYM=2 + KPAR=2` 在这个 B2 测试里可行，这只支持把它作为
`ordered_high_symmetry_dfpt_recipe`，**不能推广为全局默认**。

---

## 6. D0：修订版（job 55210 进行中）

### D0-A 成本与超胞收敛

- 体系：`BCC-Al-Ta-2`（软模验证）+ 一个稳定控制体（`BCC-Cr-Ta-2`，池内且审计为 stable）
- 超胞：2×2×2 / 3×3×3 / 4×4×4
- 方法：DFPT / FD(`sym_reduce=True`)
- 记录：wall、core-hour、内存峰值、电子迭代数、**微扰/位移数**、
  force constants、ASR、permutation、共同 q 点频率、ZPE、热自由能

**k 点策略（v2 新增）**：各超胞用 `mesh = round(10 / n)`，
使 `n × mesh ≈ 10`，即**原胞布里渊区的绝对采样密度固定**。
2×2×2 → 5×5×5（legacy 值）；3×3×3 → 3×3×3；4×4×4 → 3×3×3。

**q 点比较策略（v2 新增，评审指出的风险）**：

不同超胞的实空间截断不同，可直接解析的 commensurate q 点也不同：

| 超胞 | q=(½,½,0) 是否 commensurate |
|---|---|
| 2×2×2 | **是**（精确） |
| 3×3×3 | **否**（需插值） |
| 4×4×4 | **是**（精确） |

所以比较分两类：

1. **共同 q 点**：各超胞都可可靠表示的 q（Γ、共同高对称点、共同 commensurate 点）
2. **插值热力学**：每套力常数用**相同的倒空间采样密度**比较 ZPE / F / Cv / S / 频率分布

**并且必须单独盯住 q=(½,½,0)**，因为软模就在这里——
不能因为某个超胞不直接包含该 q 就误判为"没有软模"。

### D0-B 对称性策略

`ISYM=2 + sym_reduce=True` vs `ISYM=0 + sym_reduce=False`：
记录实际位移数、VASP 内部微扰数、MAGMOM、结构对称性、力常数、频率差异。

### D0-C Nb POTCAR 预审计（**阻塞项**）

`/nfs_ssd/documents/potpaw_PBE` 下**没有标准 `Nb`**，只有 `Nb_pv`(ZVAL 11) / `Nb_sv`(ZVAL 13)。
池内其余 6 个都有标准变体（Al 3 / Cr 6 / Hf 4 / Ta 5 / Ti 4 / Zr 4，与 legacy 实测一致）。

**在选定之前不要启动含 Nb 的批量计算**（影响 28 个 BCC 端基）。
审计内容：纯 Nb + Nb-Ta / Nb-Ti / Nb-Zr + 一个目标多元体系，
比较晶格常数、总能、k 点/ENCUT 收敛、声子稳定性、与现有参考能的一致性。
选定后全库不得混用，provenance 记录 `{element, potcar_variant, zval, potcar_sha256}`。

---

## 7. 修订记录（v1 → v2）

| # | v1 | v2 | 来源 |
|---|---|---|---|
| 1 | "DFPT 成本由原胞决定，超胞不进成本公式" | **删除**。改为"两者都随超胞变贵，谁快必须实测" | 外部评审，成立 |
| 2 | "池内 66/84 有软模" | 改为"66/84 的 legacy `thermal_properties` 出现模式丢失"；并补上 14 个可审计端基的**直接 q 点审计**（代理指标 14/14 正确） | 外部评审，成立；并补了新证据 |
| 3 | `0.1 THz` 作为门禁阈值 | 标记 `provisional`、recipe-specific、不全局冻结 | 外部评审，成立 |
| 4 | D1 直接改造 `QhaWorker` | 改为新增并行 `DFPTQhaWorker`，验证后再加开关 | 外部评审，成立 |
| 5 | `max\|Φ\| > 50` 作为普适 fail | 降为诊断指标；ASR 改为相对+绝对双条件 | 外部评审，成立 |
| 6 | `ISYM=2` / `sym_reduce=True` 作为全局默认 | 改为按结构类型与磁性条件化的策略表 | 外部评审，成立 |
| 7 | D0 只比"作业数" | D0 改为比较 wall/core-hour/内存/共同 q 点/插值热力学，并加 k 点策略与 commensurate q 讨论 | 外部评审，成立 |
| 8 | — | 新增 §3.3(b)：`SER-Ta-2` 的软模疑为 HCP 模板伪影 | 本轮回溯发现 |

**v1 中经复核仍然成立的**：方法等价性、`KPAR=2` 配方、`ISYM=0`/`sym_reduce=False` 的成本代价、
BCC-Al-Ta 的 N 点软模不是方法伪影、legacy E(V) 质量 ≠ 声子质量、
软模分级必须在 thermal_properties 之前。

---

## 8. 待拍板

1. **超胞与 q 网格**（决定方法，也决定成本）—— D0-A 已给出 2³/3³，
   4³ 待定（外推 5.5–6.0 h 一个作业）
2. **`Nb` POTCAR 变体** —— 阻塞 28 个 BCC 端基
3. **软模分级阈值**，以及 `finite_q_soft_mode` 是"只出诊断"还是"可进生产但标注"
4. **是否保留有限位移路径**作为交叉验证（评审建议保留，我同意）

## 9. 下一步顺序

```
修订设计（本文档 v2）           ← 已完成
    ↓
D0-A/B/C（超胞、成本、对称性、POTCAR）  ← 55210 进行中
    ↓
DFPTQhaWorker 最小实现
    ↓
BCC-Al-Ta E1 三方重演
    ↓
稳定控制端基验证（BCC-Cr-Ta-2）
    ↓
软模分级 + fail-closed 验收
    ↓
小批量 triage（14 个可审计端基先做）
    ↓
再决定是否扩展到 84 个
```
