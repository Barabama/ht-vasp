> ⚠️ **本文件是过程记录，部分内容已过时。**
> 当前结论以 **[docs/legacy-qha-rebuild-report.md](../../docs/legacy-qha-rebuild-report.md)** 为准。
> 已知过时点：`PARAM_FREEZE.md` 里的 `KPAR=2` **是错的**（应为 `NCORE=8`，见 D0-D）；
> `REVIEW_RESPONSE_AND_D0.md` §3.3 的成本数字全部虚高 7.6 倍；
> `E1_DFPT_vs_FD.md` 的成本对比需按 `NCORE=8` 重算。

# E1：DFPT vs 有限位移（BCC-Al-Ta-2）

**日期：** 2026-10-10 · **作业：** 55208（DFPT）、55209（有限位移）· **结果：** `/nfs_ssd/tmp/legacy-qha-pool/e1-comparison.json`

**结论：两种方法在物理上等价（差 ~2%），E1 的判据可以照用；但这次试验意外给出了一个更重要的数字——池内 79% 的端基在 legacy 里就有软模。**

---

## 1. 试验设计

几何**固定**在 legacy 弛豫后的胞（`BCC-Al-Ta-2/R3/CONTCAR` = `dfpt/QHA-04/POSCAR-unitcell`，
a = 3.2246015525754168 Å，B2 2 原子 → 2×2×2 超胞 = 16 原子），**不重新弛豫**——
这样方法差异里不会混进几何差异。其余量全部对齐：`ENCUT=520`、`MP 5×5×5`、
`ISMEAR=1/SIGMA=0.2`、`ISPIN=2`、`EDIFF=1e-5`、POTCAR `Al`(ZVAL 3)+`Ta`(ZVAL 5)。

两条路线的**信息量相同**：DFPT 对称性约化到 2 个微扰（原子 1 和 9，x 方向），
phonopy 也约化到 2 个位移。所以两者可以喂进同一套 phonopy 机器直接比。

## 2. 三个配方发现

### 2.1 legacy 的 DFPT 并行配方在本集群的 VASP 6.4.3 上直接跑不起来

第一次提交（job 55207）在第一个微扰就中止：

```
|     VASP internal routines have requested a change of the k-point set.
|     Unfortunately, this is only possible if NPAR=number of nodes.
|     Please remove the tag NPAR from the INCAR file and restart the
|     calculation.
|       ---->  I REFUSE TO CONTINUE WITH THIS SICK JOB ... BYE!!! <----
```

legacy 用 `NCORE=4`（VASP 5.4.4 允许），48 rank 下 NPAR = 12 ≠ 节点数 1，
而 `IBRION=8` 中途要改 k 点集。**BCC-Ti-V 试点踩过同一个坑**
（`archive/.../DFPT.failed-53674-*/INCAR` 里就是 `NCORE = 4`）。

**正式配方：`KPAR = 2`，不写 `NCORE`。** 已提交并成功。

### 2.2 `main_qha.py` 的全局 `ISYM=0` 会让 DFPT 贵 24 倍

`base.py:92` 默认 `ISYM: 0`，`main_qha.py` 没有覆盖。
对 DFPT 而言 `ISYM=0` 关闭对称性约化 → 48 个微扰而不是 2 个。
E1 的输入显式写 `ISYM = 2`。

### 2.3 `sym_reduce=False` 让有限位移的作业数从 2 变成 12

`htvasp/workflows/qha.py:350` 给 `PhononMaker` 传 `sym_reduce=False`，
而 atomate2 默认是 `True`。实测（同一结构、同一超胞）：

```
phonopy is_symmetry=True   (atomate2 sym_reduce=True ) ->   2 displacements
phonopy is_symmetry=False  (当前 ht-vasp)              ->  12 displacements
```

`generate_phonon_displacements` 里就是 `Phonon(..., is_symmetry=sym_reduce)`。
**6 倍**，不是 24 倍，但仍然可观。

| 路线 | 每体积点作业 | 84 端基 × 10 体积 |
|---|---|---|
| **DFPT**（VASP 内部解全部微扰） | 1 | **840** |
| 有限位移 `sym_reduce=True` | 2 + 1 静态 | ~2 520 |
| 有限位移 `sym_reduce=False`（当前） | 12 + 1 静态 | **~10 920** |

## 3. 三方对比（同一批 q 点）

把 legacy 交付的 `FORCE_CONSTANTS` 也拉进来，用同一套 phonopy 算：

| 量 | legacy DFPT | 本次 DFPT | 本次有限位移 |
|---|---|---|---|
| max \|Φ\| / eV·Å⁻² | 6.6466 | 6.7774 | 6.6759 |
| ASR 残差 / eV·Å⁻² | **0.0562** | **0.0000** | 0.0003 |
| Γ 声学 / THz | −0.0231 ×3 | **−0.0001** ×3 | −0.0154 ×3 |
| Γ 光学 / THz | 8.4872 ×3 | 8.5204 ×3 | 8.5188 ×3 |
| **N(½,½,0) 软模** / THz | **−3.644** | **−3.5588** | **−3.6327** |
| N 其它 / THz | 4.0651, 7.9604, 8.5565×2 | 4.016, 8.0452, 8.6081×2 | 4.0141, 8.0259, 8.5809×2 |
| H(½,½,½) / THz | 2.3282 ×3, 7.2361 ×3 | 2.4376 ×3, 7.389 ×3 | 2.4098 ×3, 7.1685 ×3 |
| ZPE / kJ·mol⁻¹ | 5.7462 | 5.8017 | 5.7624 |

**读法：**

- **方法差异 ≤ 2.3%**（N 软模：−3.644 / −3.5588 / −3.6327，跨度 0.085 THz）；
  Γ 光学差 0.39%；ZPE 差 0.97%。**DFPT 与有限位移在本体系等价。**
- **DFPT 的 ASR 是精确的**（0.0000），有限位移有 0.0003 的噪声，
  legacy 交付的 FC 反而有 0.0562 的 ASR 违反（它的 Γ 声学也因此是 −0.0231 而非 0）。
- **Γ 声学噪声量级 0.0001 – 0.0231 THz；真实软模 3.56 – 3.64 THz。相差 150 倍。**

## 4. 意外的主要发现：池内 79% 的端基有软模

`B2 AlTa` 在 **N 点 (½,½,0)** 有真实的软模，三条独立路线都复现（见上表）。
这不是方法伪影。

顺着这条线索查了 legacy 全部 84 个池内端基的 `thermal_properties`：

```
池内端基总数                      84
legacy 的 yaml 有模式被丢弃的       66   (79%)
  其中 10/10 个体积点都丢的        54
按相:  BCC 21 / FCC 41 / SER 4
没有任何丢弃的                    18
```

（"丢弃"= `num_modes > num_integrated_modes`，即低于 cutoff 的模式被 phonopy 静默排除；
这批数据 cutoff=0，所以就是虚频。）

**这改变了规划的性质**：软模不是边缘情况，而是池内的多数情况。
B2 铝化物（Al-Ti、Al-Ta、Al-Cr）和 Zr/Hf/Ti/Ta 基端基的 N 点或 Γ 附近软模是已知现象，
这里的数据与之相符。

**因此新库不能只做"通过/拒绝"二元门，必须带软模分级**，例如：

```
stable                         -> 进生产 QHA
gamma_acoustic_only            -> 进生产 QHA（|f| 阈值以上才算软模）
finite_q_soft_mode             -> 记录位置/幅度/本征矢；默认不进生产，出诊断产物
```

阈值这次有实测依据：Γ 声学噪声 ≤ 0.023 THz，真实软模 ≥ 3.5 THz。

## 5. 对上一版记录的更正

`E1_DFPT_vs_FD.md` 初版写 `BCC-Al-Ta-2`"是池内质量最好的 legacy 端基之一"。**这句话不完整。**

它在 **E(V)** 上确实好（`eos_r2 = 0.99999`，`G(0)` 偏差 −0.33 meV），
但它有 −3.6 THz 的 N 点软模，10/10 个体积点的 yaml 都在丢模式
（QHA-04 丢 16 344 / 82 944 个加权模式）。它是**E(V) 干净但有软模**这一类。

这对 E1 反而是好事：它让这次试验必须复现一个软模，比只比一个稳定结构严格。

## 6. 待冻结

- **Nb 没有标准 POTCAR**：`/nfs_ssd/documents/potpaw_PBE` 下只有 `Nb_pv`(ZVAL 11)
  与 `Nb_sv`(ZVAL 13)，没有 `Nb`。池内其余 6 个都有标准变体
  （Al 3 / Cr 6 / Hf 4 / Ta 5 / Ti 4 / Zr 4，与 legacy 实测一致）。
  这是池级决定，影响全部 28 个含 Nb 的 BCC 端基。
- 软模分级的具体阈值与"是否仍可用于 QHA"的判据，需要用户确认。
