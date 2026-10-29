# TB2J Wannier 参数问题：完整发现链路 + 待确认的解决方案

> **用法**：把下面「提示词正文」整段复制给一个新的 AI 会话。
> 目标仓库：`/nfs_hdd/2025/gaominliang/ht-vasp`。
> 本文档最后一节「事实附录」是给**你**（派发人）核对用的，不要复制进去。

> **关于本文档的路径**：这些提示词文件在 `master`（以及 `feat-dfpt-qha`）的 `docs/prompts/` 下。
> 你在 `master` 上开的新分支会继承它们，不会消失。
> **你手上粘贴的正文就是权威版本**；如果正文里引用了某个 `docs/prompts/*.md` 而你找不到，
> 停下来向派发人索取，不要猜内容。

---

## 提示词正文

你的任务是**查清并锁定** TB2J 的 Wannier90 参数该怎么改。

**注意你的角色：你不是来执行一个既定方案的。方案还没有定——连"哪个是主因"都还没有定（见 §2.11 的反证）。** 下面第 2 节给你完整的发现链路和全部证据，第 3 节给你 6 个候选方案（A–F）和它们的证据强度。你要做的是：**独立复核证据 → 设计判别实验 → 反复与派发人确认 → 拿到明确批准后才动手改代码**。

**在派发人书面确认之前，不要改 `htvasp/tb2j/input_set.py` 里的任何参数，不要重跑任何 TB2J 作业。**

### 0. 先读

按顺序：

1. `/nfs_hdd/2025/gaominliang/AGENTS.md`
2. `/nfs_hdd/2025/gaominliang/ht-vasp/CLAUDE.md`
3. `/nfs_hdd/2025/gaominliang/ht-vasp/docs/TB2J与SQS数据审计-20261025.md` —— §1（TB2J 覆盖与数值问题）、§1.3（结构层面的问题）
4. `/nfs_hdd/2025/gaominliang/tb2j_share/TB2J磁交换计算与磁参数构建指南.md` —— 注意 §3 的脚本模板，那个模板里的 `.win` 就是问题的源头
5. `/nfs_hdd/2025/gaominliang/sof-ace/docs/TB2J流程与工作量估算.md` —— §5、§6、§7

然后 `git status`，**不要 reset、clean、stash、覆盖或删除任何已有文件**。

### 1. 分支

**2026-10-28 分支事实（与旧版本文档不同，以本节为准）：`master` 已经包含 PR #3 的全部内容，TB2J 代码就在 master 上。**

```bash
git branch -vv
#   feat-dfpt-qha c6671c4 [origin/feat-dfpt-qha]   ← master 的祖先，已合并
# * feat-tb2j-sqs 6002f2a
#   master       6002f2a [origin/master]
```

- **不要 `git checkout feat-tb2j`。** 本地分支已删；远端 `origin/feat-tb2j` 落后 master **99** 个提交、领先 **0** 个，是 master 的真子集。checkout 会给你建 tracking 分支，`merge master` 变成 99 提交的大合并。
- **不会再"切分支就丢 `htvasp/dfpt_qha/`"**——`master` 现在有 `htvasp/dfpt_qha/`、`main_dfpt_qha.py`、`main_spin_screen.py`、`htvasp/tb2j/context.py`、`htvasp/tb2j/magnetic_contract.py`。调查要用的 DFPT-QHA 代码与 `data/endpoints/` 数据都在同一个工作区里，**不需要跨分支搬任何东西**。

**你要做的：在当前 HEAD（`feat-tb2j-sqs` @ `6002f2a`，= master）上开一个只读调查分支。**

```bash
git checkout -b feat-tb2j-wannier-params
./.conda/bin/python -m pytest tests/ -q   # 记基线（test_tb2j.py 的"12 passed"是假绿，见 #1 提示词）
```

### 2. 发现链路（按时间顺序，每一步都有可复核的文件）

#### 2.1 起点：一个看起来正常的数字

审计 44 个 TB2J 归档时发现，池内只有 2 个端基有 TC：

```
BCC-Nb-Nb   BMAGN = 0.0073 μB   TC(MFA) = 1.44 K
BCC-Ti-Ti   BMAGN = 0.0197 μB   TC(MFA) = 33.83 K
```

**0.007 μB 的磁矩配 1.44 K 的居里温度，物理上讲不通。** 这是第一个疑点。

复核路径：`data/endmembers/BCC-Nb-Nb/BCC-Nb-Nb-tb2j.json` 的 `magnetic_observables`。

#### 2.2 疑点一：J 的量级

同一个流水线的对照：

| 端基                                  | max\|J_iso\| (meV) | 折合 k_B·T(300K) |
| ------------------------------------- | -----------------: | ----------------: |
| BCC-Nb-Nb                             |             0.0347 |            0.0013 |
| BCC-Ti-Ti                             |             0.7332 |             0.028 |
| **BCC-Co-Fe**（同一流程，池外） |    **23.33** |    **0.90** |

差 2–3 个数量级。**但这本身还不能定罪**——Nb/Ti 的磁矩确实很小。

#### 2.3 疑点二：`dH/|m|` 不守恒

如果交换劈裂 ∝ 磁矩，那么 `max|H_up − H_dn| / |m|` 应该在各体系间大致相同。实测（从 `wannier90.{1,2}_hr.dat` 直接算）：

| 端基      | max\|m\| (μB) | max\|dH\| (meV) | dH/\|m\| (meV/μB) |
| --------- | -------------: | --------------: | -----------------: |
| BCC-Nb-Nb |         0.0073 |           2.334 |      **320** |
| BCC-Ti-Ti |         0.0198 |          14.314 |      **723** |
| BCC-Co-Fe |         2.8161 |           8.213 |      **2.9** |

**差 250 倍。** 这是第二个疑点，但同样不是定论——不同的态密度形状可以让这个比值不同。

#### 2.4 疑点三：TB2J 自己的警告

`exchange.out` 头部就写着：

```
Warning: Please check if the noise level of Wannier function Hamiltonian to
make sure it is much smaller than the exchange values.
```

把 `max|H_up − H_dn|` 当噪声、`max|J_iso|` 当信号：

| 端基      | 噪声 (meV) | max\|J\| (meV) |          J/噪声 |
| --------- | ---------: | -------------: | --------------: |
| BCC-Nb-Nb |      2.334 |         0.0347 | **0.015** |
| BCC-Ti-Ti |     14.314 |         0.7332 | **0.051** |
| BCC-Co-Fe |      8.213 |          23.33 |            2.84 |

**池内 2 个有 TC 的端基，J 分别只有噪声的 1/67 和 1/20。** 但注意：噪声和信号都来自同一个 `H_up − H_dn`，所以这个比值度量的是"J 有多大比例是 Wannier 拟合误差"，不是"J 是不是物理的"。

#### 2.5 决定性的一步：`.wout` 里一直写着没收敛

上面三条都还是间接证据。直接去看 Wannier90 自己的日志：

```bash
grep -l "Disentanglement convergence criteria not satisfied" \
    data/endmembers/*/tb2jflow/2-w90_static/wannier90.1.wout | wc -l
# 44
grep -l "Disentanglement convergence criteria not satisfied" \
    data/endmembers/*/tb2jflow/2-w90_static/wannier90.2.wout | wc -l
# 44
```

**44/44 个 up 通道 + 44/44 个 down 通道，全部撞了 200 次迭代上限。**

（`wannier90.wout`——不带编号的那个——是 VASP 的 setup pass，**没有 DIS 迭代表，不要拿它判**。它的大小不是常量：44 个 run 里实测 15421 / 15581 / 15741 / 16061 字节四种。按 `<-- DIS` 行识别，不要按字节数。）

#### 2.6 为什么是 200 次上限：窗口放错了位置

`htvasp/tb2j/input_set.py:265-268` 生成：

```
dis_win_min = -50      dis_win_max = 50
dis_froz_min = -10     dis_froz_max = 10
```

wannier90 官方文档（3.1.0 `doc/user_guide/parameters.tex`）：这些是**绝对能量，单位 eV**。

实测 E_F 与窗口的关系（E_F 从各端基的 `2-w90_static/OUTCAR.gz` 读）：

| 端基                |        E_F (eV) | `dis_froz` 相对 E_F                      |
| ------------------- | --------------: | ------------------------------------------ |
| BCC-Ni-Ni           |            4.53 | −14.5 … +5.5（跨过 E_F）                 |
| SER-Co              |            5.40 | −15.4 … +4.6                             |
| **BCC-Nb-Nb** | **11.48** | **−21.5 … −1.5（全在 E_F 之下）** |
| SER-Nb              |           11.25 | −21.3 … −1.3                            |
| BCC-Al-Nb           |           10.20 | −20.2 … −0.2                            |

对 `BCC-Nb-Nb`（band 跨 −42.9 … 71.2 eV）：

- `dis_win = −50…50` 只覆盖 **66.8%** 的流形（33.2% 在 50 eV 以上）
- `dis_froz = −10…10` 只覆盖 **4.0%**，而且全在 E_F 下方 1.5–21.5 eV —— **冻住的是半芯态，不是 E_F 附近的 d 流形**

**对比：如果窗口按 E_F 写（`dis_froz = E_F−10 … E_F`，`dis_win = E_F−15 … E_F+15`）**，`dis_froz` 覆盖 7.7%、`dis_win` 覆盖 23.6%——是同一个数量级的流形，但位置对。

> **但这条只对 E_F > 10 eV 的 3 个端基成立。** 其余 41 个的窗口本来就跨过 E_F，它们同样不收敛。见 §2.11——那里有一个单变量对照，把嫌疑指向 k 点数而不是窗口。

#### 2.7 窗口是哪来的：Fe 原型被当成通用参数

```bash
git log -1 --format="%h %ad %s" --date=short 063c12d
# 063c12d 2026-08-16 feat(tb2j): TB2J magnetic exchange workflow (R3→W90→TB2J)
```

`tb2j_test/fe_bcc_wannier/wannier90.win` 的时间戳是 **Aug 13 15:34**，内容和 `input_set.py` 生成的逐字相同：

```
num_wann = 9   num_bands = 16
dis_win_min = -50   dis_win_max = 50
dis_froz_min = -10  dis_froz_max = 10
```

Fe bcc 的 **E_F ≈ 1.17 eV**，所以那个窗口对它恰好跨过 E_F（`dis_froz` 覆盖 46.6% 的 band）。**这个"对"是巧合，不是设计。**

#### 2.8 为什么早期原型没暴露：三个条件同时变了

`tb2j_test/` 里的实测：

| 目录                     | num_wann/num_bands |            DIS 迭代 | 结果                                      |
| ------------------------ | ------------------ | ------------------: | ----------------------------------------- |
| `fe_bcc_wannier`       | 9/16               |                   3 | 收敛（空态只有 7 个，解纠缠几乎无自由度） |
| `fe_bcc_wannier_4x4x4` | 9/16               |                  71 | **收敛**                            |
| `fe_bcc_wannier_8x8x8` | 9/16               | 200 / **181** | **up 撞上限；down 第 181 步收敛**   |
| `SER-Co/W90`           | 18/64              |                 200 | **撞上限**                          |
| `BCC-Co-Ni/W90`        | 18/64              |                 200 | **撞上限**                          |
| `FCC-Fe-Mn/W90`        | 36/64              |                 200 | **撞上限**                          |

**同一个 Fe bcc，k 网格从 4×4×4 换到 8×8×8，up 通道就从收敛翻成不收敛（down 通道仍然收敛）。** 而 `num_wann/num_bands` 从 9/16 变成 18/64（空态从 7 个变成 46 个）之后，三个体系全部不收敛。

（`fe_bcc_wannier` 的"3 步收敛"要打折：那个 run 的空态只有 7 个，解纠缠几乎没有自由度，所以它不是"窗口对了"的证据。）

#### 2.9 为什么流程没拦住：门查错了对象

`tb2j_test/submit_tb2j.sh:228-234` 的"W90 收敛检查"：

```python
from pymatgen.io.vasp.outputs import Vasprun
v = Vasprun('vasprun.xml', parse_dos=False, parse_potcar_file=False)
ok = v.converged          # ← 这是 VASP 的 SCF 收敛
exit(0 if ok else 1)
```

Wannier90 是以**库模式**跑在 VASP 里的，所以 VASP 的 SCF 收敛了，这道门就过了。真正的判据在 `wannier90.{1,2}.wout` 里，而 `htvasp/tb2j/` 全包 grep `wout` / `disentangle` / `Omega` 是 **0 命中**。

而且 VASP 自己也不报。**注意这条的出处**：`W90.log`（113 行、第 91 行 `Running wannier90 in library mode for spin up   (check wannier90.wout)`，注意括号前是三个空格）只存在于 `tb2j_test/SER-Co/W90/`，**44 个归档 run 里一个 `W90.log` 都没有**（`find data -name 'W90.log*'` → 0）。所以"日志明确告诉你去查 wout"这句话对归档 run 不成立——归档 run 连这个提示都没有，只有 `.wout` 本身。

#### 2.10 文档知道要查，但没人实现

`sof-ace/docs/TB2J流程与工作量估算.md` 第 191 行和 218 行两处都写着：

> 重点检查 **Wannier90 拟合质量**：.wout spread（3d < 1 Å）、`assigned_basis.txt` 轨道归属、能带复现、TB2J 输出磁矩 vs DFT 磁矩一致性

**文档写了要查，代码里没有实现，所以 44 个 run 全部带着未收敛的 Wannier 基底进了数据库。**

#### 2.11 一个反证：窗口不是唯一原因（**派发人最初的说法被数据推翻**）

我把 44 个端基按 E_F 与 `dis_froz_max = 10 eV` 的关系分组，再按"冻结窗内有多少个 band"统计：

| 分组                            |         数量 | 冻结 band 数 / num_wann（中位） |
| ------------------------------- | -----------: | ------------------------------: |
| E_F < 10（窗口跨过 E_F）        | **41** |                  **0.65** |
| E_F > 10（窗口整个在 E_F 之下） |  **3** |                            0.15 |

（分组数 41/3 已独立复核：从 44 个 `2-w90_static/OUTCAR.gz` 的 `E-fermi` 取最后一个值，> 10 eV 的**恰好**是 `BCC-Al-Nb 10.20`、`BCC-Nb-Nb 11.48`、`SER-Nb 11.25`，全部 E_F 范围 4.53–11.48 eV。**但"冻结 band 数/k 中位 0.65 / 0.15"这两个数没复现**，见下。）

**41 个体系的窗口位置是"对"的，它们照样全部不收敛。** 所以窗口错位不是唯一原因，甚至可能不是主因。

我另外做了两个相关性检验，都不支持"窗口位置决定难度"：

- `correlation(冻结 band 数/k, log|ΔΩ@200|)` = **+0.056**（几乎无关）
- `correlation(E_F − dis_froz_max, log|ΔΩ@200|)` = **+0.312**（弱正相关，方向甚至与"窗口越偏越难"相反）

（这两个相关系数是派发人自己算的，**2026-10-26 的独立复核没有复现**——复核 agent 拿不到同样的"冻结 band 数"定义。你要么自己重新定义并重算，要么**不要引用这两个数**。结论本身（窗口位置不是主因）由 §2.11 的 k 点对照支撑，不依赖这两个相关系数。）

**真正强的证据是另一个变量：W90 SCF 的 k 点数。**

`tb2j_test` 里有一组**接近单变量的对照**——同一个 Fe bcc（`POSCAR` 的 a 都是 2.87）、同一个窗口、同样 `num_wann/num_bands = 9/16`、`ENCUT`/`NBANDS`/`ISPIN`/`EDIFF`/`ISMEAR`/`SIGMA` 全部相同，**KPOINTS 从 4×4×4 改成 8×8×8**。两个目录的 `POSCAR`、`INCAR`、`wannier90.win`、`wannier90.1.win`、`wannier90.2.win` 都逐字节相同（2026-10-26 复核）：

> **但严格说不是"只改了一个变量"**：`submit.sh` 里 `--ntasks` 也不同（4×4×4 用 4，8×8×8 用 16）。而且 **8×8×8 的 down 通道其实是收敛的**（spin 2 在第 181 步满足判据，只有 spin 1 撞了 200 上限）。所以"k 点多→难收敛"这条证据**比文档原先写的更强也更窄**：它只在 up 通道成立，且并行度也变了。做判别实验时必须把 `ntasks` 一起钉死。

| 目录                                                                                                                                                                                                              | k 点 | 收敛                                    |
| ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ---: | --------------------------------------- |
| `fe_bcc_wannier_4x4x4`                                                                                                                                                                                          |   64 | **是**                            |
| `fe_bcc_wannier_8x8x8`                                                                                                                                                                                          |  512 | **up 否（200）** / down 是（181） |
| `fe_bcc_wannier`（默认 8×8×8，但**这个 run 没跑完**——`wannier90.1.wout` 只有 166 行停在 DIS 表头，**连 `wannier90.2.wout` 都没有**，无 `*.hr.dat`、CHG/CHGCAR 为 0 字节、没有 TB2J 交接） |   — | 不可判                                  |

而 44 个归档 run **全部是 512 个 k 点**。

同时 `num_bands − num_wann`（解纠缠的空态数）分布是：

| num_wann / num_bands | 空态 | 端基数 |
| -------------------- | ---: | -----: |
| 9 / 32               |   23 |      5 |
| 18 / 64              |   46 |     23 |
| 36 / 64              |   28 |     16 |

**所以你至少有三个候选原因，而派发人最初只强调了窗口。** 你的复核任务里，这一条优先级最高：**先判断 k 点数、空态数、窗口位置三者哪个是主因**，再谈改什么。这也直接决定了方案 A/D/E 值不值得做。

（`num_bands = ceil((num_wann+20)/ntasks)*ntasks`，见 `input_set.py:144-146`，`num_bands_extra` 默认 20、`ntasks` 默认 32。
**修正一处**：`num_wann` 本身**不是** `ntasks` 的副作用——它是 `len(structure) × 9`（`input_set.py:102-122`，`ORBITALS_PER_ATOM=9`），由胞的原子数决定；只有 `num_bands` 依赖 `ntasks`。
而且单一默认 `ntasks=32` 就能同时给出 9→32、18→64、36→64 这三组，**不需要三个不同的 `ntasks`**。所以"空态数"= `num_bands − num_wann` 是 `num_wann`（胞大小）和 `ntasks` 的联合结果。）

### 3. 六个候选方案，以及各自的证据强度

**先说清楚：没有任何一个方案已经被证明是对的。** 下面按"证据强度"排序，你要做的是把它们变成可判别的实验。

#### 方案 A：窗口改成相对 E_F

```
dis_win_min = E_F - 15     dis_win_max = E_F + 15
dis_froz_min = E_F - 10    dis_froz_max = E_F
```

- **支持**：绝对窗口在 E_F 从 1.2 变到 11.5 eV 的 7 元素池里必然错位，这是数学事实（§2.6）。
- **反对/未知**：wannier90 的窗口是**绝对**的，VASP 的 `WANNIER90_WIN` 也接受绝对数值。要写成相对 E_F，必须在**运行时**知道 E_F——而 E_F 是 SCF 之后才有的，`wannier90.win` 是 SCF 之前写的。**这是本方案最大的实现障碍。**
  - 可能的绕法：(i) 先用较低精度 SCF 估 E_F，再写 win（两步 SCF，成本翻倍）；(ii) 用 VASP 的 `WANNIER90_WIN` 标签在 SCF 后覆盖（VASP 6.x 支持，但需要核实 6.4.3 的具体行为）；(iii) 按元素/体系预先估算 E_F 的偏移量，写成绝对数值（脆弱，且每个体系不同）。
  - **这条必须先查清 VASP 6.4.3 的 `WANNIER90_WIN` 行为，再谈可行性。**

#### 方案 B：只提高 `dis_num_iter`

```
dis_num_iter = 1000   （或 2000）
```

- **支持**：从 190→200 步的局部衰减率外推，**35/44 会在 307–12263 次迭代内收敛**（中位 1001）。这是 **spin-1 单通道**的外推，且『收敛』是乐观上界——down 通道差很多（SER-Mn 682→1923、BCC-Mn-Mn 894→7159）。
- **反对**：**9/44 看起来已平台化**（衰减率 ≥ 1，`BCC-Al-Al`、`BCC-Co-Mn`、`BCC-Fe-Mn`、`BCC-Mn-Ni`、`BCC-Nb-Ti`、`FCC-Mn-Co`、`FCC-Mn-Fe`、`SER-Co`、`SER-V`）。**这个判定的强度要说清**：它们的 q 只比 1 高 0.1–3%，换窗口（180→190 / 195→200）后只有 **7 个**在两个自旋上都稳，`BCC-Mn-Ni` 与 `BCC-Nb-Ti` 会掉出去。所以『提高迭代也救不了』对 7 个成立、对 2 个是刀锋判定。
- **反对**：`BCC-Al-Nb` 需要 12263 次——成本不可接受。
- **反对**：迭代多不等于结果对。窗口错位时，它收敛到的是"这个错误窗口下的最优子空间"。

#### 方案 C：放宽 `dis_conv_tol`

`dis_conv_tol` 默认 `1e-10`。**单位陷阱**：wannier90 手册写 "Units are Å²"，但源码
`src/disentangle.F90:1988` 里 `delta_womegai = womegai1/womegai - 1.0_dp` 是**相对变化**，
`.wout` 表头也写 `Delta (frac.)`。手册与代码不一致，以代码为准——所以 1e-10 是**相对**判据。

- **支持**：这个默认值不是为金属 d 带调的。**而且实测差距很小**：88 个通道到 200 步时 `|Delta|` 中位 **5.67e-6**、最小 1.16e-8、最大 1.32e-3——放宽到 1e-6 只能救 18/88，放宽到 1e-5 能救 56/88，放宽到 1e-3 才几乎全过。也就是说 **1e-5 ~ 1e-4 这一档是"能过门"和"还算收敛"的分界**，1e-6 基本救不了任何东西。
- **反对**：**这是最危险的方案**——它让"没收敛"变成"收敛"，而判据本身没有任何物理含义。如果最后选它，必须同时证明 Ω_I 在 1e-8 和 1e-10 之间对 J 的影响小于某个阈值。

#### 方案 D：用 VASP 原生的 `LWANNIER90_AUTO_WINDOW`（**可能是最干净的**）

VASP 6.4.3 有一个本仓库代码里从未用过的 INCAR 标签（`*.py`/`*.sh` grep 0 命中；只有这份文档提到它）：

```
LWANNIER90_AUTO_WINDOW = .TRUE.
```

源码 `/nfs_ssd/softwares/vasp/vasp.6.4.3/src/mlwf.F`（**2026-10-26 逐行复核，行号已校正**）：

- `:297` `CALL PROCESS_INCAR(..., 'LWANNIER90_AUTO_WINDOW', ...)`，`:132` 声明默认 `.FALSE.`；
- `:569` `IF (LWANNIER90_AUTO_WINDOW) CALL ESTIMATE_AUTO_WINDOW(MLWF,W)`，在 `:560` 写完本征值之后；
- `:884-909` 的 `ESTIMATE_AUTO_WINDOW`：`:894` 先赋 `CUTOFF_DES%MU(1)`（**下一行就被覆盖，是残留代码**），`:902-903` 对 k/spin/band 取 min/max，`:907` 才是生效的 `dis_froz_max = ELIM(1,NUM_WANN) - 1e-5`；
- `:1821-1823` 把 `dis_froz_max` 写进 `<seed>.win`（**不是** `:1818-1822`，那是 `dis_froz_min` 分支；四个窗口分别写在 `:1812-1814` / `:1815-1817` / `:1818-1820` / `:1821-1823`）；
- `:3325-3328` 的 `FILTER_WIN_FILE` 确实会丢掉同名的 `dis_froz_max` 行——**但它过滤的不是磁盘上的 `wannier90.win`，而是 INCAR 里 `WANNIER90_WIN` 那个字符串**（调用点 `:1787` `CALL FILTER_WIN_FILE(MLWF, WANNIER90_WIN, WANNIER90_WIN_FILTERED)`）。所以"VASP 算出的值优先"只对 `WANNIER90_WIN` 这条路径成立。

> ### ⚠ 方案 D 现在**直接加这个标签是空操作**——这是本方案最大的坑，必须先解决
>
> `mlwf.F:300-311`：
>
> ```fortran
> WANNIER90_WIN = ""
> CALL PROCESS_INCAR(LOPEN, IU5, IU5, 'WANNIER90_WIN', WANNIER90_WIN, ...)   ! :300-301
> IF (WANNIER90().AND.IERR/=0) THEN                                          ! :308
>    INQUIRE(FILE=DIR_APP(1:DIR_LEN)//'wannier90.win',EXIST=LWIN_FOUND,SIZE=IWIN_SIZE)
>    IF (LWIN_FOUND.AND.IWIN_SIZE>0) THEN
>       WANNIER90_WIN_LEGACY_MODE = .TRUE.                                    ! :311
> ```
>
> 而写 win 的两个调用点都被这个开关挡住：
>
> ```fortran
> :596   IF (IO%IU6>=0.AND..NOT.WANNIER90_WIN_LEGACY_MODE) CALL MLWF_WANNIER90_WIN(P_MLWF)
> :1046  IF (IO%IU6>=0.AND..NOT.WANNIER90_WIN_LEGACY_MODE) CALL MLWF_WANNIER90_WIN(MLWF,LSETUP=.TRUE.)
> ```
>
> `ht-vasp` 现在正好踩中：`input_set.py:330` **预写** `wannier90.win`，INCAR 里**没有** `WANNIER90_WIN`
> （`W90_INCAR_DEFAULTS` 只有 `LWANNIER90 / LWANNIER90_RUN / LWRITE_MMN_AMN / KPAR / LREAL / NSW / IBRION / NELM / LORBIT`；
> 归档的 `2-w90_static/INCAR.gz` 里 `grep -i wannier` 只出 `LWANNIER90` 和 `LWANNIER90_RUN`）。
> → **legacy mode 打开 → VASP 根本不写 win → `LWANNIER90_AUTO_WINDOW` 算完就丢掉了。**
>
> 要让方案 D 生效，必须**同时**改掉写入方式：要么不再预写 `wannier90.win`、改把整份 win 内容塞进 INCAR 的 `WANNIER90_WIN` 字符串；要么接受 legacy 模式并另找路径。**这一条会显著抬高方案 D 的成本，评估时不能只算"加一行 INCAR"。**

**2026-10-28 补的三条源码事实（都在 `mlwf.F` 里，方案 D 评估要用）：**

1. **win 会被写两次**。`MLWF_WANNIER90_SETUP:1046` 在 SCF 前写一次（此时 `DIS_FROZ_MAX` 还是 `:1022` 刚重置的 HUGE，所以没有 auto 值）；SCF 后 `:569` `ESTIMATE_AUTO_WINDOW` 设值，`:596` 再写一次，**最终落盘的是带 auto 值的版本**。所以不存在"先写的被后写的盖掉"的问题。
2. **`wannier90.win` 之后会被复制成 `wannier90.1.win` / `wannier90.2.win`，wannier90.x 跑的是复制件**（`MLWF_WANNIER90_RUN:1892` 的 `COPYFILE`）。所以 auto 值**确实会进到真正被 wannier90 读的那个文件**——这一条排除了"D 写了但没生效"的怀疑（前提是 legacy mode 已关掉）。
3. **`FILTER_WIN_FILE` 只会在 VASP 自己也算出同名值时才丢掉用户的行**。由于 `DIS_WIN_MIN/MAX`、`DIS_FROZ_MIN` 在 VASP 里**只被赋值成 HUGE、从不被算出来**（`:1021-1023` 之后只有 `:1812/:1818` 的读取判断），用户 `WANNIER90_WIN` 里的 `dis_win_*` 和 `dis_froz_min` 会被**原样保留**；只有 `dis_froz_max` 会被 auto 值替换。也就是说方案 D 的语义精确等于"只覆盖内窗上界"——和之前从源码推的一致。

**2026-10-28 实测（BCC-Nb-Nb，E_F = 11.4766 eV，5 个变体，429e，每个约 3 min）**

| 变体 | 磁盘上有 win | `WANNIER90_WIN` | `AUTO_WINDOW` | VASP 写了 win? | 写出的 `dis_froz_max` | 投影 | legacy 提示 |
|---|---|---|---|---|---|---|---|
| A（=现状） | 有（预写） | 无 | 无 | **否**，md5 不变 | 10（原文件没被碰） | 来自文件 | 有 |
| B | 有（预写） | 无 | **有** | **否**，md5 与 A 逐字节相同 | 10（原文件没被碰） | 来自文件 | 有 |
| **C** | 无 | **有** | **有** | **是** | **14.0246854（VASP 算的）** | `Nb: s,p,d` 保留 | 无 |
| D | 无 | **有** | 无 | **是** | 10（用户值保留） | 保留 | 无 |
| E | 无 | 无 | 有 | **是** | 14.0246854 | **丢失 → `use_bloch_phases = .true.`** | 有 |

四条实测结论：

1. **方案 D 在生产配置下确实是空操作**。A 与 B 落盘的 win md5 完全相同（`fdb21cd7…`），两个变体的 OUTCAR 都打出 legacy 提示各 1 次（每个自旋一次）。`LWANNIER90_AUTO_WINDOW` 算出的值被丢弃——这是**实测**，不再是源码推论。
2. **换成 `WANNIER90_WIN` 后 D 完全按设计工作**（C）：legacy 关闭（提示 0 次），VASP 重写 win，`dis_froz_max` 被替换成 14.0246854，而 `dis_win_min/max`、`dis_froz_min`、`begin projections` 全部原样保留。
3. **"干脆不预写 win" 不是可行的替代**（E）：VASP 会自行生成一份"barebones" win，但里面是 `use_bloch_phases = .true.` 且**没有投影** —— Wannier 基完全变了，这不是方案 D 想要的结果。
4. 自动值 14.0246854 与离线从 `.eig` 算出的 14.0251 差 4e-4 eV，量级一致（`.eig` 只有 6 位小数）。

**端到端确认（变体 F = C 但 `LWANNIER90_RUN = True`，让 wannier90 真的读这份 win）**

- `wannier90.win` 里 `dis_froz_max = 14.0249794`；随后复制出的 `wannier90.1.win` / `wannier90.2.win`
  **两份 md5 相同且都带这个值** → auto 值确实到达 wannier90.x 的输入。
- `wannier90.{1,2}.wout` 头部：`Outer: -50.00000 to 50.00000`、**`Inner: -10.00000 to 14.02498`**
  —— auto 窗是活的。
- 结果：**仍然不收敛**（200 次撞上限），`Final Omega_I = 14.73875157 / 14.73873761`。

对照基线（现状窗 `[-10, 10]`）的 Ω_I = 12.2939 —— auto 窗把散度从 12.29 抬到 **14.74**，
和离线 Dpure 的预测（1000 次 14.39）方向一致。**所以"让 D 生效"和"让 D 收敛"是两件独立的事。**

（顺带：F 用 24 ranks，VASP 把 NBANDS 从 64 向上取整到 72，auto 值随之从 14.0246854 变成
14.0249794，差 3e-4 eV —— 说明这个量对 SCF 细节有 ~1e-4 相对量级的敏感性，报数时别写到小数点后五位。）

**但要注意 C 的冻结窗实际是 `[-10, 14.02]`**，不是"从带底到 14.02"：用户 win 里的 `dis_froz_min = -10` 被保留，VASP 只覆盖上界。这个组合正是离线测过**不收敛**的那个（1000 次 `|Δ| = 5.65e-6`）。所以"让 D 生效"和"让 D 收敛"是两件事。

**但先说清楚它解决的是什么问题。** 窗口错位只对 **E_F > `dis_froz_max` = 10 eV** 的体系成立，44 个里只有 **3 个**：

```
BCC-Al-Nb   E_F = 10.20
BCC-Nb-Nb   E_F = 11.48
SER-Nb      E_F = 11.25
```

**其余 41 个的 E_F 都低于 10 eV，窗口确实跨过 E_F。** 所以"窗口放错位置"解释不了 44/44 全部不收敛——这是我（派发人）在 §2.6 最初的说法，**它被数据推翻了，见 §2.11**。方案 E 仍然值得测，但它不是"显然的正解"。

**这解决了方案 A 最大的实现障碍**：E_F / 带边是 SCF 之后才知道的，而 `wannier90.win` 是 SCF 之前写的；但 `LWANNIER90_AUTO_WINDOW` 让 VASP 在 SCF 之后自己算、自己覆盖，不需要两次 SCF。

实测它会给什么（从各端基的 `wannier90.1.eig` 直接算 `min_k E(band=NUM_WANN)`）：

| 端基      |               num_wann | auto `dis_froz_max` (eV) | 现在硬编码 |
| --------- | ---------------------: | -------------------------: | ---------: |
| SER-Mn    | **9**（不是 18） |                       6.70 |         10 |
| BCC-Ti-Ti |                     18 |                       7.96 |         10 |
| BCC-Nb-Nb |                     18 |            **14.03** |         10 |
| BCC-Co-Fe |                     18 |            **15.85** |         10 |
| SER-Al    |                      9 |            **25.70** |         10 |
| BCC-Al-Al |                     18 |            **25.92** |         10 |

**44 个里 36 个会得到比 10 eV 更高的冻结上限**，8 个更低（2026-10-26 独立复算一致）。低于 10 eV 的 8 个：
`SER-Mn 6.70`、`SER-Ti 7.23`、`SER-V 7.30`、`BCC-V-V 7.71`、`BCC-Ti-Ti 7.96`、`BCC-Mn-Mn 8.03`、`FCC-Mn-Mn 8.26`、`BCC-Ti-V 8.43`。
也就是说现在这个硬编码值对绝大多数体系是**偏低**的——冻住的态比应该冻的少。

（表里 `SER-Mn` 的 `num_wann` 是 **9**，不是 18；`min E(band=18)` 是 42.94 eV。6.70 来自 band 9。其余五行核对无误。）

- **支持**：VASP 原生、无需两次 SCF、按体系自适应、有源码可读。
- **未知/反对**：
  - 它**只设 `dis_froz_max`**，不设 `dis_froz_min`、也不设 `dis_win_min/max`。所以 `dis_win` 仍然要么留空（用默认 = 全部本征值），要么另想办法。
  - `ESTIMATE_AUTO_WINDOW` 的第一行 `MLWF%DIS_FROZ_MAX = CUTOFF_DES%MU(1)` 会被下一行覆盖，**看起来是残留代码**；实际生效的是 `ELIM(1, NUM_WANN) - 1e-5`。
  - 这个标签的文档极少，**必须先在 `BCC-Nb-Nb` 上实测一次**，确认它真的改了 win 且 Wannier 收敛了，再谈推广。
  - 对 `num_wann` 的选择依赖很强：`num_wann=18`（spd 全要）时 band 18 的最低能量对 Nb 是 14 eV，如果 `num_wann` 取小了，冻结上限会掉进价带里。

#### 方案 E：降 W90 SCF 的 k 点数（**§2.11 的证据指向这里**）

- **支持**：`tb2j_test` 的单变量对照——同一 Fe bcc、同一窗口、同样 9/16，**64 kpts 收敛、512 kpts 不收敛**。而 44 个归档 run 全部是 512 个 k 点。
- **反对/未知**：k 点数直接决定 `hr.dat` 的实空间插值质量，降到 4×4×4 可能让 J 的 R 空间收敛变差（TB2J 的 `--kmesh 9 9 9` 是另一层插值，两者要分开看）。而且"k 点少所以容易收敛"未必是好事——解纠缠在一个欠采样的流形上收敛，得到的 Wannier 函数可能物理上不对。
- **必须做的判别**：同一个端基（`BCC-Nb-Nb`）在 4×4×4 / 6×6×6 / 8×8×8 三个 W90 k 网格下，比较 `wannier_converged`、`Final Omega_I`、`max|J_iso|`、`BMAGN`、`TC`。如果 4×4×4 收敛但 J 与 8×8×8 差很多，说明"收敛"是欠采样的假象。

#### 方案 E 实测（2026-10-29，BCC-Nb-Nb，只改 KPOINTS，其余全部相同）

`WANNIER90_WIN` 标签路径（新 `input_set.py`），32 ranks（NBANDS 保持 64），三个 k 网格：

| W90 k 网格 | k 点 | 收敛 | 迭代 | `Final Omega_I` (s1 / s2) | 末步 \|Δ\| | E_F | max\|J\| | Σ\|J\| |
|---|---:|---|---:|---|---:|---:|---:|---:|
| 4×4×4 | 64 | **否** | 200 | 11.36487 / 11.36486 | 1.83e-5 | 11.3624 | **26.38 meV** | **1624.2 meV** |
| 6×6×6 | 216 | **否** | 200 | 11.93796 / 11.93797 | 5.78e-5 | 11.4452 | **0.0000** | **0.000** |
| 8×8×8 | 512 | **否** | 200 | **12.29390697 / 12.31226910** | 2.15e-4 | 11.4766 | 0.0026 meV | 0.119 meV |
| 归档 8×8×8（legacy 路径） | 512 | 否 | 200 | 12.29390697 / 12.31226910 | 2.15e-4 | 11.4767 | 0.0347 meV | 1.0 meV |

**三条结论：**

1. **降 k 点不解决收敛。** 4×4×4 一样撞 200 次上限。`tb2j_test` 里 Fe bcc 的"64 kpts 收敛、512 kpts 不收敛"**不能外推到这里**。
2. **8×8×8 那一行同时是代码改动的回归证明**：`WANNIER90_WIN` 路径给出的 Ω_I 与归档 legacy 路径**逐位相同**（12.29390697 / 12.31226910），E_F 也一致（11.4766 vs 11.4767）。所以换标签是行为等价的。
3. **J 完全不可用**：三个网格给出 max|J| = 26.4 / 0.0000 / 0.0026 meV、Σ|J| = 1624 / 0 / 0.119 meV —— 跨三个数量级的随机跳变。

**为什么 J 是这样，可以直接测出来：** BCC-Nb-Nb 的 SCF 在三个网格上都收敛到**非磁**（末态 `magnetization (x)` 每格点都是 0），所以 `H_up` 与 `H_dn` 本应完全相同。实际的两套 Wannier 哈密顿量之差：

| k 网格 | hr.dat 行数 | max\|H_up−H_dn\| | p99 | 中位 |
|---|---:|---:|---:|---:|
| 4×4×4 | 40 500 | **17 789 meV** | 752 meV | 9.40 meV |
| 6×6×6 | 111 132 | 18.6 meV | 0.252 meV | 0.0022 meV |
| 8×8×8 | 236 196 | 284 meV | 2.76 meV | 0.059 meV |

**J（0.003–26 meV）与这个非物理的 up/down 不对称同量级甚至更小。** 也就是说：J 是两次**独立且都未收敛**的解纠缠落在不同子空间所产生的伪信号，不是交换耦合。

**这条实测把方案 E 的判别力也削掉了**：在非磁端基上，"收敛了没有"可以判（都没收敛），"J 变没变"判不了（J 本来就没有物理意义）。要真正判 E，必须换**磁性端基**（`BCC-Co-Fe` max|J|=23 meV、`SER-Co` 71 meV 那一档）。

顺带：8×8×8 的 J（0.0026 meV）与归档（0.0347 meV）差 13×，尽管 Ω_I 逐位相同 —— 这再次说明该端基的 J 在 Wannier gauge 的复现误差之下。

#### 全池离线扫（2026-10-29，44 端基 × 2 自旋 = 88 通道，零 VASP）

用 `scripts/w90/w90sweep.py` 对归档的 `mmn/amn/eig` 重跑 wannier90（不带 `.chk`，
Ω_I 与归档逐位一致，见回归）。两个配置：

- **BC**：`dis_num_iter=2000` + `dis_conv_tol=1e-6`（窗保持现状 `[-50,50] / [-10,10]`）
- **A2**：逐端基窗 `dis_win = E_F−25..E_F+45`、`dis_froz_min = E_F−10`、
  `dis_froz_max = min_k E(band=num_wann) − 1e-5`、`dis_num_iter=1000`

| 配置 | 通道收敛 | 双自旋都收敛的端基 |
|---|---:|---:|
| BC | **81 / 88** | **39 / 44** |
| A2 | 38 / 88 | 18 / 44 |
| 现状（归档，iter=200, tol=1e-10） | 0 / 88 | 0 / 44 |

**结论 1：BC 是唯一可用的配置。** A2 的 18 个端基是 BC 的 39 个的真子集（交集 18、并集 39），
没有任何一个端基是"只有 A2 能过"的。而且**在两配置都收敛的 38 个通道上，A2 的 Ω_I
无一例外比 BC 差**（中位 +14.7%，最差 FCC-Ni-Ni +92.7%）。也就是说 A2 不只是覆盖率低，
它收敛到的是更差的子空间 —— 用 Ω_I 判，BC 完胜。

**结论 2：BC 的成本是可控的。** 81 个收敛通道的迭代数：min 93 / 中位 192 / p75 365 / max 1508。
**只有 4 个通道需要 >1000 次**，1 个需要 >1500 次。所以 `dis_num_iter=2000` 对 99% 的通道是浪费，
`dis_num_iter=600` 就能覆盖 p75 以上。

**结论 3：5 个端基两个配置都不过，且分成两类，处理方式不同。**

| 端基 | BC 失败通道 | 末步 \|Δ\| | 距 1e-6 | 性质 |
|---|---|---:|---:|---|
| `BCC-Mn-Mn` | s2 | 4.83e-6 | 5× | 差一点 |
| `BCC-Ti-V` | s1 / s2 | 5.21e-6 / 2.34e-6 | 5× / 2× | 差一点 |
| `FCC-Co-Mn` | s1 | 1.73e-4 | 173× | **真平台** |
| `FCC-Fe-Mn` | s1 | 3.09e-4 | 309× | **真平台** |
| `FCC-Ni-Mn` | s1 / s2 | 2.39e-4 / 3.11e-4 | 239× / 312× | **真平台** |

前两个（3 个通道）只要把 `dis_num_iter` 从 2000 提到 ~4000，或把容差放到 1e-5 就能过；
后三个（4 个通道）是含 Mn 的 FCC，末步 |Δ| 卡在 1e-4 量级，**加迭代和放宽容差都救不了**，
需要单独处理（换窗 / 换 `num_wann` / 判为不可用）。

**这直接对应"多层门控"：**

- **第一层（通用）**：`dis_num_iter=2000` + `dis_conv_tol=1e-6`，覆盖 **39/44** 端基。
- **第二层（少量例外）**：`BCC-Mn-Mn`、`BCC-Ti-V` 走"更多迭代 / 更宽容差"。
- **第三层（真失败）**：`FCC-Co-Mn`、`FCC-Fe-Mn`、`FCC-Ni-Mn` 需要单独决策，不能用通用门放行。

**一个必须同时做的改动：** 现在的 `dis_conv_tol=1e-10` 是**全池 0/88**。要落地 BC，
`input_set.py` 里的 `dis_num_iter` / `dis_conv_tol` 必须真的写进 win（现在两个都没写，
靠 wannier90 默认 200 / 1e-10）。这是方案 B+C 的**必要前提**，不是可选项。

#### 方案 F：改 `num_wann` / `num_bands` 的比例

现在是 `num_wann = n_atoms × 9`、`num_bands = ceil((num_wann+20)/ntasks)×ntasks`（`BCC-Nb-Nb`: 18/64，46 个空态选 18 个）。

- **支持**：空态太多会让解纠缠的自由度过大。早期 Fe 原型 9/16（7 个空态）收敛得很轻松。
- **反对**：`num_bands` 太少会截断 d 带；`NBANDS` 还受 Wannier90 接口的 NPAR 约束（`KPAR=1, NCORE=1` 时 NPAR=ntasks）。
- **未知**：`num_wann = n_atoms × 9`（spd 全要）对过渡金属是否合适？Nb 的 4d 是价态，但 `Nb_pv` 的 ZVAL=11 包含 4s4p；Wannier 投影用 `Nb: s,p,d` 是否把半芯态也拉进来了？

### 4. 你要做的事

**第一阶段：只读复核（不写代码）**

1. 复核 §2.5–2.9 的每一个数字。**自己从文件算，不要相信本文档的数字。** 特别是：
   - 44/44 的 grep 计数；
   - `BCC-Nb-Nb` 的 E_F 与 band 范围；
   - 190→200 的衰减率外推（注意这是外推，不是实测）；
   - 平台化那 9 个的判定（衰减率 ≥ 1 是否稳定，还是只是 10 步窗口的噪声？）。
2. 查清**方案 A 的实现障碍**：VASP 6.4.3 的 `WANNIER90_WIN` 标签能不能在 SCF 之后写 win？源码在 `/nfs_ssd/softwares/vasp/vasp.6.4.3/src/`，wannier90 库源码在 `/nfs_ssd/softwares/wannier90/wannier90-3.1.0/src/`。这是方案 A 能不能成立的前提。
3. 查 `assigned_basis.txt`：`data/endmembers/*/tb2jflow/3-tb2j_solve/assigned_basis.txt` 里的轨道归属对不对？`wannier90.1_centres.xyz` 的 WF 中心离原子多远？

**第二阶段：设计判别实验（仍然不跑）**

设计一个**最小**实验，能区分 A/B/C/D。要求：

- 用**一个**端基（建议 `BCC-Nb-Nb`：E_F=11.48 离窗口最远、且有完整归档可对比）；
- 只改**一个**参数（或一组必须同时改的参数），其余全部固定；
- 明确写出"如果方案 X 对，应该看到什么；如果错，应该看到什么"；
- 估算成本（VASP 侧是 R3 relax + W90 static，实测 0.06–0.14 h/端基）。

**把实验设计 + 成本 + 预期判别结果发给派发人，等书面批准。**

**第三阶段：批准后**

只跑批准的那个实验。产出报告要包含：改前/改后的 `wannier_converged`、`Final Omega_I`、`max|J_iso|`、`BMAGN`、`TC`，以及 `wannier90.{1,2}.wout` 的关键行。

### 5. 硬约束

- **在派发人书面确认前，不改 `htvasp/tb2j/input_set.py`、不跑 TB2J。**
- 不改 `data/endmembers/` 下任何已归档 run（包括不删、不覆盖、不"顺手修"）。
- 不改 `tb2j_test/`（那是历史证据，它的价值就在于保留了原始日志）。
- 不碰 `em-tdb-v2`。
- 不 commit 任何参数改动，直到方案锁定。
- 每个阶段结束都要给派发人一份结论 + 待确认项，**不要连续推进两个阶段**。

### 6. 你必须反复与派发人确认的四个问题

-1. **方案 E（降 k 点数）和方案 D（auto window）哪个先测**：§2.11 的单变量对照把嫌疑指向 k 点数，但"降 k 点让解纠缠收敛"可能只是欠采样的假象。我（派发人）倾向**两个都测，用同一个端基**，因为它们成本都很低（一次 R3+W90 = 0.06–0.14 h）。
0. **方案 D（`LWANNIER90_AUTO_WINDOW`）是不是首选**：它是 VASP 原生、不用两次 SCF、按体系自适应。但它只设 `dis_froz_max` 一个量，而且几乎没文档。要不要先只测它？

1. **方案 A 的实现路径**：如果 `WANNIER90_WIN` 不能在 SCF 后覆盖，方案 A 就要两次 SCF，成本翻倍。这值得吗？还是接受方案 B 的部分收敛 + 显式标记？
2. **9/44 平台化的那批怎么办**：如果提高迭代也救不了，它们的 J 要不要直接判为不可用？还是换窗口再试？
3. **`dis_conv_tol` 能不能动**：我（派发人）倾向**不动**——它是唯一的物理外锚，放宽它等于放弃判据。如果你认为必须动，需要给出"Omega_I 的相对变化从 1e-8 收到 1e-10 时 J 变化小于 X%"的实测。
4. **收敛了是否等于物理对**：`BCC-Co-Fe`（2.816 μB，J/噪声 2.84）的 J 看起来可信，但它同样没收敛。如果修好收敛后 J 变了 20%，那么**之前所有基于 TB2J 的结论都要重做**——这个范围有多大？

---

## 事实附录（派发人核对用，不要复制）

### 关键文件位置

| 内容                       | 路径                                                                                                                                                                                                                                                     |
| -------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `.win` 生成器            | `htvasp/tb2j/input_set.py:240-280`（`get_wannier90_win_string`），**`265-268`** 是四行窗口（`263-264` 是 `num_wann`/`num_bands`）                                                                                                      |
| VASP 自动窗口实现          | `/nfs_ssd/softwares/vasp/vasp.6.4.3/src/mlwf.F:884-909`（`ESTIMATE_AUTO_WINDOW`）、`:3325-3328`（`FILTER_WIN_FILE`）                                                                                                                             |
| win 落盘                   | `htvasp/tb2j/input_set.py:330`（`Path("wannier90.win").write_text(...)`）                                                                                                                                                                            |
| TB2J 调用                  | 在 `master` / `feat-tb2j-sqs` 上（558 行）：`run_wann2j` 在 `jobs.py:251`、`tb2j_solve` 在 `:486`。（旧文档写的 `:111`/`:263` 是 `origin/feat-tb2j` 那个 328 行旧版的编号，PR #3 之后已作废。）**按函数名 grep，不要按行号跳。** |
| 收敛判据（wannier90 源码） | `/nfs_ssd/softwares/wannier90/wannier90-3.1.0/src/disentangle.F90:2298`                                                                                                                                                                                |
| 失败打印                   | 同上 `:2041-2044`                                                                                                                                                                                                                                      |
| 参数文档                   | `/nfs_ssd/softwares/wannier90/wannier90-3.1.0/doc/user_guide/parameters.tex:832-899`                                                                                                                                                                   |
| 原始脚本模板               | `tb2j_test/submit_tb2j.sh:204-219`                                                                                                                                                                                                                     |
| 原始试验日志               | `tb2j_test/*/wannier90.{1,2}.wout`                                                                                                                                                                                                                     |
| 归档结果                   | `data/endmembers/*/tb2jflow/2-w90_static/`、`3-tb2j_solve/`                                                                                                                                                                                          |

### 迭代需求外推的原始表（190→200 局部衰减率）

35 个会收敛：`SER-Al` 307、`SER-Ni` 363、`BCC-Co-Co` 437、`BCC-Co-Ni` 439、`BCC-Fe-Fe` 440、`SER-Fe` 441、`BCC-Ni-Ni` 442、`BCC-Co-Fe` 487、`SER-Ti` 496、`BCC-Fe-Ni` 506、`BCC-Al-Ti` 570、`SER-Mn` 682、`SER-Nb` 727、`BCC-Al-V` 758、`BCC-Mn-Mn` 894、`BCC-V-V` 915、`FCC-Ni-Ni` 995、`FCC-Co-Ni` 1001、`FCC-Fe-Ni` 1015、`FCC-Ni-Co` 1078、`FCC-Mn-Mn` 1083、`FCC-Co-Co` 1115、`BCC-Ti-V` 1125、`BCC-Nb-Nb` 1176、`FCC-Fe-Co` 1216、`FCC-Ni-Fe` 1279、`FCC-Co-Fe` 1311、`FCC-Fe-Fe` 1357、`FCC-Fe-Mn` 2911、`BCC-Nb-V` 2978、`FCC-Mn-Ni` 3303、`FCC-Co-Mn` 3475、`FCC-Ni-Mn` 3934、`BCC-Ti-Ti` 5247、`BCC-Al-Nb` 12263

9 个平台化：`BCC-Al-Al`、`BCC-Co-Mn`、`BCC-Fe-Mn`、`BCC-Mn-Ni`、`BCC-Nb-Ti`、`FCC-Mn-Co`、`FCC-Mn-Fe`、`SER-Co`、`SER-V`

**外推的局限**：这是按 190→200 的局部指数衰减率外推的，如果实际是次指数衰减（很可能），上表是**乐观估计**。要得到实测值必须真跑。

### 与本次调查无关但相关的两个已知缺陷（不要顺手改）

1. `htvasp/tb2j/context.py:_bcc_context` 给纯元素 B2 胞（`BCC-Nb-Nb`）也建了两个子晶格 `{Nb1:1, Nb2:0}`、`sublattice_ratios [0.5,0.5]`。物理上纯元素只有一个子晶格。
2. 44 个 TB2J JSON 里 42 个的 `structure` 块是未弛豫模板（`BCC-Nb-Nb`: `structure.lattice.a = 2.8597` vs `cell[0][0] = 3.3031`，a 差 13.4%、体积差约 24%），而 em-tdb-v2 的 `read_tb2j_json`（`em-tdb-v2/src/emtdb_v2/magnetic/readers.py:611`，required 元组在 `:624`）**只读 `structure`，全库 grep `cell` 0 命中**。两个例外（`BCC-Fe-Ni`、`SER-Al`）只是模板恰好等于弛豫值，不是元数据正确。
