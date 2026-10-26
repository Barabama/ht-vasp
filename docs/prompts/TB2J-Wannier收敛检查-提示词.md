# 提示词：把 Wannier90 收敛检查做进 TB2J 流程

> **用法**：把下面「提示词正文」整段复制给一个新的 AI 会话。
> 目标仓库：`/nfs_hdd/2025/gaominliang/ht-vasp`，分支 `feat-tb2j`。
> 本文档最后一节「事实附录」是给**你**（派发人）核对用的，不要复制进去。

> **关于本文档的路径**：这些提示词文件在 `feat-dfpt-qha` 分支的 `docs/prompts/` 下。
> 如果你按提示切到了 `feat-tb2j`（或任何别的分支），那个目录可能不存在。
> **你手上粘贴的正文就是权威版本**；如果正文里引用了某个 `docs/prompts/*.md` 而你找不到，
> 停下来向派发人索取，不要猜内容。

---

## 提示词正文

你是 `ht-vasp` 的开发者，任务是把 **Wannier90 的收敛检查**做成 TB2J 流程里的硬门。

### 0. 先读，再动手

按顺序读：

1. `/nfs_hdd/2025/gaominliang/AGENTS.md`
2. `/nfs_hdd/2025/gaominliang/ht-vasp/CLAUDE.md`
3. `/nfs_hdd/2025/gaominliang/ht-vasp/docs/TB2J与SQS数据审计-20261025.md` —— **§1 和 §1.3 是你的任务来源**，特别是 "1.3 结构层面的问题" 之前那段关于 `.wout` 的分析
4. `/nfs_hdd/2025/gaominliang/tb2j_share/TB2J磁交换计算与磁参数构建指南.md`
5. `/nfs_hdd/2025/gaominliang/sof-ace/docs/TB2J流程与工作量估算.md` §7 第 4 条（那里写了"必查 .wout spread"）

然后 `git status`，看清工作区状态。**不要 reset、不要 clean、不要 stash、不要覆盖或删除任何已存在的文件。**

### 1. 分支

当前在 `feat-dfpt-qha`。你的工作在 **`feat-tb2j`** 上：

```bash
git checkout feat-tb2j
git merge master          # 见下方"合并前必须确认"
```

**合并前必须确认的事**（不要跳过）：

- `feat-tb2j` 是 `master` 的**祖先**（`git merge-base --is-ancestor feat-tb2j master` 为真），所以 `master` 领先它 16 个提交，合并是快进或接近快进。
- 但 `feat-dfpt-qha`（你刚离开的分支）**不在** `master` 的祖先链上，它有 89 个提交是 master 没有的，其中包括 `htvasp/dfpt_qha/`、`main_dfpt_qha.py`、`main_spin_screen.py`、`htvasp/tb2j/magnetic_contract.py`、`htvasp/tb2j/context.py`。
- 也就是说：**`feat-tb2j` 合并 `master` 之后，仓库里不会有 DFPT-QHA 和 spin screen 的代码。** 这是预期状态——你的任务不需要它们。但如果你发现某个改动需要 `htvasp/dfpt_qha/`，**停下来报告**，不要自己跨分支搬文件。
- 合并后跑一次 `./.conda/bin/python -m pytest tests/ -q`，记下基线（哪些通过、哪些本来就失败）。`tests/test_nscf.py::test_flow_makers` 已知失败，与本任务无关。

### 2. 问题是什么

`htvasp/tb2j/input_set.py` 的 `get_wannier90_win_string()` 生成：

```
dis_win_min = -50
dis_win_max = 50
dis_froz_min = -10
dis_froz_max = 10
```

这是**绝对能量（eV）窗口**（wannier90 官方文档：`dis_win_min` "Units are eV"，默认值是系统最低本征值；`dis_froz_*` 同理）。

它来自 `tb2j_test/` 里 8 月 13 日的 Fe bcc 原型（E_F ≈ 1.17 eV，窗口恰好跨过 E_F），`063c12d` 原样搬进了包。7 元素池的 E_F 从 4.5 到 11.5 eV 不等，所以这个固定窗口对每个体系的意义都不同。

**窗口位置对 3 个端基是明确错的**（`BCC-Al-Nb`、`BCC-Nb-Nb`、`SER-Nb`，E_F > 10 eV，`dis_froz` 整个落在 E_F 之下）。**但另外 41 个的窗口是跨过 E_F 的，它们同样不收敛**——所以窗口不是唯一原因，别把它当成已定的结论（详见 `docs/prompts/TB2J-Wannier参数-调查链与待定方案.md` §2.11）。

**后果已经实测**（44 个归档 run，`data/endmembers/*/tb2jflow/2-w90_static/wannier90.{1,2}.wout`）：

- **44/44 个 up 通道 + 44/44 个 down 通道全部撞了 200 次迭代上限**，`Disentanglement convergence criteria not satisfied`。
- 迭代到 200 时 `Delta` 仍有 1e-8 ~ 1e-3，而判据是 **1e-10**（`dis_conv_tol` 默认值）。
  **注意单位**：wannier90 的用户手册写 `dis_conv_tol` "Units are Å²"，但源码里实际比较的是
  `delta_womegai = womegai1/womegai - 1.0_dp`（`src/disentangle.F90:1988`），即**相对变化**，
  `.wout` 的表头也写着 `Delta (frac.)`。**手册与代码不一致，以代码为准**——这条本身就是个陷阱，
  写解析器时不要按 Å² 去换算。
- 按 190→200 步的局部衰减率外推：**35/44 会收敛（中位 1001 次迭代，最多 12263 次）**，**9/44 已经平台化（衰减率 ≥ 1，再给迭代也不会收敛）**。

**而流程里没有任何一处读这个结果。** `htvasp/tb2j/` 全包 grep `wout` / `disentangle` / `Omega` 是 0 命中；`tb2j_test/submit_tb2j.sh` 的"W90 收敛检查"只查 `Vasprun.converged`（那是 VASP 的 SCF），每次都能过。

### 3. 你要做的事

**只做检查，不改窗口。** 窗口要不要改成相对 E_F、`dis_num_iter` 要不要提高，是另一个决策（见第 5 节）。你这一步的目标是：**让"没收敛"这件事再也无法被静默通过。**

具体：

1. **加一个解析函数**，读 `wannier90.{1,2}.wout`，返回结构化结果：
   - `converged`: bool —— 以文件里 wannier90 自己写的判定为准，**不要自己重算**。成功的标志是
     `<<< Disentanglement convergence criteria satisfied >>>`
     （失败是 `<<< Warning: Maximum number of disentanglement iterations reached >>>` +
     `<<< Disentanglement convergence criteria not satisfied >>>`）。
   - `n_iterations`, `final_omega_i`（`Final Omega_I  <值> (Ang^2)`）
   - `last_delta`（最后一行 `<-- DIS` 的 Delta 列）
   - `outer_window`, `inner_window`（`Outer:` / `Inner:` 两行）
   - 若 `converged=False`，再算一个 `projected_iterations`：用 190→200 步的局部衰减率外推还要多少步才到 1e-10（与 `dis_conv_tol` 同量纲，即相对变化），衰减率 ≥ 1 时返回 `None`。**这是给读者判断"值不值得提高 dis_num_iter"用的，不是判据。**

2. **在 `tb2j_solve` 里接上**（`htvasp/tb2j/jobs.py`）。它在 `prev_dir` 里跑完 `wann2J` 之后：
   - 读两个 `.wout`；
   - 把上面那组字段放进返回 dict（连同现有字段一起，不要删任何现有 key）；
   - **不收敛时怎么办，你决定，但必须显式**。可选：(a) 抛异常；(b) 照常返回但带 `wannier_converged: false` 并在 `warnings` 里写明。**如果你选 (b)，必须在同一提交里加一条测试断言下游能看见这个标志**，否则等于没做。我倾向 (a) —— 仓库其它地方是 fail-closed 风格（`validate_dfpt_task_doc`、`check_endmembers_survive` 都是直接 raise），而且这个失败已经让 44/44 个结果不可用。
   - **但注意**：直接抛异常会让 `main_w90.py` 已有的 44 个归档 run 无法 resume。所以如果你的方案会改变"重新运行"的行为，先报告再动。

3. **把判定写进产物**。`Tb2jResult`（`htvasp/tb2j/task_doc.py`）加字段，`*-tb2j.json` 里能读到 `wannier_converged` 和 `wannier_disentanglement`。**已有的 `emtdb-magnetic-dataset/v1` 顶层字段一个都不要改**（下游 em-tdb-v2 的 `read_tb2j_json` 在读），只加新字段。

4. **测试**（`tests/test_tb2j.py` 是 pytest 可收集的，12 个用例全过，照它的风格写）：
   - 解析器对"收敛"样例返回 True —— **仓库里有现成的真样例**：`/nfs_hdd/2025/gaominliang/tb2j_test/fe_bcc_wannier_4x4x4/wannier90.1.wout`（71 步收敛，`Delta < 1.000E-10 over 3 iterations`）。
   - 解析器对"没收敛"样例返回 False —— `data/endmembers/BCC-Nb-Nb/tb2jflow/2-w90_static/wannier90.1.wout`（200 步撞上限）。
   - 缺文件、空文件、格式不认识时**不猜**：返回 `None` 或 raise，不要默认 True。
   - 两个通道分别检查：只 up 收敛、down 没收敛，必须判为不收敛。

5. **写一份短报告** `docs/TB2J-Wannier收敛门-<日期>.md`：加了什么、门在哪、对已有 44 个 run 的影响（它们现在都会被判失败）、以及**你没做的事**（没改窗口、没重跑）。

### 4. 硬约束

- 不改 `htvasp/tb2j/input_set.py` 里的 `dis_win_*` / `dis_froz_*` / `dis_num_iter`。**这是下一个会话的决策，不是你的。**
- 不改 `data/endmembers/` 下任何已归档的 run。
- 不提交任何 VASP / Slurm 作业。
- 不动 `em-tdb-v2`。
- 不改 `tb2j_test/` 里的原始试验（那是历史证据）。
- 提交信息用仓库现有的 conventional-commit 风格，带 `feat(tb2j):` 或 `fix(tb2j):` 前缀。

### 5. 你要在报告里回答的问题（不要自己拍板，列出来）

1. 9/44 平台化的那批，是窗口放错了、`num_wann`/`num_bands` 比例不合适（`num_wann=18`, `num_bands=64`，46 个空态选 18 个），还是 **W90 SCF 的 k 点数太多**（44 个全是 512 点，而 `tb2j_test` 的 4×4×4=64 点收敛、8×8×8=512 点不收敛）？**这是三个不同的原因，别混在一起。**
2. `dis_conv_tol=1e-10` 按源码是**相对变化**（手册说 Å²，与代码不符）。对 Ω_I ≈ 10–22 Å² 的体系，
   1e-10 的相对变化相当于绝对 1e-9 ~ 2e-9 Å²。这个判据是不是过严？wannier90 默认就是这个值，
   但默认值不是为金属 d 带调的。
3. 提高 `dis_num_iter` 到 1000/2000 和改成相对 E_F 的窗口，哪个更可能是正解？给出你从数据里看到的理由，不要凭直觉。
4. 有没有可能"收敛了但物理仍然错"？例如窗口对了、Ω_I 收敛，但 WF 的轨道归属（`assigned_basis.txt`）跑偏。要不要同时查 `wannier90.1_centres.xyz` 与原子位置的距离？

---

## 事实附录（派发人核对用，不要复制）

### 分支拓扑（2026-10-26 实测）

```
master            ← feat-tb2j 是它的祖先（落后 16 个提交）
                  ← feat-dfpt-qha 不是它的祖先（feat-dfpt-qha 领先 89 个）
feat-tb2j         ← 只有 5 个 tb2j 文件：__init__/input_set/jobs/maker/task_doc
feat-dfpt-qha     ← 多了 context.py 与 magnetic_contract.py，以及整个 dfpt_qha/
origin/master     ← 本地 master 领先 7 个（其中 3 个是 SQS 的 add+revert，净效果只剩 0108ef7）
```

`feat-tb2j` 合并 `master` 是安全的（`git merge-tree` 报 0 冲突）。**切到 `feat-tb2j` 时你会看到 `htvasp/dfpt_qha/`、`main_dfpt_qha.py`、`main_spin_screen.py` 从工作区消失**——只有 `feat-dfpt-qha` 有这些文件。这是分支拓扑的事实，不是丢代码；`git checkout feat-dfpt-qha` 就回来。

### 44/44 没收敛的证据

```
$ grep -l "Disentanglement convergence criteria not satisfied" \
    data/endmembers/*/tb2jflow/2-w90_static/wannier90.1.wout | wc -l
44
$ grep -l "Disentanglement convergence criteria not satisfied" \
    data/endmembers/*/tb2jflow/2-w90_static/wannier90.2.wout | wc -l
44
```

（`wannier90.wout` 是 VASP 的 setup pass，15421 字节，没有 DIS 迭代，不要拿它判。）

### 迭代需求外推（190→200 局部衰减率）

| 类 | 数量 | 说明 |
|---|---:|---|
| 会收敛 | 35 | 中位 1001 次，最多 12263（`BCC-Al-Nb`） |
| 已平台化 | 9 | 衰减率 ≥ 1，再给迭代也不收敛 |

平台化的 9 个：`BCC-Al-Al`、`BCC-Co-Mn`、`BCC-Fe-Mn`、`BCC-Mn-Ni`、`BCC-Nb-Ti`、`FCC-Mn-Co`、`FCC-Mn-Fe`、`SER-Co`、`SER-V`。

### 窗口错位的直接证据

`dis_froz = -10..10` 相对 E_F 的位置：

| 端基 | E_F (eV) | dis_froz 相对 E_F |
|---|---:|---|
| BCC-Ni-Ni | 4.53 | −14.5 … +5.5（跨过 E_F） |
| SER-Co | 5.40 | −15.4 … +4.6 |
| BCC-Nb-Nb | 11.48 | **−21.5 … −1.5（全在 E_F 之下）** |
| SER-Nb | 11.25 | **−21.3 … −1.3** |
| BCC-Al-Nb | 10.20 | **−20.2 … −0.2** |

`BCC-Nb-Nb` 的 band 跨 −42.9 … 71.2 eV，`dis_win=-50..50` 只覆盖 66.8%，`dis_froz` 只覆盖 4.0%。

### 窗口是哪来的

`063c12d`（2026-08-16）把 `tb2j_test/` 里 8 月 13–14 日 Fe bcc 原型的 `.win` 原样搬进了 `htvasp/tb2j/input_set.py`。Fe 的 E_F ≈ 1.2 eV，那个窗口对它是对的。

### 早期原型为什么没暴露

| 目录 | num_wann/num_bands | 结果 |
|---|---|---|
| `tb2j_test/fe_bcc_wannier` | 9/16 | 3 步收敛（空态只有 7 个，解纠缠几乎无自由度） |
| `tb2j_test/fe_bcc_wannier_4x4x4` | 9/16 | 71 步收敛 |
| `tb2j_test/fe_bcc_wannier_8x8x8` | 9/16 | **200 步撞上限** |
| `tb2j_test/SER-Co/W90` | 18/64 | **撞上限** |
| `tb2j_test/BCC-Co-Ni/W90` | 18/64 | **撞上限** |
| `tb2j_test/FCC-Fe-Mn/W90` | 36/64 | **撞上限** |

也就是说：**同一个 Fe bcc 体系，k 网格从 4×4×4 换到 8×8×8 就从"收敛"翻成"不收敛"**；而 `num_wann/num_bands` 从 9/16 变成 18/64 之后，三个体系全部不收敛。警告一直在 `wannier90.{1,2}.wout` 里，但当时没有人读。

### 现成的真样例（给测试用）

- **收敛**：`/nfs_hdd/2025/gaominliang/tb2j_test/fe_bcc_wannier_4x4x4/wannier90.1.wout`
  第 239–242 行：`<<< Delta < 1.000E-10 over 3 iterations >>>` / `<<< Disentanglement convergence criteria satisfied >>>` / `Final Omega_I 8.13942423 (Ang^2)`
- **不收敛**：`data/endmembers/BCC-Nb-Nb/tb2jflow/2-w90_static/wannier90.1.wout`
  第 368–372 行：`<<< Warning: Maximum number of disentanglement iterations reached >>>` / `<<< Disentanglement convergence criteria not satisfied >>>` / `Final Omega_I 12.29390697 (Ang^2)`

两个文件都在本机可读，测试可以直接引用绝对路径 + `pytest.mark.skipif(not path.exists())`。

### wannier90 官方语义（3.1.0 源码 + user_guide）

- `dis_win_min/max`：外层窗口，"Units are eV"，默认是系统最低/最高本征值。
- `dis_froz_min/max`：内层冻结窗口，"Units are eV"，不给 `dis_froz_max` 就没有冻结态。
- `dis_conv_tol`：默认 `1.0E-10`。**手册写 "Units are Å²"，但代码比较的是相对变化
  `womegai1/womegai - 1.0_dp`（`disentangle.F90:1988`），`.wout` 表头也写 `Delta (frac.)`。
  以代码为准。**
- `dis_conv_window`：默认 3，即连续 3 步都小于容差才算收敛。
- `dis_num_iter`：默认 200。
- 判据代码：`src/disentangle.F90:2298`，`dis_converged = all(abs(history) .lt. dis_conv_tol)`。
- 失败打印：`src/disentangle.F90:2041-2044`。

### 已知的旁支风险（写进报告，不要顺手改）

- `htvasp/tb2j/context.py:_bcc_context` 把 B2 的 (0,0,0) 映射到 sublattice 0、(½,½,½) 映射到 1，**对纯元素 B2（`BCC-Nb-Nb`）也这么做**，于是同核胞得到 `{Nb1:1, Nb2:0}` 两个"子晶格"、`sublattice_ratios [0.5,0.5]`。物理上纯元素只有一个子晶格。
- 44 个 TB2J JSON 里 42 个的 `structure` 块是**未弛豫的模板胞**（`BCC-Nb-Nb`: `structure.lattice.a=2.8597` vs `cell[0][0]=3.3031`），而 em-tdb-v2 的 `read_tb2j_json` 读的正是 `structure`。
