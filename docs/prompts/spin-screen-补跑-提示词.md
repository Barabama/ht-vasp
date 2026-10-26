# 提示词：补跑 66 个端基的 spin screen

> **用法**：把下面「提示词正文」整段复制给一个新的 AI 会话。
> 目标仓库：`/nfs_hdd/2025/gaominliang/ht-vasp`。
> 本文档最后一节「事实附录」是给**你**（派发人）核对用的，不要复制进去。

---

## 提示词正文

你负责给 `ht-vasp` 的 7 元素端基池补跑 **spin screen**（自旋态筛选）。

### 0. 先读

按顺序：

1. `/nfs_hdd/2025/gaominliang/AGENTS.md`
2. `/nfs_hdd/2025/gaominliang/ht-vasp/CLAUDE.md`
3. `/nfs_hdd/2025/gaominliang/ht-vasp/docs/TB2J与SQS数据审计-20261025.md` —— §3 与 §4 说明了为什么这件事必须先做
4. `/nfs_hdd/2025/gaominliang/ht-vasp/main_spin_screen.py` —— **你要用的工具，从头读到尾**
5. `/nfs_hdd/2025/gaominliang/ht-vasp/docs/legacy-qha-rebuild-report.md` §8.0a(j)（L1₂ 上 Cr 的两个方向差 23 vs 184 meV/atom）

然后 `git status`。**不要 reset、clean、stash，不要覆盖或删除任何已存在的文件。**

### 1. 分支与环境

当前在 `feat-dfpt-qha`。**这个分支是对的，不要切走**——`main_spin_screen.py` 和 `htvasp/dfpt_qha/` 只在这个分支上（`master` 和 `feat-tb2j` 都没有）。

```bash
git branch --show-current          # 应该是 feat-dfpt-qha
./.conda/bin/python main_spin_screen.py --dry-run --batch BCC-Al-Al-2
```

### 2. 现状：84 个端基里只筛过 18 个

`pool_recipes.json` 的 `magnetic_policy` 分布是 15 个 `confirmed_branch` + 69 个 `nonmagnetic`，看起来像全筛过了。**不是。**

- 真正被 `spin screen 2026-10-20` 测过的只有 **18 个**，全部含 Cr。
- 剩下 **66 个**的 `nonmagnetic` 是 **"默认预期"**，其中 57 个的 `magnetic_note` 原文就写着 default expectation，6 个连 note 都没有。
- 那 18 个的 `spin-screen.json` **在 429pro 上已经没有了**，但在 **B5070 的 `/workspace/gaominliang/ht-vasp/data/endpoints/*/spin-screen.json`** 还完整保留着（18 个）。**这是重要的历史证据，先取回来。**

### 3. 你要做的事

**第一步：把 18 个历史 spin-screen.json 从 B5070 取回来。**

B5070 是另一台机器、另一个用户、只有密码登录，**不能 `ssh B5070`**。用 `/tmp/b5070.py` 这个 paramiko 辅助脚本（如果它还在；不在就照 `ht-vasp/scripts/nas/nas.py` 的写法自己写一个，凭据在 `/nfs_hdd/2025/gaominliang/.env` 的 `SERVER_B5070_*`）：

```bash
./.conda/bin/python /tmp/b5070.py 'ls data/endpoints/*/spin-screen.json'
```

取回到 `data/endpoints/<name>/spin-screen.json`（**不要覆盖任何已有文件**，本地那 18 个目录里现在只有 `magnetic-state.json` 和 `dfpt-qhaflow/`，没有 spin-screen.json）。同时检查 B5070 上有没有 `spin-screenflow/` 目录——有就一起取回 store（那是原始 TaskDoc，`--reanalyse` 要用）。

**第二步：确认 66 个的名单。**

`docs/prompts/spin-screen-unscreened.txt` 是已经算好的清单（66 个端基、133 个自旋态、含每个端基的 `n_Cr` 和该跑哪几个态）。**用它之前先自己核对一遍**：

```python
# 期望：66 个，与 pool_recipes.json 的 84 个减去 /tmp/spin_results.txt 的 18 个一致
```

三个含 Cr 但没筛过的要特别注意：

| 端基 | 当前 policy | 依据 | 问题 |
|---|---|---|---|
| `BCC-Al-Cr-2` | `confirmed_branch` | report 4.9.3，283.7 meV/atom | 单个体积、未弛豫胞 |
| `FCC-Al-Cr-4` | `confirmed_branch` | report 8.0a(j)，22.96 meV/atom | **归档的 DFPT run 跑成了非磁** |
| `FCC-Cr-Al-4` | `confirmed_branch` | report 8.0a(j)，183.89 meV/atom | 单个体积 |

这三个的"磁性"结论来自报告，不是来自 spin screen。`FCC-Al-Cr-4` 尤其要重筛：它的 DFPT-QHA run 从 10 μB 起步、45 步 SCF 内塌到 0，而报告测 22.96 meV 时用的是**未弛豫**的 13.03 Å³/atom 胞，DFPT 跑到的却是 12.43——两者不可比。

**第三步：跑。**

工具是 `main_spin_screen.py`，它已经能做全部三件事（nm / fm / afm），并且有 `--slurm`。**机器和 Slurm 参数由你决定**，但要基于事实而不是猜：

**机器的关键约束（2026-10-25 实测）：**

| 节点 | 核数 | 内存 | 备注 |
|---|---:|---:|---|
| 429e | 48 | **192 GB** | 跑过 NELECT 43 的 DFPT 没 OOM |
| 429pro | 64 | **128 GB** | 32 核作业并发 2 个时发生过 21 次内核 OOM 杀 vasp_std |

- 两个节点都在 `partCPU`，**没有独立的分区可选**。`--nodelist 429e` 和 `--nodelist 429pro` 都是投到 `partCPU` 再钉节点。
- 429pro 上长期有别人的作业（`fx-dpgen-test` 已跑数天）。spin screen 是 2/4 原子胞的小作业，内存压力远小于 DFPT 的 16 原子 DFPT 步骤，但仍建议**钉 429e**，或至少在 429pro 上限制并发。
- `SlurmJobManager.get_cpu_config()` 支持 `nodelist`、`memory`、`time_limit`、`ntasks` 全部作为关键字参数；`submit_screen()` 目前**没有**把 `nodelist` 传下去（它写死 `memory="96G"`、`time_limit="24:00:00"`），需要的话你改这个函数。

**成本（实测外推）：**

- 一次 spin screen = 每个态一次 `ISIF=3` 短弛豫（`DfptRelaxMaker`，`NSW` 上限 80）。
- DFPT-QHA 的 EOS relax 实测 0.5–3.5 min CPU（已收敛的续跑，2–3 步）；SQS 从 45% 过密的模板出发是 38–55 min（33–40 步）。spin screen 的起点是 `pool_structures.json` 里的有序胞，介于两者之间。
- 66 端基 × 2–3 态 = **133 个弛豫**，粗估 20–45 CPU 小时；32 核作业并行 3–4 个，wall clock 约 1–2 h。

**跑之前必须做的检查：**

1. `squeue -p partCPU` 看有没有别人的大作业。
2. `--dry-run` 一次，确认 flow 里真的是 `spin nm` / `spin fm` / `spin afm` 三个 job，且 afm 只出现在 ≥2 个 Cr 的胞里。
3. 确认 `flow_dir` 是 `/nfs_ssd/tmp`（脚本里写死），磁盘够（`df -h /nfs_ssd`）。

**第四步：验收。** 每个端基要产出 `data/endpoints/<name>/spin-screen.json`，里面 `decision` ∈ {`confirmed_branch`, `nonmagnetic`}，`magnetic_gain_meV_per_atom` 是数值，`excluded_states` 列出来。**失败的端基要单独列出来**（`ok: false` 带 `error`），不要当成"跑过了"。

**第五步：不要动 `pool_recipes.json`。** 脚本自己的 docstring 写着："Updating `pool_recipes.json` from that file is a separate, explicit step." 把 66 个的新决策汇总成一张表放在报告里，**由人决定要不要写回**。

### 4. 硬约束

- **不改** `main_spin_screen.py` 的判定逻辑（`DECISION_THRESHOLD_MEV_PER_ATOM = 5.0`、`STRUCTURE_TOLERANCE = 0.02`）。如果你认为阈值要改，写进报告，不要动。
- **不改** `MAGMOM` 表、不改 `spin_states()` 的初猜值。
- **不改** 任何已归档的 `dfpt-qhaflow/` 结果。
- **不提交** 任何 DFPT-QHA / SQS / TB2J 作业——本任务只做 spin screen。
- **不碰** `em-tdb-v2`。
- 提交信息用 `feat(spin-screen):` 或 `chore(spin-screen):`。

### 5. 报告要回答的问题

1. 66 个跑完后，`confirmed_branch` 有多少个？与现有的 15 个相比增删了哪些？
2. 有没有端基的 `magnetic_gain` 落在 5 meV/atom 阈值附近（比如 2–10 meV/atom）？这些是最需要人工看的。
3. `FCC-Al-Cr-4` 在**弛豫后的体积**上还是非磁吗？
4. `--constrained-afm` 模式（在 fm 弛豫几何上固定体积跑 nm 和 afm）要不要对哪些端基跑？`pool_recipes.json` 里 `FCC-Hf-Cr-4`、`FCC-Zr-Cr-4` 的 note 提到"free afm collapsed the cell"，`FCC-Nb-Cr-4`、`FCC-Ta-Cr-4` 的 note 说 constrained afm 比 fm 高 0.3–1.0 meV/atom。这些说明**自由 AFM 会改结构**，constrained 模式可能更适合它们。
5. 66 个里有没有端基的胞在某个自旋态下 `structure_changed`（`lattice_deviation_vs_nm > 0.02`）？被 `excluded_states` 排除的有哪些？

---

## 事实附录（派发人核对用，不要复制）

### 18 个已筛的端基

```
BCC-Cr-Cr-2  BCC-Cr-Hf-2  BCC-Cr-Nb-2  BCC-Cr-Ta-2  BCC-Cr-Ti-2  BCC-Cr-Zr-2
FCC-Cr-Cr-4  FCC-Cr-Hf-4  FCC-Cr-Nb-4  FCC-Cr-Ta-4  FCC-Cr-Ti-4  FCC-Cr-Zr-4
FCC-Hf-Cr-4  FCC-Nb-Cr-4  FCC-Ta-Cr-4  FCC-Ti-Cr-4  FCC-Zr-Cr-4  SER-Cr-1
```

`/tmp/spin_results.txt` 里有全部 18 个的 `decision | gain`。**这 18 个的 `spin-screen.json` 只在 B5070 上还有**（429pro 上已无）。

### 归档的 DFPT-QHA run 与 screen 判定的一致性（已核，全对）

8 个有归档 run 的磁性端基，其收敛磁态与 screen 判定**全部一致**：BCC-Cr-Cr-2 和 SER-Cr-1 是 AFM（±1.134），其余 FM，`FCC-Al-Cr-4` 非磁（而 screen 说 confirmed_branch —— 这就是第 3 步要重查的那个）。

### 66 个的分布

- 含 Cr 3 个：`BCC-Al-Cr-2`、`FCC-Al-Cr-4`、`FCC-Cr-Al-4`
- 不含 Cr 63 个
- 状态数：65 个跑 nm+fmn，1 个跑 nm+fm+afm（`FCC-Al-Cr-4`，3 个 Cr）
- 合计 **133 个自旋态**

### `unresolved_endpoints()` 现在选不出任何端基

`pool_recipes.json` 里 `magnetic_policy == "unresolved"` 的数量是 **0**。所以 `main_spin_screen.py --batch`（不带名字）会说"nothing to do"。**必须显式传名字**：

```bash
./.conda/bin/python main_spin_screen.py --slurm \
    --batch $(awk '!/^#/{print $1}' docs/prompts/spin-screen-unscreened.txt)
```

`--batch` 是 `nargs="*"`（空格分隔的名字，不是文件路径）。

### 429e vs 429pro 的 OOM 证据

`docs/搬回429pro候选清单-20261023.md` §6.1：429pro 上累计 21 次内核 OOM 杀 `vasp_std`，6 个 FCC 端基全部失败；`journalctl -k` 显示每次截断前 6–11 秒都有 OOM。同一批作业在 429e（192 GB）上 7/7 成功。Slurm 的 `ConstrainRAMSpace=no`，所以 `--mem` 只是声明，不强制。
