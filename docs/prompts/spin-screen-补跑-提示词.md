# 提示词：补跑 66 个端基的 spin screen

> **用法**：把下面「提示词正文」整段复制给一个新的 AI 会话。
> 目标仓库：`/nfs_hdd/2025/gaominliang/ht-vasp`。
> 本文档最后一节「事实附录」是给**你**（派发人）核对用的，不要复制进去。

> **关于本文档的路径**：`docs/prompts/` 现在在 `master` 上（PR #3 已合入），你从 master 开的分支会继承它。
> **你手上粘贴的正文就是权威版本**；如果正文里引用了某个 `docs/prompts/*.md` 而你找不到，
> 停下来向派发人索取，不要猜内容。

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

当前在 `feat-dfpt-qha`。**这个分支可以用，不要切走**——`main_spin_screen.py` 和 `htvasp/dfpt_qha/` 都在这里。

（**2026-10-28 更新**：PR #3 合入后 `master` **也有**这两样东西了，从 `master` 开分支同样可行。仍然成立的是：**不要在 `origin/feat-tb2j` 上工作**——那是 master 的落后子集，没有 spin screen 的代码。）

```bash
git branch --show-current          # 应该是 feat-dfpt-qha
./.conda/bin/python main_spin_screen.py --dry-run --batch BCC-Al-Al-2
```

### 2. 现状：84 个端基里只筛过 18 个

`pool_recipes.json` 的真实路径是 **`/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_recipes.json`**（`main_spin_screen.py:48` 的 `POOL_DIR`；不在仓库里，`data/endpoints/` 下没有）。它的 `magnetic_policy` 分布是 15 个 `confirmed_branch` + 69 个 `nonmagnetic`，看起来像全筛过了。**不是。**

- 真正被 `spin screen 2026-10-20` 测过的只有 **18 个**，全部含 Cr。
- 剩下 **66 个**的 `nonmagnetic` 里，真正的 `magnetic_note` 分三类（注意 `magnetic_note` 是 `pool_recipes.json` 的**顶层**键，不在 `recipe` 里）：
  - **57 个**原文是 `无 Cr；记录为 nonmagnetic 但未做自旋对照 —— 这是默认预期，不是验证过的结论`；
  - **6 个** `magnetic_note` 是 `null`：`SER-Al-1 SER-Hf-1 SER-Nb-1 SER-Ta-1 SER-Ti-1 SER-Zr-1`；
  - **3 个**含 Cr 的引的是报告（见下表）。
- 那 18 个的 `spin-screen.json` **在本机（429e/429pro 共享的 `/nfs_hdd`）一个都没有**，但在 **B5070 的 `/workspace/gaominliang/ht-vasp/data/endpoints/*/spin-screen.json`** 还完整保留着（18 个），`spin-screenflow/` 也 18 个。**这是重要的历史证据，先取回来。**

### 3. 你要做的事

**第一步：把 18 个历史 spin-screen.json 从 B5070 取回来。**

B5070 是另一台机器、另一个用户、只有密码登录，**不能 `ssh B5070`**。用 `/tmp/b5070.py` 这个 paramiko 辅助脚本（如果它还在；不在就照 `ht-vasp/scripts/nas/nas.py` 的写法自己写一个，凭据在 `/nfs_hdd/2025/gaominliang/.env` 的 `SERVER_B5070_*`）：

```bash
# /tmp/b5070.py 存在（1207 B，2026-10-26 00:11）。远端 cwd 是 /home/mcmf507，
# 不是仓库，所以**必须写绝对路径**：不带路径的 `ls data/endpoints/*` 返回 0 个。
./.conda/bin/python /tmp/b5070.py 'ls -d /workspace/gaominliang/ht-vasp/data/endpoints/*/spin-screen.json | wc -l'
# -> 18
```

取回到 `data/endpoints/<name>/spin-screen.json`（**不要覆盖任何已有文件**）。

本地现状（实测）：18 个里**只有 9 个有 `data/endpoints/<name>/` 目录**（含 `magnetic-state.json` 和/或 `dfpt-qhaflow/`），另外 9 个连目录都没有：

```
缺目录的 9 个：FCC-Cr-Hf-4 FCC-Cr-Ta-4 FCC-Cr-Ti-4 FCC-Cr-Zr-4 FCC-Hf-Cr-4
                FCC-Nb-Cr-4 FCC-Ta-Cr-4 FCC-Ti-Cr-4 FCC-Zr-Cr-4
```

取回时要 `mkdir -p`。同时检查 B5070 上有没有 `spin-screenflow/` 目录——有就一起取回（那是原始 TaskDoc 和 store，`--reanalyse` 要用；本地一个都没有）。

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

这三个的"磁性"结论来自报告，不是来自 spin screen。`FCC-Al-Cr-4` 尤其要重筛，而且**原因比"体积不可比"更硬**（2026-10-26 复核过归档 OUTCAR）：

- INCAR `MAGMOM = 1*1.0 3*3.0`、ISPIN=2，OUTCAR 磁矩从 10.0 起，4 个离子步、45 个电子步内塌到 ~0：`10.0 → 4.113 → 3.140 → 0.096`（第 8 步）→ 末尾 3e-7，且 `reached required accuracy`。
- 弛豫后的胞 `a=3.6768 Å → 12.4266 Å³/atom`；EOS json `V0 = 49.72329 Å³ → 12.4308 Å³/atom`。**报告那个 13.0301 Å³/atom 是 Stage 0D 的模板胞（a=3.73539），DFPT run 根本没跑过那个体积。**
- 更要命的是：报告里 22.96 meV/atom 的"磁性分支"其实是**两个非磁态之间的能量差**——cr0（非磁，a=3.73539，13.0301 Å³/atom，E=−31.52065342）与 cr3 起步弛豫终态（**同样非磁**，a=3.66733，12.3308 Å³/atom，E=−31.61250462）差 22.9628 meV/atom。在**匹配体积**下（a=3.66733 的 EOS 点，MAGMOM 1.0 3.0 3.0 3.0）也是非磁（净矩 0.00107）。也就是说这 22.96 meV/atom **主要是体积效应，不是磁性效应**。
- 所以"`FCC-Al-Cr-4` 在弛豫后的体积上还是非磁吗"（第 5 节问题 3）这一问，答案大概率是"是"；**真正要问的是：在池子自己的弛豫体积上，fm/afm 到底有没有比 nm 低 5 meV/atom。** 这就是这次重筛的意义。

**第三步：跑。**

工具是 `main_spin_screen.py`，它已经能做全部三件事（nm / fm / afm），并且有 `--slurm`。

**机器不是你的选项：必须 `--nodelist=429e`。禁止投 429pro。**（派发人 2026-10-26 的硬要求。）
`ntasks` / `memory` / `time_limit` / 并发仍然由你定，但要基于下面的事实而不是猜。

**节点的关键事实（2026-10-26 21:52 实测）：**

| 节点 | 核数 | 内存 | 备注 |
|---|---:|---:|---|
| 429e | 48 | **192 GB** | **本任务唯一允许的节点**；跑过 NELECT 43 的 DFPT 没 OOM |
| 429pro | 64 | **128 GB** | **不许用**；32 核作业并发 2 个时发生过 21 次内核 OOM 杀 vasp_std |

- 两个节点都在 `partCPU`，**没有独立的分区可选**。`--nodelist 429e` 是投到 `partCPU` 再钉节点。
- **429e 当前状态（2026-10-28 13:11）**：`57507 FCC-Ti-Ta-4-dfpt-qha` 已于 **10-27 04:00 结束**（7/7 volume 全部验证通过，但 7 个体积都非简谐 → 该端基被判 `finite_q_soft_mode`，不进生产 TDB）。现在 429e 上只有两个 4 核小作业，`CPUAlloc = 8/48`，**空闲 40 核**，`FreeMem` ≈ 128 GB。**`ntasks=32` 可以直接跑，不用排队。**
- 只要钉在 429e，作业 pending 就是预期行为，**不要因为 pending 就换节点**。现在 429e 有 40 核空闲，`ntasks<=32` 的作业应该立刻起跑；如果被 pending，先看 `scontrol show job <id>` 的 `Reason`，再决定缩 `ntasks`。
- `SlurmJobManager.get_cpu_config()` 支持 `nodelist`、`memory`、`time_limit`、`ntasks` 全部作为关键字参数；**`submit_screen()` 目前没有把 `nodelist` 传下去**（它写死 `memory="96G"`、`time_limit="24:00:00"`）。**你必须改这个函数把 `nodelist="429e"` 传下去**，否则作业会落到 429pro。提交后立刻用 `scontrol show job <jobid> | grep ReqNodeList` 核验真的是 `429e`。

**成本（实测，不是外推）：**

- 一次 spin screen = 每个态一次 `ISIF=3` 短弛豫（`DfptRelaxMaker`，`NSW` 上限 80）。
- **历史实测**：2026-10-19 那次 18 个端基 / 44 个自旋态，从 14:14 跑到 20:09，**一个 Slurm 作业、wall clock 5.93 h**（B5070 上 18 个 `spin-screenflow/store.json` 的 `completed_at` 时间戳；单态中位 126 s、最大 3578 s）。
- **`run_screen()` 是顺序 for 循环**（`for index, name in enumerate(names, 1)` → `run_locally_custom`），一个作业内**端基之间没有并行**。所以"32 核作业并行 3–4 个"不成立。
- 133 个态按上面的量级外推：**单作业 wall clock 约 18 h**。想更快只能自己拆成多个 `sbatch`（每个作业一批端基），但 429e 现在只剩 8 核，拆也拆不出多少并发。

**跑之前必须做的检查：**

1. `squeue -p partCPU` 与 `scontrol show node 429e` 看清 429e 还剩多少核。
2. **`--dry-run` 不要和 `--slurm` 一起用。** `submit_screen(dry_run=True)` 在 `build_flow()` 之前就 `return`，只打印一条命令 JSON，**根本不显示 job 名**——看起来"通过"其实什么都没验。要验 flow 结构必须用**不带 `--slurm`** 的 dry-run：
   ```bash
   ./.conda/bin/python main_spin_screen.py --dry-run --batch BCC-Al-Al-2
   # -> states: spin nm, spin fm
   ./.conda/bin/python main_spin_screen.py --dry-run --batch FCC-Al-Cr-4
   # -> states: spin nm, spin fm, spin afm
   ```
3. 真正落盘的位置是 `data/endpoints/<name>/spin-screenflow/`（在 `/nfs_hdd`，还有 5.6 T 空闲），**不是** `/nfs_ssd/tmp`：`main_spin_screen.py:45` 的 `flow_dir = /nfs_ssd/tmp` 是死代码，全文件只出现这一次，没有任何地方消费它。`df -h /nfs_hdd` 才是要看的。

**第四步：验收。** 每个端基要产出 `data/endpoints/<name>/spin-screen.json`，里面 `decision` ∈ {`confirmed_branch`, `nonmagnetic`}，`magnetic_gain_meV_per_atom` 是数值，`excluded_states` 列出来。**失败的端基要单独列出来**（`ok: false` 带 `error`），不要当成"跑过了"。

**第五步：不要动 `pool_recipes.json`。** 脚本自己的 docstring 写着："Updating `pool_recipes.json` from that file is a separate, explicit step." 把 66 个的新决策汇总成一张表放在报告里，**由人决定要不要写回**。

### 4. 硬约束

- **不改** `main_spin_screen.py` 的判定逻辑（`DECISION_THRESHOLD_MEV_PER_ATOM = 5.0`、`STRUCTURE_TOLERANCE = 0.02`）。如果你认为阈值要改，写进报告，不要动。
- **不改** `MAGMOM` 表、不改 `spin_states()` 的初猜值。
- **不改** 任何已归档的 `dfpt-qhaflow/` 结果。
- **不提交** 任何 DFPT-QHA / SQS / TB2J 作业——本任务只做 spin screen。
- **不碰** `em-tdb-v2`。
- 提交信息用 `feat(spin-screen):` 或 `chore(spin-screen):`。

### 4.1 你会撞上的一个真 bug（已复核，请按此处理）

`main_spin_screen.py:375` 和 `:441` 都写着：

```python
output = store.get_output(uuid=flow.output.uuid, which=flow.output.index, load=True)
```

`flow.output` 是 `OutputReference`，它的 `.index` **又是一个 `OutputReference`**（`OutputReference(uuid, .index)`），不是 int。`JobStore.get_output` 只有在 `isinstance(which, int) or which in {"last","first"}` 时才走单条分支；一个 OutputReference 落到 else 分支 → **返回 list**。

（2026-10-26 用一个真 `JSONStore` 实测确认：同一 uuid 存 2 条索引时，`which="last"` → `{'a': 2}`，`which=OutputReference(uuid).index` → `[{'a': 1}, {'a': 2}]`。）

- `run_screen()` 里 `if isinstance(output, list): output = output[0] if len(output)==1 else output`，只对"恰好 1 条"有效；**同一个 flow 有 2 条索引时（`--rerun` 之后就是）会 `raise ValueError("the summariser returned a dict... expected")`，一个成功的弛豫被记成 `ok:false`。**
- `run_constrained()`（`:377-378`）也有同一个 `len==1` 的三元，但**没有**后面那个 `isinstance(output, dict)` 检查，紧接着 `:384` 就 `output.get('decision')` → list 上调用 `.get` 会 `AttributeError`（`--constrained-afm` 那条路会以 `ok:false` 记下来，不是崩掉整个 run）。

**这是本任务要跑的代码路径，你可以修，但必须：** (a) 改成 `which=flow.output.index` 之外的写法（例如 `which="last"`，或先从 `flow.output.attributes` 取 `("i", idx)`）；(b) 单独一个提交；(c) 加一条测试覆盖"同一 flow 有两条索引"的情形。**不要顺手重写判定逻辑。**

### 4.2 一个会静默骗过你的测试夹具

`tests/fixtures/pool/pool_recipes.json` 与 `pool_structures.json` 和 `/nfs_ssd/tmp/legacy-qha-pool/stage0/` 下的活文件**字节相同**（不是裁剪过的小夹具）。所以 `tests/test_spin_screen.py`（7 个用例，现在全过）读的是真池；池子一改，夹具就和活数据悄悄分叉。**跑测试前先 `md5sum` 对一下。**

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

`/tmp/spin_results.txt` 里有全部 18 个的 `decision | gain`。**这 18 个的 `spin-screen.json` 只在 B5070 上还有**（本机 `/nfs_hdd` 上一个都没有；`data/remote/B5070/` 那个本地镜像里也没有，只有 `dfpt-qhaflow/`）。

### 归档的 DFPT-QHA run 与 screen 判定的一致性（已核）

全库 `data/endpoints/*/magnetic-state.json` 共 **56 个**，其中非 `nonmagnetic` 的只有 **7 个**（不是 8）：

| 端基 | 归档 kind | 在 18 里？ | screen 判定 | 一致 |
|---|---|---|---|---|
| BCC-Cr-Cr-2 | antiferromagnetic ±1.134 | 是 | confirmed_branch | ✅ |
| BCC-Cr-Hf-2 | ferromagnetic | 是 | confirmed_branch | ✅ |
| BCC-Cr-Ti-2 | ferromagnetic | 是 | confirmed_branch | ✅ |
| BCC-Cr-Zr-2 | ferromagnetic | 是 | confirmed_branch | ✅ |
| SER-Cr-1 | antiferromagnetic ±1.134 | 是 | confirmed_branch | ✅ |
| BCC-Al-Cr-2 | ferromagnetic (Cr 2.535 μB) | **否**（在 66 里） | 报告 4.9.3 | 参考 |
| FCC-Cr-Al-4 | ferromagnetic (Cr 2.498 μB) | **否**（在 66 里） | 报告 8.0a(j) | 参考 |

另外 18 个已筛端基里有 8 个存了 `magnetic-state.json`，其中 BCC-Cr-Nb-2 / BCC-Cr-Ta-2 / FCC-Cr-Cr-4 三个是非磁，与 screen 的 `nonmagnetic` 一致；**8/8 一致**。

`FCC-Al-Cr-4` 的归档 `magnetic-state.json` 是 `nonmagnetic`（初猜 `1*1.0 3*3.0`，45 步内塌到 0），而报告说 22.96 meV/atom —— 这就是第 3 步要重查的那个。注意它**不在** 18 个已筛名单里，别把它算进"已核"。

### 66 个的分布

- 含 Cr 3 个：`BCC-Al-Cr-2`、`FCC-Al-Cr-4`、`FCC-Cr-Al-4`
- 不含 Cr 63 个
- 状态数：65 个跑 `nm,fm`，1 个跑 `nm,fm,afm`（`FCC-Al-Cr-4`，3 个 Cr，afm 向量 `[0.0,-3.0,3.0,-3.0]` —— 按位序交替，3 个 Cr 会得到净 −3 μB，不是补偿的）
- 合计 **133 个自旋态**

### `unresolved_endpoints()` 现在选不出任何端基

`pool_recipes.json` 里 `magnetic_policy == "unresolved"` 的数量是 **0**。所以 `main_spin_screen.py --batch`（不带名字）会说"nothing to do"。**必须显式传名字**：

```bash
./.conda/bin/python main_spin_screen.py --slurm \
    --batch $(awk '!/^#/{print $1}' docs/prompts/spin-screen-unscreened.txt)
```

`--batch` 是 `nargs="*"`（空格分隔的名字，不是文件路径）。

### 429e vs 429pro 的 OOM 证据

`docs/搬回429pro候选清单-20261023.md` §6.1（L127 `### 6.1 已确诊：OOM，不是算错`、L144 `累计 **21 次 OOM 杀掉 vasp_std**`）：429pro 上累计 21 次内核 OOM 杀 `vasp_std`，6 个 FCC 端基全部失败；`journalctl -k` 显示每次截断前 6–11 秒都有 OOM。同一批作业在 429e（192 GB）上 7/7 成功。

补充（2026-10-26 现场核过，**不是**那份文档的内容）：`scontrol show config | grep ConstrainRAMSpace` → `ConstrainRAMSpace = no`，`SelectType = select/cons_tres` / `CR_CORE`。所以 `--mem` 只是调度声明、**cgroup 不强制**；`scontrol show node` 的 `AllocMem=0` 也是同一个原因（内存不按作业记账，`FreeMem` 是节点实测值）。写 `--mem` 主要是让调度器别把两个大内存作业塞一起，真正的保护是**钉 429e**。
