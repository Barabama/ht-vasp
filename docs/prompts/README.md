# 待派发的提示词

四份提示词，各自独立、互不阻塞。按建议顺序。

| # | 文件 | 目标分支 | 能自己决定什么 | 必须先问的 |
|---|---|---|---|---|
| 1 | [TB2J-Wannier收敛检查-提示词.md](TB2J-Wannier收敛检查-提示词.md) | `feat-tb2j` | 检查函数怎么写、门放在哪、要不要 raise | 不改窗口参数、不重跑 |
| 2 | [TB2J-Wannier参数-调查链与待定方案.md](TB2J-Wannier参数-调查链与待定方案.md) | `feat-tb2j` | 只读复核 + 设计实验 | **改任何参数前必须书面批准**；每阶段结束停下来 |
| 3 | [spin-screen-补跑-提示词.md](spin-screen-补跑-提示词.md) | `feat-dfpt-qha` | `ntasks`/`memory`/`time_limit`/并发 | **节点固定 429e（`--nodelist=429e`），禁止 429pro**；不改判定阈值、不写回 `pool_recipes.json` |
| 4 | [ISIF8方案-科学审查-提示词.md](ISIF8方案-科学审查-提示词.md) | `feat-dfpt-qha` | 判定标签、替代方案、要不要反对 | **只读**；不改代码/数据，不跑 VASP，不 commit |

## 为什么是这三份，为什么这个顺序

**#1 和 #2 不冲突，但 #2 依赖 #1 的产物。** #1 只做检查（不改参数），所以它能立刻做、立刻合入；做完之后，"这个 run 有没有收敛"从"要人肉 grep wout"变成"JSON 里一个字段"。#2 要做的参数实验，判据就是 #1 加的那个字段。

**#3 与 #1/#2 完全独立**：它在 `feat-dfpt-qha` 上，不需要 TB2J 的代码；#1/#2 在 `feat-tb2j` 上，不需要 DFPT-QHA 的代码。两者可以并行派发。

**#4 是审查，不是执行。** 它审的是"ISIF=3 → ISIF=8"这份**还没执行**的论证，唯一交付物是判定表和替代方案。它和 #1/#2/#3 都不冲突，但**必须在改 `htvasp/sqs/maker.py` 之前跑完**——改了就变成事后追认。

## 分支事实（2026-10-26 实测）

```
master          ← feat-tb2j 是它的祖先（master 领先 16 个提交）
                ← feat-dfpt-qha 不是它的祖先（feat-dfpt-qha 领先 93 个）
feat-tb2j       ← 5 个 tb2j 文件；无 dfpt_qha/、无 main_spin_screen.py
feat-dfpt-qha   ← 有 dfpt_qha/、main_spin_screen.py、magnetic_contract.py、context.py
```

`feat-tb2j` 合并 `master` 是安全的（`git merge-tree` 报 0 冲突）。合并后 `htvasp/dfpt_qha/` 和 `main_spin_screen.py` **不在工作区**——这是拓扑事实，不是丢代码。**派发 #3 的 agent 不要切到 `feat-tb2j`；派发 #1/#2 的 agent 不要切到 `feat-dfpt-qha`。**

## 派发人 2026-10-26 追加的硬约束（已写进 #3）

**所有 Slurm 作业必须 `--nodelist=429e`，不许投 429pro。** 派发时的原话：

```
squeue 57507   partCPU  FCC-Ti-T  mcmf429  R  8:53:08  1  429e
你的这些任务排在这个 job 它后面，要在 nodelist 429e 提交，不要提交到 429pro
```

实测（2026-10-26 21:52）：

| 事实 | 值 |
|---|---|
| `57507 FCC-Ti-Ta-4-dfpt-qha` | 429e，`NumCPUs=32`，`TimeLimit=4-04:00`，11:29 起跑，日志已到 DFPT volume 6/7 |
| 429e | 48 核，`CPUAlloc=40` → **只剩 8 核**；192 GB，`FreeMem` ≈ 122 GB |
| 429pro | 64 核，`CPUAlloc=57` → 只剩 7 核；128 GB；57502 占 48 核，57503/57504 还排队 |
| `submit_screen()` | **不传 `nodelist`** → 不改就会落到 429pro |

代价要说清楚：`ntasks=32` 的作业钉在 429e 上会一直 pending 到 57507 释放（EndTime 2026-10-30T15:29，实际应在 volume 7/7 之后几小时内结束）；`ntasks<=8` 可以 backfill 立刻进去。

## 2026-10-26 复核：这批提示词改了哪些

八路独立复核（每条论断都从文件重算，不信文档里的数字）跑完之后，改了三份文档里**会误导执行 agent** 的地方。改的都是事实，没有改判定标准、没有改方案排序。

**#1 收敛检查**（`TB2J-Wannier收敛检查-提示词.md`）
- "35/44 会收敛（中位 1001）"补上：这是 **spin-1 单通道**外推，down 通道差很多（`SER-Mn` 682→1923、`BCC-Mn-Mn` 894→7159、`BCC-Ti-Ti` 5247→637）；下限是 **307**（不是 296）。
- "9/44 平台化"改成"看起来已平台化"：衰减率只比 1 高 0.1–3%，换窗口后**只有 7 个**在两个自旋上都稳。
- `wannier90.wout`（无编号）的大小**不是 15421**：44 个 run 有 15421/15581/15741/16061 四种，改按 `<-- DIS` 行识别。
- `tests/test_tb2j.py` 的"12 个全过"是**假绿**：全文件 `assert` 出现 0 次，`return False` 也算 passed（`test_structure_names` 期望 30、实际 44，照样过）。新增用例必须用 `assert`。
- `feat-dfpt-qha` 领先 master 是 **93** 个提交（原文 89 已过期）。
- wannier90 源码补了 Gamma-only 重复实现（`:2791/:2826-2830/:3075`）、`.lt.` 严格比较、`dis_conv_window` 门（`:2297`）。

**#2 参数调查**（`TB2J-Wannier参数-调查链与待定方案.md`）
- **方案 D 有一个致命前提没写**：`ht-vasp` 预写 `wannier90.win` 且 INCAR 无 `WANNIER90_WIN` → `mlwf.F:308-311` 打开 `WANNIER90_WIN_LEGACY_MODE` → `:596`/`:1046` 两处 win 写入**都被跳过**。**现在直接加 `LWANNIER90_AUTO_WINDOW` 是空操作**，必须先改写入方式。这一条把方案 D 从"加一行 INCAR"抬成"要改 input_set 的落盘路径"。
- `FILTER_WIN_FILE` 过滤的是 INCAR 的 `WANNIER90_WIN` 字符串，**不是磁盘上的 win**（调用点 `:1787`）。
- 行号校正：`dis_froz_max` 写入在 `:1821-1823`（不是 `:1818-1822`）；`jobs.py` 的行号**分支相关**（feat-tb2j `:111/:263`，feat-dfpt-qha `:251/:486`），改成按函数名 grep。
- `SER-Mn` 的 `num_wann` 是 **9**，不是 18（6.70 来自 band 9；band 18 是 42.94 eV）。
- §2.11 的 k 点对照不是严格单变量：`--ntasks` 也不同（4 vs 16）；而且 **8×8×8 的 down 通道第 181 步就收敛了**，只有 up 撞上限。结论方向不变，但证据要按这个说。
- `num_wann = len(structure)×9`（`input_set.py:102-122`），**不依赖 `ntasks`**；只有 `num_bands` 依赖。单一 `ntasks=32` 就能同时给出 9→32/18→64/36→64。
- §2.11 那两个相关系数（+0.056 / +0.312）**没被独立复现**，已标注"不要引用"。

**#3 spin screen**（`spin-screen-补跑-提示词.md`）
- 派发人的硬要求写进去了：**必须 `--nodelist=429e`，禁止 429pro**，并给出实测的节点状态与 `57507` 的 EndTime。
- 成本改写：`run_screen()` 是**顺序 for 循环**，"32 核并行 3–4 个"不成立；历史实测（2026-10-19，18 端基 / 44 态）是**单作业 wall clock 5.93 h**，133 个态外推 **≈18 h**。
- `pool_recipes.json` 的真实路径是 `/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_recipes.json`（原文没给路径）。
- `magnetic_note` 是**顶层**键：57 个是那句中文默认预期、**6 个是 `null`**（SER-Al/Hf/Nb/Ta/Ti/Zr-1）、3 个含 Cr 的引报告。
- B5070 helper 的远端 cwd 是 `/home/mcmf507`，**必须用绝对路径** `/workspace/gaominliang/ht-vasp/...`。
- `--dry-run --slurm` **不打印 job 名**（在 `build_flow()` 之前就 return）；要验 flow 结构必须用不带 `--slurm` 的 dry-run。
- `flow_dir = /nfs_ssd/tmp` 是**死代码**；实际落盘在 `/nfs_hdd`（5.6 T 空闲），`df -h /nfs_ssd` 查错了盘。
- `main_spin_screen.py:375/:441` 的 `which=flow.output.index` 是真 bug：`OutputReference.index` 又返回 `OutputReference`，`get_output` 落到 `"all"` 分支返回 list；`--rerun` 后会把成功记成 `ok:false`。
- 归档磁性端基是 **7 个不是 8 个**；`FCC-Al-Cr-4` 的 22.96 meV/atom 其实是**两个非磁态之间的体积差**，不是磁性分支。
- `tests/fixtures/pool/*.json` 与活池**字节相同**，不是裁剪夹具。

## 派发人手上还留着的决策

1. **SQS 的 ISIF 方案**——已从"ISIF=2 + 体积来源"推进到"ISIF=8 + Vegard 起始体积"，论证见 [#4](ISIF8方案-科学审查-提示词.md) §1。**仍未定**，等 #4 的审查结论。
2. **Wannier 参数最终怎么改**（方案 A–F）——留给 #2 的 agent 调查后再定，但**批准权在你**。
3. **spin screen 的 66 个结果要不要写回 `pool_recipes.json`**——#3 只出表，不写回。

## 我在写这批提示词时被数据推翻的一个说法

我最初判断"44/44 没收敛是因为 `dis_froz` 被写成了绝对值"。复核时发现：**44 个里只有 3 个的 E_F 高于 10 eV**（`BCC-Al-Nb`、`BCC-Nb-Nb`、`SER-Nb`）。其余 41 个的窗口本来就跨过 E_F，照样不收敛。

而 `tb2j_test` 里有一组单变量对照把嫌疑指向 **W90 SCF 的 k 点数**：同一个 Fe bcc、同一个窗口、同样 9/16，**64 kpts 收敛、512 kpts 不收敛**；44 个归档 run 全部是 512 kpts。

所以 #2 的提示词里我加了 §2.11（反证）和方案 E（降 k 点数），并要求 agent **先判断主因再谈改什么**。这条如果你有别的看法，直接改那份文档。

## 数据文件

- [spin-screen-unscreened.txt](spin-screen-unscreened.txt) —— 66 个未筛端基 + 每个该跑哪几个自旋态（133 个态）
