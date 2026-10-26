# 待派发的提示词

三份提示词，各自独立、互不阻塞。按建议顺序。

| # | 文件 | 目标分支 | 能自己决定什么 | 必须先问的 |
|---|---|---|---|---|
| 1 | [TB2J-Wannier收敛检查-提示词.md](TB2J-Wannier收敛检查-提示词.md) | `feat-tb2j` | 检查函数怎么写、门放在哪、要不要 raise | 不改窗口参数、不重跑 |
| 2 | [TB2J-Wannier参数-调查链与待定方案.md](TB2J-Wannier参数-调查链与待定方案.md) | `feat-tb2j` | 只读复核 + 设计实验 | **改任何参数前必须书面批准**；每阶段结束停下来 |
| 3 | [spin-screen-补跑-提示词.md](spin-screen-补跑-提示词.md) | `feat-dfpt-qha` | 提交机器（429e/429pro）、Slurm 参数、并发 | 不改判定阈值、不写回 `pool_recipes.json` |

## 为什么是这三份，为什么这个顺序

**#1 和 #2 不冲突，但 #2 依赖 #1 的产物。** #1 只做检查（不改参数），所以它能立刻做、立刻合入；做完之后，"这个 run 有没有收敛"从"要人肉 grep wout"变成"JSON 里一个字段"。#2 要做的参数实验，判据就是 #1 加的那个字段。

**#3 与 #1/#2 完全独立**：它在 `feat-dfpt-qha` 上，不需要 TB2J 的代码；#1/#2 在 `feat-tb2j` 上，不需要 DFPT-QHA 的代码。两者可以并行派发。

## 分支事实（2026-10-26 实测）

```
master          ← feat-tb2j 是它的祖先（master 领先 16 个提交）
                ← feat-dfpt-qha 不是它的祖先（feat-dfpt-qha 领先 89 个）
feat-tb2j       ← 5 个 tb2j 文件；无 dfpt_qha/、无 main_spin_screen.py
feat-dfpt-qha   ← 有 dfpt_qha/、main_spin_screen.py、magnetic_contract.py、context.py
```

`feat-tb2j` 合并 `master` 是安全的（`git merge-tree` 报 0 冲突）。合并后 `htvasp/dfpt_qha/` 和 `main_spin_screen.py` **不在工作区**——这是拓扑事实，不是丢代码。**派发 #3 的 agent 不要切到 `feat-tb2j`；派发 #1/#2 的 agent 不要切到 `feat-dfpt-qha`。**

## 派发人手上还留着的决策

1. **SQS 的 ISIF 方案**（B 方案：ISIF=2 + 体积来源）——已口头讨论，尚未定。**不在这一批提示词里。**
2. **Wannier 参数最终怎么改**（方案 A–F）——留给 #2 的 agent 调查后再定，但**批准权在你**。
3. **spin screen 的 66 个结果要不要写回 `pool_recipes.json`**——#3 只出表，不写回。

## 我在写这批提示词时被数据推翻的一个说法

我最初判断"44/44 没收敛是因为 `dis_froz` 被写成了绝对值"。复核时发现：**44 个里只有 3 个的 E_F 高于 10 eV**（`BCC-Al-Nb`、`BCC-Nb-Nb`、`SER-Nb`）。其余 41 个的窗口本来就跨过 E_F，照样不收敛。

而 `tb2j_test` 里有一组单变量对照把嫌疑指向 **W90 SCF 的 k 点数**：同一个 Fe bcc、同一个窗口、同样 9/16，**64 kpts 收敛、512 kpts 不收敛**；44 个归档 run 全部是 512 kpts。

所以 #2 的提示词里我加了 §2.11（反证）和方案 E（降 k 点数），并要求 agent **先判断主因再谈改什么**。这条如果你有别的看法，直接改那份文档。

## 数据文件

- [spin-screen-unscreened.txt](spin-screen-unscreened.txt) —— 66 个未筛端基 + 每个该跑哪几个自旋态（133 个态）
