# 提示词：审查「SQS 弛豫改用 ISIF=8」这一方案的科学论证

> **用法**：把下面「提示词正文」整段复制给一个新的 AI 会话。
> 本文档最后一节「派发人自查记录」是给**你**（派发人）核对用的，不要复制进去。

> **关于本文档的路径**：这些提示词在 `feat-dfpt-qha` 分支的 `docs/prompts/` 下。
> 如果你切到了别的分支，那个目录可能不存在。**你手上粘贴的正文就是权威版本**；
> 正文引用的数据文件如果找不到，停下来报告，不要猜内容、不要用别的数据顶替。

---

## 提示词正文

你是一名**对抗性科学审查员**。任务不是复述或赞同，而是**找出下面这份论证里的错误**。
结论是"我错了"和"我对了"同样有价值；但"看起来合理"不算通过——每一条都要落到
数据、源码或文献上。

审查对象：一份主张「SQS 混合焓的 relax 阶段应该从 `ISIF=3` 改成 `ISIF=8`」的
内部论证。它**尚未执行**，只存在于讨论中。你的意见会直接决定是否改代码。

### 0. 先读，不要先判断

按顺序读：

1. `/nfs_hdd/2025/gaominliang/AGENTS.md`
2. `/nfs_hdd/2025/gaominliang/ht-vasp/CLAUDE.md`
3. `/nfs_hdd/2025/gaominliang/ht-vasp/docs/TB2J与SQS数据审计-20261025.md` §2（SQS 部分，L186 起）
4. `/nfs_hdd/2025/gaominliang/ht-vasp/docs/SQS-mixing固化进度.md` §3.2（L140–160）与 §6.11（L420 起）
5. `/nfs_hdd/2025/gaominliang/ht-vasp/docs/legacy-qha-rebuild-report.md` §4.12.5（L1380 起，模板体积误差）、
   §8.0a(l)（L3283 起，形状弛豫判据写反的那次事故）
6. 然后 `git status`。**不要 reset、clean、stash，不要覆盖或删除任何文件。**

**VASP 源码在本机**：`/nfs_ssd/softwares/vasp/vasp.6.4.3/src/`。
凡涉及 ISIF 语义，以 `main.F` 为准，不要凭记忆或 wiki。

### 1. 被审查的论证（原文摘要）

> **主张**：SQS 混合焓的 relax 应改 `ISIF=3 → ISIF=8`，并把起始胞放到 Vegard 体积。
> **理由链**：
> 1. `ISIF=8` = 离子自由 + 体积自由 + **形状锁死**（源码依据见下）。
> 2. 现有 `ISIF=3` 的 18 个成分点里 **8 个** 弛豫成 `a/b > 1.05`，最大 **1.485**；
>    50/50 无序 bcc 固溶体没有理由这样。
> 3. 同一体系同一成分的两次独立实现（BCC-Al-Nb x=0.5）**体积一致到 0.01%、
>    形状差 26%、H_mix 差 32 meV/atom**，所以形状项是"实现间随机量"，
>    往 H_mix 里注入 32 meV/atom 量级的伪散差。
> 4. 退到 `ISIF=2` 会把体积钉死，代价是 `½·B0·V·(ΔV/V)²`：
>    Al-Nb 2.0/6.8/14.4、Al-Ti 6.2/11.4/14.7 meV/atom（x=0.25/0.50/0.75），
>    单调、只加不减 → 被 L0/L1 吸收成人造正过剩项。
> 5. `ISIF=7` 冻结离子，丢掉 11.7/56.5/70.0 meV/atom 的离子弛豫能，比 2 和 8 都差。
> 6. CALPHAD 的 L 参数是 p=0 下的过剩 Gibbs 能，要求每个成分在**自己的**平衡体积上；
>    文献推荐 `ISIF=2` 的理由（Nguyen et al., PRB 96, 014107 (2017)）是 CE 训练集
>    需要共享参考晶格，与混合焓的目标相反，不能直接搬。
> 7. 另有独立缺陷：relax 用 `ENCUT=400`、static 用 `ENCUT=520`，两个能量面；
>    x=0.5 的 static 残余力 0.5851 eV/Å，估约 17 meV/atom。只改 ISIF 不改 ENCUT，
>    这一项仍在。

**你要审查的就是这 7 条，以及它们的组合结论。**

### 2. 你必须自己验证的事实（不要采信上面的数字）

#### 2.1 ISIF=8 的语义（读源码，不要读 wiki）

```bash
sed -n '4096,4115p' /nfs_ssd/softwares/vasp/vasp.6.4.3/src/main.F   # ISIF 表
sed -n '4155,4200p' /nfs_ssd/softwares/vasp/vasp.6.4.3/src/main.F   # FACT / FACTSI / D2SIF
cat /nfs_ssd/softwares/vasp/vasp.6.4.3/src/constr_cell_relax.F      # 是不是空实现
```

要确认的：
- `ISIF=8` 下 `D2SIF` 是否**真的只剩对角**（因此剪切恒为 0，形状不可能变）。
- `ISIF=8` 下 `FACT` 是否非零（因此离子真的会动）。注意 `main.F:4160` 的条件是
  `ISIF<5 .OR. ISIF==8`——**`ISIF=7` 不满足，所以 7 确实冻离子；`ISIF=6` 也不满足**。
- `ISIF=8` 的**收敛判据**（`main.F:4201-4207`，`LSTOP2` 用
  `ABS(D2SIF(J,I))/FACTSI/NIONS > ABS(EDIFFG)`）。这个判据在 `ISIF=8` 下测的是什么量、
  量纲是什么、`EDIFFG=-0.05` 对它意味着多严。**这一条原文没审，是已知空白。**
- `ISIF=8` 与 `ISIF=2` 的差别是否**仅**为体积。原文这么断言，请独立判断
  （注意 `main.F:2261` 的 XDATCAR header 分支、`main.F:4478`/`4731`/`4753` 的
  `ISIF>=3` 分支——`ISIF=8` 会触发 `ISIF>=3` 的那些重算路径，`ISIF=2` 不会。
  这对结果有没有影响？）

#### 2.2 形状统计（18 个点）

```bash
cd /nfs_hdd/2025/gaominliang/ht-vasp
# 六个体系 × 三点：data/endmembers/BCC-{Al-Nb,Al-Ti,Al-V,Nb-Ti,Nb-V,Ti-V}/*-sqs.json
# 每个 points[i].final_structure.lattice.matrix 求三边长与夹角
```

原文的表（最长/最短边）需要你重算。注意：原文**同时**声称"角度偏离 90° 全部 < 1.5°，
所以是近四方不是真三斜"——请核这个说法，因为它与
`docs/TB2J与SQS数据审计-20261025.md` §2.3 的 "P1 / 三斜" 措辞冲突。
**两份文档里必有一份是错的，请指出是哪一份。**

#### 2.3 两次实现的可比性（这条最可能推翻原文）

```bash
cd /nfs_hdd/2025/gaominliang/ht-vasp
# 08-21 实现：data/sqs_mixing/BCC-Al-Nb/BCC-Al-Nb-sqs-mixing.json
# 09-04 实现：data/endmembers/BCC-Al-Nb/BCC-Al-Nb-sqs.json
# 原始目录：data/sqs_mixing/BCC-Al-Nb/runs/x0.5/{1-sqs_mixing_relax,2-sqs_mixing_static}
#           data/endmembers/BCC-Al-Nb/sqsflow/runs/point/BCC-Al-Nb-x0.5/{1-sqs_relax,2-sqs_static}
```

**原文只核到"体积相同、形状不同、能量差 32 meV/atom"，没有核两者的 SQS 构型是否同一
个。** 派发人后来补测的结果是：`x=0.5` 的两套构型**在立方点群下不等价**（`x=0.25`、
`x=0.75` 等价）。这直接动摇原文第 3 条——如果两个点连构型都不同，那 32 meV/atom
里有多少来自形状、多少来自构型，原文没有分开。

请独立重算这三件事，并给出你自己的分解：
1. 三个成分点上，两套构型在立方点群（48 个操作）+ bcc 平移下是否等价。
2. 若不等价，`x=0.5` 的 32 meV/atom 有多少是"构型不同"（可以用两个纯元素的
   同构型对照、或用两套各自的 x=0.25/0.75 差来定标）。
3. **在构型不同的前提下，"形状是实现间随机量"这个论断还成立吗？**
   换一种问法：一个 16 原子、特定构型的 SQS，在 ISIF=3 下弛豫成 a/b=1.371，
   是"数值噪声"，还是"这个构型真实的对称性破缺"？
   **这是本次审查最重要的一问。** 若答案是后者，则 ISIF=8 的正当理由必须改写
   （不是"消噪声"，而是"我们不要构型特定的对称性破缺，我们要随机合金的平均能量"），
   而且必须承认：锁死形状会**丢掉**一部分真实能量，只是丢掉的这部分不该进 H_mix。

#### 2.4 弹性惩罚的算式与参数

原文用 `ΔE = ½·B0·V·(ΔV/V)²`，`B0` 取自 DFPT-QHA 的**有序端基**
（`BCC-Al-Nb-2` 141.2 GPa、`BCC-Al-Ti-2` 110.9 GPa，
见 `data/endpoints/*/dfpt-qhaflow/runs/*/*-dfpt-qha.json` 的 `eos.B0_GPa`）。

要审的：
- 这个式子对**各向同性**压缩是对的，但钉死体积后离子仍会弛豫、且**形状可能变**
  （`ISIF=2` 固定形状；`ISIF=8` 放开体积锁形状）。两种钉法的弹性项是否同式？
- 用**有序**端基的 `B0` 代替**无序**合金的 `B0`，误差多大？方向如何？
- `Vegard` 体积取自同一次 run 自己的两个纯元素参考。**请先核这些参考是否自洽**：
  `Al` 在六个体系里都是 `17.000 Å³/atom`，但 `Nb` 在 `Al-Nb` 里是 `18.066`、
  在 `Nb-Ti`/`Nb-V` 里是 `17.990`（差 0.42%）。同一个纯元素参考在不同体系里不一致，
  这会直接进 Vegard。**这条原文完全没提，请判定它是否足以改变结论。**

#### 2.5 ENCUT 不一致的代价

原文的 17 meV/atom 是用 `k ≈ 10 eV/Å²` 手估的。请判断这个估算是否有依据，
以及残余力 0.5851 eV/Å 到能量的换算该怎么做才站得住。
数据：`data/sqs_mixing/BCC-Al-Nb/runs/x0.5/2-sqs_mixing_static/` 的 OUTCAR，
以及同目录 `1-sqs_mixing_relax/`。

#### 2.6 一个原文没提、但会影响所有结论的事实

已归档的 08-21 run 与**当前提交的代码**产出的 INCAR **不一样**：

| | 归档 08-21 x=0.5 relax | 当前提交代码会产出 |
|---|---|---|
| `ALGO` | `Normal` | **`Fast`** |
| `MAGMOM` | `8*1.0 8*3.0` | **`16*0.6`**（pymatgen 默认） |

（`main_sqs.py` 的 `DEFAULT_GLOBAL_INCAR` 有 `"ALGO": "Fast"`，relax 阶段不覆盖它；
`MAGMOM` 键根本不存在。归档值来自 `cf93be6` 时代的另一份配置。）

请判断：**当前代码能否复现已归档的结果？** 如果不能，那么所有基于归档数据的
ISIF 论证都只在"旧配方"上成立，改 ISIF 时必须同时说明配方是否一起变。
另外 `ALGO=Fast`（RMM-DIIS）在 16 原子金属胞上是否可靠——仓库里有一条记忆说
`BCC-Ti-V` 的 DFPT 在 16 site 胞上 `Fast` 不收敛、正式配方是 `ALGO=Normal`。
请去核实那条说法的出处，并判断它对 SQS 的 relax 是否适用。

### 3. 你要给出的判定

对**上面第 1 节的 7 条**逐条给一个标签：

- `CONFIRMED` —— 独立复算通过，且措辞没有夸大
- `OVERSTATED` —— 事实对，但结论/措辞超出证据（请给出你认为站得住的措辞）
- `WRONG` —— 事实错（请给出正确值和你算的方式）
- `UNVERIFIABLE` —— 现有数据不足以判断（请说明缺什么）

然后单独回答：

1. **第 2.3 节第 3 问**（构型特定对称性破缺 vs 数值噪声）——这是决定性的一问。
2. `ISIF=8` 会不会引入**新的**偏差？特别是：形状锁死在一个非平衡的立方胞上，
   对**有序化倾向强**的体系（Al-Nb、Al-Ti 的 H_mix 到 −427 meV/atom）
   是不是把"该畸变的能量"人为抬高？有没有更合适的选择（例如 `ISIF=3` 但加
   `ISYM`/对称性约束，或对同一个 SQS 做多种子取平均）？
3. 如果目标是"每个成分在自己平衡体积上的 H_mix"，`ISIF=8` 是不是**唯一**或**最好**的
   做法？请至少给出一个替代方案并比较。
4. 六个体系里 `BCC-Nb-Ti`/`Nb-V`/`Ti-V` 的 `a/b` 都在 1.001–1.035，
   `Al-Nb`/`Al-Ti`/`Al-V` 才畸变。**畸变与 Al 含量强相关**——
   这是化学效应（Al 的原子半径/电子结构）还是体积起点太远（模板 a=2.8597 Å，
   bcc Al 是 3.2396、bcc Nb 是 3.3059）？原文把起点问题放在附注里，
   请判断它是不是**真正的主因**。

### 4. 硬约束

- **只读**。不改任何代码、不改任何已归档数据、不跑 VASP、不提交 Slurm、不碰 `em-tdb-v2`。
- 不改 git 状态（不 commit、不 branch、不 stash）。
- 结论要能落地成"改 / 不改 / 改成什么"，不要停在"需要更多研究"。
- 你算出来的每个数都要给出**可复现的路径或一行命令**，否则标 `UNVERIFIABLE`。

### 5. 交付格式

写一份 `docs/ISIF8方案-科学审查-<日期>.md`，结构：

1. **一句话结论**：这个方案该不该执行。
2. **逐条判定表**（7 条 × 标签 + 一句话理由）。
3. **决定性问题的回答**（第 3 节第 1 问）。
4. **你发现的、原文完全没想到的问题**（这一节最重要，允许为空，但要说明你找过什么）。
5. **如果你反对**：给出你认为正确的方案，以及它的代价（重跑量、CPU 秒）。
6. **如果你赞成**：列出执行前必须一起改的东西（含你判断会被测试挡住的那些）。
7. **你没做的事**。

---

## 派发人自查记录（不要复制）

### 我自己已知的薄弱点

- **第 3 条（32 meV/atom 是形状造成的）**：我只核到体积相同、形状不同、能量差 32。
  后来才发现两套构型在 x=0.5 **不等价**（x=0.25/0.75 等价）。所以"32 全归形状"是
  **未分解的**，这条大概率会被判 `OVERSTATED`。
- **"形状是噪声"这个措辞**：严格说，一个 16 原子特定构型的 SQS 在 ISIF=3 下畸变，
  是**这个构型真实的对称性破缺**，不是数值噪声。我把它叫"噪声"是为了表达
  "它不该进 CALPHAD 的 H_mix"，但这个措辞会误导。预期被判 `OVERSTATED`。
- **弹性惩罚的 B0**：用的是**有序**端基的 B0，不是无序合金的。未量化误差。
- **17 meV/atom**：`k ≈ 10 eV/Å²` 是手估，没有依据。
- **Nb 参考体积不一致**（18.066 vs 17.990）：我在正文里**没提**，是写这份提示词时
  才发现的。它会进 Vegard，进而进第 4 条的算式。这一条我自己没有结论。
- **ALGO/MAGMOM 与归档不符**：也是写提示词时才实测到的。当前提交代码产出
  `ALGO=Fast` + `MAGMOM=16*0.6`，归档是 `ALGO=Normal` + `8*1.0 8*3.0`。
  这意味着"当前代码 → 重跑"本身就不复现归档。这一条足以独立地要求一次重跑。
- **ISIF=8 的收敛判据**：`LSTOP2` 在 ISIF=8 下只看 `D2SIF` 对角（非对角被清零），
  所以判据实际是在测体积收敛。`EDIFFG=-0.05` 对这个量意味着多严，我没算。
  原文完全没审，是已知空白。

### 事实附录：ISIF 决策的溯源（已核）

`git blame` 三处引入点（`13fb500` 2026-08-19、`cf93be6` 2026-08-22、
`3832fc4` 2026-10-03）的 **commit body 全为空**；代码里唯一注释是
`htvasp/sqs/maker.py:53` 的 `# SQS-mix uses a full cell relaxation in this milestone.`
（陈述句，非理由）。`docs/SQS-mixing固化进度.md` L147/L428 也只是转述。

两个模板来源：
- `/nfs_ssd/softwares/atat/examples/vasp.wrap` 与 `glue/vasp/vasp.wrap`：
  `NSW=41 / IBRION = 2 / ISIF = 3 / ... / DOSTATIC`（`ezvasp:66` 同为 `ISIF = 3`）
- `htvasp/workflows/static.py:83`（`322e96b` 2026-03-25）的 `DoubleRelaxMaker`

`main_sqs.py:104` 的 `DEFAULT_GLOBAL_INCAR` 里写的是 `"ISIF": 2`，
被 `:199`/`:237` 和 `maker.py:54` 三处覆盖成 3 —— 同一文件里两个相反的值并存。

`tests/test_sqs_workflow.py:54-59`（`test_relax_maker_uses_full_cell_relaxation`）
**显式传 `relax_incar={"ISIF": 2}` 再断言结果 == 3**，所以任何改动都会先撞这个用例。

### 已知会变红的测试

- `tests/test_sqs_workflow.py::test_relax_maker_uses_full_cell_relaxation`（断言 == 3）
- 若改 `main_sqs.py` 的 `_calculation_settings()`，`settings_fingerprint` 会变，
  旧 `batch.json` 会被 fingerprint 检查拒（`main_sqs.py:298-306` 的逻辑）。
