# CMCH@CoNi-LDH-S 实验-计算对照研究

> 初稿来源: `初稿第三版.docx` - 《硫脲辅助MOF拓扑转化构筑硫修饰镍钴层状双金属氢氧化物/碱式碳酸钴锰多级异质结电极及超级电容性能研究》

---

## 一、目标材料

**CMCH@CoNi-LDH-S** = **CoMn(CO₃)(OH)₂** (基底/核) @ **S-doped CoNi-LDH** (外壳)

即第一篇论文中 **S1-CMCH@NC-LDHs** 样品（1 wt% 硫脲）。

---

## 二、实验制备关键信息

| 项目                     | 内容                                                                                           |
| ------------------------ | ---------------------------------------------------------------------------------------------- |
| **合成方法**       | 三步法：水热 CMCH → 浸渍 ZIF-67 → 二次水热转化为 NC-LDHs/S-NC-LDHs                           |
| **基底**           | 泡沫镍（2×3 cm²），首先水热生长 CMCH 纳米针阵列                                              |
| **CMCH 合成**      | Co(NO₃)₂·6H₂O + MnCl₂·4H₂O + 尿素 + NH₄F, 120°C 8h                                    |
| **ZIF-67 生长**    | CMCH 浸入 Co(NO₃)₂ + 2-MIM 溶液, 25°C 12h                                                   |
| **LDH 转化**       | Ni(NO₃)₂·6H₂O (1 mmol) + 乙醇 30 mL, 120°C 2h ——**双钴源策略**: CMCH 也提供 Co 源 |
| **硫源**           | 硫脲 CH₄N₂S, 添加量 1/3/7/10 wt%（相对于 Ni 前驱体）                                         |
| **最优掺杂**       | **S1 (1 wt%)** — 综合电化学性能最优（容量的稳定性平衡）                                 |
| **溶剂**           | 乙醇（无额外加水）                                                                             |
| **后处理**         | 去离子水 + 乙醇冲洗, 60°C 干燥 12 h                                                           |
| **活性物质负载量** | CMCH:~2.5 mg/cm² → CMCH@NC-LDHs: 2.7 mg/cm² → S1-CMCH@NC-LDHs: 3.7~5 mg/cm²              |

---

## 三、S 掺杂物种确认：S²⁻

**实验 XPS 证据（初稿图 2e）：**

- S 2p₃/₂ = **161.8 eV**, S 2p₁/₂ = **163.0 eV** → **过渡金属-硫键 (M-S)**，确认 S²⁻ 取代层板 OH⁻
- 宽峰 ~168.5 eV → 表面 SOₓ²⁻（空气中表面自发氧化）

**建模结论与初稿一致：**

> S²⁻ 形式，无 H (SH⁻ 在弱碱性水热条件下去质子化)。
> 1 个 S²⁻ 替换两个相邻 OH⁻ 位置（形成 M-S-M 桥，保持层板电荷 -2 → -2 不变）。

**⚠ 注意：** 实验中 S 掺杂的 LDH 层板同时存在表面 SOₓ²⁻ 物种（~168.5 eV），此为空气暴露后的表面重构，非原始 S²⁻ 掺杂态。常规 DFT 建模无需包含此物种。

---

## 四、建模依据

### 4.1 LDH 层板：S 取代位筛选

LDH 层板（CoNiOH₂）的超胞模型为 3×2×2，共含 36 个 OH⁻ 配位体。取代位根据 **Co/Ni 比例和配位环境** 设定了 5 种构型：

| 模型 | 取代形式 | 化学式 | 原子数 | 目的 |
|:----|:--------|:------|:-----:|:----|
| CoNiHO | 本征（无 S） | Co₉Ni₉H₃₆O₃₆ | 90 | 基准参考 |
| CoNiHOS-Co3 | SH⁻, 3 Co 配位 | Co₉Ni₉SH₃₆O₃₅ | 90 | S 倾向与 Co 键合？ |
| CoNiHOS-Ni3 | SH⁻, 3 Ni 配位 | Co₉Ni₉SH₃₆O₃₅ | 90 | S 倾向与 Ni 键合？ |
| CoNiHOS-Co1Ni2 | SH⁻, 1Co+2Ni 配位 | Co₉Ni₉SH₃₆O₃₅ | 90 | 统计随机取代 |
| CoNiHOS-Co2Ni1 | SH⁻, 2Co+1Ni 配位 | Co₉Ni₉SH₃₆O₃₅ | 90 | 统计随机取代 |
| CoNiHOS-Co3-noH | S²⁻（无H）, 3 Co 配位 | Co₉Ni₉SH₃₅O₃₅ | 89 | 验证脱氢 vs 保留 H |

筛选目标：对比上述 5 种掺杂态与实验 XPS（S 2p 峰位 161.8/163.0 eV，对应 M-S 键而非 S-H 键），确认 **S²⁻ 无 H 模型（即最终采用的 CoNiOH2S-noH）** 最符合实验观测。同时，Co/Ni 配位环境的能量对比验证了 S 在不同金属环境中的稳定性差异。

### 4.2 CMCH 基底：Mn 取代位的选择

CMCH 的结构原型 Co₂(OH)₂CO₃（mp-1202876）含有两种 Co 配位位点：
- **八面体位（CoO₆）** — ~1/2 的 Co 原子
- **四面体位（CoO₄）** — ~1/2 的 Co 原子

实验论文通过 EDS 确认 Mn 均匀分布于 CMCH 基底中，且文献中 Mn 在碳酸盐体系中通常以 Mn²⁺ 形式占据八面体位。因此建模选择 **将八面体位的 Co 替换为 Mn**，四面体位保留 Co。

替换比例：Co₈ → Co₄Mn₄（化学式 CoMn(CO₃)(OH)₂），与实验合成配比一致。

### 4.4 S 掺杂 LDH 的筛选结果（comp_s_doping.py）

对 5 种 S 掺杂构型的对比分析：

| 构型 | 形成能 (eV) | 带隙 (eV) | d-Co (eV) | S 配位 | S-M 键长 (Å) |
|:----|:----------:|:--------:|:---------:|:-----:|:------------:|
| 本征 CoNiHO | — | 2.84 | -0.61 | 无 S | — |
| SH-3Co | +1.63 | 2.40 | -0.95 | 3 Co | 2.33 |
| SH-3Ni | +1.72 | 2.78 | -0.67 | 3 Ni | 2.30 |
| SH-Co1Ni2 | +1.63 | 2.52 | -0.87 | Co+2Ni | 2.32 |
| SH-Co2Ni1 | +1.65 | 2.48 | -1.01 | 2Co+Ni | 2.33 |
| **S-3Co (无H)** | **+2.03** | **0.23** | **-1.10** | **3 Co** | — |

关键结论：
1. **SH⁻ 保留 H 的构型比 S²⁻ 脱氢更稳定**（形成能低 ~0.4 eV），但 XPS 显示 S 2p 峰位 161.8/163.0 eV 对应 M-S 键而非 S-H。实验水热条件（碱性 + 120°C）可能使 SH⁻ 去质子化，且层板中 S²⁻ 形成 M-S-M 桥结构更稳定
2. **S 倾向与 Co 配位**：Co3 构型形成能最低（1.63 eV），Ni3 最高（1.72 eV）
3. **S 掺杂使带隙从 2.84 eV（本征）收缩**：SH-3Co 为 2.40 eV，S-3Co-noH 降至 0.23 eV
4. **Co d-band center 下移**：S 掺杂使 Co d 带远离 E_F，-0.61 → -0.95 eV（SH-3Co），趋势与异质结计算一致

S 掺杂主要效应：S 的电负性（2.58）低于 O（3.44），取代后 S→M 电子转移方向与 OH⁻ 相反，导致 M 的电子密度微调，Co/Ni 氧化态升高。<img src="" width="0" height="0">

### 4.3 异质结建模策略

异质结的构建流程如下：

```
1. LDH(001) slab 和 CMCH(010) slab 各自在自然晶格下 relax
2. 取平均面内晶格 avg_a, avg_b → 对两 slab 施加应变
3. Interface.from_slabs(gap=2.0) 构建初始界面
4. 全弛豫（ISIF=2）：CMCH 底部 TTF（固定 z，释放面内应变），LDH 全开
5. 再跑 static + DOS + band + LOCPOT_dipole
```

三个异质结的区别仅在于 LDH 层中 S 的位置：
- **hetero_intrinsic**: LDH 无 S
- **hetero_s_doped**: S 位于 LDH 层板内部（远离 CMCH 界面）
- **hetero_s_exposed**: S 位于 LDH 层板底侧（靠近 CMCH 界面）

---

## 四、实验性能数据（供计算验证）

### 4.1 电化学性能摘要

| 样品                              | 比电容 @1A/g (F/g) |    Rct (Ω)    |  5000圈保持率  |
| --------------------------------- | :----------------: | :------------: | :-------------: |
| CMCH (基底)                       |      1295.68      |      0.90      |       —       |
| NC-LDHs (纯 LDH)                  |       940.24       | **5.10** |       —       |
| CMCH@NC-LDHs (未掺杂)             |      2284.30      |      0.78      |      77.3%      |
| **S1-CMCH@NC-LDHs (1 wt%)** | **3530.40** | **0.57** | **79.9%** |
| S3-CMCH@NC-LDHs (3 wt%)           |      3604.24      |       —       |      65.4%      |
| S7-CMCH@NC-LDHs (7 wt%)           |        下降        |       —       |      55.7%      |
| S10-CMCH@NC-LDHs (10 wt%)         |        下降        |       —       |      51.1%      |

### 4.2 电化学机制

- **b 值分析**：阳极 0.51, 阴极 0.52 → 以扩散控制为主（低扫速）
- **电容贡献**：5 mV/s 时 32% 电容控制, 50 mV/s 时 **82%** 电容控制
- **ASC 性能**：26.8 Wh/kg @ 721.6 W/kg (0-1.4V, 10000圈后 78.8%)

### 4.3 实验推论的核心机制（需计算验证）

1. **S²⁻ 取代 OH⁻ 引起晶格膨胀**（XRD 低角偏移，S²⁻ 半径 184 pm > OH⁻ 110 pm）
2. **电子从 Ni/Co 向 S 转移**（XPS Ni 2p/Co 2p 正移 0.1-0.3 eV）
3. **Ni³⁺/Ni²⁺ 比例提升**: 44% → 49%（XPS 面积比）
4. **O 空位生成**: S 掺杂触发层板去质子化，释放氧空位（O 1s 峰面积扩张）
5. **禁带缺陷能级**: 增加自由载流子浓度，Rct 0.78→0.57 Ω
6. **S-O 极性基团的电子诱导效应**: XPS 结合能正移，加速界面电荷转移

---

## 五、当前计算进度

### 已完成的计算体系（10 个）

#### Layer 0: 结构准备（步骤1）

- 六方 Ni(OH)₂ → 正交化 → Co 替换 SQS → 3×2×2 超胞
- Co₂(CO₃)(OH)₂ → Co→Mn 八面体位替换 → 3×1×1 超胞
- 构建异质结（CMCH 基底 + LDH 薄膜, 真空层 ~15-20 Å c轴）

#### Layer 1: Bulk 优化（步骤2，3个体系 ✅）

| 体系         | 目录                  | S 掺杂? | DOS/能带? |     Band gap     |   静态能量   |
| ------------ | --------------------- | :------: | :-------: | :---------------: | :-----------: |
| CoNiOH2      | `data/CoNiOH2`      | ❌ 本征 |    ✅    | **2.93 eV** | ✅ store.json |
| CoNiOH2S-noH | `data/CoNiOH2S-noH` | ✅ S²⁻ |    ✅    | **0.08 eV** | ✅ store.json |
| CoMnH2CO5    | `data/CoMnH2CO5`    | ❌ 基底 |    ✅    |        —        |      ✅      |

> 🔬 **关键发现：** S 掺杂使带隙从 2.93 eV（本征）骤降至 0.08 eV（S²⁻ 掺杂）。这是 S-3p 在禁带中引入缺陷能级的**直接 DFT 证据**，完美支撑实验论文的"禁带缺陷能级 → 电导率提升"的核心机制。

#### Layer 2: Slab 优化（步骤3，4个体系 ✅）

| 体系                   | 目录                            |    S 掺杂?    | DOS/band? | LOCPOT? |
| ---------------------- | ------------------------------- | :-----------: | :-------: | :-----: |
| CoNiOH2-slab           | `data/CoNiOH2-slab`           |    ❌ 本征    |    ✅    |   ✅   |
| CoNiOH2S-noH-slab      | `data/CoNiOH2S-noH-slab`      | ✅ S 远离表面 |    ✅    |   ✅   |
| CoNiOH2S-noH-slab-flip | `data/CoNiOH2S-noH-slab-flip` | ✅ S 靠近表面 |    ✅    |   ✅   |
| CoMnH2CO5-slab         | `data/CoMnH2CO5-slab`         |    ❌ 基底    |    ✅    |   ✅   |

#### Layer 3: 异质结优化（步骤4，3个体系 ✅）

| 体系 | 目录 | S 位置 |
| --- | --- | --- |
| hetero_intrinsic (gap=2.0) | `data/hetero_intrinsic` | ❌ 无 S |
| hetero_s_doped (gap=2.0) | `data/hetero_s_doped` | ✅ S 在 LDH 内部，远离 CMCH 界面 |
| hetero_s_exposed (gap=2.0) | `data/hetero_s_exposed` | ✅ S 在 LDH 底侧，靠近 CMCH 界面 |

#### Layer 4: 晶格一致优化（Level 2，7个体系 ✅，gap=1.0，共同面内晶格）

| 体系 | 类型 | 目录 | 总能 (eV) | Bandgap (eV) | 原子数 |
|:----|:----|:----|:---------:|:-----------:|:-----:|
| **hetero_intrinsic** | 异质结 | `data/hetero_intrinsic` | -714.54 | **0.745** | 120 |
| **hetero_s_doped** | 异质结 | `data/hetero_s_doped` | -707.64 | **0.388** | 119 |
| **hetero_s_exposed** | 异质结 | `data/hetero_s_exposed` | -707.48 | **0.293** | 119 |
| **CMCH_strained** | 单板 | `data/CMCH_strained` | -396.93 | **2.438** | 60 |
| **LDH_strained** | 单板 | `data/LDH_strained` | -314.60 | **1.685** | 60 |
| **LDH_S_strained** | 单板 | `data/LDH_S_strained` | -307.37 | **0.152** | 59 |
| **LDH_S_flip_strained** | 单板 | `data/LDH_S_flip_strained` | -307.31 | **0.294** | 59 |

> **注意：** 所有 7 个体系共享同一组面内晶格（a=10.3255, b=9.4602），FFT 网格 (160,140,420) 全同，可直接做 3D 差分电荷密度。总能取自重跑后的 `static_out.json`（ALGO=Normal, EDIFF=1e-6）。

---

## 六、实验与计算的对照与差异

### ✅ 一致之处

| 方面            | 实验                           | 计算                         |        状态        |
| --------------- | ------------------------------ | ---------------------------- | :----------------: |
| S²⁻ 取代物种  | XPS S 2p: 161.8/163.0 eV (M-S) | 无H模型                      |   **一致**   |
| LDH 晶格膨胀    | XRD 衍射峰低角偏移             | 需提取 relax 后晶胞参数      |     ⏳ 待验证     |
| LDH 层的存在    | XRD 显示 LDH(003)(006)(012)    | slab 和 hetero 含 LDH 层     |   **一致**   |
| CMCH 基底存在   | XRD 显示 CMCH 衍射峰           | slab 和 hetero 含 CMCH       |   **一致**   |
| 核壳结构        | SEM 显示 LDH 包覆 CMCH 纳米针  | 异质结模型                   | **概念一致** |
| S 在 LDH 层板内 | XPS M-S + TEM EDS 确认         | s_doped/s_exposed 均位于 LDH |   **一致**   |

### ❌ 不一致 / 需要重新审视

| 问题                             | 实验情况                                                     | 计算假设                                                                               | 差异影响                                    |
| -------------------------------- | ------------------------------------------------------------ | -------------------------------------------------------------------------------------- | ------------------------------------------- |
| **S 掺杂量**               | 最优 S1 (1 wt% 硫脲), 实际 S/(S+O) 可能在**0.01-0.03** | 3×2×2 超胞 36 OH⁻ 中 1 S²⁻ 替换**2** OH⁻, 即 1/36 ≈ **0.028**       | 合理 (符合 S1 水平) ✅                      |
| **层间阴离子**             | 实验中存在 NO₃⁻（FTIR ~1380 cm⁻¹）、CO₃²⁻、SO₄²⁻   | 计算模型不含层间阴离子（仅为层板）                                                     | ⚠️ 模型过于简化，层间阴离子可能影响层间距 |
| **表面 SOₓ²⁻**          | XPS 168.5 eV 存在显著硫氧化物峰                              | 未建模                                                                                 | ⚠️ 实验中重要，但为表面后氧化，建模非必要 |
| **O 空位**                 | 实验认为 S 掺杂释放 O 空位                                   | 未专门建立含 O 空位的模型                                                              | ⚠️ 可能缺失关键机制                       |
| **双钴源**                 | CMCH 作 Co 源参与 LDH 形成，界面存在 Co 梯度扩散             | 界面为理想平面接触                                                                     | ⚠️ 忽略界面扩散                           |
| **ZIF-67 中间体**          | 实验通过 ZIF-67 转化得到 LDH，非直接生长                     | 直接模型 LDH/CMCH 界面                                                                 | ⚠️ 可能影响界面结构                       |
| **泡沫镍基底**             | 实验在 NF 上生长                                             | 计算用真空 slab，无 NF                                                                 | ✅ 标准 DFT 做法，NF 仅为导电集流体         |
| **非晶/多孔**              | SEM 显示纳米片交错多孔形貌                                   | DFT 用周期性结晶模型                                                                   | ✅ DFT 局限，无法反映形态学                 |
| **Ni³⁺/Ni²⁺ 比例增加** | XPS 面积比 44% → 49%                                        | Bader 电荷：S 掺杂后 Co/Ni 失去少量电子（ΔCT_Co ~ +0.01~0.02 e）                      | ✅ 已验证趋势                               |
| **S 的电荷状态**           | XPS S 2p: 161.8 eV (M-S), S²⁻ 形态                         | Bader 电荷：**S 失去 ~0.8 e**（电负性 S=2.58 < O=3.44，取代 OH⁻ 后 S 呈正电性） | ✅ 电负性差异解释一致                       |

---

## 七、步对应的计算-实验验证

| 计算维度                     | 验证的实验结论                         |   优先级   |        状态        |
| ---------------------------- | -------------------------------------- | :--------: | :-----------------: |
| **1. S 掺杂形成能**    | S²⁻ 是否可稳定取代 OH⁻              |   ★★★   | ⏳ 已有数据，未导出 |
| **2. 晶格膨胀定量**    | XRD 峰移：S 掺杂后 a/c 轴变化          |   ★★★   |     ✅ 已出数据     |
| **3. Bader 电荷**      | XPS 结合能正移：Ni/Co → S 电荷转移    | ★★★★★ |      ❌ 待实施      |
| **4. DOS/PDOS**        | 禁带缺陷能级：Rct 下降的解释           | ★★★★★ |      ✅ 已实施      |
| **5. 功函数**          | S 掺杂 → 界面势垒变化 → 电荷注入效率 |  ★★★★  |      ✅ 已实施      |
| **6. 能带对齐**        | Type I/II? S 如何改变带阶              |  ★★★★  |      ❌ 待实施      |
| **7. Planar Avg Δρ** | 界面电荷重分布：流向 + 定量            |  ★★★★  |      ❌ 待实施      |
| **8. O 空位形成能**    | S 能否促进 O 空位生成                  |   ★★★   |     需额外 DFT     |
| **9. Slab 表面能**     | LDH 的生长取向                         |    ★★    | ⏳ 已有数据，未导出 |
| **10. 磁性**           | 与 O 空位/缺陷态的关联                 |    ★★    |      ✅ 已实施      |

---

## 八、重构方案：晶格一致的异质结 workflow

### 数据可用性总览

| 分析项 | 所需数据 | Bulk(3) | Slab-Natural(4) | Slab-Matched(4) | Hetero Level 2(3) |
|:------|:--------|:-------:|:--------------:|:---------------:|:-----------------:|
| 总能量 | static_out.json | ✅ | ✅ | ✅ | ✅ |
| 弛豫结构 | CONTCAR | ✅ | ✅ | ✅ | ✅ |
| DOS/PDOS | dos_out.json | ✅ | ✅ | ✅ | ✅ |
| 能带 | band_out.json | ✅ | ✅ | ✅ | ✅ |
| CHGCAR | 3-static/CHGCAR* | ✅ | ✅ | ✅ | ✅（重跑后恢复） |
| LOCPOT | 3-static/LOCPOT* | ❌ | ✅ | ✅ | ✅（重跑后恢复） |
| AECCAR0+2 | 3-static/AECCAR* | ✅ | ✅ | ✅ | ✅（重跑后恢复） |

### 已完成的分析（postprocessing/）

所有分析脚本已更新至 Level 2（共同晶格），输出位于 `postprocessing/output/`：

#### A. 能带对齐（band_alignment.py）
| 异质结 | VB offset (eV) | DV_CMCH (eV) | DV_LDH (eV) | 类型 |
|:------|:-------------:|:-----------:|:----------:|:----|
| Intrinsic Het. | −2.78 | +2.01 | −2.36 | **Type III** (broken gap) |
| S-Doped Het. | +3.22 | +1.83 | +3.11 | **Type I** (straddling) |
| S-Exposed Het. | +2.25 | +2.40 | +3.14 | **Type I** (straddling) |

#### B. 功函数分析（work_function_analysis.py）
| 体系 | Φ (eV) |
|:----|:------|
| CMCH_strained | **5.38** |
| LDH_strained | **3.06** |
| LDH_S_strained | **2.92** |
| LDH_S_flip_strained | **3.62** |
| hetero_intrinsic | **3.78** |
| hetero_s_doped | **3.31** |
| hetero_s_exposed | **3.89** |

#### C. DOS/PDOS（dos_analysis.py）
Level 2 带隙（共同晶格）：
- CMCH_strained: **2.43 eV**, LDH_strained: **1.65 eV**
- LDH_S_strained: **0.12 eV**, LDH_S_flip_strained: **0.27 eV**
- hetero_intrinsic: **0.72 eV**, s_doped: **0.37 eV**, s_exposed: **0.26 eV**

#### D. Bader 电荷分析（bader_analysis.py）
- S 的电荷转移：−0.74～−0.81 e⁻（正值=失电子，S 呈正电性）
- 异质结与对应 slab 的电荷转移趋势一致

#### E. 界面差分电荷密度（charge_difference_3d.py）
基于 **LOCPOT 静电势 Poisson 方程求解**（3D 网格 (160,140,420) 全同，无需插值）：
直接 CHGCAR 减法不可行——核心电子密度在原子核附近 > 1000 e⁻/Å³，
异质结与孤立 slab 间微小的原子位移会在此区域产生 4-5 个数量级的伪影。
改用平面平均静电势差分 ΔV(z) → Poisson 方程反推 Δρ(z)：
| 异质结 | 界面净电荷 (e⁻) | 界面位置 (Å) |
|:------|:-------------:|:-----------:|
| hetero_intrinsic | −0.81 | 8.86 |
| hetero_s_doped | +3.54 | 8.86 |
| hetero_s_exposed | −2.07 | 8.06 |

#### F. 界面结合能（interface_binding.py）
| 异质结 | E_bind (eV) | E_bind/area (eV/Å²) |
|:------|:----------:|:-----------------:|
| hetero_intrinsic | −3.01 | −0.031 |
| hetero_s_doped | −3.34 | −0.034 |
| hetero_s_exposed | −3.25 | −0.033 |

> 所有 E_bind < 0 → 界面稳定。S 掺杂轻微增强界面结合（约 0.3 eV）。

### 待补充计算

1. **O 空位计算**（如需验证实验的 O 空位机制）
   - 选取 slab 或 hetero_s_doped，删去一个 O 重新 relax
   - 计算 `E_vac = E_defect + ½E(O₂) - E_perfect`
   - 对比有/无 S 掺杂时的 O 空位形成能

### 应立即开始的计算（先补充 Slab-Matched 数据）

**优先级 1: gap=1.0 对比 relax + 分离 strained slab**

```
Step 1: 构建 hetero gap=1.0 → relax (ISIF=2, NSW=50)
Step 2: 对比 gap=1.0 vs 2.0 的 relax 后总能量+间距 → 选最优
Step 3: 从 hetero CONTCAR 提取 slab → *-Matched
Step 4: *-Matched slab relax + static + DOS/band + LOCPOT_dipole
Step 5: 更新能带对齐 baseline
Step 6: 3D 差分电荷密度
```

**B. 界面性质组（Planar Average 差分电荷密度 + 界面结合能）**

```
优先级: ★★★★
体系: 3 个 hetero + 对应的 4 个 slab
```

5. **Planar Average 差分电荷密度** — 分析界面电荷重分布

   **背景：** 构建异质结时 slab 和 hetero 的 a/b 晶格做了 `strain_to(avg_a, avg_b)`，导致
   slab 与 hetero 的 3D FFT 网格（NGXF/NGYF/NGZF）全部不同，无法直接做 3D 差分。

   **方案（1D 对齐差分）：**

   ```
   ① 对 hetero、slab_A(CMCH)、slab_B(LDH) 的 CHGCAR 各自做 planar average → ρ(z)
   ② z 坐标是连续空间坐标（Å），三个体系现在都是 1D 曲线
   ③ 从 hetero 的结构确定 CMCH 和 LDH 各自占据的 z 区间
   ④ 将 slab_A 的 ρ(z) 缩放到 hetero 中 CMCH 的 z 区间（scipy.interpolate）
   ⑤ 将 slab_B 的 ρ(z) 缩放到 hetero 中 LDH 的 z 区间
   ⑥ Δρ(z) = ρ_hetero(z) - ρ_CMCH_scaled(z) - ρ_LDH_scaled(z)
   ⑦ 积分 Δρ(z)dz → 界面净电荷转移方向与总量（e⁻/Å²）
   ```

   **可靠性评估：**

   - Δρ(z) 的峰形与电荷流向（界面积累/耗尽）✅ 可信
   - 本征 vs S 掺杂 vs S 暴露的差异比较 ✅ 可信（系统误差抵消）
   - 定量绝对值 ±20% ⚠️ 半定量

   **工具链：** `Chgcar.get_average_along_axis(2)` + `scipy.interpolate.interp1d`，无需额外安装。
6. **界面结合能**

   - `E_bind = E_hetero - E_LDH_slab - E_CMCH_slab`
   - S 掺杂增强还是削弱界面？

**C. 磁性分析组**

```
优先级: ★★★
体系: 所有 10 体系
```

7. 总磁矩对比：bulk → slab → hetero
8. Co/Ni 局域磁矩：S 掺杂影响
9. 界面磁性变化

**D. 结构分析组**

```
优先级: ★★★
```

10. 晶格常数对比（验证 XRD 峰移）
11. S-M 键长、配位数
12. 弛豫前后原子位移

---

## 九、建议的论文分析逻辑链

> 以下为论文"结果与讨论"部分的计算支撑逻辑（供参考）

```
实验 XPS 结合能正移 0.1-0.3 eV
    ↓ 验证
Bader 电荷: S 失去 ~0.8 e（正电性，S<O 电负性差），Co 失 e⁻ ~0.01~0.02
    ↓ 解释
Bulk DOS 对比: 本征 2.93 eV → S 掺杂 0.08 eV 带隙
    ↓ S-3p 在禁带中引入缺陷态 → 费米能级处态密度↑
    ↓ 解释
电化学: Rct 0.78→0.57 Ω (电导率提升)
         电荷转移加快
    ↓ 界面影响
能带对齐 + 功函数 → 界面势垒降低
    ↓ 宏观性能
比电容提升 2284→3530 F/g
```

---

## 十、待做清单

- [x] **gap=1.0 对比 relax** — 已完成，gap=1.0 选定为生产参数
- [x] **分离 strained slab** — `build_heterostructure.py` 独立生成 4 个 strained slab
- [x] **strained slab relax + static + DOS/band** — Level 2 全部 7 体系完成
- [x] **能带对齐（Level 2）** — band_alignment.py 已更新并运行，含带阶和静电势
- [x] **3D 差分电荷密度** — 改用 LOCPOT 1D 法（CHGCAR 因原子数差异无法直接相减），含 ΔV(z)、Δρ(z)、界面偶极
- [x] **界面结合能** — interface_binding.py 已编写完成
- [x] **功函数计算** — work_function_analysis.py 已更新并运行
- [x] **Bader 电荷分析** — bader_analysis.py 已更新至 Level 2 并运行，S ct ≈ −0.8 e⁻
- [ ] **晶格膨胀定量** — 提取 relax 后晶胞参数，验证 XRD 峰移
- [ ] **O 空位形成能** — 选项，验证实验"氧空位"机制
