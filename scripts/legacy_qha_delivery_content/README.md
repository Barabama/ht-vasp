# legacy DFPT-QHA 交付（A 线 → B 线）

**给 B 线的交接说明。先读本文件，再读 `LAYOUT_PROPOSAL.md`，再读各端基的 `FINDINGS.md`。**

---

## 1. 这是什么

`ht-vasp`（A 线）从 NAS 逐字节镜像了两个 legacy DFPT-QHA 端基，做了 fail-closed 审计，
并把结果整理成 B 线可以直接消费的形态。**本目录是交付层，不是镜像，也不是标准。**

| 端基 | NAS 来源 | 有 VASP 层？ | 审计状态 | applicability |
|---|---|---|---|---|
| `BCC-Nb-Ti-2` | `…\16symbols-cxy-hmd\BCC-Nb-Ti-2` | **无** | `unverifiable` | `unverifiable` |
| `BCC-Co-Ti-2` | `…\16symbols-cxy-hmd\BCC-Co-Ti-2` | 有（完整） | `fail`（37 条） | `degraded_defective` |

镜像本体（不可变，A 线未改动一个字节）：

```
ht-vasp/data/legacy-dfpt-qha-pilot/BCC-Nb-Ti/legacy-as-received/   30 文件 / 5.65 MiB
ht-vasp/data/legacy-dfpt-qha-pilot/BCC-Co-Ti/legacy-as-received/  487 文件 / 31.01 MiB
```

两者的 `source-manifest.json` 记录了每个文件的路径 / 大小 / mtime / sha256，
以及端基级内容摘要（`ca13b516…` / `f284934c…`）。

## 2. 两个端基一句话

- **`BCC-Co-Ti-2`：数据齐全、可完整重建，但重建出来的曲线带确定性缺陷。**
  E(V) 与声子层都能从目录内文件重算，而且逐点与随附 yaml 吻合到 1e-6 kJ/mol——
  但全部 10 个体积点的静态与 DFPT 都在 `ICHARG=1` 且无 CHGCAR 的条件下
  退回原子叠加初猜，且 QHA-08 的 DFPT 基态 SCF 发散（E = +4359 eV）。
  报告的 Gibbs 曲线在 1920→1930 K 跳 0.515 eV/cell。**不能直接当生产观测。**

- **`BCC-Nb-Ti-2`：只有曲线，没有 VASP。**
  13 个 `.dat` + 8 个 thermal_properties yaml，没有任何 `INCAR` / `POSCAR` /
  `FORCE_CONSTANTS` / `OUTCAR` / `phonopy.yaml`。曲线能读，
  但组成、结构、方法、生产者、参考态**五项全部 unresolved**。
  0 条 fail-closed 失败不代表通过，只代表没有东西可查。

## 3. 端基目录里有什么

```
<endpoint>/
├── DATASET.json        逐字段 claim + evidence + status（给机器）
├── FINDINGS.md         结论与缺陷（给人）
├── INVENTORY.json      本目录所有文件的路径 / 大小 / sha256，含与镜像 manifest 的核对结果
├── audit/audit.json    fail-closed 审计报告原文（legacy-qha-audit/v1）
├── curves/             *.dat 逐字节副本
├── phonons/            thermal_properties-*.yaml 逐字节副本
├── recipe/             INCAR* / KPOINTS* / mesh.conf 副本 + POTCAR_IDENTITY.txt
├── raw/RAW_INDEX.json  原始 VASP 层指针（路径 + sha256），不复制字节
└── structure/          胞；逐体积胞在 structure/volumes/
```

**每个副本都在生成时与镜像 manifest 对过 sha256**
（`INVENTORY.json.files[].matches_mirror_manifest`，全部为 `true`）。

### 字段约定

- `INVENTORY.json.files[].delivery_path` 是**相对交付根**的路径
  （形如 `BCC-Co-Ti-2/curves/v-e.dat`），所以整棵树搬到别处也仍然有效。
- `matches_mirror_manifest` 为 `null` 的条目是**非镜像来源**：
  `BCC-Nb-Ti-2/structure/from-ht-vasp/*`（从 ht-vasp 导入）与
  `recipe/POTCAR_IDENTITY.txt`（生成物）。这三个没有镜像 sha256 可对。
- `raw/RAW_INDEX.json` 里的 `mirror_path` 是**相对镜像根**的路径，
  镜像根记在同文件的 `mirror_root` 字段里。

## 4. B 线最需要的三个字段怎么用

`DATASET.json.claims` 的每个字段都是 `{value, status, evidence}` 三元组，
`status` 取值：

| status | 含义 | 对应 §B.5.6 的建议处置 |
|---|---|---|
| `resolved_with_evidence` | A 线有文件级证据 | 可作为 `metadata_source` 的证据来源 |
| `partially_resolved` | 部分有证据，其余不可恢复 | 降级；缺失部分必须另找来源 |
| `plausible_corroborated_not_proven` | 与证据一致但未被证明 | **不得**当作 resolved，需 profile 或用户确认 |
| `unresolved` | 无依据 | 必须由 profile 或用户给定，**禁止推定** |

按端基看：

| claim | BCC-Co-Ti-2 | BCC-Nb-Ti-2 |
|---|---|---|
| `composition_reduced` | `resolved_with_evidence` | `plausible_corroborated_not_proven` |
| `structure_prototype` | `resolved_with_evidence` | `unresolved` |
| `anchor_cell` | `resolved_with_evidence` | `unresolved` |
| `volume_grid` | `resolved_with_evidence` | `partially_resolved` |
| `calculation_method` | `resolved_with_evidence`（DFPT） | `unresolved` |
| `producer` | `resolved_with_evidence` | `unresolved` |
| `units` | `resolved_with_evidence` | `partially_resolved` |
| `temperature_grid` | `resolved_with_evidence` | `resolved_with_evidence` |
| `reference_role` | `unresolved` | `unresolved` |
| `energy_reference_state` | `unresolved` | `unresolved` |

**注意 `reference_role` 与 `energy_reference_state` 在两个端基里都是 `unresolved`。**
这不是 A 线偷懒：两个端基里确实没有任何 SER / GHSER / 元素参考的声明。
按 §B.5.6 的硬门槛，缺这两项就不能把 eV/cell 换成 J/mol-atom，
也就不能进 alignment / fit / TDB——这一项需要用户拍板。

## 5. A 线顺手发现的一条对 B 线也重要的教训

VASP 5.4.4 在电子步耗尽 `NELM` 时**不打印任何 `unconverged` 标记**。
`BCC-Co-Ti-2` 的 QHA-08 就是这样：基态 SCF 从 RMM 13 起发散到 +4359 eV，
但 OUTCAR 里 3 个 `aborting loop because EDIFF is reached` 全来自线性响应阶段，
"未收敛"一个字都没有。**任何"数标记"的 fail-closed 门都看不见它。**

本审计因此加了两道物理门：

1. DFPT 基态能量 vs 同体积 R3f 静态（QHA-08：280 265 meV/atom；健康点 4.6–28.7）；
2. 基态 SCF 迭代数是否撞 `NELM`（并同时报出末段能量跨度，以区分
   "发散"（4.9e6 meV）与"只是没到 EDIFF"（0.5–1.5 meV））。

B 线在做 applicability 判定时，**不要只依赖标记计数**。

## 6. 怎么重建本交付

生成部分（曲线 / 声子 / 配方 / 胞 / 审计副本 / 指针 / INVENTORY）由脚本产生：

```bash
cd /nfs_hdd/2025/gaominliang/ht-vasp

# 0) 清空交付目录（脚本会拒绝覆盖已存在的生成文件，这是刻意的保护）
rm -rf data/legacy-dfpt-qha-delivery

# 1) 审计（只读，不改镜像）
.conda/bin/python scripts/legacy_qha_audit.py \
    --root data/legacy-dfpt-qha-pilot/BCC-Nb-Ti --mesh 24
.conda/bin/python scripts/legacy_qha_audit.py \
    --root data/legacy-dfpt-qha-pilot/BCC-Co-Ti --mesh 24

# 2) 交付层
C=scripts/legacy_qha_delivery_content
D=data/legacy-dfpt-qha-delivery
.conda/bin/python scripts/build_legacy_qha_delivery.py \
    --root data/legacy-dfpt-qha-pilot/BCC-Nb-Ti --outdir $D --label BCC-Nb-Ti-2 \
    --external-structure-json data/endmembers/BCC-Nb-Ti/BCC-Nb-Ti-static.json \
    --external-structure-template data/poscars/BCC-Nb-Ti.vasp \
    --content-dir $C
.conda/bin/python scripts/build_legacy_qha_delivery.py \
    --root data/legacy-dfpt-qha-pilot/BCC-Co-Ti --outdir $D --label BCC-Co-Ti-2 \
    --content-dir $C

# 3) 根索引
.conda/bin/python scripts/build_legacy_qha_delivery.py --outdir $D --index-only
```

人工撰写的文档（`README.md` / `LAYOUT_PROPOSAL.md` / 各端基的
`DATASET.json` / `FINDINGS.md` / `structure/SOURCE.md`）放在
`ht-vasp/scripts/legacy_qha_delivery_content/`，由 `--content-dir` 覆盖进来——
所以整棵树是可重建的，`rm -rf` 之后重跑上面三条命令就能复原。

`build_legacy_qha_delivery.py` 拒绝覆盖已存在的文件，**不会**误伤人工文档。

## 7. 明确没有被做的事

- 没有提交任何 VASP / Slurm 作业。
- 没有修改任何科学输入（`ICHARG` / `ALGO` / `ENCUT` / k 网格 / POTCAR / 磁矩）。
- 没有改动 `legacy-as-received/` 下的任何字节，也没有删除任何东西。
- 没有把缺陷"修好"再交付：缺陷照原样交付并量化。
- 没有替用户冻结 §B.5.7 的任何一项。

`BCC-Nb-Ti-2` 的复现需要新跑 VASP，前置条件见其 `FINDINGS.md` §7，等用户决定。

## 8. A 线自己的待做（**与 B 线无关，B 线不必处理**）

用户已同意、但**列为待做、尚未执行**的方案：

> **模式跟随（mode following）**：拿软模的本征矢当位移方向 → 构造畸变结构 →
> 做静态 / 受约束弛豫 → 比较能量与对称性，以判断该软模对应的是真实结构失稳
> 还是数值伪影。

背景（BCC-Ti-V 试点）：R 点 (½,½,½) 的 −3.796883 THz 是 3 重简并、
且**只有 V 子晶格动**（Ti 的质量加权振幅精确为 0.0000）；候选畸变空间群为
V 沿 [100]→P4mm(99)、[110]→Amm2(38)、[111]→R3m(160)。
软模的本征矢分析见 `data/legacy-dfpt-qha-pilot/BCC-Ti-V/OPEN_QUESTIONS_20261008.md` §1.7。

此项与 B 线的交付无关，放在这里只是为了让 A 线的待办有据可查。
