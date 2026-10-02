# SQS-mixing 固化到 HT-VASP：进度与设计记录

> 状态：基础模型层、结构级 `Relax -> Static` 计算层、生成/编排层、CLI 和 `v1` 数据契约已落地；真实 endpoint 验证部分成功，points/collect 尚未完成
>
> 生产者：`ht-vasp`；后续唯一消费者：`em-tdb`
>
> 最后核对：2026-08-24

本文记录 SQS-mixing 为什么出现、它和 `sof-ace`/`em-tdb` 的关系、在
`ht-vasp` 中经历过的原型实现、已经完成的计算验证，以及下一版 OJ 风格工作流
的落地边界。文中严格区分“已经运行过的实验代码/数据”、“当前已落地的基础模型”
和“后续正式接口”。

## 1. SQS-mixing 的由来

### 1.1 CALPHAD 参数拟合需要什么数据

`em-tdb` 的目标不是保存某个 VASP 目录，而是为 CALPHAD 相模型提供可以审计、重算
和拟合的热力学观测。对于二元置换型无序相，最直接的 0 K 输入是若干内部成分的
DFT 混合焓：

```text
H_mix(x) = [E_SQS(x) - n_A E_A - n_B E_B] / (n_A + n_B)
```

其中 `E_SQS` 是含有 A/B 的无序近似超胞总能，`E_A`、`E_B` 是相同计算设置下的
纯元素端点每原子能量。对本阶段的 BCC 单亚晶格模型，混合焓被转换为
`J/mol-atom`，再由 `em-tdb` 拟合 Redlich-Kister 过剩 Gibbs 能：

```text
G_excess = x_A (1 - x_A) Σν L^(ν) (2 x_A - 1)^ν
```

这里的 `L^(ν)` 是无序相模型的交互参数。它不是 B2、L1_0 等有序相的双亚晶格
参数，也不是磁交换 `J_ij`、`TC` 或 `BMAGN`。

### 1.2 为什么使用 SQS

真实无序固溶体需要统计平均，直接用一个有限周期结构无法完全表示随机占位。
Special Quasirandom Structure（SQS，特殊准随机结构）通过在有限超胞中匹配目标的
近邻相关函数，使一个可计算的周期结构近似随机合金。它提供了两个必要条件：

1. 每个成分有明确的整数原子计数，能直接交给 VASP；
2. 结构、生成器、迭代次数和质量指标可以随结果一起保存，形成可追溯数据。

因此，SQS-mixing 不是独立的磁性工作流，而是 `DFT/SQS 结构 -> 混合焓 -> CALPHAD G/L 参数` 这条数据链中的 DFT 采样阶段。

### 1.3 它从哪里进入当前工作区

最初的工作流在 `sof-ace` 中形成。历史脚本
[`sof-ace/scripts/fit_l_parameters_sqs.py`](/nfs_hdd/2025/gaominliang/sof-ace/scripts/fit_l_parameters_sqs.py)
会读取 SQS 目录和 VASP 能量，计算混合焓并拟合 `L0/L1`。这条链路证明了科学需求，
但不适合作为正式生产边界，原因包括：

- 结构生成、VASP 结果搜索、能量解析和 TDB 写入耦合在脚本中；
- 依赖目录命名和历史 `out*.tdb`，缺少统一的计算设置 fingerprint；
- 端点会随每个二元体系重复计算；
- `em-tdb` 若直接读取目录，无法区分 relax/static、POTCAR 或旧参数；
- SQS 生成质量和结构 checksum 没有稳定地进入下游数据契约。

因此职责逐步收敛为：

```text
ht-vasp：生成 SQS、执行 Relax -> Static、输出规范 JSON
em-tdb：只接受规范 JSON，重算 Hmix、拟合 G/L、导出 TDB
sof-ace：保留 CE/MC/SQS 探索和 ESPEI 研究流程，不作为 em-tdb 的生产输入
```

这也是把 SQS-mixing 固化到 `ht-vasp` 的根本原因：VASP 输入、任务状态、结构和能量
必须在同一个计算生产者中形成可验证的边界，而不是让 TDB 工具反向猜测计算目录。

## 2. 设计演进时间线

### 2.1 初始数据契约：五点方案

提交 `7f83d8d`（2026-08-19）首次定义了
`htvasp-sqs-mixing/v1` 数据契约。最初设想是 BCC 16-site 超胞的五个内部点，
例如 `x_A = 0.125, 0.25, 0.50, 0.75, 0.875`，并要求：

- 两个纯元素参考和所有 SQS 点使用同一套 static 设置；
- 只接受成功 task document；
- 明确 `BCC_A2`、单亚晶格和 `parameter_target = G`；
- 在 JSON 中保存结构、能量、单位、SQS provenance 和计算 fingerprint。

该契约还明确了一个重要边界：SQS 数据拟合的是无序 BCC_A2 的 `G` 参数，不能把
它写成有序 B2 的 reciprocal `L` 参数。

### 2.2 第一轮工作流原型

提交 `d7ab160` 和 `13fb500` 先后尝试加入 Worker 与入口脚本。原型采用一个脚本
同时生成结构、运行端点、运行 SQS 点和收集 JSON，能够验证 atomate2 的
`RelaxMaker -> StaticMaker` 串接，但阶段职责还没有完全分离。

随后主线对该实验功能做过回退（`4bcf542`、`cf0b2ab`、`d5698bb`）。这些回退说明
该功能当时还不是稳定主线 API；当前正式设计不能把回退前的目录或接口当成兼容承诺。

### 2.3 `feat-sqs-mixing` 实验实现

提交 `cf93be6`（2026-08-22）在独立分支加入了较完整的原型，主要内容为：

- `htvasp/model/sqs.py`：`SqsSystem`、`SqsPlan`、成分点、SQS 生成适配、
  `SqsMixingData` Pydantic schema；
- `htvasp/workflows/sqs_mixing.py`：一个结构对应的 `Relax -> Static` Worker；
- `main_sqs.py`：生成、运行、收集、检查和 Slurm 入口；
- `tests/test_sqs_model.py`：离线模型和 JSON round-trip 测试；
- `scripts/backfill_tb2j_structures.py`：与 SQS 无关的 TB2J 辅助脚本。

这个原型已把“有序端基”和“无序 SQS”分开：`Endmember` 保持方向敏感的
`phase-A-B` 语义，`SqsSystem` 将 `phase-A-B` 与 `phase-B-A` 规范化为同一无序二元
体系。它仍把纯元素端点放进每个体系的运行目录，因此暴露了端点重复计算问题。

### 2.4 从五点改为三点

实际生产成本评估后，当前正式 recipe 改为三个内部点：

```text
x_A = 0.25, 0.50, 0.75
```

原因是当前第一阶段只拟合 BCC_A2 的低阶（`L0 + L1`）模型，三个点足以支撑低阶
smoke/生产数据，同时显著减少每个二元体系的 VASP 数量。旧五点设计保留为历史
记录，不再作为新分支的输入要求，也不需要向后兼容旧结果。

## 3. 当前实验参数和物理边界

### 3.1 默认原型和超胞

三种端基原型沿用 Materials Project 模板：

| 相  | MP 模板      |      默认超胞 | 原型基元数 | 默认站点数 |
| --- | ------------ | ------------: | ---------: | ---------: |
| BCC | `mp-1487`  | `2 x 2 x 2` |          2 |         16 |
| FCC | `mp-2258`  | `2 x 2 x 2` |          4 |         32 |
| HCP | `mp-20112` | `2 x 2 x 1` |          8 |         32 |

BCC、FCC、HCP 的默认规模分别是 16、32、32 个站点。文档和 CLI 不再使用笼统的
“16-site”说法。正式实现应在 `generate` 阶段从结构对象计算 `n_sites`，并把它写入
manifest。模板只应从本地快照读取，计算阶段不调用 Materials Project API；更新模板
必须是显式、可审阅的操作。

### 3.2 VASP 参数

当前实验已经使用用户确认的三点生产参数：

| 阶段     | 关键设置                                                                         |
| -------- | -------------------------------------------------------------------------------- |
| Relax    | `ENCUT=400 eV`，`PREC=Normal`，`EDIFF=1e-5`，`EDIFFG=-0.05`，`ISIF=3`  |
| Static   | `ENCUT=520 eV`，`PREC=Accurate`，`EDIFF=1e-6`，`NSW=0`，`IBRION=-1`    |
| 两阶段   | `ISPIN=2`，保留初始 `MAGMOM`，PBE，atomate reciprocal density `64`         |
| 输出     | `LCHARG=False`，`LAECHG=False`，`LORBIT=10`                                |
| 文件传递 | Relax 写`WAVECAR=True`；Static 读取 Relax 的 `WAVECAR`，自身 `LWAVE=False` |

这套设置的目的不是追求所有电子结构输出，而是得到可比较的结构和 static 总能。
SQS-mixing 第一阶段不消费 DOS、Bader、TB2J 或磁交换文件。

## 4. 已经完成的实验验证

### 4.1 BCC-Al-Nb

日志 [`logs/BCC-Al-Nb-sqs.log`](/nfs_hdd/2025/gaominliang/ht-vasp/logs/BCC-Al-Nb-sqs.log)
显示旧原型已完成：

- 纯元素端点：`BCC-Al-Al`、`BCC-Nb-Nb`；
- 混合点：`x0.25`、`x0.5`、`x0.75`；
- 每个点均为一个 Relax 后接一个 Static；
- 最终写出旧原型文件
  `data/sqs_mixing/BCC-Al-Nb/BCC-Al-Nb-sqs-mixing.json`。
  正式流程改为在 `data/endmembers/BCC-Al-Nb/sqsflow/` 保存阶段数据，并在其上级目录
  写出 `data/endmembers/BCC-Al-Nb/BCC-Al-Nb-sqs.json`。

已生成的数据还包括 manifest、SQS POSCAR、生成 provenance、static task document
和运行目录。该结果证明了参数、文件 handoff 和 Hmix 计算路径可以在本地
atomate2/jobflow runner 中跑通，但它属于旧原型的实验产物；由于 `data/` 和运行日志
被 `.gitignore` 忽略，不能视为已提交的稳定 fixture。

### 4.2 BCC-Al-Ti

manifest 已生成三个内部点，`BCC-Al-Al` 和 `BCC-Ti-Ti` 端点也已完成。`x0.25`
的 Relax 失败，日志中出现：

- `ZBRENT: fatal error in bracketing`；
- `VasprunXMLValidator: Premature end of data`；
- 后续因为没有成功的 Relax 输出而找不到 `sqs mixing static`。

因此 `x0.5`、`x0.75` 没有形成完整的 mixing JSON。该失败暴露出旧入口把“运行失败
后的收集”混在一起：它移动了部分目录后才在收集阶段报错。新流程必须在阶段边界
先验证 task document 和端点 fingerprint，缺失时直接阻止下一阶段。

### 4.3 当前证据的正确解读

上述运行说明旧原型已经验证了：

1. 3 点 recipe 能构造合法整数成分；
2. Relax 到 Static 的 WAVECAR handoff 可用；
3. 规范 JSON 可以保存结构、能量、单位和 provenance；
4. 端点重复计算确实存在；
5. 失败任务需要结构化的阶段状态，而不是依赖目录是否存在。

它们不说明以下内容已经完成：

- 新 `feat-sqs-mix` 分支的 OJ 风格实现；
- FCC/HCP 的 TDB 映射；
- `em-tdb` 的 JSON 校验、拟合和 TDB 导出；
- 生产 Slurm 批量依赖关系；
- 多个 SQS 构型的统计平均或不确定度分析。

## 5. 为什么要改成 OJ 风格编排

旧 `main_sqs.py` 用一个 `run_tick` 混合了四种语义：生成、端点计算、成分点计算和
收集。这样做在小规模 smoke 中简单，但会产生三个实际问题：

1. 每个 `phase-A-B` 目录重新计算 A/A、B/B，跨二元体系无法去重；
2. 端点不完整时仍可能启动 SQS 点，失败会延迟到收集阶段才暴露；
3. 目录、task document 和最终 JSON 之间没有清晰的阶段契约。

新设计借鉴 OJ/TB2J 的边界：入口脚本或阶段编排器决定“本阶段应该运行哪些任务”，
Worker 只负责一个结构的 `Relax -> Static` Flow，Maker 只负责构造 Flow，中间 JSON
负责跨阶段和跨 Slurm 任务传递可靠结果。

### 5.1 四个显式阶段

```text
generate  ->  endpoints  ->  points  ->  collect
```

| 阶段          | 负责内容                                      | 禁止内容                           |
| ------------- | --------------------------------------------- | ---------------------------------- |
| `generate`  | 模板超胞、SQS 三个 POSCAR、manifest           | 启动 VASP、读取旧运行结果          |
| `endpoints` | 对所有体系的唯一纯元素端点运行 Relax→Static  | 运行 SQS 点、拟合 Hmix             |
| `points`    | 在端点完整且 fingerprint 匹配后运行三个内部点 | 重算端点、猜测缺失结果             |
| `collect`   | 读取规范结果，重算 Hmix，写最终 JSON          | 解析 OUTCAR、扫描旧 TDB/ESPEI 输出 |

正式入口按阶段提交，便于检查端点完成情况和重试单点。

### 5.2 端点去重

例如声明：

```text
BCC-Al-Nb
BCC-Al-Ti
BCC-Nb-Ti
```

端点阶段只生成：

```text
BCC-Al-Al
BCC-Nb-Nb
BCC-Ti-Ti
```

`phase-A-B` 和 `phase-B-A` 在无序 SQS 模型中是同一个体系；元素必须规范化排序。
端点缓存键不能只有目录名，还必须包含：phase、element、supercell、模板 checksum、
计算设置 fingerprint。修改模板、POTCAR 或 INCAR 后自动得到新 cache key，避免误用
旧结果。

## 6. 新分支的目标模块边界

当前 `feat-sqs-mix` 的 HEAD 为 `0108ef7`，只包含基础 ht-vasp/TB2J 维护；SQS
源码仍在历史分支 `feat-sqs-mixing` 的 `cf93be6`，没有被新分支继承。该历史提交
可以作为数据模型和 `Relax -> Static` 原型的参考，但不能直接迁移：它把模型、SQS
引擎适配、结果解析和最终 JSON schema 集中在一个 743 行的 `sqs.py` 中，且引用了
提交中实际不存在的 `htvasp/model/templates.py`。

新实现按以下边界落地，不继续扩展历史 `run_tick`：

```text
htvasp/model/             领域对象、模板快照和跨阶段数据契约
htvasp/sqs/               SQS 生成、atomate2 jobs、Maker 和阶段编排
htvasp/workflows/         Worker 用户接口
main_sqs.py               CLI 和 Slurm 调度入口
```

### 6.1 `htvasp/model/sqs.py`

只放不执行计算的领域模型：

- `SqsSystem`：规范化 `phase-A-B`/`phase-B-A`；
- `SqsCompositionPoint`：整数计数、`x_A` 和站点数；
- `SqsPlan`：模板引用、超胞、成分点和生成参数；
- `SqsBatchPlan`：多个体系及其去重后的端点集合；
- `SqsEndpointSpec`、`SqsPointSpec`：跨阶段计算任务描述和稳定 ID；
- `SqsStructureArtifact`：结构路径、组成和 structure checksum；
- `SqsCalculationSpec`：计算设置 fingerprint、任务类型和输入结构引用。

该模块不运行 VASP，不提交 Slurm，不调用 `sqsgenerator`，不解析 VASP 输出，也不
计算最终 Hmix。endpoint registry 属于阶段编排层，不由单个模型对象持有。

### 6.2 `htvasp/model/templates.py`

只管理本地模板快照：读取结构、验证文件/结构 checksum，并返回不可变的
`TemplateSnapshot`。模板记录至少包含 phase、源结构信息、source species 和 checksum。
Materials Project 网络访问只能由显式模板更新命令触发，不能隐藏在生产计算中。

### 6.3 `htvasp/model/sqs_schema.py`

从计划模型中分离最终数据契约：

- `StaticObservation`：Static task document 标准化后的总能量、每原子能量、站点数、
  state 和 checksum；
- `PureReference`：端点最终结构和 Static observation；
- `SqsMixingPoint`：成分、SQS provenance、最终结构和 Hmix；
- `CalculationProvenance`：Relax/Static 设置、POTCAR、k 点和 fingerprint；
- `CalphadMapping`、`SqsMixingData`：`em-tdb` 唯一输入，schema 为
  `htvasp-sqs-mixing/v1`。

这些对象可以从标准化 JSON 读写，但不负责扫描 VASP 目录。最终 Hmix 只由 Static
总能量和端点每原子能量重算。

### 6.4 `htvasp/sqs/generator.py`

定义 SQS 生成器接口和 `sqsgenerator` 适配器。输入为模板结构和
`SqsCompositionPoint`，输出为 `SqsStructureArtifact` 及生成 provenance：engine/version、
iterations、seed、配置、objective 和结构 checksum。生成阶段已经知道所有成分点，
不需要使用 OJ 的动态 `Response(detour=...)` 机制。

### 6.5 `htvasp/sqs/jobs.py` 和 `htvasp/sqs/task_doc.py`

`jobs.py` 只放轻量的 `@job` 函数，例如把 Static task document 转换成
`StaticObservation`。它只能消费上游 `static_job.output`，不能解析 OUTCAR、扫描目录
或寻找其他任务的结果。

`task_doc.py` 定义计算层输出模型，负责成功状态、能量单位、站点数和 checksum 校验。

### 6.6 `htvasp/sqs/maker.py`

定义通用的 `SqsCalculationMaker`，不区分 endpoint 和 point 的 Flow 结构：

```python
class SqsCalculationMaker(Maker):
    def make(
        self,
        structure: Structure,
        spec: SqsCalculationSpec,
        prev_dir: str | Path | None = None,
    ) -> Flow: ...
```

每个 Flow 固定为：

```text
sqs relax -> sqs static -> sqs observation
```

其中 Static 必须使用 `relax_job.output.structure` 和
`relax_job.output.dir_name` 完成结构及 WAVECAR 交接。第一版使用 atomate2 的
`RelaxMaker` 和 `StaticMaker`，后续不引入 `DoubleRelaxMaker`；因此这里的
`Relax -> Static` 是稳定的工作流语义，而不是对其它 Relax 策略的概括。

endpoint 和 point 的差异只体现在 `SqsCalculationSpec.kind` 和元数据，不复制两套
Maker 实现。

### 6.7 `htvasp/workflows/sqs_mixing.py`

`SqsMixingWorker` 负责一个结构对应的一次计算：

- 构造并持有 `SqsCalculationMaker`；
- 生成一个 `Relax -> Static -> Observation` Flow；
- 调用基类 `run_flow`，使用既有 JSONStore、flow_dir/store_dir 和 resume 机制；
- 获取最终 observation。

Worker 不负责生成多个 SQS 结构、endpoint 去重、阶段前置检查、批量 Hmix 或 Slurm
批处理。这样它与 `OJWorker`、`Tb2jWorker` 的用户接口一致，同时保持计算边界清晰。

### 6.8 `htvasp/sqs/orchestrator.py`

阶段编排器负责跨任务状态，而不是构造 VASP 细节：

- `generate`：读取模板、生成 POSCAR、manifest 和唯一 endpoint registry；
- `endpoints`：按 endpoint cache key 去重并运行端点任务；
- `points`：检查端点完整且 fingerprint 一致后，运行三个独立成分点；
- `collect`：读取 observation JSON，重算 Hmix，生成最终 `SqsMixingData`；
- `status`/`validate`：检查 manifest、任务状态和数据契约。

端点 cache key 至少包括 phase、element、supercell、template checksum、POTCAR
functional 和计算设置 fingerprint。缺失端点或 fingerprint 不匹配时，`points`
阶段必须在启动 VASP 前失败。

### 6.9 `main_sqs.py`

入口脚本只负责阶段调用和 Slurm 提交，不散落 atomate2 task document 处理逻辑。日常
接口采用易记的 positional command；`--name/--batch` 和动作旗标继续兼容：

```bash
python main_sqs.py generate
python main_sqs.py endpoints --slurm
python main_sqs.py points --slurm
python main_sqs.py collect
python main_sqs.py status
python main_sqs.py validate

# 单个体系：命令后追加 --name NAME
python main_sqs.py status --name BCC-Al-Nb

# 旧 flag 形式仍可用
python main_sqs.py --name BCC-Al-Nb --status
```

`main_sqs.py` 不把四个阶段合并为一个 `run_tick`。正式 Slurm 生产按阶段提交，便于
在端点阶段完成后检查、重试单个任务，再开放 points 阶段。

### 6.10 第一阶段工程实现状态

第一阶段已在当前 `feat-sqs-mix` 工作树落地，范围限定为不执行计算的基础模型层：

- `htvasp/model/templates.py`：本地模板记录、快照、结构 checksum 和 JSON registry；
- `htvasp/model/sqs.py`：`SqsSystem`、`SqsCompositionPoint`、`SqsPlan`、
  `SqsStructureArtifact` 和 `SqsCalculationSpec`；
- `tests/test_sqs_plan.py`：名称规范化、整数成分、模板校验、manifest round-trip 和
  endpoint/point 稳定 ID 测试；
- `htvasp/model/__init__.py`：导出基础模型，并延迟加载历史 `Endmember`，使模型层
  的离线检查不要求 Materials Project API key；
- `.gitignore`：为 `htvasp/model/templates.py` 增加精确的 unignore 规则，避免被历史
  `*temp*` 规则误忽略。

该阶段明确不包含 `sqsgenerator` 适配、VASP、atomate2/jobflow、Worker、Maker、
Slurm、Static task document 转换和最终 Hmix schema。模板 registry 当前要求显式
注册记录或 JSON registry；从现有 `data/poscars` 自动发现模板不属于本阶段承诺。

### 6.11 第二阶段工程实现状态

第二阶段已在当前工作树落地结构级计算层，范围是构造并验证一个可被既有
`Worker.run_flow` 执行的 atomate2/jobflow Flow；不包含 SQS 结构生成、批量阶段编排、
Slurm 提交或 Hmix 收集：

- `htvasp/sqs/maker.py`：`SqsCalculationMaker` 组合
  `RelaxMaker -> StaticMaker -> sqs observation`。Relax 使用 `ISIF=3`；Static 使用
  `ISTART=0`、`ICHARG=2`、`ALGO=Normal`、`NELM=200` 做干净的自洽计算，不复制
  Relax 的 `WAVECAR`。Static 的输入结构和阶段目录引用仍分别来自 Relax 的
  `output.structure` 与 `output.dir_name`；
- `htvasp/sqs/jobs.py`：`make_static_observation` 是结果标准化 job，只消费
  `static_job.output`，兼容字典、Pydantic model 和嵌套 task document；它拒绝失败状态、
  缺失最终结构、缺失能量及非有限能量，并生成稳定的 task-document SHA-256；
- `htvasp/sqs/task_doc.py`：`SqsStaticObservation` 校验成功状态、总能量、每原子能量、
  站点数、最终结构和 checksum，并检查总能量与每原子能量的一致性；
- `htvasp/workflows/sqs_mixing.py`：`SqsMixingWorker` 负责计算层用户接口，支持以
  `SqsCalculationSpec` 或 calculation ID 标识任务，复用基类的 JSONStore、resume 和
  `get_result` 机制；
- `tests/test_sqs_workflow.py`：覆盖 Flow job 顺序、`OutputReference` 结构/目录交接、
  WAVECAR handoff、observation 标准化、task-document 形状 smoke、失败输入拒绝、
  calculation spec 传递和结果 job 名称。

本阶段的验收结果：SQS 专项测试共 24 项通过；目标模块 `compileall` 和
`git diff --check` 通过。
测试未运行 VASP 或 Slurm，因此这里只证明 Flow 构造、输入引用和结果契约，不宣称已有
新的生产计算结果。

### 6.12 第三阶段工程实现状态

第三阶段已落地 SQS 生成与批次编排的第一部分，边界限定为离线准备和阶段门禁，不启动
VASP、不提交 Slurm：

- `htvasp/sqs/generator.py`：定义 `SqsGenerator` 协议和 `GeneratedSqs` 结果；
  `SqsgeneratorAdapter` 使用项目 Conda 环境中的公开 `sqsgenerator.optimize` API，输入
  必须是已展开超胞，输出必须通过站点数和元素计数校验；记录 engine/version、iterations、
  seed、seed 是否实际生效、配置、objective 和结构 checksum。当前版本没有文档化的随机
  seed 配置项，因此只记录 seed，不宣称可复现；
- `htvasp/sqs/orchestrator.py`：定义文件化的 `SqsBatchManifest`、`SqsSystemRecord`、
  `SqsEndpointRecord` 和 `SqsOrchestrator`。`generate` 生成 endpoint/point POSCAR、系统
  manifest 和批次 manifest；endpoint cache key 包含 phase、element、supercell、模板
  checksum 和计算设置 fingerprint；共享端点在批次内去重；已有 batch 默认拒绝重建，只有
  显式 `overwrite=True` 才允许覆盖；
- points 门禁要求所有 endpoint 状态为 successful，并检查 endpoint observation 文件存在、
  task document schema、站点数、元素组成、artifact checksum 和 settings fingerprint，
  任一条件不满足时在启动 points 前失败；
- `tests/test_sqs_orchestrator.py`：覆盖两个体系共享 endpoint 的去重、manifest round-trip、
  endpoint cache key 变化、缺失 artifact/observation 检测、overwrite 保护、endpoint 未完成
  时阻止 points，以及全部 endpoint 成功后开放 points。

本阶段验收结果：SQS 三组测试共 20 项通过；另用 `sqsgenerator 0.5.6` 对 16-site BCC、
`Al4Nb12`、2 iterations 做真实生成 smoke，输出站点数、组成和 checksum 校验通过。
本阶段的离线边界不启动 VASP 或提交 Slurm；CLI、Worker 批量执行和最终收集在后续阶段
单独验收。

### 6.13 第四阶段工程实现与真实 endpoint 验证状态

第四阶段已完成软件入口和数据契约的落地，但真实生产闭环尚未完成：

- `main_sqs.py` 已统一为 positional command 和兼容 flag 的 portal 入口：
  `generate`、`endpoints`、`points`、`collect`、`status`、`validate`；旧的
  `--generate`、`--endpoints`、`--points`、`--collect`、`--status`、`--validate` 仍可用。
  endpoint/point 阶段支持 `--slurm` 和 `--dry-run`。无 `--name` 时默认批量处理六个配置体系；
  计算 recipe 和 Slurm 资源全部在脚本内部维护；
- `htvasp/model/sqs_schema.py` 已实现 `htvasp-sqs-mixing/v1`，包括 BCC_A2 单亚晶格
  映射、Static observation、Hmix 重算和 eV/atom 到 J/mol-atom 的换算；
- 历史实验目录 `data/sqs_mixing/BCC-Al-Nb-v2` 已完成一次真实 Slurm endpoint smoke
  （作业 `53188`）；正式流程默认使用 `data/endmembers/{structure_name}/sqsflow/`。
  `BCC-Nb` 的 Relax、Static 和 observation 成功；`BCC-Al` 的 Relax 成功，但 Static
  因 `BRMIX: very serious problems` 及不完整 `vasprun.xml` 失败。失败 JSON、
  `store.json` 和运行目录均保留，未覆盖旧证据；
- `v2` 失败分析确认 Relax 的 `WAVECAR` 不适合作为本 Static 的默认 restart：Relax 和
  Static 的 ENCUT、晶格及 KPOINTS 不一致，且 Nb 曾出现 `while reading WAVECAR`。因此
  当前实现已改为 Static `ISTART=0, ICHARG=2, ALGO=Normal, NELM=200`，不复制 WAVECAR，
  并通过新的 settings fingerprint 隔离后续重跑；
- 修复后的 SQS 专项回归为 30 项通过；全量回归为 72 项通过、1 项已有的 NSCF 基线失败
  （`tests/test_nscf.py::test_flow_makers` 期待 4 个 maker，而当前实现为 2 个）。

白天只完成新 batch 的离线生成、`status/validate` 和配置审阅，不进行长期轮询、
不提交新的 Slurm 作业。endpoint 全部成功前，`points` 阶段保持门禁关闭。

## 7. JSON 数据边界

### 7.1 端点 JSON

端点 JSON 是 `endpoints` 到 `points/collect` 的轻量边界，至少保存：

- `endpoint_id`、phase、element、supercell；
- template checksum、settings fingerprint；
- Relax 前后结构；
- static 总能量、每原子能量、站点数和收敛状态；
- jobflow/task document 路径或 checksum；
- 生成时间和计算 provenance。

### 7.2 最终 SQS mixing JSON

最终文件使用：

```text
schema_version = "htvasp-sqs-mixing/v1"
```

核心内容为：

- 明确 `calphad_mapping`；第一阶段为 `BCC_A2`、单亚晶格、`parameter_target=G`；
- 两个纯元素 reference；
- `x_A=0.25, 0.50, 0.75` 三个成功点；
- 每点 `site_counts`、初始/最终结构、static observation；
- `h_mix_eV_per_atom` 和 `h_mix_J_per_mol`；
- SQS engine/version、迭代预算、seed（若工具不真正支持 seed，必须明确注明）、
  objective 和 structure checksum；
- relax/static 参数和整体 settings fingerprint。

能量换算固定为：

```text
1 eV/atom = 96485.33212331002 J/mol-atom
```

`em-tdb` 只能消费该 JSON，不读取 POSCAR、`store.json`、OUTCAR、历史 `out*.tdb`
或 ESPEI 生产 TDB。它应从 JSON 重算 Hmix，数值不一致就拒绝拟合。

## 8. 目前的风险和未完成项

### 8.1 已知风险

- `sqsgenerator` 的公开 API 不保证 seed 真正控制随机性；必须保存生成结构、工具
  版本、配置和 checksum，不能仅声称 seed 可复现；
- SQS 是有限超胞的一次无序实现，正式科学研究可能需要多个构型平均和误差估计；
- `ISPIN=2` 和初始 MAGMOM 已保留，但弱磁/非磁端点的磁态稳定性仍需单独审阅；
- VASP 的 XML 可能因作业中断而不完整，不能只依赖 XML validator 的报错文本；
- FCC/HCP 可以生成结构和运行 VASP，但当前 `em-tdb` P2 mapping 只定义 BCC_A2；
- 低阶 `L0+L1` 只适合当前第一阶段，不能把三点拟合结果解释成完整热力学模型。

### 8.2 尚未完成

1. ~~补齐本地模板快照和基础 SQS 领域模型~~（已完成，见 6.10）；
2. ~~实现 `SqsCalculationMaker`、observation job 和 `SqsMixingWorker`~~（已完成，见 6.11）；
3. ~~建立 endpoint 跨体系去重和 fingerprint 校验~~（generate 与 points 门禁已完成，见 6.12）；
4. ~~建立 fake task document 的离线测试，不依赖 VASP~~（计算层和阶段门禁已完成）；
5. ~~实现 `generate` 阶段及 `sqsgenerator` 适配器~~（已完成，见 6.12）；
6. ~~建立 `endpoints`、`points`、`collect`、`status` 和 `validate` 的正式 CLI/执行入口~~（已完成，见 6.13）；
7. ~~对 `BCC-Al-Nb` 做新阶段式入口的 dry-run~~（已完成；真实 endpoint smoke 见 6.13）；
8. 用新 fingerprint 重跑两个 endpoint；成功后再提交三个 points；
9. 实现 `collect` 的真实成功验收，并确认 `htvasp-sqs-mixing/v1` 可被 `em-tdb` 消费；
10. 加入 TDB `G/L` 写出和 pycalphad 回读测试；
11. 再决定是否扩展 FCC/HCP、更多成分点和多 SQS 统计平均。

### 8.3 数据目录约定

SQS-mixing 与现有端基工作流统一使用 `data/endmembers/{structure_name}/` 作为
structure-level 根目录：

```text
{structure_name}/
├── qhaflow/
├── staticflow/
├── tb2jflow/
├── sqsflow/
├── {structure_name}-qha.json
├── {structure_name}-static.json
├── {structure_name}-tb2j.json
└── {structure_name}-sqs.json
```

其中 `sqsflow/` 保存 `batch.json`、系统 manifest、SQS POSCAR、observation、失败记录
和运行目录；`collect` 的最终输出保存为同级的 `{structure_name}-sqs.json`。正式入口只
使用上述固定布局。历史 `data/sqs_mixing/` 仅作为旧实验和失败证据，不再被新任务读取。

### 8.4 推荐入口

`main_sqs.py` 与 `main_qha.py`、`main_w90.py` 一样，将计算 recipe、Slurm 资源和生产
结构列表保存在脚本内部。日常命令不需要再传递 INCAR、POTCAR、KPOINTS、SQS engine
或 scratch 目录参数。

当前默认批量列表是已经有 `BCC_A2` CALPHAD 映射的 Al-Nb-Ti-V 六个二元体系：

```text
BCC-Al-Nb  BCC-Al-Ti  BCC-Al-V
BCC-Nb-Ti  BCC-Nb-V   BCC-Ti-V
```

推荐命令：

```bash
# 生成全部默认体系的 SQS 结构；单个体系追加 --name NAME
python main_sqs.py generate
python main_sqs.py generate --name BCC-Al-Nb

# 先做离线检查，再批量提交 endpoint 阶段
python main_sqs.py validate
python main_sqs.py endpoints --slurm

# endpoint 全部成功后，再批量提交三个 SQS 成分点
python main_sqs.py status
python main_sqs.py points --slurm

# points 全部成功后，批量写出各体系的 {name}-sqs.json
python main_sqs.py collect
```

单个体系可以用 `--name` 选择；不提供 `--name` 时默认批量。提交前可使用
`--slurm --endpoints --dry-run` 或 `--slurm --points --dry-run` 检查提交对象；dry-run
不会调用 `sbatch`。`--rerun` 只用于明确的重建/重跑，不应在已有生产证据上随意使用。

`--root`、模板/超胞/INCAR/POTCAR/KPOINTS 和 Slurm 资源参数均不接受；旧动作 flag
与 positional command 的语义相同。旧 manifest 若缺少当前 `recipe_version` 或 fingerprint 不一致，会被拒绝；必须使用
`--rerun` 在正式目录重新生成后才能继续计算。旧的 `data/sqs_mixing/` 数据不会被新入口
读取。

### 8.5 `sof-ace` 历史逻辑提取与职责边界确认

本轮对 [`sof-ace/scripts/fit_l_parameters_sqs.py`](/nfs_hdd/2025/gaominliang/sof-ace/scripts/fit_l_parameters_sqs.py)、
[`sof-ace/scripts/fit_l_parameters.py`](/nfs_hdd/2025/gaominliang/sof-ace/scripts/fit_l_parameters.py)
及相关 CALPHAD 文档进行了只读核对。历史实现包含以下科学步骤：

1. 从 SQS 结构目录和 VASP 结果目录建立映射。旧脚本通过 POSCAR 的结构 fingerprint
   找到 `result.json`，并从结果中读取 total energy；端点能量使用脚本内硬编码的
   `REF_EV_PER_ATOM`，超胞站点数也固定为 16。
2. 计算静态混合焓。对二元 A-B 体系，历史公式为

   ```text
   H_mix(eV/atom) = [E_SQS - n_A E_A - n_B E_B] / N
   H_mix(J/mol-atom) = H_mix(eV/atom) × 96485.332...
   ```

3. 对每个二元体系使用加权最小二乘拟合低阶参数。历史实现的形式为

   ```text
   H_mix(x) = x²(1-x)² [L0 + L1(2x-1)]
   ```

   权重为 `x²(1-x)²`，并输出 RMSE、标准误和拟合点信息。
4. 将拟合结果写入 TDB，并用 pycalphad 读取验证。历史脚本还需要删除旧的 `VVxxxx`
   参数和函数，避免旧 CE/MC 参数与 SQS 参数叠加。

这些内容证明了从 SQS DFT 能量到 BCC 无序相 Redlich-Kister 参数的科学需求和公式路径，
但历史脚本把目录扫描、能量来源、拟合、TDB 文本改写和 pycalphad 验证耦合在一起，不能
作为 `ht-vasp` 的生产接口。尤其不能把硬编码端点能量、固定 16-site 假设、旧 `VVxxxx`
清理逻辑或 TDB 写出逻辑迁移到 VASP 生产者中。

职责边界最终确认如下：

| 项目 | 当前/目标职责 | 明确不负责的内容 |
| --- | --- | --- |
| `ht-vasp` | 生成 SQS；运行 `RelaxMaker -> StaticMaker`；校验 task document、结构、站点数和 fingerprint；输出 `htvasp-sqs-mixing/v1` DFT 数据 JSON | L0/L1 拟合、TDB 写出、pycalphad 相图和 TDB 回读 |
| `em-tdb` | 负责读取 `htvasp-sqs-mixing/v1`，重算/检查 Hmix，拟合 Redlich-Kister 参数，写入并验证 TDB | 扫描 `ht-vasp` VASP 目录、猜测 POSCAR/OUTCAR 能量、消费历史 `sof-ace` 目录 |
| `sof-ace` | 保留历史 CE/MC/SQS 探索和方法参考 | 不作为新 `em-tdb` 生产 JSON 输入 |

截至本次核对，`em-tdb` 已有 `emtdb/sqs_mixing.py` 和 `main.py sqs` 的 P2-A 部分实现：
它能够读取 `htvasp-sqs-mixing/v1`、校验端点/Hmix，并按标准无序 RK 形式
`x(1-x)[L0+L1(2x-1)]` 计算拟合结果和独立 TDB fragment。可是该实现当前硬编码要求
exactly five 个 SQS composition points，而 `ht-vasp` 正式 recipe 是
`x_A=0.25, 0.50, 0.75` 三点；因此它还不能消费当前正式生产 JSON，不能视为已完成的
生产闭环。尚待 `em-tdb` 将输入点数约定与 v1 对齐，并补齐与最终 TDB/回归验证的正式
集成。当前 SQS 固化的有效性结论必须限定为：

> `ht-vasp` 的 DFT 数据生产边界已经有效固化；`em-tdb` 已有但尚未与三点 v1 契约对齐的
> P2-A 拟合草稿，完整的 SQS → CALPHAD 生产闭环尚未完成，后续实现归属 `em-tdb`。

### 8.6 `ht-vasp` DFT JSON 生产边界验收

当前实现满足以下生产边界：

- `generate` 只生成本地模板超胞、SQS POSCAR、provenance、checksum 和 manifest，不启动
  VASP；
- endpoint/point 的 Worker 只执行 `RelaxMaker -> StaticMaker -> observation`，Static
  总能量来自成功的 Static task document；
- `collect` 只读取已校验的 Static observation 和初始/最终结构，重算并保存 Hmix 数据，
  不执行参数拟合；
- `htvasp-sqs-mixing/v1` 保存端点能量、SQS 点能量、结构、站点计数、生成器 provenance、
  计算 fingerprint、Hmix 和明确的 BCC_A2 单亚晶格语义元数据；
- schema 校验可以从点能量和端点每原子能量复算 `h_mix_eV_per_atom`，并检查
  `J/mol-atom` 换算、站点数和 fingerprint 一致性；
- 没有 `L0`、`L1` 拟合结果字段，也没有 TDB 输出副作用。

因此，`ht-vasp` 当前可验收为“可审计 DFT JSON 生产者”，但 `{structure_name}-sqs.json`
不是已经拟合完成的 TDB 参数文件；它是供 `em-tdb` 后续消费的输入数据集。

## 9. 验收标准

新的 SQS-mixing 固化完成前，至少应满足：

- 反向名称不会重复生成体系；
- 同一元素端点在一个批次中只计算一次；
- 端点缺失或 fingerprint 不匹配时，`points` 阶段在启动 VASP 前失败；
- 每个任务都有明确的 Relax/Static 成功状态和结构/能量 checksum；
- 只能从 static 总能计算 Hmix，不能混用 relax 能量；
- JSON round-trip 后能重算 Hmix，单位和站点数一致；
- `em-tdb` 不需要扫描任何 ht-vasp 目录；
- 旧五点目录、旧参数和 ESPEI/历史 TDB 不会被新流程静默消费；
- 至少完成 `BCC-Al-Nb` 的离线 generate/dry-run/collect smoke，再考虑生产 Slurm。

第一阶段基础模型工程的额外验收标准：

- 本地模板读取不联网，结构 checksum 稳定且可校验；
- `SqsSystem` 能规范化正反向体系名称；
- `SqsCompositionPoint` 拒绝非整数、非正数和不匹配的站点数；
- manifest JSON round-trip 后计划、模板引用和成分点不变；
- endpoint/point calculation spec 具有稳定、可审计的 ID；
- 基础模型模块不依赖 VASP、Slurm、`sqsgenerator` 或 VASP 输出文件；
- 离线测试可在没有 VASP、Materials Project API 和 SQS 引擎的环境中运行。

## 10. 相关入口

- 当前新分支：`ht-vasp` `feat-sqs-mix`，HEAD `0108ef7`；
- 历史原型分支：`feat-sqs-mixing`，SQS 实现提交 `cf93be6`；
- 原始 SQS 契约（历史）：提交 `7f83d8d` 的 `docs/SQS混合焓数据契约-v1.md`；
- SQS 参数说明：[`docs/sqsgenerator 参数配置指南.md`](</nfs_hdd/2025/gaominliang/ht-vasp/docs/sqsgenerator%20参数配置指南.md>)；
- CALPHAD/L 参数历史实现：[`sof-ace/scripts/fit_l_parameters_sqs.py`](/nfs_hdd/2025/gaominliang/sof-ace/scripts/fit_l_parameters_sqs.py)；
- 后续消费方：[`em-tdb/CLAUDE.md`](/nfs_hdd/2025/gaominliang/em-tdb/CLAUDE.md)。
