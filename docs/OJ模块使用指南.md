# OJ 模块使用指南

## 概述

`htvasp.workflows.oj` 模块提供了 OstravaJ 磁性交换相互作用计算的工作流支持。该模块遵循 atomate2 的设计模式，支持灵活的工作流组合。

## 架构

```
htvasp/oj/
├── __init__.py      # 模块入口
├── config.py        # OJConfig 配置类
├── input_set.py     # OJ 输入文件生成
├── jobs.py          # @job 任务函数
├── maker.py         # OJMaker, OJSimpleMaker, OJRelaxMaker
└── task_doc.py      # 输出结果模型

htvasp/workflows/
└── oj.py            # OJWorker 用户接口
```

## 快速开始

### 方式1: 使用 OJWorker (推荐用于简单场景)

```python
from pymatgen.core import Structure
from htvasp.workflows import OJWorker
from htvasp.oj import OJConfig

structure = Structure.from_file("POSCAR")
config = OJConfig(j_count=3, magnetic_ion_types=["Fe", "Co"])

worker = OJWorker(
    worker_name="fe_co_exchange",
    config=config,
    global_incar={"ENCUT": 520, "ISPIN": 2},
    vasp_args={"vasp_cmd": "mpirun -np 48 vasp_std"},
)

result = worker.run_flow(
    name="FeCo",
    structure=structure,
    flow_dir="./oj_work",
    resume=True,  # 支持续算
)
```

### 方式2: 使用 OJMaker (推荐用于 Flow 组合)

```python
from pymatgen.core import Structure
from htvasp.oj import OJMaker, OJConfig
from jobflow import run_locally

structure = Structure.from_file("POSCAR")
config = OJConfig(j_count=3, magnetic_ion_types=["Fe", "Co"])

maker = OJMaker(config=config, vasp_cmd="mpirun -np 48 vasp_std")
flow = maker.make(structure)

# 使用 jobflow 执行
run_locally(flow, root_dir="./oj_work")
```

### 方式3: 使用 @job 函数自定义 Flow

```python
from htvasp.oj import oj_generate, oj_vasp, oj_solve
from jobflow import Flow, run_locally

gen_job = oj_generate(structure, config)
vasp_job = oj_vasp(gen_job.output["flip_dirs"], vasp_cmd="vasp_std")
solve_job = oj_solve(gen_job.output["run_dir"], vasp_job.output["results"])

flow = Flow([gen_job, vasp_job, solve_job])
run_locally(flow, root_dir="./oj_work")
```

## 与 atomate2 Maker 组合

### Relax + OJ 组合

```python
from atomate2.vasp.flows.core import DoubleRelaxMaker
from htvasp.oj import OJRelaxMaker, OJConfig

config = OJConfig(j_count=3, magnetic_ion_types=["Fe"])

maker = OJRelaxMaker(
    relax_maker=DoubleRelaxMaker(),
    config=config,
    vasp_cmd="mpirun -np 48 vasp_std",
)

flow = maker.make(structure)
```

### 手动组合多个 Maker

```python
from atomate2.vasp.jobs.core import StaticMaker
from htvasp.oj import OJMaker, OJConfig
from jobflow import Flow

static_maker = StaticMaker()
oj_maker = OJMaker(config=config)

static_job = static_maker.make(structure)
oj_flow = oj_maker.make(static_job.output.structure)

combined_flow = Flow([static_job, *oj_flow.jobs])
```

## 配置参数

### OJConfig

| 参数 | 类型 | 说明 |
|------|------|------|
| `j_count` | int | J 对数量 (与 dist_cutoff 二选一) |
| `dist_cutoff` | float | 距离截断 (Å) (与 j_count 二选一) |
| `magnetic_ion_types` | list[str] | 磁性离子类型 |
| `noncollinear` | bool | 是否使用非共线磁性 |
| `base_spin` | float / list | 基础自旋值 |
| `extend_poscar` | tuple | 超胞扩展因子 |
| `kppa` | float | K点密度 (默认为 5000) |
| `incar` | dict | INCAR 参数设置 |

### INCAR 参数自定义

通过 `global_incar` 参数覆盖默认 INCAR 设置：

```python
worker = OJWorker(
    worker_name="fe",
    config=OJConfig(j_count=3, magnetic_ion_types=["Fe"]),
    global_incar={
        "ENCUT": 520,
        "ISPIN": 2,
        "MAGMOM": "4*5.0 8*0.0",  # 自定义磁矩
        "ISMEAR": 0,  # 设置为高斯展宽
    },
)
```

## 输出

### OJTaskDocument

```python
from htvasp.oj.task_doc import OJTaskDocument

doc = OJTaskDocument.from_oj_solution(solution)
print(doc.J_reprs)  # ["Fe-Fe-1", "Fe-Co-1", ...]
print(doc.Js)       # [10.5, -5.2, ...]
print(doc.Tc_MFA)   # 300.0
print(doc.Tc_RPA)   # 280.0
```

### 输出结果格式

典型输出结果包含以下字段：

```python
{
    "J_reprs": ["Fe-Fe-1", "Fe-Fe-2"],  # J 对表示
    "Js": [10.5, -2.3],  # 交换参数值 (meV)
    "Tc_MFA": 350.2,  # 居里温度 (平均场近似)
    "Tc_RPA": 320.5,  # 居里温度 (RPA近似)
    "run_dir": "/path/to/run/dir",  # 运行目录
    "num_configs": 10,  # 生成的磁性构型数量
    "vasp_success_rate": 1.0,  # VASP 计算成功率
}
```

## Slurm 支持

通过在 `vasp_args` 中设置 Slurm 相关参数：

```python
worker = OJWorker(
    worker_name="fe_co",
    config=config,
    vasp_args={
        "vasp_cmd": "srun vasp_std",  # 使用 srun 启动 VASP
        # Slurm 作业参数可在外部脚本中设置
    },
)

result = worker.run_flow(
    name="FeCo",
    structure=structure,
    flow_dir="./oj_work",
    resume=True,
)
```

**注意**：实际的 Slurm 作业提交通常通过外部脚本或集群管理系统进行，`vasp_cmd` 应设置为适合集群环境的命令。

## 典型目录结构

OJ 工作流执行后的目录结构如下：

```
./oj_work/
├── store.json          # Jobflow 存储文件
├── 1-oj-Fe_oj_generate/  # 第1步：生成磁性构型
│   ├── flip000/        # 第1个磁性构型
│   ├── flip001/        # 第2个磁性构型
│   └── ...
├── 2-oj-Fe_oj_vasp/    # 第2步：VASP 计算
└── 3-oj-Fe_oj_solve/   # 第3步：求解交换参数
    └── OJ_solution.json  # 最终结果文件
```

## 工作流程

1. **生成磁性构型**：调用 `oj_generate` 生成指定数量的磁性构型（flip 目录）
2. **VASP 计算**：调用 `oj_vasp` 串行执行每个磁性构型的 VASP 计算
3. **求解交换参数**：调用 `oj_solve` 计算 J 参数和居里温度 Tc

## 与旧版对比

| 特性 | 旧版 oj.py | 新版 oj 模块 |
|------|-----------|-------------|
| 架构 | 单文件混合 | 模块化分离 (htvasp/oj/) |
| Maker | 无 | OJMaker, OJSimpleMaker, OJRelaxMaker |
| Job | 定义但未用 | 完整集成 (oj_generate, oj_vasp, oj_solve) |
| Flow 组合 | 不支持 | 支持 |
| 续算功能 | 独立实现 | 统一使用 run_locally_custom(resume=True) |
| INCAR 自定义 | 不支持 | 支持 global_incar 参数 |
| 目录结构 | 冗余复杂 | 按 Job 拆分，结构清晰 |
| J_reprs 处理 | 列表作为 dict keys | 转换为字符串格式 |
| 测试 | 无 | 完整的单元测试 |
