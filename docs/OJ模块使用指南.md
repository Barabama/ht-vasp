# OJ 模块使用指南

## 概述

`htvasp.workflows.oj` 模块提供了 OstravaJ 磁性交换相互作用计算的工作流支持。该模块遵循 atomate2 的设计模式，支持灵活的工作流组合。

## 架构

```
htvasp/workflows/oj/
├── __init__.py      # 模块入口
├── config.py        # OJConfig, OJIncarSettings
├── input_set.py     # OJInputSet, OJInputSetGenerator
├── jobs.py          # @job 任务函数
├── maker.py         # OJExchangeMaker, OJRelaxExchangeMaker
├── task_doc.py      # OJTaskDocument, OJSolutionSummary
└── worker.py        # OJWorker (兼容层)
```

## 快速开始

### 方式1: 使用 OJWorker (推荐用于简单场景)

```python
from pymatgen.core import Structure
from htvasp.workflows import OJWorker, OJConfig

structure = Structure.from_file("POSCAR")
config = OJConfig(j_count=3, magnetic_ion_types=["Fe", "Co"])

worker = OJWorker(
    worker_name="fe_co_exchange",
    oj_config=config,
    vasp_args={"vasp_cmd": "mpirun -np 48 vasp_std"},
)

result = worker.run_flow(
    name="FeCo",
    structure=structure,
    flowdir="./oj_work",
)
```

### 方式2: 使用 OJExchangeMaker (推荐用于 Flow 组合)

```python
from pymatgen.core import Structure
from htvasp.workflows.oj import OJExchangeMaker, OJConfig
from jobflow import run_locally

structure = Structure.from_file("POSCAR")
config = OJConfig(j_count=3, magnetic_ion_types=["Fe", "Co"])

maker = OJExchangeMaker(oj_config=config, vasp_cmd="mpirun -np 48 vasp_std")
flow = maker.make(structure, workdir="./oj_work")

# 使用 jobflow 执行
run_locally(flow)
```

### 方式3: 使用 @job 函数自定义 Flow

```python
from htvasp.workflows.oj import oj_generate_job, oj_vasp_batch_job, oj_solve_job
from jobflow import Flow, run_locally

gen_job = oj_generate_job(structure, config, workdir="./oj_work")
vasp_job = oj_vasp_batch_job(gen_job.output["flip_dirs"], vasp_cmd="vasp_std")
solve_job = oj_solve_job(gen_job.output["run_dir"])

flow = Flow([gen_job, vasp_job, solve_job])
run_locally(flow)
```

## 与 atomate2 Maker 组合

### Relax + OJ 组合

```python
from atomate2.vasp.flows.core import DoubleRelaxMaker
from htvasp.workflows.oj import OJRelaxExchangeMaker, OJConfig

config = OJConfig(j_count=3, magnetic_ion_types=["Fe"])

maker = OJRelaxExchangeMaker(
    relax_maker=DoubleRelaxMaker(),
    oj_config=config,
    vasp_cmd="mpirun -np 48 vasp_std",
)

flow = maker.make(structure, workdir="./relax_oj_work")
```

### 手动组合多个 Maker

```python
from atomate2.vasp.jobs.core import StaticMaker
from htvasp.workflows.oj import OJExchangeMaker, OJConfig
from jobflow import Flow

static_maker = StaticMaker()
oj_maker = OJExchangeMaker(oj_config=config)

static_job = static_maker.make(structure)
oj_flow = oj_maker.make(static_job.output.structure, workdir="./oj_work")

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

### OJIncarSettings

默认 INCAR 设置，可通过 `user_incar_settings` 覆盖。

## 输出

### OJTaskDocument

```python
from htvasp.workflows.oj import OJTaskDocument

doc = OJTaskDocument.from_oj_solution(solution)
print(doc.J_reprs)  # ["Fe-Fe-1", "Fe-Co-1", ...]
print(doc.Js)       # [10.5, -5.2, ...]
print(doc.Tc_MFA)   # 300.0
print(doc.Tc_RPA)   # 280.0
```

### OJSolutionSummary

轻量级输出摘要：

```python
from htvasp.workflows.oj import OJSolutionSummary

summary = OJSolutionSummary.from_dict(solution)
print(summary.to_dict())
```

## Slurm 支持

```python
worker = OJWorker(
    worker_name="fe_co",
    oj_config=config,
    vasp_args={
        "use_slurm": True,
        "vasp_cmd": "srun vasp_std",
        "ntasks": 48,
        "memory": "20G",
    },
)

# 提交后监控
result = worker.run_flow(
    name="FeCo",
    structure=structure,
    flowdir="./oj_work",
    monitor=True,  # 等待 Slurm 任务完成
)
```

## 与旧版对比

| 特性 | 旧版 oj.py | 新版 oj/ |
|------|-----------|----------|
| 架构 | 单文件混合 | 模块化分离 |
| Maker | 无 | OJExchangeMaker |
| Job | 定义但未用 | 完整集成 |
| Flow 组合 | 不支持 | 支持 |
| 输出 | dict | OJTaskDocument |
| 测试 | 无 | test_oj_refactor.py |
