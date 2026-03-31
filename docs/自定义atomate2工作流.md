# 如何开发一个新的 atomate2 工作流

> 一个 atomate2 计算工作流的解剖结构（即，你需要写什么？）

## 概述

每个 atomate2 工作流都是 jobflow 的 `Flow` 类的实例，它是一个 `Job` 和/或其他 `Flow` 对象的集合。因此，**你的最终目标是生成一个 `Flow`**。

## 核心概念

### 1. Maker 类

`Maker` 是创建工作流的主要入口点，它包含一个工厂方法 `make()`，根据输入生成 `Flow`。

```python
class ExampleMaker(Maker):
    def make(self, structure: Structure) -> Flow:
        # 根据输入坐标返回一个 Flow
        return Flow([job1, job2, ...])
```

通常，`Maker` 类包含大部分计算参数和其他设置，这些参数用于正确设置计算。

### 2. Job 与 @job 装饰器

大部分逻辑可以像普通 Python 函数一样编写，然后通过 `@job` 装饰器转换为 `Job`：

```python
from jobflow import job

@job
def my_calculation(structure, params):
    # 计算逻辑
    return result
```

### 3. InputGenerator / InputSet

编写计算输入文件通常通过 `pymatgen InputSet` 类完成。`InputSet` 是一个类似字典的容器，指定需要写入的文件及其内容。

通过 `InputGenerator` 类创建 `InputSet`：

```python
class ExampleInputGenerator(InputGenerator):
    def get_input_set(self, structure: Structure) -> InputSet:
        return InputSet(...)
```

在 atomate2 中，将 `InputGenerator` 作为 Maker 的类参数使用：

```python
class ExampleMaker(Maker):
    input_set_generator: ExampleInputGenerator = field(
        default_factory=ExampleInputGenerator
    )

    def make(self, structure: Structure) -> Flow:
        input_set = self.input_set_generator.get_input_set(structure)
        input_set.write_inputs()  # 写入输入文件
        return Flow([...])
```

### 4. TaskDocument

大多数 atomate2 工作流以**任务文档**的形式返回结构化输出。任务文档是 emmet 的 `BaseTaskDocument` 类实例，定义存储计算输出的模式。

已有模式包括：
- `MaterialsDoc`: 固体材料计算数据
- VASP、Q-Chem、FEFF 等的专用模式

也可以通过 cclib 解析其他代码的输出。

---

## 新建 atomate2 工作流的组成

| 组件 | 说明 | 必需 |
|------|------|------|
| `Maker` | 实际生成工作流 | ✅ |
| `Job` / `Flow` | 定义工作流中的离散步骤 | ✅ |
| `InputGenerator` | 生成输入文件 | 可选 |
| `TaskDocument` | 定义输出数据模式 | 可选 |

---

## 代码组织结构

由于 MP 软件生态系统的分布式设计，编写完整的新工作流可能涉及多个 GitHub 仓库：

### 工作流代码 (Job, Flow, Maker)
→ 放在 **atomate2** 仓库

### InputSet / InputGenerator
→ 放在 **pymatgen** 仓库

如果需要从头创建（比如处理 pymatgen 尚未支持的代码），可以先放在 atomate2 中快速迭代，成熟后再迁移。

### TaskDocument 模式
→ 优先检查 **emmet** 是否有匹配的模式
→ 其次检查 **cclib**（可通过 `atomate2.common.schemas.TaskDocument` 导入）
→都没有则新建于 **atomate2** 或 **cclib**

---

## 最小示例

```python
from atomate2 import Maker
from jobflow import Flow, job
from pymatgen.core import Structure

class MyWorkflowMaker(Maker):
    name = "my_workflow"

    def make(self, structure: Structure) -> Flow:
        # 步骤1: 准备输入
        prep_job = self._prepare_inputs(structure)

        # 步骤2: 执行计算
        calc_job = self._run_calculation(prep_job.output)

        # 步骤3: 后处理
        post_job = self._post_process(calc_job.output)

        return Flow([prep_job, calc_job, post_job])

    @job
    def _prepare_inputs(self, structure):
        # 生成输入文件
        return {"ready": True}

    @job
    def _run_calculation(self, inputs):
        # 执行计算
        return {"energy": -100.0}

    @job
    def _post_process(self, results):
        # 处理结果
        return {"final_result": results["energy"]}
```

---

## 与 htvasp OJWorker 的对比

| 概念 | atomate2 方式 | htvasp OJWorker 方式 |
|------|---------------|---------------------|
| 工作流组装 | `Flow` + `Maker` | `OJWorker.run_flow()` |
| Job 定义 | `@job` 装饰器 | 普通方法 |
| 输入生成 | `InputGenerator` | `OJWorker._write_incar()` |
| 执行模式 | 分布式 Job | 单进程同步 (`subprocess`) |

OJWorker 不完全符合 atomate2 的设计哲学，但实现了类似的功能结构，且更轻量、解耦更好。
