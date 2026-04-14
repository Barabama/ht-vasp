# 自定义 atomate2 工作流

## 核心概念

### 1. Maker

```python
class MyMaker(Maker):
    def make(self, structure: Structure) -> Flow:
        return Flow([job1, job2])
```

### 2. @job 装饰器

```python
from jobflow import job

@job
def my_job(structure):
    return result
```

### 3. InputGenerator

```python
class MyInputGenerator(VaspInputGenerator):
    @property
    def incar_updates(self):
        return {"ENCUT": 500}

    def get_input_set(self, structure):
        return super().get_input_set(structure)
```

## 最小示例

```python
from atomate2 import Maker
from jobflow import Flow, job

class MyMaker(Maker):
    name = "my_flow"

    def make(self, structure):
        job1 = self._step1(structure)
        job2 = self._step2(job1.output)
        return Flow([job1, job2])

    @job
    def _step1(self, s):
        return {"ready": True}

    @job
    def _step2(self, inputs):
        return {"done": True}
```

## 组件对照

| 组件 | 说明 |
|------|------|
| Maker | Flow 工厂 |
| @job | 可执行步骤 |
| VaspInputGenerator | 输入文件生成 |
| TaskDocument | 输出模式 |