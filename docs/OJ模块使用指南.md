# OJ 模块使用指南

## 快速开始

```python
from pymatgen.core import Structure
from htvasp.workflows import OJWorker

structure = Structure.from_file("POSCAR")
worker = OJWorker(
    worker_name="fe_co",
    j_count=3,
    magnetic_ion_types=["Fe", "Co"],
    vasp_args={"vasp_cmd": "mpirun -np 48 vasp_std"},
)
worker.run_flow("FeCo", structure, "./oj_work")
```

## 核心组件

| 组件 | 位置 | 功能 |
|------|------|------|
| OJWorker | workflows/oj.py | 用户接口 |
| OJMaker | oj/maker.py | Flow 生成 |
| OJInputSetGenerator | oj/input_set.py | 输入文件生成 |
| oj_generate | oj/jobs.py | 生成磁性构型 |
| create_flip_jobs | oj/jobs.py | 并行执行 VASP |
| oj_solve | oj/jobs.py | 解算 J 参数 |
| OJResult | oj/task_doc.py | 结果模型 |

## OJResult 字段

- `J_reprs`: J 对表示 ["Fe-Fe-1", ...]
- `Js`: 交换参数 (eV)
- `Tc_MFA`, `Tc_RPA`: 居里温度 (K)
- `E_DLM`: DLM 能量
- `D_stiff`: 旋波刚度
- `condition_number`: 条件数
- `is_complete`: 解是否完整
- `configs_data`: 各配置详细数据