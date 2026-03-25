# HT-VASP 项目扩展规划 - CALPHAD 端基数据库构建

## 项目背景

基于第三代 CALPHAD 热力学数据库构建需求，扩展 HT-VASP 项目以支持：
1. **磁性参数计算**（高优先级）
2. **SQS 结构生成和混合焓计算**
3. **端基热力学性质计算**（基于现有 QHA 流程）

## 现有代码架构分析

### 核心模块理解

#### 1. `htvasp/workflows/qha.py` - QHA 计算流程
**核心功能**：
- 使用 `atomate2.vasp.flows.qha.QhaMaker` 进行准谐近似计算
- 包含 R3 结构优化（`initial_relax_maker`）
- 包含 EOS 静态计算（`eos_relax_maker`）
- 包含声子计算（`phonon_maker`）

**数据提取机制**：
```python
# 从 JobStore 查询数据
for doc in self.store.query(
    criteria={"name": {"$regex": r"phonon static eos deformation \d+"}},
    properties=["uuid", "index", "name"],
):
    output = self.store.get_output(uuid=uuid, which="last", load=True)
    energy = output["output"]["energy"]  # 提取能量
```

**关键发现**：
- QHA 流程已经包含 R3 优化步骤
- 可以从 `output["output"]` 中提取能量和其他属性
- 需要探索 atomate2 的 TaskDocument 结构以提取磁矩

#### 2. `htvasp/workflows/relax.py` - 结构优化流程
**核心功能**：
- 使用 `DoubleRelaxMaker` 进行 R7 + R3 双重优化
- **适用场景**：超大包结构的优化（不适合端元计算）
- **优化策略**：先 R7（体积优化）再 R3（全优化）

**关键理解**：
- relax.py 不是为端元设计，端元的 R3 优化已在 qha.py 中实现
- relax.py 适合需要体积优化的复杂结构

## 扩展规划（优先级排序）

### Priority 1: 磁性参数计算模块

#### 1.1 `htvasp/workflows/magnetic.py` - 磁性参数提取

**设计思路**：
1. **输入**：已优化的端基结构（来自 QHA 流程的 R3 结果）
2. **数据来源**：使用 atomate2 框架从 JobStore 提取
3. **输出**：BMAGN 和 TC 参数

**技术实现**：

**方案 A：从现有 QHA 流程提取（推荐）**
```python
class MagneticWorker(Worker):
    def extract_magnetic_from_qha(self, store_path: Path):
        """
        从 QHA 计算结果中提取磁性参数
        """
        self.store = JobStore(JSONStore(store_path))
        self.store.connect()
        
        # 查询 R3 优化结果
        relax_job = self.store.query_one(
            criteria={"name": {"$regex": "r3_relax"}},
            properties=["uuid", "index", "name"],
            sort={"index": -1},
        )
        
        output = self.store.get_output(uuid=relax_job["uuid"], which="last", load=True)
        
        # 提取磁矩（需要确认 atomate2 TaskDocument 的结构）
        # 可能的字段：output["output"]["total_magnetization"]
        # 或：output["output"]["magnetization"]
        
        return {
            "total_magnetization": output["output"].get("total_magnetization"),
            "energy": output["output"]["energy"],
        }
```

**方案 B：独立静态计算（备用）**
```python
class MagneticWorker(Worker):
    def run_static_calculation(self, structure: Structure):
        """
        对已优化结构进行静态计算以提取磁性参数
        """
        # 使用 atomate2 的 StaticMaker
        from atomate2.vasp.jobs.core import StaticMaker
        from atomate2.vasp.sets.core import StaticSetGenerator
        
        static_maker = StaticMaker(
            input_set_generator=StaticSetGenerator(
                user_incar_settings={
                    "ISPIN": 2,
                    "ISMEAR": 0,
                    "SIGMA": 0.05,
                    "EDIFF": 1e-6,
                    "LORBIT": 11,  # 输出磁矩信息
                }
            )
        )
        
        # 运行计算并提取结果
        # ...
```

**TC 计算**：
- 使用 OstravaJ 软件包（预印本论文）
- 需要安装和配置 OstravaJ
- 输入：磁性结构、磁矩配置
- 输出：居里温度估算值

**数据结构设计**：
```python
class MagneticData(TypedDict):
    name: str
    structure: dict[str, Any]
    total_magnetization: float  # μB
    bmag: float  # BMAGN 参数
    tc: float | None  # 居里温度 (K)，可能为 None
    energy: float  # eV
```

### Priority 2: SQS 结构生成和计算模块

#### 2.1 `htvasp/workflows/sqs.py` - SQS 结构生成和 VASP 计算

**技术选型**：
- **SQS 生成工具**：`sqsgenerator` Python 库（而非 ATAT）
- **理由**：
  - atomate2 0.0.23 版本中**没有内置的 SQS 功能**
  - `sqsgenerator` 更快、更适配、Python 原生支持
  - 可以与 atomate2 的 Job 和 Flow 框架无缝集成
- **调研结果**：
  - atomate2 当前版本（0.0.23）不包含 `atomate2.common.jobs.transform.SQS`
  - 该功能可能在开发分支或未来版本中
  - 需要自行实现 SQS 生成和计算流程

**设计思路**：
1. **SQS 生成**：使用 sqsgenerator 生成无序固溶体结构
2. **VASP 计算**：R3 优化 → R2 静态计算
3. **数据提取**：从 JobStore 提取能量和磁矩

**核心流程**：
```python
class SqsWorker(Worker):
    def generate_sqs_structure(
        self,
        base_structure: Structure,
        composition: dict[str, float],  # {"A": 0.5, "B": 0.5}
        supercell_size: int = 16,
    ) -> Structure:
        """
        使用 sqsgenerator 生成 SQS 结构
        """
        from sqsgenerator import sqs_optimize
        
        # 生成 SQS 结构
        sqs_structure = sqs_optimize(
            structure=base_structure,
            composition=composition,
            supercell_size=supercell_size,
        )
        
        return sqs_structure
    
    def run_sqs_calculation(
        self,
        sqs_structure: Structure,
        store_path: Path,
    ):
        """
        对 SQS 结构进行 R3 优化和 R2 静态计算
        """
        # R3 优化
        relax_r3_maker = DoubleRelaxMaker.from_relax_maker(
            RelaxMaker(
                name="sqs_r3_relax",
                input_set_generator=RelaxSetGenerator(
                    user_incar_settings={
                        "ISIF": 3,  # 全优化
                        "ISPIN": 2,  # 磁性
                        # ... 其他参数
                    }
                )
            )
        )
        
        # R2 静态计算
        static_maker = StaticMaker(
            name="sqs_r2_static",
            input_set_generator=StaticSetGenerator(
                user_incar_settings={
                    "ISIF": 2,  # 固定体积
                    "ISPIN": 2,
                    # ... 其他参数
                }
            )
        )
        
        # 构建工作流并运行
        # ...
        
        # 从 JobStore 提取结果
        self.store.connect()
        
        # 查询 R2 静态计算结果
        static_job = self.store.query_one(
            criteria={"name": {"$regex": "sqs_r2_static"}},
            properties=["uuid"],
            sort={"index": -1},
        )
        
        output = self.store.get_output(uuid=static_job["uuid"], which="last", load=True)
        
        return {
            "energy": output["output"]["energy"],
            "total_magnetization": output["output"].get("total_magnetization"),
            "n_atoms": len(sqs_structure),
        }
```

**数据结构设计**：
```python
class SqsData(TypedDict):
    name: str
    composition: dict[str, float]  # {"A": 0.5, "B": 0.5}
    structure: dict[str, Any]
    energy: float  # eV
    n_atoms: int
    total_magnetization: float | None  # μB
    reference_energies: dict[str, float]  # {"A": E_FCC_A, "B": E_FCC_B}
    mixing_enthalpy: float  # eV/atom
```

### Priority 3: 参考态能量管理

#### 3.1 参考态能量提取和存储

**设计思路**：
- 从现有 QHA 计算结果中提取纯元素参考态能量
- 确保参考态与 SQS 结构同结构（如都是 FCC）

**实现方案**：
```python
class ReferenceEnergyManager:
    def extract_reference_energies(
        self,
        element: str,
        structure_type: str,  # "FCC", "BCC", etc.
        qha_store_path: Path,
    ):
        """
        从 QHA 计算结果中提取参考态能量
        """
        store = JobStore(JSONStore(qha_store_path))
        store.connect()
        
        # 查询 EOS 静态计算结果
        eos_job = store.query_one(
            criteria={"name": {"$regex": "phonon static eos deformation 0"}},
            properties=["uuid"],
        )
        
        output = store.get_output(uuid=eos_job["uuid"], which="last", load=True)
        
        return {
            "element": element,
            "structure_type": structure_type,
            "energy": output["output"]["energy"],
        }
```

### Priority 4: 整合流程

#### 4.1 `htvasp/workflows/calphad_preprocess.py` - CALPHAD 前处理主流程

**流程设计**：

**Phase 1: 端基热力学性质（基于现有 QHA）**
```python
def run_endmember_calculations(
    self,
    elements: list[str],
    structure_type: str = "FCC",
):
    """
    计算端基热力学性质
    """
    # 1. 运行 QHA 计算（使用现有 QhaWorker）
    qha_worker = QhaWorker(...)
    
    for element in elements:
        structure = self.get_structure(element, structure_type)
        qha_data = qha_worker.run_flow(
            struct_name=f"{element}_{structure_type}",
            structure=structure,
            flow_dir=self.flow_dir / f"{element}_{structure_type}",
        )
        
        # 2. 提取磁性参数（使用 MagneticWorker）
        magnetic_worker = MagneticWorker(...)
        magnetic_data = magnetic_worker.extract_magnetic_from_qha(
            store_path=qha_worker.store.store.path
        )
        
        # 3. 存储结果
        self.results["endmembers"][element] = {
            "qha": qha_data,
            "magnetic": magnetic_data,
        }
```

**Phase 2: SQS 混合焓计算**
```python
def run_sqs_calculations(
    self,
    element_a: str,
    element_b: str,
    compositions: list[float] = [0.25, 0.5, 0.75],
):
    """
    计算 SQS 混合焓
    """
    sqs_worker = SqsWorker(...)
    ref_manager = ReferenceEnergyManager()
    
    # 1. 获取参考态能量
    ref_energy_a = ref_manager.extract_reference_energies(
        element=element_a,
        structure_type="FCC",
        qha_store_path=self.results["endmembers"][element_a]["store_path"],
    )
    ref_energy_b = ref_manager.extract_reference_energies(
        element=element_b,
        structure_type="FCC",
        qha_store_path=self.results["endmembers"][element_b]["store_path"],
    )
    
    # 2. 对每个成分点生成 SQS 并计算
    for x_b in compositions:
        # 生成 SQS 结构
        sqs_structure = sqs_worker.generate_sqs_structure(
            base_structure=self.get_structure(element_a, "FCC"),
            composition={element_a: 1-x_b, element_b: x_b},
        )
        
        # 运行 VASP 计算
        sqs_data = sqs_worker.run_sqs_calculation(
            sqs_structure=sqs_structure,
            store_path=self.flow_dir / f"sqs_{element_a}_{element_b}_{x_b}",
        )
        
        # 3. 计算混合焓
        mixing_enthalpy = self.calculate_mixing_enthalpy(
            sqs_energy=sqs_data["energy"],
            ref_energy_a=ref_energy_a["energy"],
            ref_energy_b=ref_energy_b["energy"],
            x_b=x_b,
            n_atoms=sqs_data["n_atoms"],
        )
        
        sqs_data["mixing_enthalpy"] = mixing_enthalpy
        self.results["sqs"][f"{element_a}_{element_b}_{x_b}"] = sqs_data
```

**混合焓计算公式**：
```python
def calculate_mixing_enthalpy(
    self,
    sqs_energy: float,
    ref_energy_a: float,
    ref_energy_b: float,
    x_b: float,
    n_atoms: int,
) -> float:
    """
    计算混合焓
    
    ΔH_mix(x) = [E_SQS - ((1-x)*E_A + x*E_B)] / N_atoms
    
    Args:
        sqs_energy: SQS 结构的总能 (eV)
        ref_energy_a: 纯元素 A 的参考态能量 (eV/atom)
        ref_energy_b: 纯元素 B 的参考态能量 (eV/atom)
        x_b: 元素 B 的成分
        n_atoms: SQS 结构的原子数
    
    Returns:
        混合焓 (eV/atom)
    """
    reference_energy = (1 - x_b) * ref_energy_a + x_b * ref_energy_b
    mixing_enthalpy = (sqs_energy - reference_energy) / n_atoms
    return mixing_enthalpy
```

## 技术实现细节

### 1. atomate2 SQS 功能调研

**调研结果**：
- **atomate2 版本**：0.0.23（最新稳定版）
- **SQS 功能状态**：**不存在**
  - 检查了 `atomate2.common.jobs` 模块
  - 没有找到 `SQS` 或 `transform.SQS` 类
  - 没有找到相关的 SQS 转换功能

**解决方案**：
- 使用 `sqsgenerator` Python 库
- 自定义实现与 atomate2 兼容的 SQS Job
- 可以与 atomate2 的其他 Job 和 Flow 无缝集成

**未来展望**：
- atomate2 可能在未来版本中添加 SQS 功能
- 可以关注 atomate2 的 GitHub 仓库更新
- 当前方案可以平滑迁移到未来的官方实现

### 2. atomate2 数据提取机制

**TaskDocument 结构探索**：
```python
# atomate2 的 TaskDocument 包含以下关键字段（需要验证）：
{
    "output": {
        "energy": float,  # 总能量
        "structure": dict,  # 优化后的结构
        "total_magnetization": float,  # 总磁矩
        "forces": list,  # 原子受力
        "stress": list,  # 应力张量
        # ... 其他字段
    },
    "input": {
        "structure": dict,  # 输入结构
        "parameters": dict,  # VASP 参数
        # ...
    }
}
```

**验证方法**：
```python
# 在测试中验证 TaskDocument 结构
output = store.get_output(uuid=job_uuid, which="last", load=True)
print(json.dumps(output, indent=2))  # 查看完整结构
```

### 2. sqsgenerator 使用

**安装**：
```bash
pip install sqsgenerator
```

**基本用法**：
```python
from sqsgenerator import sqs_optimize
from pymatgen.core import Structure

# 加载基础结构
base_structure = Structure.from_file("POSCAR_FCC")

# 生成 SQS
sqs_structure = sqs_optimize(
    structure=base_structure,
    composition={"Ni": 0.5, "Al": 0.5},
    supercell_size=16,
    iterations=1000,
)

# 保存结果
sqs_structure.to(filename="POSCAR_SQS")
```

**与 atomate2 集成**：
```python
from atomate2.common.jobs import Job
from jobflow import Flow
from sqsgenerator import sqs_optimize

class SQSGenerationJob(Job):
    """自定义 SQS 生成 Job，与 atomate2 兼容"""
    
    def __init__(
        self,
        base_structure: Structure,
        composition: dict[str, float],
        supercell_size: int = 16,
        name: str = "SQS Generation",
    ):
        super().__init__(name=name)
        self.base_structure = base_structure
        self.composition = composition
        self.supercell_size = supercell_size
    
    def run(self):
        """执行 SQS 生成"""
        sqs_structure = sqs_optimize(
            structure=self.base_structure,
            composition=self.composition,
            supercell_size=self.supercell_size,
        )
        return {"sqs_structure": sqs_structure}

# 使用示例
sqs_job = SQSGenerationJob(
    base_structure=fcc_structure,
    composition={"Ni": 0.5, "Al": 0.5},
)

# 可以与 atomate2 的其他 Job 组合成 Flow
flow = Flow([sqs_job, relax_job, static_job])
```

### 3. OstravaJ 集成

**安装和配置**：
```bash
# OstravaJ 是新软件，需要从源码安装
git clone https://github.com/.../OstravaJ.git
cd OstravaJ
pip install -e .
```

**使用方法**：
```python
from ostravaj import calculate_tc

# 计算 TC
tc = calculate_tc(
    structure=magnetic_structure,
    magnetic_moments=magnetic_moments,
    method="heisenberg",  # 或其他方法
)

print(f"Estimated Curie temperature: {tc} K")
```

## 依赖管理

### 更新 `pyproject.toml`

```toml
dependencies = [
    "pymatgen>=2022.0,<2026.0",
    "atomate2>=0.0.20,<0.0.23",
    "mp-api>=0.30",
    "seekpath>=2.0",
    "phonopy>=2.13",
    "sqsgenerator>=0.2",  # 新增：SQS 生成
]

[project.optional-dependencies]
dev = ["pytest>=6.0", "black"]
magnetic = ["ostravaj>=0.1"]  # 可选：TC 计算
```

### 创建 `requirements.txt`

```txt
# Core dependencies with version constraints
pymatgen>=2022.0,<2026.0
atomate2>=0.0.20,<0.0.23
mp-api>=0.30
seekpath>=2.0
phonopy>=2.13

# SQS generation
sqsgenerator>=0.2

# Optional: Magnetic TC calculation
# ostravaj>=0.1  # Uncomment if available
```

## 测试计划

### 1. 磁性参数提取测试

```python
# tests/test_magnetic.py
def test_extract_magnetic_from_qha():
    """测试从 QHA 结果提取磁性参数"""
    # 1. 运行一个简单的 QHA 计算
    # 2. 使用 MagneticWorker 提取磁性参数
    # 3. 验证提取的磁矩是否正确
    pass

def test_magnetic_static_calculation():
    """测试独立的磁性静态计算"""
    # 1. 创建一个磁性结构
    # 2. 运行静态计算
    # 3. 验证磁矩提取
    pass
```

### 2. SQS 生成和计算测试

```python
# tests/test_sqs.py
def test_sqs_generation():
    """测试 SQS 结构生成"""
    # 1. 使用 sqsgenerator 生成 SQS
    # 2. 验证结构的成分和对称性
    pass

def test_sqs_calculation():
    """测试 SQS VASP 计算"""
    # 1. 生成 SQS 结构
    # 2. 运行 R3 优化和 R2 静态计算
    # 3. 验证能量提取
    pass

def test_mixing_enthalpy():
    """测试混合焓计算"""
    # 1. 准备参考态能量
    # 2. 计算 SQS 能量
    # 3. 计算并验证混合焓
    pass
```

### 3. 整合流程测试

```python
# tests/test_calphad_preprocess.py
def test_endmember_workflow():
    """测试端基计算流程"""
    pass

def test_sqs_workflow():
    """测试 SQS 计算流程"""
    pass

def test_complete_workflow():
    """测试完整的 CALPHAD 前处理流程"""
    pass
```

## 实施计划

### 阶段 1：磁性参数模块（1-2 周）
- [ ] 探索 atomate2 TaskDocument 结构
- [ ] 实现 `MagneticWorker.extract_magnetic_from_qha()`
- [ ] 实现 OstravaJ 集成（如果可用）
- [ ] 编写测试用例

### 阶段 2：SQS 模块（2-3 周）
- [ ] 安装和测试 sqsgenerator
- [ ] 实现 `SqsWorker.generate_sqs_structure()`
- [ ] 实现 R3→R2 计算流程
- [ ] 实现混合焓计算
- [ ] 编写测试用例

### 阶段 3：参考态管理（1 周）
- [ ] 实现参考态能量提取
- [ ] 实现参考态数据库管理
- [ ] 编写测试用例

### 阶段 4：整合流程（2 周）
- [ ] 实现 `CalphadPreprocessWorker`
- [ ] 整合端基和 SQS 计算
- [ ] 实现完整的数据输出
- [ ] 编写端到端测试

### 阶段 5：文档和优化（1 周）
- [ ] 编写详细的 API 文档
- [ ] 编写使用示例
- [ ] 性能优化
- [ ] 错误处理和边缘情况

## 预期成果

1. **完整的 CALPHAD 前处理流程**
   - 端基热力学性质（QHA + 磁性）
   - SQS 混合焓计算
   - 自动化的数据提取和存储

2. **高质量的第一性原理数据**
   - 磁性参数：BMAGN、TC
   - 混合焓：ΔH_mix(x)
   - 参考态能量：E_FCC(A)、E_FCC(B)

3. **与 atomate2 框架的无缝集成**
   - 使用 JobStore 进行数据管理
   - 使用 atomate2 的 Maker 和 Flow
   - 遵循 atomate2 的最佳实践

4. **可扩展的架构**
   - 模块化设计
   - 支持不同的晶体结构
   - 支持批量处理

## 风险和挑战

1. **atomate2 TaskDocument 结构**
   - 需要验证磁矩字段的确切名称
   - 可能需要查看 atomate2 源码

2. **OstravaJ 可用性**
   - 软件较新，可能不稳定
   - 需要备用方案（如简化模型）

3. **SQS 生成质量**
   - sqsgenerator 的性能需要验证
   - 可能需要调整参数

4. **计算资源**
   - SQS 计算需要较多资源
   - 需要合理规划计算任务

## 总结

本扩展规划基于对现有 HT-VASP 代码的深入理解，采用 atomate2 框架的最佳实践，实现了 CALPHAD 端基数据库构建所需的前处理和 VASP 计算功能。通过模块化设计和优先级排序，确保了项目的可实施性和可扩展性。
