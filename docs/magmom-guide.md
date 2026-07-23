# 磁矩（Magnetic Moment）获取方式

## 方式一：从 JSON 获取（最常用）

输出结果存储在 `*-static.json` 或 `*-qha.json`，解析即可。

**关键路径：**

| 字段 | JSON 路径 | 示例值 |
|------|----------|--------|
| 初始磁矩设置 | `data["input"]["incar"]["MAGMOM"]` | `[3.0, 3.0]` |
| 每个原子收敛磁矩 | `data["output"]["structure"]["sites"][i]["properties"]["magmom"]` | `0.137` |
| 总磁矩 | `data["calcs_reversed"][0]["output"]["outcar"]["total_magnetization"]` | `1.5246` |
| 磁矩密度 | `data["calcs_reversed"][0]["output"]["mag_density"]` | `0.052` |
| 轨道分解磁矩 | `data["calcs_reversed"][0]["output"]["outcar"]["magnetization"]` | `[{"s","p","d","tot"}, ...]` |
| 是否自旋极化 | `data["input"]["incar"]["ISPIN"]` | `2` |
| SCF 步数 | `data["calcs_reversed"][0]["output"]["num_electronic_steps"]` | `[17]` |

**Python 代码片段：**

```python
import json

with open("SER-Al-static.json") as f:
    data = json.load(f)

cr = data["calcs_reversed"][0]["output"]

# 每个原子的磁矩
site_magmoms = [s["properties"]["magmom"] for s in cr["structure"]["sites"]]

# 总磁矩
total_mag = cr["outcar"]["total_magnetization"]

# 轨道分解（s/p/d）
orbital_mag = cr["outcar"]["magnetization"]  # list[dict]
```

**注意事项：**
- `LORBIT = 11` 时才有轨道分解磁矩
- `ISPIN = 1`（非自旋极化）时磁矩始终为 0
- JSON 中 `is_metal` 字段可能为 `None`（非必需）
- 此方式不包含每个 SCF 步的磁矩演化轨迹（需读原始 OUTCAR）
- **结构注意**：`data["output"]` 与 `cr["structure"]` 中都有 `sites`，两者 magmom 一致

---

## 方式二：从原始 VASP 文件 + pymatgen 获取

原始计算结果存储在 `store_dir`（例如 `staticflow/3-static/`），文件均经 gzip 压缩。

### 所需文件

| 文件 | 用途 |
|------|------|
| `OUTCAR.gz` | 完整的 VASP 输出，含每个 SCF 步的磁矩 |
| `vasprun.xml.gz` | XML 格式的结构化输出 |
| `CONTCAR.gz` | 最终结构的 POSCAR |

### 方法 A：从 OUTCAR 读取（推荐）

```python
from pymatgen.io.vasp import Outcar

outcar = Outcar.from_file("OUTCAR.gz")

# 最终磁矩
total_mag = outcar.total_magnetization  # float
site_mag = outcar.magnetization         # list[dict], 轨道分解
spin_mag = outcar.spin_mag              # list[float], 每个原子的自旋磁矩
charges = outcar.charge                 # list[dict], 每个原子的电荷

# 每个 SCF 步的磁矩演化（若有）
if hasattr(outcar, "all_magnetizations"):
    for step, mag in enumerate(outcar.all_magnetizations):
        print(f"SCF step {step}: mag = {mag}")
```

### 方法 B：从 vasprun.xml 读取

```python
from pymatgen.io.vasp import Vasprun

vasprun = Vasprun("vasprun.xml.gz")

# 最终结构（含磁矩）
final_struct = vasprun.final_structure
site_magmoms = [site.properties.get("magmom") for site in final_struct]

# 每个离子步的结构（含磁矩轨迹）
for step, struct in enumerate(vasprun.structures):
    mags = [s.properties.get("magmom", 0) for s in struct]
    print(f"Ionic step {step}: {mags}")

# 直接读取 OUTCAR 中的磁矩
outcar = vasprun.outcar
```

### 方法 C：从 CONTCAR + OUTCAR 组合

```python
from pymatgen.core import Structure
from pymatgen.io.vasp import Outcar

struct = Structure.from_file("CONTCAR.gz")
outcar = Outcar.from_file("OUTCAR.gz")

# 将磁矩注入结构对象
for i, site in enumerate(struct):
    site.properties["magmom"] = outcar.magnetization[i]["tot"]
```

### gzip 文件处理

存储目录下的文件均为 `.gz` 格式，pymatgen 的 `from_file()` 直接支持 gzip 读取，无需手动解压。

### 使用示例（批量处理）

```python
from pathlib import Path
from pymatgen.io.vasp import Outcar

for outcar_path in Path("data/endmembers").rglob("*/staticflow/*/OUTCAR.gz"):
    outcar = Outcar.from_file(outcar_path)
    print(f"{outcar_path.parent.name}: total_mag={outcar.total_magnetization}")
```

---

## 两种方式对比

| 维度 | JSON 方式 | OUTCAR 方式 |
|------|----------|------------|
| 速度 | 极快，无需解压 | 较慢，需解压 gzip |
| 易用性 | 直接读取，结构简单 | 需安装 pymatgen |
| SCF 步演化 | ❌ 无 | ✅ 有（`all_magnetizations`） |
| 轨道分解 | ✅ 有 | ✅ 有 |
| 结构-磁矩关联 | ✅ 直接对应 | ✅ 可注入 |
| 依赖 | 无 | `pymatgen` |
