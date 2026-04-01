# sqsgenerator 参数配置指南（Markdown 格式）

> 适用于 Coding AI 参考：结构化参数说明 + 配置示例 + 注意事项

---

## 📋 概述

`sqsgenerator` 使用类 `dict` 的配置方式（默认 YAML 格式）来生成特殊准随机结构（SQS）。以下参数均为 YAML/JSON 中的 key，可直接用于 Python dict 配置。

---

## 🔧 核心参数速查表

| 参数名 | 必填 | 默认值 | 类型 | 作用简述 |
|--------|------|--------|------|----------|
| `composition` | ✅ | - | `dict` | 定义输出结构的元素组成与原子数 |
| `structure` | ✅ | - | `str`/`dict` | 输入晶体结构（文件路径或内联定义） |
| `iteration_mode` | ❌ | `random` | `str` | 结构生成模式：`random` / `systematic` |
| `sublattice_mode` | ❌ | `interact` | `str` | 亚晶格处理方式：`interact` / `split` |
| `iterations` | ❌ | `1000` | `int` | 随机模式下的迭代次数 |
| `shell_radii` | ❌ | 自动计算 | `list[float]` | 配位壳层半径边界（Å） |
| `shell_weights` | ❌ | `1/i` | `dict[int,float]` | 各壳层权重（距离越远权重越低） |
| `pair_weights` | ❌ | 空心矩阵 | `ndarray` | 原子对权重矩阵 $p_{cn}$ |
| `target_objective` | ❌ | 零矩阵 | `ndarray` | 目标 SRO 参数 $\alpha_n$ |
| `prefactors` | ❌ | 零矩阵 | `ndarray` | 键合预因子 $f$ |
| `bin_width` | ❌ | `0.05` | `float` | 直方图分箱宽度（用于自动识别壳层） |
| `peak_isolation` | ❌ | `0.25` | `float` | 峰隔离阈值 [0.0, 1.0] |
| `thread_config` | ❌ | 所有物理核 | `list[int]` | 每 rank 使用的线程数 |

---

## 📦 参数详解

### 1. `composition`（必填）
定义输出结构中各元素的原子数量。

```yaml
# 示例1：三元合金 Ti18Al18Mo18
composition:
  Ti: 18
  Al: 18
  Mo: 18

# 示例2：含空位的 fcc-Al（56 Al + 8 空位）
composition:
  Al: 56
  Va: 8   # 用 "Va" 或 "0" 表示空位

# 示例3：亚晶格约束（split 模式）
# (Ti0.25Al0.25)(B0.25N0.25)，Ti/Al 占 Ti 位，B/N 占 N 位
composition:
  - sites: "Ti"
    Ti: 16
    Al: 16
  - sites: "N"
    N: 16
    B: 16
```

> ⚠️ 原子总数必须严格匹配晶格位置数；若使用 `which` 指定位点，数量需与选中位点数一致。

---

### 2. `structure`（必填）
输入晶体结构，支持文件路径或内联定义。

```yaml
# 方式1：内联定义 CsCl (B2) 结构
structure:
  lattice:
    - [4.123, 0.0, 0.0]
    - [0.0, 4.123, 0.0]
    - [0.0, 0.0, 4.123]
  coords:  # 分数坐标（非笛卡尔）
    - [0.0, 0.0, 0.0]
    - [0.5, 0.5, 0.5]
  species:
    - Cs
    - Cl

# 方式2：读取文件（自动检测 ase/pymatgen）
structure:
  file: "cs-cl.vasp"  # 支持 POSCAR/CONTCAR/CIF 等

# 方式3：指定读取器 + 参数
structure:
  file: "md.traj"
  reader: "ase"
  args:
    index: -1  # 读取轨迹最后一帧
```

#### `structure.supercell`（可选）
创建超胞：
```yaml
structure:
  supercell: [3, 3, 3]  # 3×3×3 超胞
  file: "cs-cl.cif"
```

---

### 3. 迭代与亚晶格控制

#### `iteration_mode`
```yaml
iteration_mode: "random"   # 默认：随机打乱
# 或
iteration_mode: "systematic"  # 字典序遍历全空间（忽略 iterations 参数）
```
> ⚠️ `systematic` 仅支持 `sublattice_mode: interact`；配置前请用 `sqsgenerator compute total-permutations` 评估计算量。

#### `sublattice_mode`
```yaml
sublattice_mode: "interact"  # 默认：整体优化，原子种类固定在亚晶格上
# 或
sublattice_mode: "split"     # 各亚晶格独立优化（仅支持 random 模式）
```

#### `iterations`
```yaml
iterations: 1000  # random 模式下的采样次数（默认 10³）
```

---

### 4. 配位壳层参数

#### `shell_radii`
```yaml
# 手动指定壳层半径边界（单位：Å）
shell_radii: [2.5, 4.0, 6.0]  # 壳层1: (0,2.5], 壳层2: (2.5,4.0], ...

# split 模式下可为每个亚晶格单独指定
shell_radii:
  - [2.0, 3.5]   # sublattice 0
  - [2.2, 3.8]   # sublattice 1
```

#### 自动壳层检测参数
```yaml
bin_width: 0.05        # 直方图分箱宽度（Å）
peak_isolation: 0.25   # 峰隔离阈值：邻箱高度 < 25% 则视为独立壳层
```

#### `shell_weights`
```yaml
# 壳层索引 → 权重（用于目标函数中的 wⁱ）
shell_weights:
  1: 1.0   # 第一壳层（默认）
  2: 0.5   # 第二壳层权重减半
  # 未指定的壳层将被忽略
```

---

### 5. 目标函数参数（多维数组）

> 📌 原子种类按**原子序数升序**内部重排！  
> 例：CsCl → `[Cl(Z=17), Cs(Z=55)]`；TiAlMo → `[Al(13), Ti(22), Mo(42)]`

#### `pair_weights`（$p_{cn}$）
原子对区分权重，默认空心矩阵：$p = \mathbb{1}_N - I_N$

```yaml
# 2D 输入：(N_species, N_species)，自动沿壳层维度堆叠
pair_weights:
  - [0, 1, 1]
  - [1, 0, 1]
  - [1, 1, 0]

# 3D 输入：(N_shells, N_species, N_species)，直接指定每壳层权重
```

#### `target_objective`（$\alpha_n$，目标 SRO）
```yaml
# 标量：全填充
target_objective: 0  # 默认：完全随机

# 2D 矩阵：(N_species, N_species)
target_objective:
  - [1, 1]
  - [1, 1]           # CsCl 聚类结构目标

# 3D 矩阵：(N_shells, N_species, N_species)
target_objective:
  - [[1, -1, 0],
     [-1, 1, 0],
     [0, 0, 1]]      # 三元合金自定义目标
```
> ⚠️ $\alpha_{cn} = \alpha_{nc}$，输入矩阵必须对称，否则抛出 `BadSettings`。

#### `prefactors`（$f$，键合预因子）
输入格式同 `target_objective`，默认零矩阵。可理解为期望键合数的倒数。

---

### 6. 并行配置

#### `thread_config`
```yaml
# 单进程：使用所有物理核（默认）
thread_config: 8

# MPI 多 rank：为每个 rank 指定线程数
thread_config: [4, 4, 4]  # 3 ranks，每 rank 4 线程
```

---

## 🧩 数组输入解释规则（重要！）

以下参数支持 1D/2D/3D 输入，AI 生成代码时需严格校验维度：

| 参数 | 1D 输入 | 2D 输入 | 3D 输入 |
|------|---------|---------|---------|
| `shell_radii` | 全局壳层边界 | `split` 模式下：`[sublattice_0, sublattice_1, ...]` | - |
| `pair_weights` | - | `(N_sp, N_sp)` → 沿壳层堆叠 | `(N_shells, N_sp, N_sp)` 直接使用 |
| `target_objective` | 标量 → 全填充 | `(N_sp, N_sp)` → 沿壳层堆叠 | `(N_shells, N_sp, N_sp)` 直接使用 |
| `prefactors` | 同 `target_objective` | 同 `target_objective` | 同 `target_objective` |

> ✅ 所有矩阵输入必须对称（$\alpha_{cn} = \alpha_{nc}$）  
> ✅ 原子种类顺序 = 按原子序数升序排列

---

## 🚀 常用配置模板（AI 可直接复用）

### 模板1：基础随机 SQS（三元合金）
```yaml
structure:
  file: "fcc-supercell.vasp"
composition:
  Ti: 18
  Al: 18
  Mo: 18
iteration_mode: random
iterations: 5000
shell_weights:
  1: 1.0
  2: 0.5
```

### 模板2：亚晶格约束（岩盐结构掺杂）
```yaml
structure:
  file: "TiN-64.vasp"
sublattice_mode: split
composition:
  - sites: "Ti"
    Ti: 16
    Al: 8
    V: 8
  - sites: "N"
    N: 32
iteration_mode: random
iterations: 2000
```

### 模板3：目标导向有序化（CsCl 聚类）
```yaml
structure:
  file: "CsCl.vasp"
composition:
  Cs: 1
  Cl: 1
target_objective:
  - [1, 1]
  - [1, 1]
shell_weights:
  1: 1.0
```

---

## ⚠️ 注意事项（AI 生成代码时必须校验）

1. **原子数守恒**：`composition` 中原子总数 == 晶格位置数（或 `which` 选定位点数）
2. **矩阵对称性**：`target_objective`/`pair_weights` 必须对称，否则报错
3. **模式兼容性**：`systematic` 仅支持 `interact` 亚晶格模式
4. **壳层索引**：`shell_weights` 中未指定的壳层将被忽略；指定不存在的索引会报错
5. **原子序排序**：多维数组的物种维度按原子序数升序排列，非输入顺序
6. **单位统一**：所有距离参数单位为 **Å**，坐标必须为**分数坐标**

---

## 🔍 调试与验证命令

```bash
# 查看解析后的参数（如自动计算的 shell_weights）
sqsgenerator params show input.yaml -p shell_weights

# 预估系统排列总数（systematic 模式前必查）
sqsgenerator compute total-permutations -i input.yaml

# 预估计算时间
sqsgenerator compute estimated-time -i input.yaml --iterations 1000
```

---

> 📌 **Coding AI 提示**：生成配置时优先使用 `composition` + `structure.file` 简化输入；多维数组参数建议封装为 `numpy` 数组校验函数，确保维度与对称性符合要求。

*文档基于 sqsgenerator v9+，参考：https://sqsgenerator.readthedocs.io/en/latest/parameters.html*