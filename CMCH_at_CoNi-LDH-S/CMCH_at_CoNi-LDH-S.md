

```markdown

## 一、实验制备关键信息总结

| 项目 | 内容 |
|------|------|
| **目标材料** |  CoMn(CO₃)(OH)₂ @ S掺杂的 CoNi-LDH 复合电极（简称 CMCH@CoNi-LDH-S） |
| **基底** | 泡沫镍（2×3 cm²），预先生长有 CMCH@ZIF-67 |
| **硫源** | 硫脲（CH₄N₂S） |
| **硫加入量** | 1、3、7、10 wt%（相对于溶液或前驱体质量）<br>实验人员实际掺杂浓度区间：**0.01 ~ 0.03**（推测为原子比 S/(S+O) 或摩尔分数） |
| **金属前驱体** | Ni(NO₃)₂·6H₂O（1 mmol，溶于30 mL乙醇） |
| **反应环境** | 水热法，120 °C，2 小时 |
| **溶剂** | 乙醇（30 mL）—— **无额外加水** |
| **后处理** | 去离子水 + 乙醇冲洗，60 °C 干燥 12 h |

---

## 二、关于取代物种的判断：SH⁻ 还是 S²⁻？

### 1. 反应化学环境分析
- **硫脲在热/水热条件下分解**：  
  CH₄N₂S + 2H₂O → CO₂ + 2NH₃ + H₂S（或 HS⁻/S²⁻，取决于pH）
- **溶剂为乙醇**：含水量低，但仍可能由硫脲分解产生少量H₂S。
- **无额外酸碱调节**：体系中存在微量NH₃（来自硫脲分解），呈弱碱性。
- **LDH层板**：CoNi-LDH层板含有大量 OH⁻ 基团，可与硫物种发生**阴离子交换**或**配体取代**。

### 2. 热力学与结构兼容性
| 取代物种 | 电荷 | 半径（Å） | 能否稳定存在于层板中？ | 典型掺杂形式 |
|----------|------|-----------|------------------------|----------------|
| **OH⁻**  | -1   | ~1.4      | 是（本征）              | — |
| **SH⁻**  | -1   | ~1.8      | 可能，但S-H键易断裂     | 少见，不稳定 |
| **S²⁻**  | -2   | ~1.84     | 会改变层板电荷平衡，但文献常见 | **M-S-M** 桥联形式，稳定 |

### 3. 文献与常识
- **大多数硫掺杂LDH（如CoNi-LDH、NiFe-LDH）使用硫脲或硫代乙酰胺合成时**，最终嵌入层板的是 **S²⁻**，取代 OH⁻ 位置，形成 **M-S-M** 结构，而不是保留 SH⁻ 末端。
- **SH⁻ 极易去质子化**：在弱碱性或加热条件下，SH⁻ → S²⁻ + H⁺，H⁺被NH₃中和。
- 若为 SH⁻ 取代，红外谱图应有明显的 S-H 伸缩振动（~2550 cm⁻¹），此类文献未见报道。

### 4. 结论用于建模
> **应使用 S²⁻ 替换模型中的 OH⁻**，而不是 SH⁻。  
> 每个 S²⁻ 取代两个相邻的 OH⁻ 位置（维持电荷平衡），或同时引入一个金属空位 / 额外阳离子。

---

## 三、对建模的具体建议（3×3×2 超胞）

- **实验掺杂浓度**：0.01 ~ 0.03（原子比 S / (S+O)）  
- **超胞中总OH⁻数目估算**：  
  以常见CoNi-LDH（M₂(OH)₆ 层板单元）为例，3×3×2 超胞约含 36 个 OH⁻。  
- **替换个数**：  
  - 若浓度 0.01 → 替换 0.36 ≈ 0 或 1 个 S²⁻（但0.01偏低，可能取0.02）  
  - 更合理：浓度 0.02 → 替换 ~0.72 ≈ **1 个 S²⁻**（对应约 2 个 OH⁻ 被 1 个 S²⁻ 替换）  
  - 浓度 0.03 → 替换 ~1.08 ≈ **1~2 个 S²⁻**（对应 2~4 个 OH⁻）

- **电荷补偿机制**：  
  每个 S²⁻ 取代两个 OH⁻ 时，层板电荷不变（-2 → -2）。若只取代一个 OH⁻ 位置，则需额外引入阴离子（如层间 CO₃²⁻）或金属空位。

建议优先尝试 **1个S²⁻ 替换两个相邻的OH⁻**（形成 M-S-M 桥），并与实验XPS（S 2p 峰位 ~161-163 eV，典型S²⁻）对照验证。
```

CoNi-LDH所使用的原型结构：Ni(OH)2, mp-27912.json

```cif
# generated using pymatgen
data_Ni(HO)2
_symmetry_space_group_name_H-M   P-3m1
_cell_length_a   3.12512265
_cell_length_b   3.12512265
_cell_length_c   4.47344830
_cell_angle_alpha   90.00000000
_cell_angle_beta   90.00000000
_cell_angle_gamma   120.00000000
_symmetry_Int_Tables_number   164
_chemical_formula_structural   Ni(HO)2
_chemical_formula_sum   'Ni1 H2 O2'
_cell_volume   37.83617158
_cell_formula_units_Z   1
loop_
 _symmetry_equiv_pos_site_id
 _symmetry_equiv_pos_as_xyz
  1  'x, y, z'
  2  '-x, -y, -z'
  3  '-y, x-y, z'
  4  'y, -x+y, -z'
  5  '-x+y, -x, z'
  6  'x-y, x, -z'
  7  'y, x, -z'
  8  '-y, -x, z'
  9  'x-y, -y, -z'
  10  '-x+y, y, z'
  11  '-x, -x+y, -z'
  12  'x, x-y, z'
loop_
 _atom_type_symbol
 _atom_type_oxidation_number
  Ni2+  2.0
  H+  1.0
  O2-  -2.0
loop_
 _atom_site_type_symbol
 _atom_site_label
 _atom_site_symmetry_multiplicity
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
  Ni2+  Ni0  1  0.00000000  0.00000000  0.00000000  1
  H+  H1  2  0.33333333  0.66666667  0.55907427  1
  O2-  O2  2  0.33333333  0.66666667  0.77512013  1
```


```json
{
  "crystal_structure": {
    "atomic_positions": [
      {
        "wyckoff": "1a",
        "element": "Ni",
        "a": "0",
        "b": "0",
        "c": "0"
      },
      {
        "wyckoff": "2d",
        "element": "H",
        "a": "1/3",
        "b": "2/3",
        "c": "0.559074"
      },
      {
        "wyckoff": "2d",
        "element": "O",
        "a": "1/3",
        "b": "2/3",
        "c": "0.77512"
      }
    ],
    "symmetry": {
      "crystal_system": "Trigonal",
      "lattice_system": "Hexagonal",
      "hall_number": "-P 3 2\"",
      "international_number": 164,
      "symbol": "P3̅m1",
      "point_group": "3̅m"
    },
    "structure_meta": {
      "number_of_atoms": 5,
      "density": "4.07 g·cm⁻³",
      "dimensionality": "",
      "possible_oxidation_states": "O²⁻, H⁺, Ni²⁺"
    }
  }
}
```

需要将该六方晶胞正交化，得到：

```
Ni1 H2 O2
1.0
        3.1251225471         0.0000000000         0.0000000000
        0.0000000000         5.4128708839         0.0000000000
        0.0000000000         0.0000000000         4.4734482765
   Ni    H    O
    2    4    4
Direct
     0.000000000         0.000000000         0.000000000
     0.500000000         0.500000000         0.000000000
     0.000000000         0.333333343         0.559074270
     0.000000000         0.666666657         0.440925730
     0.499999985         0.833333328         0.559074270
     0.500000015         0.166666672         0.440925730
     0.000000000         0.333333343         0.775120130
     0.000000000         0.666666657         0.224879870
     0.499999985         0.833333328         0.775120130
     0.500000015         0.166666672         0.224879870
```

```
Title               Ni1 H2 O2

Lattice type        P
Space group name    P 1
Space group number  1
Setting number      1

Lattice parameters

   a        b        c       alpha    beta     gamma
 3.12512  5.41287  4.47345  90.0000  90.0000  90.0000

Unit-cell volume = 75.672336 Å^3
```

接下来需要做CoMn(OH)2CO3的建模，它与rosasite（玫瑰砷酸锌矿） (Cu,Zn)₂(CO₃)(OH)₂，是完全同构体，由于这两种个晶体MP数据库都没有，找到了Co2(OH)2CO3，mp-1202876


```cif
# generated using pymatgen
data_Co2H2CO5
_symmetry_space_group_name_H-M   P2_1/c
_cell_length_a   3.19535090
_cell_length_b   12.30609000
_cell_length_c   9.67280607
_cell_angle_alpha   90.00000000
_cell_angle_beta   93.05194829
_cell_angle_gamma   90.00000000
_symmetry_Int_Tables_number   14
_chemical_formula_structural   Co2H2CO5
_chemical_formula_sum   'Co8 H8 C4 O20'
_cell_volume   379.81727776
_cell_formula_units_Z   4
loop_
 _symmetry_equiv_pos_site_id
 _symmetry_equiv_pos_as_xyz
  1  'x, y, z'
  2  '-x, -y, -z'
  3  '-x, y+1/2, -z+1/2'
  4  'x, -y+1/2, z+1/2'
loop_
 _atom_type_symbol
 _atom_type_oxidation_number
  Co2+  2.0
  H+  1.0
  C4+  4.0
  O2-  -2.0
loop_
 _atom_site_type_symbol
 _atom_site_label
 _atom_site_symmetry_multiplicity
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
  Co2+  Co0  4  0.18224700  0.21350700  0.98270900  1
  Co2+  Co1  4  0.49057900  0.11283400  0.71467500  1
  H+  H2  4  0.03165800  0.01858300  0.89184700  1
  H+  H3  4  0.06570400  0.60074500  0.00054500  1
  C4+  C4  4  0.49639900  0.13191800  0.25883600  1
  O2-  O5  4  0.01482100  0.09261900  0.85615500  1
  O2-  O6  4  0.01693100  0.65376900  0.92503900  1
  O2-  O7  4  0.47271500  0.23077900  0.31021700  1
  O2-  O8  4  0.47842800  0.55116900  0.15829600  1
  O2-  O9  4  0.49946500  0.12019400  0.12516400  1
```

```json
{
  "crystal_structure": {
    "lattice": {
      "lattice_kind": "Conventional",
      "a": "3.20 Å",
      "b": "12.31 Å",
      "c": "9.67 Å",
      "α": "90.00 º",
      "β": "93.05 º",
      "ɣ": "90.00 º",
      "Volume": "379.82 Å³"
    },
    "atomic_positions": [
      {
        "wyckoff": "4e",
        "element": "Co",
        "a": "0.817753",
        "b": "0.786493",
        "c": "0.017291"
      },
      {
        "wyckoff": "4e",
        "element": "Co",
        "a": "0.509421",
        "b": "0.887166",
        "c": "0.285325"
      },
      {
        "wyckoff": "4e",
        "element": "H",
        "a": "0.031658",
        "b": "0.018583",
        "c": "0.891847"
      },
      {
        "wyckoff": "4e",
        "element": "H",
        "a": "0.934296",
        "b": "0.100745",
        "c": "0.499455"
      },
      {
        "wyckoff": "4e",
        "element": "C",
        "a": "0.503601",
        "b": "0.631918",
        "c": "0.241164"
      },
      {
        "wyckoff": "4e",
        "element": "O",
        "a": "0.500535",
        "b": "0.620194",
        "c": "0.374836"
      },
      {
        "wyckoff": "4e",
        "element": "O",
        "a": "0.527285",
        "b": "0.730779",
        "c": "0.189783"
      },
      {
        "wyckoff": "4e",
        "element": "O",
        "a": "0.478428",
        "b": "0.551169",
        "c": "0.158296"
      },
      {
        "wyckoff": "4e",
        "element": "O",
        "a": "0.016931",
        "b": "0.846231",
        "c": "0.425039"
      },
      {
        "wyckoff": "4e",
        "element": "O",
        "a": "0.985179",
        "b": "0.907381",
        "c": "0.143845"
      }
    ],
    "symmetry": {
      "crystal_system": "Monoclinic",
      "lattice_system": "Monoclinic",
      "hall_number": "-P 2ybc",
      "international_number": 14,
      "symbol": "P2₁/c",
      "point_group": "2/m"
    },
    "structure_meta": {
      "number_of_atoms": 40,
      "density": "3.71 g·cm⁻³",
      "dimensionality": "",
      "possible_oxidation_states": "O²⁻, H⁺, C⁴⁺, Co²⁺"
    }
  }
}
```

Title               Co8 H8 C4 O20

Lattice type        P
Space group name    P 1
Space group number  1
Setting number      1

Lattice parameters

   a        b        c       alpha    beta     gamma
 3.19535 12.30609  9.67281  90.0000  93.0519  90.0000

Unit-cell volume = 379.817275 Å^3

```vasp
Co8 H8 C4 O20
1.0
   0.0000000000000000    3.1953509024188560    0.0000000000000000
  12.3060899999999993    0.0000000000000000    0.0000000000000000
   0.0000000000000000   -0.5149933576206951   -9.6590868701706771
Co H C O
8 8 4 20
direct
```

该晶胞的基矢向量近似正交，同样的需要转为正交晶胞

注意到经常出现Mn-O配位八面体的，而结构中Co-O配位四面体和八面体比例正好1:1，
所以计划将Co-O配位八面体位置的Co替换为Mn，得到CoMn(OH)2CO3

CMCH@CoNi-LDH-S这个目标算异质结，
拆解步骤：
1. 确定晶面，对于Ni(OH)2,层状结构暴露(001)面;对于Co(OH)2CO3,解离面是(010)面
2. 晶胞的正交转换,对于Ni(OH)2,直接六方转为正交即可，对于Co(OH)2CO3,需先转为正交晶胞,再将(010)面转为(001)面
3. 晶胞的扩包匹配晶格,对于Ni(OH)2,扩包需要考虑S掺杂浓度在1-3%之间(预计对正交晶胞扩包3*2*2)，对于Co(OH)2CO3,预计是正交晶胞的1*3*1,两者的扩包系数需要计算失配率来调整



a 失配	b 失配	面积失配
当前 3×2 vs 1×3	3.13%	12.15%	9.03%
优化 3×3 vs 1×5	3.12%	1.63%	1.50%


建议3*2*1和1*3*1，然后新的晶格常数取平均值，
2. 然后对扩包的3*2*1的Ni(OH)2，需要将一半的Ni替换为Co并且SQS（可以mcsqs相关代码适当参考tmp_gen_conioh.py），生成CoNiOH2-321.vasp，
3. 对扩包的1*3*1的Co(OH)2CO3，需要用户手动使用vesta将Co八面体替换为Mn不需要sqs，生成CoMn(OH)2CO3-131.vasp
4. 生成的模型是周期性的，要用户使用ms构建异质结