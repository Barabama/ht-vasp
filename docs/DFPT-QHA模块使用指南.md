# DFPT-QHA 模块使用指南

## 1. 它做什么

从 VASP 的 **DFPT**（`IBRION=8`）动力学矩阵出发，算每个体积的简谐声子自由能
$F(V,T)$，再做准谐分析得到 $G(T)$、$V(T)$ 和 EOS。

与 `QhaWorker`（有限位移 + `PhononMaker`）的区别只有一个来源：

| | QhaWorker | DfptQhaWorker |
|---|---|---|
| 力常数来源 | 有限位移（phonopy 生成位移，逐个静态计算） | VASP DFPT 的 Hessian（`vasprun.xml`） |
| 超胞上的计算次数 | 3N 个位移（对称约化后） | 每个体积 1 次 |
| 超胞 | atomate2 `PhononMaker` 生成 | **phonopy 自己的超胞**（顺序必须一致，见 §6） |

## 2. 流程图

```
EOS Maker (relax at each volume)
  -> spawn_dfpt_volume_jobs       动态展开，每个体积一个 VASP 作业
       -> collect_dfpt_task_doc   vasprun <dynmat><hessian> -> DfptTaskDoc
       -> validate_dfpt_task_doc   形状/原子数/配方一致性，失败即停
  -> collect_validated_volumes
  -> classify_dfpt_stability      逐体积精确 q 点判定；**不停流**
  -> dfpt_thermal_properties      F(V,T)
  -> analyze_dfpt_qha             G(T)、V(T)、EOS
```

## 3. 入口

```bash
# 单个端基
python main_dfpt_qha.py -n BCC-Ta-Ta-2 --dry-run
python main_dfpt_qha.py -n BCC-Ta-Ta-2 --slurm

# 批量（每行一个端基，# 注释）
python main_dfpt_qha.py --batch config/dfpt_qha_batches/batch-a-429pro.txt --slurm --ntasks 32

# 失败/超时后续跑（同一 run id，复用已完成体积）
python main_dfpt_qha.py -n BCC-Ta-Ta-2 --resume

# 查看结果
python main_dfpt_qha.py -n BCC-Ta-Ta-2 --check
```

端基的胞、配方、k 网格、POTCAR、展宽全部来自
`/nfs_ssd/tmp/legacy-qha-pool/stage0/pool_recipes.json`，入口脚本不重复定义。
超胞矩阵、strain 窗口、温度网格、DFPT 热积分 mesh 和 32 原子任务是入口参数。

`--slurm` 提交的是**整个 flow**（含后处理），不是单个 VASP 作业；分配内由
`run_locally_custom` 顺序执行。

## 4. 产物

`data/endpoints/<name>/dfpt-qhaflow/runs/<run_id>/<name>-dfpt-qha.json`

关键字段：

| 字段 | 含义 |
|---|---|
| `eos` | BM3 拟合：`V0_angstrom3`、`B0_GPa`、`B0prime`、`E0_eV` |
| `free_energies_eV` | 每个体积、每个温度的振动自由能（**每个配方胞**） |
| `gibbs_eV` / `equilibrium_volume_angstrom3` | 每个温度的最小 Gibbs 与对应体积 |
| `qha_applicability` | 端点级稳定性判定（见 §5） |
| `non_harmonic_volume_indices` | 哪些体积不满足简谐适用条件 |
| `dynamically_unstable_conditional` | 是否属于"条件产物" |
| `mesh_diagnostics` | 稠密 mesh 的最低频（仅诊断，不参与判定） |

单位：自由能 eV/配方胞，体积 Å³/配方胞，B0 GPa。

## 5. 稳定性判定与放行策略

判定在**超胞能严格表示的 q 点**上做（由 phonopy 的
`DynmatToForceConstants(primitive, supercell)` 给出，不是网格采样）。分档：

```
1. numerical_unresolved                结构/数值检查失败（形状、ASR、ordering、倒格矢）
2. gamma_soft_mode                     Γ 点光学支虚频
3. finite_q_soft_mode                  非 Γ 的精确 q 点有虚频
4. incommensurate_soft_mode_suspected  稠密 mesh 最低 < -0.5 THz，且无精确 q 点覆盖
5. harmonic_QHA_applicable             以上都没有
```

**判定不停流**：`dfpt_thermal_properties` 与 `analyze_dfpt_qha` 照常运行，
`qha_applicability` 与 `non_harmonic_volume_indices` 写进产物。原因：虚频是
物理结果，不是计算失败；停流会让每个软端基只剩诊断、没有 $G(T)$，而补算要
重新分析全部 store。

下游（em-tdb-v2）必须读 `qha_applicability` 再决定能否进入 production TDB。
非 `harmonic_QHA_applicable` 的产物是**条件产物**：简谐求和中被剔除的模式
仍然影响熵，不能当作无标记的 $G(T)$ 使用。

阈值：`SOFT_FREQUENCY_TOL_THz = 1e-4`，Γ 声学容差 `0.1 THz`（只给前三个
声学支），mesh 可疑阈值 `0.5 THz`，边缘软模确认阈值 `1.0 THz`。

## 6. 数据完整性

- **力常数来源**：优先 `vasprun.xml` 的 `<dynmat><hessian>`（完整、精确对称、
  与 phonopy 的约化位移列表无关）；缺失时回退 DYNMAT 文本回放。两者在归档
  算例上对称化后一致到 2.4e-9 eV/Å²。
- **求和规则**：原始 Hessian 的 ASR 残差可达 1e-2（健康）到 0.36（损坏）。
  collector 在残差 ≤ 0.1 时施加平移不变 + 置换对称修正，并把修正前后残差写进
  provenance；超过 0.1 不修正、原样交给分类作业判 `numerical_unresolved`。
- **原子顺序**：DFPT 作业必须跑在 `phonopy.supercell` 上（`supercell_structure`），
  分类作业再做一次 ordering 断言。pymatgen 的 `make_supercell` 顺序不同，
  跨顺序挂载力常数仍会给出"看起来合理"的频率。
- **输入审计**：`audit_vasp_inputs` 比较 `INCAR.orig/KPOINTS.orig/POSCAR.orig`
  与实际运行的文件，记录 custodian 改了什么；缺 `.orig` 只记 provenance 缺口，
  不判失败。
- **电子/物理门**：`TaskDoc.state == successful`；DFPT 基态能量与同体积静态
  能量相差 > 100 meV/atom 即拒绝（NELM 截断的发散 SCF 不会自己报错）。

## 7. 已知限制

- 2×2×2 超胞只能精确表示有限个 q 点；通过判定**不证明**非公度 q 稳定。
  产物里 `supercell_convergence: not_checked` 明说了这一点。
- 边缘软模只记录 `confirmation_required`，**不会**自动发起更大超胞的确认。
- `analyze_dfpt_qha` 输出的 `equilibrium_volume_angstrom3` 是采样体积上的
  argmin，**阶梯状**。连续曲线由 `continuous_thermal_eos` 给出：每个温度对
  F(V,T) 做 BM3 拟合，得到连续的 V(T)、G(T)、B(T)，不覆盖原 JSON。
  入口：`python main_dfpt_qha.py -n <name> --refit-vt`（或 `--batch`）。
  拟合残差是可用性的判据：本池里 `harmonic_QHA_applicable` 的端基残差
  0.3–2.1 meV/atom，`finite_q_soft_mode` 的端基升到 5–21 meV/atom，
  BCC-Ti-Ti-2 甚至出现负的体模量——两套独立的诊断（力常数判定与热力学
  积分）指向同一批端基。
- fcc A=A 端基的原胞是 1 原子、配方胞 4 原子，自由能换算靠
  `n_atoms_cell / n_atoms_primitive`；改这一项会让所有 fcc 端基差 4 倍。

## 8. 相关测试

```bash
.conda/bin/python -m pytest tests/test_dfpt_qha_*.py tests/test_input_audit.py \
    tests/test_local_resume.py tests/test_spin_screen.py tests/test_vasp_handlers.py -q
```

`tests/fixtures/BCC-Al-Ta-2-dfpt-volume-1/` 是一个真实的 DFPT 作业目录
（vasprun + DYNMAT + 输入），collector/ordering/stability 的测试不依赖 scratch。
