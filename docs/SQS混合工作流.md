# SQS 混合工作流

`main_sqs.py` 将有序端基和无序 SQS 分成两条路径。`Endmember` 继续保留有序、方向敏感的 `phase-A-B` 结构；`SqsSystem` 将 `phase-A-B` 与 `phase-B-A` canonicalize 为同一个无序二元体系。

## 模板

`data/templates/registry.json` 固化三个原型快照：

| phase | Materials Project ID | source species roles |
| --- | --- | --- |
| BCC | `mp-1487` | `Ni`, `Al` |
| FCC | `mp-2258` | `Au`, `Cu` |
| HCP | `mp-20112` | `Sn`, `Ni` |

普通端基和 SQS 生成只读本地快照，并验证 SHA256；不会在计算任务中访问网络。
`source_species` 是有序端基替换的角色映射，不是 MP 返回的 site 首次出现顺序。
模板更新必须显式调用 `TemplateRegistry.update_from_mp`，更新后应审阅 registry 和结构 diff。

也可以显式执行 `python main_sqs.py --sync-template BCC`（`FCC`/`HCP` 同理）；该命令是唯一会访问 Materials Project 的入口。

## SQS 数据流

1. `--generate` 使用版本化的默认超胞，并生成 `x=1/4, 1/2, 3/4` 三个精确可表示成分点：BCC 模板为 `2×2×2`（2-site 原胞，共 16 sites），FCC 模板为 `2×2×2`（4-site 原胞，共 32 sites），HCP 模板为 `2×2×1`（8-site 原胞，共 32 sites）。也可以用 `--supercell A B C` 显式指定；显式指定后仍必须使所有成分点能在该超胞中精确表示。
2. 每个纯元素参考和 SQS 点单独执行 `RelaxMaker -> StaticMaker`。这里不复用 `StaticWorker`，但 Worker 只负责 flow 组合，不负责 SQS 命名或热混合计算。
   Relax 使用 `ENCUT=400 eV`, `PREC=Normal`, `EDIFF=1e-5`, `EDIFFG=-0.05`, `ISIF=3`；Static 使用 `ENCUT=520 eV`, `PREC=Accurate`, `EDIFF=1e-6`。两阶段均保持 `ISPIN=2` 和初始 `MAGMOM`，Relax 写出 WAVECAR 供 Static 读取，Static 不再写 WAVECAR/CHGCAR/AECCAR，并保留 `LORBIT=10`。
3. `--collect` 从静态总能量重新计算

   `Hmix = (E_SQS - n_A E_A/atom - n_B E_B/atom) / n_sites`。

   输出能量单位为 `eV`，`hmix_eV_per_atom` 为 `eV/atom`，并按
   `96485.33212331002 J/(mol*eV)` 写出 `h_mix_J_per_mol`。结果使用
   `htvasp-sqs-mixing/v1` 契约，显式包含 `calphad_mapping`、`CalculationProvenance`、
   两份 `references`、每个点的 `StaticObservation`、SQS provenance 和 checksum；
   `em-tdb` 的 P2-A 只接受 BCC 的该 JSON 且当前生产拟合使用三个点；FCC/HCP
   目前可以生成和完成 VASP flow，但其 TDB mapping 需后续 schema 扩展，收集时会明确失败。

## 入口示例

```bash
python main_sqs.py --name BCC-Al-Nb --generate
python main_sqs.py --name BCC-Al-Nb --dry-run
python main_sqs.py --name BCC-Al-Nb
python main_sqs.py --name BCC-Al-Nb --collect
python main_sqs.py --name BCC-Al-Nb --check
```

SQS 生成依赖可选包 `sqsgenerator`（`pip install -e '.[sqs]'`）。缺少该依赖时直接失败，不使用随机回退结构。VASP、Slurm 和大规模生成不属于离线 smoke 验证范围。

当前 `sqsgenerator` API 没有公开的随机种子配置项，因此 `seed` 会写入
provenance，但不会被宣称为可复现性保证；若需要严格复现，应固定工具版本、配置、模板
checksum，并保留生成出的结构及其 `structure_sha256`。
