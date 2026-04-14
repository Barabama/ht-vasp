# run_flow() 目录管理优化设计

## 背景

`Worker.run_flow()` 方法需要优化目录管理策略，以提高计算速度并支持从持久化存储恢复任务。

## 目标

1. 使用快速存储（如 `/tmp`、SSD）执行 VASP 计算，提高 I/O 性能
2. 计算完成后将结果移动到持久化存储（如 NFS）
3. 支持 resume 功能：从持久化存储恢复已完成的任务

## 设计决策

| 决策点 | 选择 | 理由 |
|--------|------|------|
| flow_dir 命名 | `{name}-{flow.uuid[:8]}` | 兼顾可读性和唯一性 |
| resume 恢复源 | 从 store_dir 恢复 | store_dir 是持久化存储 |
| 异常处理 | 无论成功失败都移动到 store_dir | 失败数据也有调试价值 |
| store_dir 结构 | `store_dir/{name}/` | 简洁直观，同名覆盖 |

## 修改范围

仅修改 `htvasp/workflows/base.py` 中的 `run_flow()` 方法。

## 流程设计

```
┌─────────────────────────────────────────────────────────────┐
│  1. 生成 flow 和目录路径                                      │
│     - self.flow = self._make_flow(structure)                │
│     - flow_dir = Path(flow_dir) / f"{name}-{self.flow.uuid[:8]}" │
│     - final_store_dir = Path(store_dir) / name              │
├─────────────────────────────────────────────────────────────┤
│  2. Resume 处理                                              │
│     if resume and final_store_dir.exists():                 │
│         shutil.copytree(final_store_dir, flow_dir)          │
├─────────────────────────────────────────────────────────────┤
│  3. 执行 Flow                                                │
│     run_locally_custom(                                      │
│         self.flow,                                           │
│         store=self.store,                                    │
│         root_dir=flow_dir,                                   │
│         resume=resume,                                       │
│         ...                                                  │
│     )                                                        │
├─────────────────────────────────────────────────────────────┤
│  4. 移动结果到 store_dir                                      │
│     if final_store_dir.exists():                             │
│         shutil.rmtree(final_store_dir)                       │
│     shutil.move(flow_dir, final_store_dir)                   │
└─────────────────────────────────────────────────────────────┘
```

## 代码变更

### run_flow() 方法

```python
def run_flow(
    self,
    name: str,
    structure: Structure,
    flow_dir: Path | str = "/tmp",
    store_dir: Path | str = ".",
    ensure_success: bool = True,
    raise_immediately: bool = False,
    resume: bool = True,
):
    """
    Run the flow for the structure.

    Args:
        name: Name of the system
        structure: Structure to run the worker on
        flow_dir: Flow directory for running the flow (fast storage, e.g., /tmp)
        store_dir: Directory to store the results (persistent storage)
        ensure_success: Raise an error if the flow did not finish successfully
        raise_immediately: Raise an error immediately if a job fails
        resume: Resume from previous completed jobs
    """
    import shutil
    import uuid

    self.flow = self._make_flow(structure)

    # 生成目录路径
    flow_dir = Path(flow_dir) / f"{name}-{self.flow.uuid[:8]}"
    final_store_dir = Path(store_dir).resolve() / name

    # Resume: 从 store_dir 复制到 flow_dir
    if resume and final_store_dir.exists():
        log.info(f"Resuming from {final_store_dir}, copying to {flow_dir}")
        shutil.copytree(final_store_dir, flow_dir)

    # 初始化 store
    self.store = JobStore(
        JSONStore(str(flow_dir / "store.json"), read_only=False),
        additional_stores={"data": MemoryStore()},
    )

    log.info(f"Running flow in {flow_dir}")
    try:
        run_locally_custom(
            self.flow,
            store=self.store,
            root_dir=flow_dir,
            ensure_success=ensure_success,
            raise_immediately=raise_immediately,
            resume=resume,
        )
    except Exception as e:
        log.critical(f"Flow execution failed: {e}")
        raise
    finally:
        # 无论成功失败，移动 flow_dir 到 store_dir
        try:
            if final_store_dir.exists():
                shutil.rmtree(final_store_dir)
            shutil.move(str(flow_dir), str(final_store_dir))
            log.info(f"Moved results to {final_store_dir}")
        except Exception as move_error:
            log.error(f"Failed to move results to store_dir: {move_error}")
```

## 异常处理

1. **执行失败**：捕获异常后重新抛出，但 `finally` 块确保目录移动
2. **移动失败**：记录错误日志但不中断流程（避免掩盖原始异常）

## 测试用例

1. **首次运行**：flow_dir 不存在，store_dir 不存在 → 正常执行并移动
2. **Resume 成功**：store_dir 存在 → 复制到 flow_dir → 继续执行
3. **执行失败**：flow_dir 仍移动到 store_dir
4. **同名覆盖**：store_dir 已存在同名目录 → 删除后移动

## 兼容性

- 现有调用代码无需修改（参数默认值保持不变）
- `run_locally_custom()` 无需修改
