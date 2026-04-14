# run_flow() 目录管理优化实现计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 优化 `run_flow()` 方法的目录管理，实现快速存储执行、持久化存储保存、以及从持久化存储恢复任务的功能。

**Architecture:** 在 `Worker.run_flow()` 中实现完整的目录生命周期管理：生成 flow_dir 路径 → resume 时从 store_dir 复制 → 执行 flow → 移动到 store_dir。修改范围仅限 `htvasp/workflows/base.py`。

**Tech Stack:** Python 3.10+, pathlib, shutil, jobflow

---

## 文件结构

| 文件 | 职责 |
|------|------|
| `htvasp/workflows/base.py` | 修改 `run_flow()` 方法，实现目录管理逻辑 |
| `tests/test_relax.py` | 更新测试用例以验证新行为 |

---

### Task 1: 实现 run_flow() 目录管理逻辑

**Files:**
- Modify: `htvasp/workflows/base.py:105-156`

- [ ] **Step 1: 替换 run_flow() 方法实现**

将现有的 `run_flow()` 方法替换为新实现：

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

- [ ] **Step 2: 验证语法正确**

Run: `python -c "from htvasp.workflows.base import Worker; print('Import OK')"`
Expected: "Import OK"

- [ ] **Step 3: 提交更改**

```bash
git add htvasp/workflows/base.py
git commit -m "feat: optimize run_flow() directory management

- Use flow_dir with name-uuid[:8] subdirectory naming
- Copy from store_dir to flow_dir on resume
- Move flow_dir to store_dir after execution (success or failure)
- Ensure results are always persisted to store_dir"
```

---

### Task 2: 更新测试用例

**Files:**
- Modify: `tests/test_relax.py:151-171`

- [ ] **Step 1: 更新 test_relax_locally() 测试函数**

修改测试函数以适配新的目录结构：

```python
def test_relax_locally():
    """Run Relax workflow locally (requires VASP)"""
    structure = get_al_bcc_structure()
    flow_dir = Path("/nfs_ssd/tmp")
    store_dir = Path("./temp")
    json_path = store_dir / "relax-Al" / "relax_Al.json"

    # 清理旧数据
    if (store_dir / "relax-Al").exists():
        shutil.rmtree(store_dir / "relax-Al")

    worker = RelaxWorker(
        vasp_args={
            "vasp_cmd": ["vasp_std"],
        },
    )

    result = worker.run_flow(
        name="relax-Al",
        structure=structure,
        flow_dir=flow_dir,
        store_dir=store_dir,
    )

    # 验证结果保存在 store_dir 下
    assert (store_dir / "relax-Al").exists(), "Results should be in store_dir"
    print(f"Results saved to: {store_dir / 'relax-Al'}")
```

- [ ] **Step 2: 提交测试更新**

```bash
git add tests/test_relax.py
git commit -m "test: update test_relax_locally for new directory structure"
```

---

## 自检清单

- [ ] spec 中所有需求都有对应任务实现
- [ ] 没有使用占位符（TBD、TODO 等）
- [ ] 类型和方法名称在所有任务中保持一致
- [ ] 每个步骤都有完整的代码或命令
