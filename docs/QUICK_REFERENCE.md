# FireWorks + MongoDB 快速参考卡片

## 核心概念

| 组件 | 作用 | 配置文件 |
|------|------|---------|
| **LaunchPad** | 工作流编排引擎（连接 MongoDB） | `~/.fireworks/mcmf_launchpad.yaml` |
| **FWorker** | 计算节点身份标识 | `~/.fireworks/mcmf_fworker.yaml` |
| **QueueAdapter** | Slurm 调度器适配器 | `~/.fireworks/mcmf_qadapter.yaml` |
| **JobStore** | atomate2 数据存储抽象 | `configs/jobflow.yaml` |
| **Rocket** | 任务执行器 | `rlaunch rapidfire` |

---

## 常用命令速查

### LaunchPad 管理

```bash
# 查看工作流状态
lpad get_wflows
lpad get_wflows -s READY      # 仅查看就绪任务
lpad get_wflows -s RUNNING    # 仅查看运行中任务
lpad get_wflows -s COMPLETED  # 仅查看已完成任务
lpad get_wflows -s FIZZLED    # 仅查看失败任务

# 查看详细任务信息
lpad get_fws -i <fw_id> -d all

# 添加工作流
lpad add workflow.yaml

# 重置 LaunchPad（清空所有数据，谨慎！）
lpad reset

# 重新运行失败任务
lpad rerun_fws -s FIZZLED
lpad rerun_fws -i <fw_id>

# 检测丢失的运行任务
lpad detect_lost_runs
lpad revive
```

### Rocket 启动

```bash
# 前台运行（测试用）
rlaunch rapidfire

# 指定配置文件目录
rlaunch -c /path/to/configs rapidfire

# 限制启动次数
rlaunch rapidfire --nlaunches 10

# 调整轮询间隔（秒）
rlaunch rapidfire --sleep 5

# 使用特定 FWorker
rlaunch rapidfire --fworker node2_worker

# 后台运行（不推荐，优先使用 systemd）
nohup rlaunch rapidfire > worker.log 2>&1 &
```

### QueueAdapter 提交

```bash
# 手动提交单个任务到队列
qlaunch singleshot

# 批量提交多个任务
qlaunch multi_shot --nlaunches 5

# 查看待提交任务
qlaunch check
```

### MongoDB 管理

```bash
# 连接 MongoDB Shell
mongosh "mongodb://user:pass@192.168.1.101:27017/fireworks_db"

# 查看数据库
show dbs
use fireworks_db
show collections

# 查看集合大小
db.fireworks.countDocuments()
db.launches.countDocuments()

# 备份数据库
mongodump --host 192.168.1.101 --username admin --password pass \
          --authenticationDatabase admin --db fireworks_db \
          --out /backup/$(date +%Y%m%d)

# 恢复数据库
mongorestore --host 192.168.1.101 --username admin --password pass \
             --authenticationDatabase admin /backup/20260407
```

### systemd 服务管理

```bash
# 启动服务
sudo systemctl start fireworks-worker

# 停止服务
sudo systemctl stop fireworks-worker

# 重启服务
sudo systemctl restart fireworks-worker

# 查看状态
sudo systemctl status fireworks-worker

# 查看日志
journalctl -u fireworks-worker -f
journalctl -u fireworks-worker --since "1 hour ago"

# 开机自启
sudo systemctl enable fireworks-worker

# 禁用开机自启
sudo systemctl disable fireworks-worker
```

---

## 典型工作流程

### 1. 提交 VASP 计算任务

```python
from htvasp.slurm.fireworks_adapter import FireWorksAdapter
from pathlib import Path

# 创建适配器
adapter = FireWorksAdapter()

# 创建 VASP 任务
fw = adapter.create_vasp_firework(
    workdir=Path("/path/to/calculation"),
    vasp_cmd="srun vasp_std",
    name="my-vasp-job",
    conda_env="htvasp",
)

# 提交到 LaunchPad
wf_id = adapter.submit_workflow([fw])
print(f"Submitted workflow ID: {wf_id}")
```

### 2. 创建工作流（带依赖关系）

```python
from fireworks import Firework, ScriptTask, Workflow

# 创建任务
fw1 = Firework([ScriptTask.from_str("echo Step 1")], name="step1")
fw2 = Firework([ScriptTask.from_str("echo Step 2")], name="step2")
fw3 = Firework([ScriptTask.from_str("echo Step 3")], name="step3")

# 定义依赖：fw2 依赖 fw1，fw3 依赖 fw2
links = {
    fw2.fw_id: [fw1.fw_id],
    fw3.fw_id: [fw2.fw_id],
}

# 创建工作流
wf = Workflow([fw1, fw2, fw3], links_dict=links, name="my-workflow")

# 提交
from fireworks import LaunchPad
lp = LaunchPad.auto_load()
wf_id = lp.add_wf(wf)
```

### 3. 监控任务状态

```python
from htvasp.slurm.fireworks_adapter import FireWorksAdapter

adapter = FireWorksAdapter()

# 获取工作流状态
status = adapter.get_workflow_status(wf_id=123)
print(status)

# 输出示例：
# {
#     'wf_id': 123,
#     'name': 'my-workflow',
#     'states': {
#         1: {'name': 'step1', 'state': 'COMPLETED'},
#         2: {'name': 'step2', 'state': 'RUNNING'},
#         3: {'name': 'step3', 'state': 'WAITING'},
#     }
# }
```

### 4. 处理失败任务

```bash
# 查看所有失败任务
lpad get_wflows -s FIZZLED

# 查看失败原因
lpad get_fws -i <fw_id> -d all

# 重新运行单个任务
lpad rerun_fws -i <fw_id>

# 重新运行所有失败任务
lpad rerun_fws -s FIZZLED

# Python API
adapter.rerun_fizzled_workflows()
```

---

## 故障排查流程图

```
任务未执行？
├─ 检查 LaunchPad 连接
│  └─ lpad get_wflows
│     ├─ 成功 → 继续
│     └─ 失败 → 检查 mcmf_launchpad.yaml 和 MongoDB 状态
│
├─ 检查是否有 READY 任务
│  └─ lpad get_wflows -s READY
│     ├─ 有任务 → 继续
│     └─ 无任务 → 提交新任务或检查依赖
│
├─ 检查 Rocket 是否运行
│  └─ ps aux | grep rlaunch
│     ├─ 正在运行 → 继续
│     └─ 未运行 → 启动 rlaunch 或 systemd 服务
│
├─ 检查 QueueAdapter 配置
│  └─ cat ~/.fireworks/mcmf_qadapter.yaml
│     ├─ 配置正确 → 继续
│     └─ 配置错误 → 修正并重试
│
└─ 检查 Slurm 作业状态
   └─ squeue -u $(whoami)
      ├─ 作业在队列 → 等待调度
      ├─ 作业失败 → 查看错误日志
      └─ 无作业 → 检查 qlaunch 日志
```

---

## 常见错误及解决方案

| 错误信息 | 原因 | 解决方案 |
|---------|------|---------|
| `ServerSelectionTimeoutError` | MongoDB 无法连接 | 检查 MongoDB 服务、防火墙、mcmf_launchpad.yaml |
| `Authentication failed` | 用户名或密码错误 | 验证 mcmf_launchpad.yaml 中的凭据 |
| `Invalid partition name` | Slurm 分区名称错误 | 运行 `sinfo` 查看可用分区，更新 mcmf_qadapter.yaml |
| `No FireWorks are ready to run` | 没有 READY 状态的任务 | 提交新任务或检查任务依赖 |
| `sbatch: command not found` | Slurm 未安装或 PATH 错误 | 确认 Slurm 已安装且 PATH 正确 |
| `Permission denied` | 文件权限问题 | 检查工作目录和脚本的执行权限 |

---

## 性能调优参数

### MongoDB 优化 (`/etc/mongod.conf`)

```yaml
storage:
  wiredTiger:
    engineConfig:
      cacheSizeGB: 4  # 根据服务器内存调整（建议 25-50%）

operationProfiling:
  mode: slowOp
  slowOpThresholdMs: 100
```

### FireWorks 优化 (`mcmf_qadapter.yaml`)

```yaml
# 增加并发度
max_jobs: 10

# 减少轮询间隔
rapidfire:
  sleep_time: 5  # 默认 10 秒
```

### Slurm 优化 (`/etc/slurm/slurm.conf`)

```ini
SchedulerType=sched/backfill  # 使用后填算法提高利用率
SlurmctldParameters=idle_on_node_suspend
```

---

## 环境变量清单

```bash
# 必需的环境变量
export JOBFLOW_CONFIG_FILE=/path/to/jobflow.yaml
export ATOMATE2_CONFIG_FILE=/path/to/atomate2.yaml

# 可选的环境变量
export FIREWORKS_CONFIG_DIR=/path/to/configs  # FireWorks 配置目录
export OMP_NUM_THREADS=1                      # VASP 线程数
```

---

## 端口清单

| 服务 | 端口 | 协议 | 说明 |
|------|------|------|------|
| MongoDB | 27017 | TCP | 数据库访问 |
| slurmctld | 6817 | TCP | Slurm 控制器 |
| slurmd | 6818 | TCP | Slurm 计算节点守护进程 |
| munge | 6820 | TCP | 认证服务 |

---

## 有用的一行命令

```bash
# 实时监控 FireWorks 状态
watch -n 5 'lpad get_wflows'

# 清理 30 天前的日志
find logs/fireworks -name "*.log" -mtime +30 -delete

# 批量取消所有 Slurm 作业
scancel -u $(whoami)

# 查看 MongoDB 数据库大小
mongosh --quiet --eval "db.stats().dataSize" fireworks_db

# 导出所有完成的工作流
lpad dump_wflows -s COMPLETED -o completed.json

# 统计各状态任务数量
lpad get_wflows -s ALL | grep -oP 'state=\K\w+' | sort | uniq -c
```

---

## 联系与支持

- **项目文档**: `docs/fireworks-mongodb-deployment-guide.md`
- **配置示例**: `configs/README.md`
- **官方文档**:
  - [FireWorks](https://materialsproject.github.io/fireworks/)
  - [atomate2](https://materialsproject.github.io/atomate2/)
  - [jobflow](https://materialsproject.github.io/jobflow/)
