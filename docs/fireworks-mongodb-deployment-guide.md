# FireWorks + MongoDB + Slurm 集成部署指南

## 目录

1. [架构与协作机制分析](#1-架构与协作机制分析)
2. [Ubuntu 24.04 双节点部署方案](#2-ubuntu-2404-双节点部署方案)
3. [安全性与网络配置](#3-安全性与网络配置)
4. [与 htvasp 模块集成](#4-与-htvasp-模块集成)
5. [故障排查与最佳实践](#5-故障排查与最佳实践)

---

## 1. 架构与协作机制分析

### 1.1 FireWorks LaunchPad 与 MongoDB 的协同工作

FireWorks 采用 **LaunchPad** 作为工作流编排引擎，其底层完全依赖 MongoDB 进行状态管理和任务调度。

#### 核心组件关系图

```
┌─────────────────────────────────────────────────────────────┐
│                     atomate2 Workflows                       │
│                  (JobFlow-based workflows)                   │
└──────────────────┬──────────────────────────────────────────┘
                   │
                   ▼
┌─────────────────────────────────────────────────────────────┐
│                    JobStore Abstraction                      │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐  │
│  │  DocsStore   │    │  DataStore   │    │  Additional  │  │
│  │ (MongoDB or  │◄──►│ (GridFS or   │◄──►│   Stores     │  │
│  │  JSONStore)  │    │  MongoDB)    │    │              │  │
│  └──────────────┘    └──────────────┘    └──────────────┘  │
└──────────────────┬──────────────────────────────────────────┘
                   │
                   ▼
┌─────────────────────────────────────────────────────────────┐
│                   FireWorks LaunchPad                        │
│  ┌──────────────────────────────────────────────────────┐   │
│  │          MongoDB Collections (Workflow State)         │   │
│  │  ┌──────────┐ ┌──────────┐ ┌──────────┐            │   │
│  │  │ fireworks│ │ launches │ │  locks   │            │   │
│  │  │(Tasks)   │◄►(History) │◄►(Mutex)   │            │   │
│  │  └──────────┘ └──────────┘ └──────────┘            │   │
│  │  ┌──────────┐ ┌──────────┐                          │   │
│  │  │fireworks_│ │fireworks_│                          │   │
│  │  │  OFFLINE │ │  ARCHIVE │                          │   │
│  │  └──────────┘ └──────────┘                          │   │
│  └──────────────────────────────────────────────────────┘   │
└──────────────────┬──────────────────────────────────────────┘
                   │
                   ▼
┌─────────────────────────────────────────────────────────────┐
│                    FWorker + QueueAdapter                    │
│  ┌──────────────┐    ┌──────────────┐                       │
│  │   FWorker    │    │QueueAdapter  │                       │
│  │  (Identity)  │◄──►│ (Slurm/PBS)  │                       │
│  └──────────────┘    └──────────────┘                       │
└──────────────────┬──────────────────────────────────────────┘
                   │
                   ▼
┌─────────────────────────────────────────────────────────────┐
│                    Slurm Scheduler                           │
│  ┌──────────┐  ┌──────────┐  ┌──────────┐                 │
│  │ Node 1   │  │ Node 2   │  │   ...    │                 │
│  │(Master)  │  │(Compute) │  │          │                 │
│  └──────────┘  └──────────┘  └──────────┘                 │
└─────────────────────────────────────────────────────────────┘
```

#### MongoDB 数据库集合详解

LaunchPad 在 MongoDB 中创建以下核心集合（collections）：

| 集合名称 | 用途 | 关键字段 |
|---------|------|---------|
| `fireworks` | 存储所有待执行的任务（Firework）状态 | `fw_id`, `state` (WAITING/READY/RUNNING/COMPLETED/FIZZLED), `spec` (任务参数), `created_on`, `updated_on` |
| `launches` | 记录每次任务执行的历史和日志 | `launch_id`, `fw_id`, `state`, `launch_dir`, `action` (执行的命令), `stdout`, `stderr` |
| `offline_launches` | 离线模式下的任务队列 | 同 `launches` |
| `fireworks_OFFLINE` | 离线模式下的 Firework 状态 | 同 `fireworks` |
| `fireworks_ARCHIVED` | 已完成任务的归档 | 历史数据 |
| `locks` | 分布式锁，防止多 Worker 竞争同一任务 | `name`, `lock_id`, `expiration` |

**状态流转机制：**

```
WAITING → READY → RUNNING → COMPLETED
                    ↓
                 FIZZLED (失败)
                    ↓
                 DEFUSED (手动取消)
```

1. **任务提交阶段**：通过 `lpad add <workflow.yaml>` 将 Firework 写入 `fireworks` 集合，初始状态为 `WAITING` 或 `READY`
2. **Rocket 拉取阶段**：`rlaunch rapidfire` 从 LaunchPad 查询 `state: READY` 的任务，原子性地将其状态改为 `RUNNING`（通过 MongoDB 的 `findAndModify` 操作确保不重复领取）
3. **任务执行阶段**：Rocket 在工作目录执行任务脚本，实时将 stdout/stderr 写回 `launches` 集合
4. **完成阶段**：任务成功则状态更新为 `COMPLETED`，失败则为 `FIZZLED`（可配置自动重试）

### 1.2 atomate2 JobStore 抽象层

atomate2 使用 **maggma** 库的 `JobStore` 作为数据持久化抽象层，支持多种后端存储：

#### JobStore 架构

```python
from maggma.stores import MongoStore, GridFSStore, JSONStore, MemoryStore
from jobflow import JobStore

# 方案 1: 纯 MongoDB 存储（推荐用于集群环境）
store = JobStore(
    docs_store=MongoStore(database="atomate_db", collection_name="docs"),
    additional_stores={
        "data": GridFSStore(database="atomate_db", collection_name="blobs"),
    }
)

# 方案 2: 混合存储（当前项目使用的轻量级方案）
store = JobStore(
    JSONStore("store.json"),           # 元数据存本地 JSON
    additional_stores={"data": MemoryStore()}  # 大数据存内存
)
```

**关键区别：**

| 特性 | JSONStore + MemoryStore | MongoStore + GridFSStore |
|------|------------------------|--------------------------|
| 适用场景 | 单机开发/测试 | 生产集群/多用户协作 |
| 数据持久化 | 单文件 JSON | 分布式 MongoDB |
| 断点续算 | ✅ 支持 | ✅ 支持 |
| 多节点访问 | ❌ 不支持 | ✅ 支持 |
| 大数据存储 | ❌ 内存限制 | ✅ GridFS 分片存储 |
| 性能 | 快（无网络开销） | 中等（依赖网络） |

**与 FireWorks 的交互流程：**

1. atomate2 workflow 生成 VASP 输入文件
2. FireWorks Rocket 执行 VASP 计算
3. 计算结果通过 `JobStore` 写入 MongoDB：
   - 小数据（结构、能量等）→ `docs_store` (MongoDB 文档)
   - 大数据（波函数、电荷密度）→ `additional_stores["data"]` (GridFS)
4. 后续任务从 MongoDB 读取前序任务结果

### 1.3 FireWorks 如何调度任务到 Slurm 节点

FireWorks 通过 **QueueAdapter** 实现与 Slurm 的集成：

```
┌─────────────────────────────────────────────────────────┐
│  rlaunch rapidfire                                      │
│  ├── 从 LaunchPad 获取 READY 任务                         │
│  ├── 通过 QueueAdapter 生成 Slurm 脚本                    │
│  ├── 执行 sbatch 提交作业                                 │
│  └── 轮询 squeue 检查作业状态                             │
└─────────────────────────────────────────────────────────┘
                    │
                    ▼
┌─────────────────────────────────────────────────────────┐
│  Slurm Controller (slurmctld)                           │
│  ├── 接收 sbatch 请求                                    │
│  ├── 根据资源需求调度到可用节点                            │
│  └── 启动 slurmd 执行任务                                │
└─────────────────────────────────────────────────────────┘
                    │
                    ▼
┌─────────────────────────────────────────────────────────┐
│  Compute Node (slurmd)                                  │
│  ├── 创建工作目录                                        │
│  ├── 执行 VASP 计算                                      │
│  └── 返回退出码                                          │
└─────────────────────────────────────────────────────────┘
```

---

## 2. Ubuntu 24.04 双节点部署方案

### 2.1 环境假设

```
节点规划：
├── 429pro (master / 192.168.9.144)
│   ├── Slurm Controller (slurmctld)
│   ├── MongoDB Server
│   └── FireWorks LaunchPad + Worker (429pro_worker)
│
└── 429e (compute / 192.168.9.143)
    ├── Slurm Worker (slurmd)
    └── FireWorks Worker (429e_worker)

共享存储：
└── NFS: /home/mcmf429/nfs_hdd (已挂载到两个节点)
```

### 2.2 MongoDB 部署

#### 步骤 1: 选择部署节点

**推荐：在主节点 (node1) 部署 MongoDB**

理由：
- MongoDB 需要稳定的网络和磁盘 I/O
- 主节点通常负载较低，适合运行数据库服务
- 避免计算节点重启影响数据库可用性
- 双节点环境下无需独立存储节点

#### 步骤 2: 安装 MongoDB 7.0 (Ubuntu 24.04)

```bash
# 在 node1 上执行
sudo apt update
sudo apt install -y gnupg curl

# 导入 MongoDB GPG 密钥
curl -fsSL https://www.mongodb.org/static/pgp/server-7.0.asc | \
   sudo gpg -o /usr/share/keyrings/mongodb-server-7.0.gpg --dearmor

# 添加 MongoDB 源
echo "deb [ signed-by=/usr/share/keyrings/mongodb-server-7.0.gpg ] \
http://repo.mongodb.org/apt/ubuntu jammy/mongodb-org/7.0 multiverse" | \
sudo tee /etc/apt/sources.list.d/mongodb-org-7.0.list

# 安装 MongoDB
sudo apt update
sudo apt install -y mongodb-org

# 启动并设置开机自启
sudo systemctl enable mongod
sudo systemctl start mongod

# 验证安装
mongosh --eval "db.version()"
```

#### 步骤 3: 配置 MongoDB 用户认证

```bash
# 进入 MongoDB Shell
mongosh

# 创建管理员用户
use admin
db.createUser({
  user: "admin",
  pwd: "your_secure_admin_password",
  roles: [ { role: "userAdminAnyDatabase", db: "admin" } ]
})

# 创建 atomate2 专用数据库和用户
use atomate_db
db.createUser({
  user: "atomate_user",
  pwd: "your_atomate_password",
  roles: [
    { role: "readWrite", db: "atomate_db" },
    { role: "dbAdmin", db: "atomate_db" }
  ]
})

# 创建 FireWorks 专用数据库和用户
use fireworks_db
db.createUser({
  user: "fireworks_user",
  pwd: "your_fireworks_password",
  roles: [
    { role: "readWrite", db: "fireworks_db" },
    { role: "dbAdmin", db: "fireworks_db" }
  ]
})

exit
```

#### 步骤 4: 配置 MongoDB 网络访问

编辑 `/etc/mongod.conf`:

```yaml
# network interfaces
net:
  port: 27017
  bindIp: 127.0.0.1,192.168.1.101  # 监听本地和局域网 IP

# security
security:
  authorization: enabled  # 启用认证
```

重启 MongoDB：

```bash
sudo systemctl restart mongod
```

#### 步骤 5: 测试远程连接

```bash
# 在 node2 上测试连接
mongosh "mongodb://atomate_user:your_atomate_password@192.168.1.101:27017/atomate_db"

# 如果连接成功，应该看到：
# Current Mongosh Log ID: ...
# Connecting to: mongodb://...
# Using MongoDB: 7.0.x
```

### 2.3 FireWorks 配置

#### 步骤 1: 安装 FireWorks

```bash
# 在所有节点上执行
conda activate htvasp
pip install fireworks
```

#### 步骤 2: 创建 mcmf_launchpad.yaml

在 **所有节点** 的 `~/.fireworks/` 目录下创建配置文件：

```bash
mkdir -p ~/.fireworks
```

文件内容 (`~/.fireworks/mcmf_launchpad.yaml`):

```yaml
# LaunchPad 配置 - 连接到 node1 的 MongoDB
host: 192.168.1.101
port: 27017
name: fireworks_db
username: fireworks_user
password: your_fireworks_password
logdir: /home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/logs/fireworks
strm_lvl: INFO
user_indices: []
wf_user_indices: []
```

**字段说明：**

| 字段 | 说明 | 示例值 |
|------|------|--------|
| `host` | MongoDB 服务器 IP | `192.168.1.101` |
| `port` | MongoDB 端口 | `27017` |
| `name` | 数据库名称 | `fireworks_db` |
| `username` | 认证用户名 | `fireworks_user` |
| `password` | 认证密码 | `your_fireworks_password` |
| `logdir` | 日志目录 | 建议使用 NFS 共享路径 |
| `strm_lvl` | 日志级别 | `INFO` / `DEBUG` / `WARNING` |

#### 步骤 3: 验证 LaunchPad 连接

```bash
# 测试连接
lpad get_wflows

# 应该输出：
# Found 0 workflows

# 重置 LaunchPad（清空所有任务，谨慎使用）
lpad reset
```

#### 步骤 4: 创建 mcmf_fworker.yaml

在每个计算节点上创建 FWorker 配置 (`~/.fireworks/mcmf_fworker.yaml`):

```yaml
# FWorker 配置 - 标识计算节点
name: node1_worker  # 每个节点使用不同的名称
category: ''
query: '{}'
env:
  vasp_cmd: "srun vasp_std"
  vasp_gam_cmd: "srun vasp_gam"
  atomate2_config_path: "/home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/atomate2.yaml"
```

**字段说明：**

| 字段 | 说明 | 示例值 |
|------|------|--------|
| `name` | Worker 唯一标识 | `node1_worker`, `node2_worker` |
| `category` | 任务分类过滤 | `''` (接受所有任务) 或 `'vasp'` |
| `query` | MongoDB 查询过滤器 | `'{}'` (无过滤) |
| `env` | 环境变量 | VASP 命令路径等 |

#### 步骤 5: 创建 mcmf_qadapter.yaml (Slurm 集成)

这是**最关键**的配置文件，定义如何将 FireWorks 任务转换为 Slurm 作业。

文件位置：`~/.fireworks/mcmf_qadapter.yaml`

```yaml
_fw_name: CommonAdapter
rocket_launch: rlaunch -c /home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/configs/fireworks rapidfire
submit_cmd: sbatch
status_cmd: squeue -j $JOB_ID
jobs_file: jobs.json
remove_submit_file: true
logdir: /home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/logs/fireworks/qlaunch
remove_launch_file: true

# Slurm 脚本模板
template: |
  #!/bin/bash

  #SBATCH -J {{name}}
  #SBATCH -o {{logdir}}/{{name}}.out
  #SBATCH -e {{logdir}}/{{name}}.err
  #SBATCH -N {{nodes}}
  #SBATCH -n {{ntasks}}
  #SBATCH -c {{cpus_per_task}}
  #SBATCH -p {{partition}}
  #SBATCH -t {{time_limit}}
  #SBATCH --mem={{memory}}

  # 加载环境
  . /opt/miniconda3/etc/profile.d/conda.sh
  conda activate htvasp

  # 切换到工作目录
  cd {{launch_dir}}

  # 执行 Rocket
  {{rocket_launch}}
```

**模板变量说明：**

| 变量 | 说明 | 默认值 |
|------|------|--------|
| `{{name}}` | 作业名称 | Firework ID |
| `{{logdir}}` | 日志目录 | 配置中的 `logdir` |
| `{{launch_dir}}` | 任务工作目录 | 自动生成 |
| `{{nodes}}` | 节点数 | 见下方自定义 |
| `{{ntasks}}` | 任务数 | 见下方自定义 |
| `{{cpus_per_task}}` | 每任务 CPU 数 | 见下方自定义 |
| `{{partition}}` | Slurm 分区 | `partCPU` |
| `{{time_limit}}` | 时间限制 | `100:00:00` |
| `{{memory}}` | 内存限制 | `20G` |

#### 步骤 6: 自定义队列适配器参数

如果需要更精细的控制，可以创建 `mcmf_qadapter.yaml` 的高级版本：

```yaml
_fw_name: CommonAdapter
rocket_launch: rlaunch -c /home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/configs/fireworks rapidfire
submit_cmd: sbatch
status_cmd: squeue -j $JOB_ID
jobs_file: /home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/configs/fireworks/jobs.json
remove_submit_file: true
logdir: /home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/logs/fireworks/qlaunch
remove_launch_file: true

# 自定义提交参数
nodes: 1
ntasks: 48
cpus_per_task: 1
partition: partCPU
time_limit: '100:00:00'
memory: 20G

# 预执行命令
pre_rocket: |
  . /opt/miniconda3/etc/profile.d/conda.sh
  conda activate htvasp
  module purge
  module load vasp/6.4.3

template: |
  #!/bin/bash

  #SBATCH -J {{name}}
  #SBATCH -o {{logdir}}/{{name}}.out
  #SBATCH -e {{logdir}}/{{name}}.err
  #SBATCH -N {{nodes}}
  #SBATCH -n {{ntasks}}
  #SBATCH -c {{cpus_per_task}}
  #SBATCH -p {{partition}}
  #SBATCH -t {{time_limit}}
  #SBATCH --mem={{memory}}

  {% if pre_rocket %}
  {{pre_rocket}}
  {% endif %}

  cd {{launch_dir}}
  {{rocket_launch}}
```

### 2.4 配置 atomate2 JobStore 使用 MongoDB

在项目根目录创建 `jobflow.yaml`:

```yaml
# jobflow.yaml - atomate2 JobStore 配置
STORE:
  docs_store:
    type: MongoStore
    database: atomate_db
    collection_name: docs
    host: 192.168.1.101
    port: 27017
    username: atomate_user
    password: your_atomate_password
    key: uuid

  additional_stores:
    data:
      type: GridFSStore
      database: atomate_db
      collection_name: blobs
      host: 192.168.1.101
      port: 27017
      username: atomate_user
      password: your_atomate_password
      key: job_uuid

JOB_STORE: "${STORE}"
```

同时创建 `atomate2.yaml`:

```yaml
# atomate2.yaml - atomate2 全局配置
VASP_CMD: "srun vasp_std"
VASP_GAM_CMD: "srun vasp_gam"
VASP_NCL_CMD: "srun vasp_ncl"

# 其他配置...
```

### 2.5 测试 FireWorks + Slurm 集成

#### 测试 1: 创建简单测试任务

创建 `test_workflow.py`:

```python
"""Test FireWorks + Slurm integration."""

from fireworks import Firework, FWorker, LaunchPad, ScriptTask, Workflow

# 初始化 LaunchPad
lp = LaunchPad.auto_load()

# 创建简单任务
fw = Firework(
    [ScriptTask.from_str('echo "Hello from FireWorks on Slurm!"')],
    name="test-fireworks-slurm",
    spec={"_allow_fizzled_parents": True}
)

# 添加到 LaunchPad
lp.add_wf(fw)
print(f"Added Firework with fw_id: {fw.fw_id}")
```

运行：

```bash
cd /home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp
python test_workflow.py
```

#### 测试 2: 通过 QueueAdapter 提交到 Slurm

```bash
# 使用 qlaunch 提交任务到队列
qlaunch -c /home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/configs/fireworks singleshot

# 查看 Slurm 作业状态
squeue -u $(whoami)

# 查看 FireWorks 任务状态
lpad get_wflows -d all
```

#### 测试 3: 启动 rapidfire 监听服务

```bash
# 前台运行（测试用）
rlaunch -c /home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/configs/fireworks rapidfire

# 这会持续监听新任务并自动提交到 Slurm
```

---

## 3. 安全性与网络配置

### 3.1 防火墙配置 (UFW)

在 **node1 (MongoDB 服务器)** 上配置：

```bash
# 启用 UFW
sudo ufw enable

# 允许 SSH
sudo ufw allow ssh

# 允许 Slurm 通信 (假设使用默认端口)
sudo ufw allow 6817/tcp  # slurmctld
sudo ufw allow 6818/tcp  # slurmd

# 仅允许受信任的计算节点访问 MongoDB
sudo ufw allow from 192.168.1.102 to any port 27017 proto tcp

# 拒绝其他所有访问 MongoDB 的请求
sudo ufw deny 27017/tcp

# 查看规则
sudo ufw status verbose
```

### 3.2 MongoDB 安全加固

#### 1. 禁用 HTTP 接口

在 `/etc/mongod.conf` 中确认：

```yaml
net:
  http:
    enabled: false
```

#### 2. 启用 TLS/SSL (可选，高安全要求)

```yaml
net:
  ssl:
    mode: requireSSL
    PEMKeyFile: /etc/ssl/mongodb.pem
    CAFile: /etc/ssl/ca.pem
```

#### 3. 限制 MongoDB 绑定 IP

```yaml
net:
  bindIp: 127.0.0.1,192.168.1.101  # 不要使用 0.0.0.0
```

### 3.3 Slurm 安全配置

在 `/etc/slurm/slurm.conf` 中：

```ini
# 仅允许特定用户提交作业
AuthType=auth/munge
AuthInfo=/var/run/munge/munge.socket.2

# 限制资源使用
DefMemPerCPU=2048
MaxMemPerCPU=4096

# 会计收集（可选）
AccountingStorageType=accounting_storage/filetxt
AccountingStorageLoc=/var/log/slurm/accounting.log
```

---

## 4. 与 htvasp 模块集成

### 4.1 现有代码分析

您当前的 `htvasp/slurm/manager.py` 实现了直接的 `sbatch --wrap` 提交方式：

```python
# 当前实现（简化）
def submit_command(self, command: str, config: SlurmConfig, ...):
    cmd_parts = [
        "sbatch",
        f"--job-name={config.job_name}",
        f'--wrap="{command}"',
        ...
    ]
    subprocess.run(cmd, shell=True, ...)
```

**优点：**
- 简单直接，无需额外配置
- 适合单次任务提交

**缺点：**
- 无法利用 FireWorks 的工作流编排能力
- 不支持自动重试、断点续算
- 难以管理复杂依赖关系

### 4.2 集成方案 A: 保留现有逻辑 + FireWorks 并行

**适用场景：** 逐步迁移，保持向后兼容

#### 步骤 1: 创建 FireWorks 适配器类

在 `htvasp/slurm/` 下创建 `fireworks_adapter.py`:

```python
"""FireWorks integration adapter for htvasp."""

from pathlib import Path
from typing import Any

from fireworks import Firework, LaunchPad, ScriptTask, Workflow
from fireworks.core.rocket_launcher import rapidfire

from htvasp.slurm.manager import SlurmConfig, SlurmJobManager


class FireWorksAdapter:
    """Adapter to integrate htvasp with FireWorks."""

    def __init__(self, config_dir: Path | str | None = None):
        """Initialize FireWorks adapter.

        Args:
            config_dir: Directory containing FireWorks config files.
                       If None, uses ~/.fireworks
        """
        self.config_dir = Path(config_dir) if config_dir else None
        self.lp = LaunchPad.auto_load()
        self.slurm_manager = SlurmJobManager()

    def create_vasp_firework(
        self,
        workdir: Path | str,
        vasp_cmd: str = "srun vasp_std",
        name: str = "vasp-job",
        conda_env: str = "htvasp",
        slurm_config: SlurmConfig | None = None,
    ) -> Firework:
        """Create a Firework for VASP calculation.

        Args:
            workdir: Working directory for VASP calculation.
            vasp_cmd: VASP command to execute.
            name: Firework name.
            conda_env: Conda environment name.
            slurm_config: Slurm configuration.

        Returns:
            Firework object.
        """
        workdir = Path(workdir).resolve()

        # Create execution script
        script = f"""#!/bin/bash
. /opt/miniconda3/etc/profile.d/conda.sh
conda activate {conda_env}
cd {workdir}
{vasp_cmd}
"""

        # Write script to workdir
        script_path = workdir / "run_vasp.sh"
        script_path.write_text(script)
        script_path.chmod(0o755)

        # Create Firework
        fw = Firework(
            [ScriptTask.from_str(f"bash {script_path}")],
            name=name,
            spec={
                "_allow_fizzled_parents": True,
                "metadata": {
                    "workdir": str(workdir),
                    "slurm_config": slurm_config.__dict__ if slurm_config else {},
                }
            }
        )

        return fw

    def submit_workflow(
        self,
        fireworks: list[Firework],
        dependencies: dict[int, list[int]] | None = None,
    ) -> int:
        """Submit workflow to LaunchPad.

        Args:
            fireworks: List of Firework objects.
            dependencies: Dict mapping fw_id to list of parent fw_ids.

        Returns:
            Workflow ID.
        """
        if len(fireworks) == 1 and not dependencies:
            # Single Firework
            wf_id = self.lp.add_wf(fireworks[0])
        else:
            # Create workflow with dependencies
            wf = Workflow(fireworks, links_dict=dependencies or {})
            wf_id = self.lp.add_wf(wf)

        return wf_id

    def launch_rapidfire(
        self,
        nlaunches: int = -1,  # -1 means infinite
        sleep_time: int = 10,
        timeout: int | None = None,
    ):
        """Launch rapidfire listener.

        Args:
            nlaunches: Number of launches (-1 for infinite).
            sleep_time: Sleep time between checks (seconds).
            timeout: Timeout in seconds (None for no timeout).
        """
        rapidfire(
            self.lp,
            fworker=None,  # Use default FWorker
            nlaunches=nlaunches,
            sleep_time=sleep_time,
            timeout=timeout,
        )
```

#### 步骤 2: 修改现有工作流以支持 FireWorks

在 `htvasp/workflows/base.py` 中添加 FireWorks 支持：

```python
"""Base worker class with FireWorks support."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Literal

from htvasp.slurm.fireworks_adapter import FireWorksAdapter


class Worker(ABC):
    """Abstract base class for workers."""

    @abstractmethod
    def __init__(
        self,
        worker_name: str,
        vasp_args: dict[str, Any],
        potcar_functional: Literal["PBE", "PBE_54", "PBE_64"] = "PBE_64",
        global_incar: dict[str, Any] | None = None,
        use_fireworks: bool = False,  # New parameter
        fireworks_config_dir: Path | str | None = None,
        **kwargs,
    ):
        self.worker_name = worker_name
        self.use_fireworks = use_fireworks

        if use_fireworks:
            self.fw_adapter = FireWorksAdapter(fireworks_config_dir)
        else:
            from htvasp.slurm.manager import SlurmJobManager
            self.slurm_manager = SlurmJobManager()

    @abstractmethod
    def run_flow(
        self,
        name: str,
        structure: Structure,
        flowdir: Path | str,
        dir_format: str = "{name}",
        store_path: str | None = None,
    ) -> dict[str, Any] | None:
        pass

    def _submit_to_fireworks(
        self,
        workdir: Path,
        command: str,
        name: str,
        **kwargs,
    ) -> int:
        """Submit job via FireWorks.

        Args:
            workdir: Working directory.
            command: Command to execute.
            name: Job name.

        Returns:
            Firework ID.
        """
        fw = self.fw_adapter.create_vasp_firework(
            workdir=workdir,
            vasp_cmd=command,
            name=name,
            **kwargs,
        )
        return self.fw_adapter.submit_workflow([fw])
```

### 4.3 集成方案 B: 完全迁移到 FireWorks (推荐)

**适用场景：** 新项目或愿意重构的项目

#### 步骤 1: 使用 atomate2 原生 FireWorks 支持

atomate2 已经内置了 FireWorks 集成。修改您的工作流以使用 `atomate2.vasp.jobs` 中的 Maker 类，它们天然支持 FireWorks。

示例（OJ 工作流改造）：

```python
"""OJ workflow with FireWorks support."""

from pathlib import Path
from typing import Any

from atomate2.vasp.jobs.base import BaseVaspMaker
from jobflow import Flow, Job
from maggma.stores import MongoStore
from jobflow import JobStore

from htvasp.oj.maker import OJMaker


class OJWorkerFireWorks:
    """OJ Worker using FireWorks backend."""

    def __init__(
        self,
        worker_name: str,
        vasp_args: dict[str, Any],
        oj_incar: dict[str, Any] | None = None,
        j_count: int = 4,
        extend_poscar: tuple[int, int, int] = (2, 2, 2),
    ):
        self.worker_name = worker_name
        self.oj_maker = OJMaker(
            vasp_args=vasp_args,
            oj_incar=oj_incar,
            j_count=j_count,
            extend_poscar=extend_poscar,
        )

        # Use MongoDB-backed JobStore
        self.store = JobStore(
            docs_store=MongoStore(
                database="atomate_db",
                collection_name="oj_docs",
                host="192.168.1.101",
                port=27017,
                username="atomate_user",
                password="your_atomate_password",
            ),
        )

    def create_flow(self, name: str, structure, flowdir: Path | str) -> Flow:
        """Create OJ workflow as a Flow.

        Args:
            name: Calculation name.
            structure: Input structure.
            flowdir: Output directory.

        Returns:
            JobFlow Flow object.
        """
        flowdir = Path(flowdir) / name
        flowdir.mkdir(parents=True, exist_ok=True)

        # Create OJ job
        oj_job = self.oj_maker.make(structure)
        oj_job.name = f"{name}-oj"

        # Create Flow
        flow = Flow([oj_job], output=oj_job.output)
        flow.name = f"{name}-oj-flow"

        return flow

    def submit_to_fireworks(self, flow: Flow) -> int:
        """Submit Flow to FireWorks.

        Args:
            flow: JobFlow Flow object.

        Returns:
            Workflow ID.
        """
        from fireworks import Workflow
        from atomate2.common.flows import flow_to_workflow

        # Convert JobFlow to FireWorks Workflow
        wf = flow_to_workflow(flow, self.store)

        # Submit to LaunchPad
        from fireworks import LaunchPad
        lp = LaunchPad.auto_load()
        wf_id = lp.add_wf(wf)

        return wf_id
```

#### 步骤 2: 创建便捷的提交脚本

创建 `scripts/submit_fireworks.py`:

```python
#!/usr/bin/env python
"""Submit htvasp workflows to FireWorks."""

import argparse
from pathlib import Path

from pymatgen.core import Structure

from htvasp.workflows.oj import OJWorkerFireWorks


def main():
    parser = argparse.ArgumentParser(description="Submit workflows to FireWorks")
    parser.add_argument("poscar", type=Path, help="POSCAR file")
    parser.add_argument("--name", required=True, help="Calculation name")
    parser.add_argument("--flowdir", type=Path, default="./flows", help="Flow directory")
    parser.add_argument("--worker", choices=["oj", "relax", "static"], required=True,
                        help="Worker type")

    args = parser.parse_args()

    # Load structure
    struct = Structure.from_file(args.poscar)

    # Create worker based on type
    if args.worker == "oj":
        worker = OJWorkerFireWorks(
            worker_name=args.name,
            vasp_args={"encut": 520},
            j_count=4,
        )
    else:
        raise ValueError(f"Unsupported worker type: {args.worker}")

    # Create and submit flow
    flow = worker.create_flow(args.name, struct, args.flowdir)
    wf_id = worker.submit_to_fireworks(flow)

    print(f"Submitted workflow with ID: {wf_id}")


if __name__ == "__main__":
    main()
```

使用：

```bash
python scripts/submit_fireworks.py POSCAR --name test-oj --worker oj
```

---

## 5. 原子服务守护：后台运行 rlaunch

### 5.1 方案 A: systemd 服务（推荐）

#### 步骤 1: 创建 systemd 服务文件

创建 `/etc/systemd/system/fireworks-worker.service`:

```ini
[Unit]
Description=FireWorks Rocket Launcher Worker
After=network.target mongod.service slurmctld.service
Wants=mongod.service

[Service]
Type=simple
User=mcmf429
Group=mcmf429
WorkingDirectory=/home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp

# 加载 conda 环境
Environment="PATH=/opt/miniconda3/envs/htvasp/bin:/opt/miniconda3/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin"
Environment="CONDA_DEFAULT_ENV=htvasp"

# 启动命令
ExecStart=/opt/miniconda3/envs/htvasp/bin/rlaunch -c /home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/configs/fireworks rapidfire --nlaunches infinite --sleep 10

# 重启策略
Restart=always
RestartSec=10

# 日志
StandardOutput=append:/home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/logs/fireworks/worker.log
StandardError=append:/home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/logs/fireworks/worker.err

# 资源限制
LimitNOFILE=65536
MemoryMax=4G

[Install]
WantedBy=multi-user.target
```

#### 步骤 2: 启用并启动服务

```bash
# 重新加载 systemd
sudo systemctl daemon-reload

# 启用开机自启
sudo systemctl enable fireworks-worker

# 启动服务
sudo systemctl start fireworks-worker

# 查看状态
sudo systemctl status fireworks-worker

# 查看日志
journalctl -u fireworks-worker -f
```

#### 步骤 3: 为多个节点创建独立服务

在 **node2** 上创建 `/etc/systemd/system/fireworks-worker-node2.service`:

```ini
[Unit]
Description=FireWorks Rocket Launcher Worker (Node 2)
After=network.target slurmd.service

[Service]
Type=simple
User=mcmf429
Group=mcmf429
WorkingDirectory=/home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp

Environment="PATH=/opt/miniconda3/envs/htvasp/bin:/opt/miniconda3/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin"
Environment="CONDA_DEFAULT_ENV=htvasp"

# 使用 node2 特定的 FWorker 配置
ExecStart=/opt/miniconda3/envs/htvasp/bin/rlaunch -c /home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/configs/fireworks rapidfire --fworker node2_worker --nlaunches infinite --sleep 10

Restart=always
RestartSec=10

StandardOutput=append:/home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/logs/fireworks/worker-node2.log
StandardError=append:/home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/logs/fireworks/worker-node2.err

LimitNOFILE=65536
MemoryMax=4G

[Install]
WantedBy=multi-user.target
```

### 5.2 方案 B: tmux/screen 会话（开发调试用）

```bash
# 使用 tmux
tmux new -s fireworks-worker
conda activate htvasp
rlaunch -c configs/fireworks rapidfire --nlaunches infinite --sleep 10
# 按 Ctrl+B 然后 D 分离会话

# 重新附加
tmux attach -t fireworks-worker

# 使用 screen
screen -S fireworks-worker
conda activate htvasp
rlaunch -c configs/fireworks rapidfire --nlaunches infinite --sleep 10
# 按 Ctrl+A 然后 D 分离

# 重新附加
screen -r fireworks-worker
```

### 5.3 方案 C: nohup（最简单但不推荐用于生产）

```bash
nohup rlaunch -c configs/fireworks rapidfire --nlaunches infinite --sleep 10 > logs/fireworks/worker.log 2>&1 &

# 记录 PID
echo $! > logs/fireworks/worker.pid

# 停止
kill $(cat logs/fireworks/worker.pid)
```

---

## 6. 完整部署检查清单

### 6.1 前置条件

- [ ] Ubuntu 24.04 双节点系统已安装
- [ ] 基本系统工具已安装：`curl`, `gnupg`（脚本会自动检测并安装）
- [ ] Slurm 集群已配置并正常运行
- [ ] NFS 共享存储已挂载
- [ ] Conda 环境 `htvasp` 已创建
- [ ] VASP 已安装并可正常执行

### 6.2 MongoDB 部署

- [ ] MongoDB 7.0 已安装在 node1
- [ ] 管理员用户和数据库用户已创建
- [ ] `/etc/mongod.conf` 已配置（bindIp, authorization）
- [ ] MongoDB 服务正在运行：`systemctl status mongod`
- [ ] 远程连接测试成功

### 6.3 FireWorks 配置

- [ ] FireWorks 已安装：`pip install fireworks`
- [ ] `~/.fireworks/mcmf_launchpad.yaml` 已创建
- [ ] `~/.fireworks/mcmf_fworker.yaml` 已创建（每个节点）
- [ ] `~/.fireworks/mcmf_qadapter.yaml` 已创建
- [ ] `lpad get_wflows` 测试通过
- [ ] `qlaunch singleshot` 测试通过

### 6.4 atomate2 配置

- [ ] `jobflow.yaml` 已创建（MongoDB JobStore）
- [ ] `atomate2.yaml` 已创建（VASP 命令路径）
- [ ] 环境变量 `JOBFLOW_CONFIG_FILE` 指向 `jobflow.yaml`

### 6.5 安全配置

- [ ] UFW 防火墙已启用
- [ ] MongoDB 端口仅对受信任 IP 开放
- [ ] MongoDB 认证已启用
- [ ] Slurm Munge 密钥已正确配置

### 6.6 守护进程

- [ ] systemd 服务文件已创建
- [ ] 服务已启用并启动：`systemctl status fireworks-worker`
- [ ] 日志正常输出

### 6.7 端到端测试

- [ ] 提交测试任务：`python test_workflow.py`
- [ ] 任务被 Slurm 调度执行
- [ ] 结果正确写入 MongoDB
- [ ] 任务状态可在 LaunchPad 中查询

---

## 7. 常见问题与故障排查

### 问题 1: MongoDB 连接超时

**症状：**
```
pymongo.errors.ServerSelectionTimeoutError: 192.168.1.101:27017: timed out
```

**解决：**
```bash
# 检查 MongoDB 是否运行
sudo systemctl status mongod

# 检查防火墙
sudo ufw status

# 检查 MongoDB 绑定的 IP
grep bindIp /etc/mongod.conf

# 测试网络连接
ping 192.168.1.101
nc -zv 192.168.1.101 27017
```

### 问题 2: FireWorks 认证失败

**症状：**
```
pymongo.errors.OperationFailure: Authentication failed
```

**解决：**
```bash
# 验证凭据
mongosh "mongodb://fireworks_user:your_password@192.168.1.101:27017/fireworks_db"

# 检查 mcmf_launchpad.yaml 中的密码是否正确
cat ~/.fireworks/mcmf_launchpad.yaml
```

### 问题 3: Slurm 作业提交失败

**症状：**
```
sbatch: error: Batch job submission failed: Invalid partition name specified
```

**解决：**
```bash
# 查看可用分区
sinfo

# 检查 mcmf_qadapter.yaml 中的 partition 配置
grep partition ~/.fireworks/mcmf_qadapter.yaml

# 确保分区名称匹配（区分大小写）
```

### 问题 4: rlaunch 找不到任务

**症状：**
```
No FireWorks are ready to run
```

**解决：**
```bash
# 检查是否有 WAITING/READY 任务
lpad get_wflows -s WAITING
lpad get_wflows -s READY

# 检查 FWorker category 是否匹配
cat ~/.fireworks/mcmf_fworker.yaml

# 手动触发一个任务
lpad detect_lost_runs
lpad revive
```

### 问题 5: 任务一直处于 FIZZLED 状态

**症状：**
```
lpad get_wflows -s FIZZLED 显示大量失败任务
```

**解决：**
```bash
# 查看详细错误信息
lpad get_fws -i <fw_id> -d all

# 查看执行日志
cat logs/fireworks/qlaunch/<job_name>.err

# 重置并重试
lpad rerun_fws -i <fw_id>

# 批量重置所有 FIZZLED 任务
lpad rerun_fws -s FIZZLED
```

---

## 8. 性能优化建议

### 8.1 MongoDB 优化

```yaml
# /etc/mongod.conf
storage:
  dbPath: /var/lib/mongodb
  journal:
    enabled: true
  wiredTiger:
    engineConfig:
      cacheSizeGB: 4  # 根据服务器内存调整

operationProfiling:
  mode: slowOp
  slowOpThresholdMs: 100
```

### 8.2 FireWorks 优化

```yaml
# mcmf_qadapter.yaml
# 增加并发度
max_jobs: 10  # 同时运行的最大作业数

# 调整轮询间隔
rapidfire:
  sleep_time: 5  # 减少等待时间
  nlaunches: -1  # 无限运行
```

### 8.3 Slurm 优化

```ini
# /etc/slurm/slurm.conf
# 增加调度器线程数
SchedulerType=sched/backfill

# 优化作业启动
SlurmctldParameters=idle_on_node_suspend
```

---

## 9. 监控与维护

### 9.1 监控脚本

创建 `scripts/monitor_fireworks.sh`:

```bash
#!/bin/bash
# Monitor FireWorks status

echo "=== FireWorks Status ==="
echo "Total workflows: $(lpad get_wflows -c all | wc -l)"
echo "Running: $(lpad get_wflows -s RUNNING | wc -l)"
echo "Completed: $(lpad get_wflows -s COMPLETED | wc -l)"
echo "FIZZLED: $(lpad get_wflows -s FIZZLED | wc -l)"

echo ""
echo "=== Slurm Jobs ==="
squeue -u $(whoami)

echo ""
echo "=== MongoDB Status ==="
mongosh --quiet --eval "db.stats()" fireworks_db

echo ""
echo "=== Disk Usage ==="
du -sh /home/mcmf429/nfs_hdd/2025/gaominliang/ht-vasp/logs/fireworks/
```

### 9.2 定期维护任务

```bash
# 每周清理旧日志
find logs/fireworks -name "*.out" -mtime +30 -delete
find logs/fireworks -name "*.err" -mtime +30 -delete

# 每月归档完成的 FireWorks
lpad archive_wflows -s COMPLETED -o archived_workflows.json

# 备份 MongoDB
mongodump --host 192.168.1.101 --username admin --password your_password \
          --authenticationDatabase admin --db fireworks_db \
          --out /backup/mongodb/$(date +%Y%m%d)
```

---

## 10. 总结

本指南详细说明了如何在 Ubuntu 24.04 双节点 Slurm 集群上部署和集成 FireWorks + MongoDB + atomate2。关键要点：

1. **架构理解**：FireWorks 通过 LaunchPad 管理 MongoDB 中的任务状态，通过 QueueAdapter 与 Slurm 集成
2. **MongoDB 部署**：推荐在主节点部署，配置认证和网络访问控制
3. **FireWorks 配置**：三个核心配置文件（mcmf_launchpad.yaml, mcmf_fworker.yaml, mcmf_qadapter.yaml）
4. **与 htvasp 集成**：提供两种方案，推荐逐步迁移到 atomate2 原生 FireWorks 支持
5. **守护进程**：使用 systemd 服务确保高可用性
6. **安全加固**：配置 UFW 防火墙和 MongoDB 认证

按照本指南操作，您可以构建一个可扩展、可靠的分布式材料计算工作流平台。
