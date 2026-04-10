# FireWorks + MongoDB Configuration Files

This directory contains configuration templates and examples for integrating FireWorks with MongoDB and Slurm.

## Directory Structure

```
configs/
├── fireworks/
│   ├── mcmf_launchpad.yaml.template          # LaunchPad (MongoDB) connection
│   ├── mcmf_fworker_node1.yaml.template      # FWorker config for node1
│   ├── mcmf_fworker_node2.yaml.template      # FWorker config for node2
│   └── mcmf_qadapter.yaml.template           # Slurm QueueAdapter config
├── systemd/
│   ├── fireworks-worker.service.template       # systemd service for node1
│   └── fireworks-worker-node2.service.template # systemd service for node2
├── jobflow.yaml.template                 # atomate2 JobStore (MongoDB)
├── atomate2.yaml.template                # atomate2 global settings
└── README.md                             # This file
```

## Quick Start

### 1. Copy Templates to Active Configurations

```bash
# On node1 (master with MongoDB)
cp configs/fireworks/mcmf_launchpad.yaml.template ~/.fireworks/mcmf_launchpad.yaml
cp configs/fireworks/mcmf_fworker_node1.yaml.template ~/.fireworks/mcmf_fworker.yaml
cp configs/fireworks/mcmf_qadapter.yaml.template ~/.fireworks/mcmf_qadapter.yaml

# On node2 (compute node)
cp configs/fireworks/mcmf_launchpad.yaml.template ~/.fireworks/mcmf_launchpad.yaml
cp configs/fireworks/mcmf_fworker_node2.yaml.template ~/.fireworks/mcmf_fworker.yaml
cp configs/fireworks/mcmf_qadapter.yaml.template ~/.fireworks/mcmf_qadapter.yaml
```

### 2. Update Configuration Values

Edit `~/.fireworks/mcmf_launchpad.yaml` and update:
- `host`: MongoDB server IP (e.g., `192.168.1.101`)
- `username`: Your MongoDB username
- `password`: Your MongoDB password (**IMPORTANT: Change default!**)

Edit `~/.fireworks/mcmf_fworker.yaml` and update:
- `name`: Unique worker name per node
- `env`: Environment variables for your system

### 3. Setup MongoDB (Node1 Only)

See the deployment guide for detailed instructions:
```bash
docs/fireworks-mongodb-deployment-guide.md
```

Or use the automated setup script:
```bash
scripts/setup_fireworks.sh
```

### 4. Test Connection

```bash
# Test LaunchPad connection
lpad get_wflows

# Should output: Found 0 workflows
```

### 5. Submit a Test Workflow

```python
from fireworks import Firework, ScriptTask

fw = Firework([ScriptTask.from_str('echo "Hello FireWorks!"')])
from fireworks import LaunchPad
lp = LaunchPad.auto_load()
lp.add_wf(fw)
print(f"Submitted Firework with ID: {fw.fw_id}")
```

### 6. Start the Worker

Option A: Using systemd (recommended for production)
```bash
sudo systemctl enable fireworks-worker
sudo systemctl start fireworks-worker
sudo systemctl status fireworks-worker
```

Option B: Manual (for testing)
```bash
rlaunch rapidfire
```

## Configuration File Reference

### mcmf_launchpad.yaml

Connects FireWorks to MongoDB.

**Key fields:**
- `host`: MongoDB server IP
- `port`: MongoDB port (default: 27017)
- `name`: Database name
- `username`: MongoDB username
- `password`: MongoDB password
- `logdir`: Log directory path

### mcmf_fworker.yaml

Identifies the compute node to FireWorks.

**Key fields:**
- `name`: Unique worker identifier
- `category`: Task category filter ('' for all)
- `query`: MongoDB query filter ('{}' for none)
- `env`: Environment variables

### mcmf_qadapter.yaml

Defines how FireWorks submits jobs to Slurm.

**Key fields:**
- `_fw_name`: Adapter type (CommonAdapter)
- `rocket_launch`: Command to execute tasks
- `submit_cmd`: Slurm submission command (sbatch)
- `template`: Slurm job script template

### jobflow.yaml

Configures atomate2 JobStore with MongoDB backend.

**Key sections:**
- `STORE.docs_store`: MongoDB collection for metadata
- `STORE.additional_stores`: GridFS for large data

## Troubleshooting

### MongoDB Connection Failed

```bash
# Check MongoDB is running
sudo systemctl status mongod

# Check firewall
sudo ufw status

# Test connection manually
mongosh "mongodb://user:pass@192.168.1.101:27017/fireworks_db"
```

### FireWorks Cannot Connect

```bash
# Verify credentials
cat ~/.fireworks/mcmf_launchpad.yaml

# Test with lpad command
lpad get_wflows -d all
```

### Slurm Jobs Not Starting

```bash
# Check QueueAdapter config
cat ~/.fireworks/mcmf_qadapter.yaml

# Check Slurm partition names
sinfo

# View qlaunch logs
tail -f logs/fireworks/qlaunch/*.err
```

## Advanced Usage

### Multiple Workers Per Node

Create multiple FWorker configs with different categories:

```yaml
# mcmf_fworker_vasp.yaml
name: node1_vasp
category: vasp
```

Submit tasks with matching category:
```python
fw = Firework([...], spec={"_category": "vasp"})
```

### Custom Slurm Resources

Modify `mcmf_qadapter.yaml` template section:

```yaml
template: |
  #!/bin/bash
  #SBATCH --gres=gpu:1        # Request GPU
  #SBATCH --constraint=highmem # Node constraint
  ...
```

### Monitoring

```bash
# Watch FireWorks status
watch -n 5 'lpad get_wflows'

# Monitor Slurm queue
watch -n 5 'squeue -u $(whoami)'

# View worker logs
tail -f logs/fireworks/worker.log
```

## Security Notes

1. **Never commit passwords** to version control
2. Use strong passwords for MongoDB users
3. Restrict MongoDB access with UFW firewall
4. Enable MongoDB authentication
5. Use TLS/SSL for remote connections (optional)

## Further Reading

- [FireWorks Documentation](https://materialsproject.github.io/fireworks/)
- [atomate2 Documentation](https://materialsproject.github.io/atomate2/)
- [Deployment Guide](../docs/fireworks-mongodb-deployment-guide.md)
