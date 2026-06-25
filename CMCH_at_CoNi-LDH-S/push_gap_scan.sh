#!/bin/bash
# push_gap_scan.sh — 推送 gap=0.5, gap=1.5 到 B5070 集群并提交
#
# 说明:
#   gap=1.0 → 本地 job 47620 (429Pro)
#   gap=2.0 → 本地 job 47621 (429Pro)
#
# 用法:
#   ./push_gap_scan.sh          # 推送 + 提交
#   ./push_gap_scan.sh --status # 查看状态

set -e

REMOTE="mcmf507@10.144.144.11"
REMOTE_BASE="/workspace/gaominliang/ht-vasp"
LOCAL_BASE="/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S"
RSYNC_OPTS="-avz --progress"

echo "============================================"
echo "  推送 gap=0.5, gap=1.5 → B5070"
echo "============================================"

# ===== 1. 推送 =====
echo ""
echo "=== 推送数据 ==="

for g in 0.5 1.5; do
    echo "  → gap=${g}"
    rsync $RSYNC_OPTS "$LOCAL_BASE/data/gap_scan/intrinsic_gap_${g}/" \
        "$REMOTE:$REMOTE_BASE/data/gap_scan/intrinsic_gap_${g}/"
done

echo "=== 推送完成 ==="

# ===== 2. 远程提交 =====
echo ""
echo "=== 在 B5070 上提交 ==="

ssh "$REMOTE" bash << EOSSH
set -e
BASE="$REMOTE_BASE"
LOG_DIR="\$BASE/CMCH_at_CoNi-LDH-S/logs"
mkdir -p "\$LOG_DIR"

for g in 0.5 1.5; do
    WORKDIR="\$BASE/data/gap_scan/intrinsic_gap_\${g}"
    JOB_NAME="gap\${g}_fresh"

    if [ ! -f "\$WORKDIR/POSCAR" ]; then
        echo "  ⚠ gap=\${g}: POSCAR 不存在, 跳过"
        continue
    fi

    # 确认 INCAR 正确
    echo "  gap=\${g}: \$(grep 'EDIFF' \$WORKDIR/INCAR), \$(grep 'NSW' \$WORKDIR/INCAR)"

    sbatch --job-name="\$JOB_NAME" \
           --output="\$LOG_DIR/\${JOB_NAME}.log" \
           --mem=20G --nodes=1 --ntasks=1 --partition=partGPU \
           --time=100:00:00 --gpus-per-task=1 \
           --wrap=". /etc/profile.d/modules.sh && ulimit -s unlimited && module purge && module load vasp-gpu && cd \$WORKDIR && srun vasp_std"
    echo "  ✅ gap=\${g} 已提交"
done

echo "=== B5070 提交完成 ==="
EOSSH

echo ""
echo "============================================"
echo "  全部完成"
echo "  本地: gap=1.0 (47620), gap=2.0 (47621)"
echo "  B5070: gap=0.5, gap=1.5"
echo "============================================"
