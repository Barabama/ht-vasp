#!/bin/bash
# push_continue.sh — 将 gap=1.0 / gap=2.0 的续跑任务推送到 B5070 集群并提交
#
# 用法:
#   ./push_continue.sh              # 推送 + 提交全部
#   ./push_continue.sh --submit     # 同默认 (推送 + 提交)
#   ./push_continue.sh --push-only  # 仅推送数据, 不提交
#   ./push_continue.sh --status     # 查看 B5070 上作业状态

set -e

REMOTE="mcmf507@10.144.144.11"
REMOTE_BASE="/workspace/gaominliang/ht-vasp"
LOCAL_BASE="/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S"
RSYNC_OPTS="-avz --progress --exclude=WAVECAR --exclude=CHGCAR --exclude=CHG --exclude=AECCAR*"

# ===== 需要推送的体系和续跑说明 =====
# gap=1.0: data/gap_scan/intrinsic_gap_1.0/ → 取 CONTCAR 作为 POSCAR 续跑 150 步
# gap=2.0: data/hetero_intrinsic/2-relax_2/  → 取 CONTCAR 作为 POSCAR 续跑 150 步

declare -A JOBS
JOBS["gap1.0"]="data/gap_scan/intrinsic_gap_1.0"
JOBS["gap2.0"]="data/hetero_intrinsic/2-relax_2"

# ===== 推送 =====
push_data() {
    echo "=== [$(date +%H:%M:%S)] 推送数据到 B5070 ==="

    for job_name in "${!JOBS[@]}"; do
        local_dir="${JOBS[$job_name]}"
        echo "  → $job_name ($local_dir)"

        # 目标: B5070 上以工作目录存储
        # gap1.0 → data/gap_scan/intrinsic_gap_1.0_continue/
        # gap2.0 → data/hetero_intrinsic/2-relax_2_continue/
        if [ "$job_name" = "gap1.0" ]; then
            remote_dir="$REMOTE_BASE/data/gap_scan/intrinsic_gap_1.0_continue"
        else
            remote_dir="$REMOTE_BASE/data/hetero_intrinsic/2-relax_2_continue"
        fi

        rsync $RSYNC_OPTS "$LOCAL_BASE/$local_dir/" "$REMOTE:$remote_dir/"
    done

    # 推送脚本自身
    rsync -avz "$0" "$REMOTE:$REMOTE_BASE/CMCH_at_CoNi-LDH-S/"
    echo "=== 推送完成 ==="
}

# ===== 在 B5070 上准备和提交 =====
submit_jobs() {
    echo "=== [$(date +%H:%M:%S)] 在 B5070 上准备续跑 ==="

    ssh "$REMOTE" bash << 'EOF'
set -e
BASE="/workspace/gaominliang/ht-vasp"
CONDA="$BASE/.conda"

# ---- gap=1.0 续跑 ----
echo "--- gap=1.0: 准备续跑 ---"
cd "$BASE/data/gap_scan/intrinsic_gap_1.0_continue"

# CONTCAR → POSCAR (续跑)
if [ -f "CONTCAR" ]; then
    cp CONTCAR POSCAR
    echo "  CONTCAR → POSCAR"
fi

# 修改 INCAR: NSW=150, EDIFF=1E-5, ISTART=0 (无WAVECAR, 从头开始)
if [ -f "INCAR" ]; then
    sed -i 's/^NSW\s*=.*/NSW = 150/' INCAR
    sed -i 's/^EDIFF\s*=.*/EDIFF = 1E-5/' INCAR
    # 确保有 ISTART
    if ! grep -q "^ISTART" INCAR; then
        echo "ISTART = 0" >> INCAR
    fi
    # 确保有 LDIPOL (slab 体系需要)
    if ! grep -q "^LDIPOL" INCAR; then
        echo "LDIPOL = .TRUE." >> INCAR
        echo "IDIPOL = 3" >> INCAR
    fi
    echo "  INCAR 已更新: NSW=150, EDIFF=1E-5"
fi

# 更新 KPOINTS
if [ -f "KPOINTS" ]; then
    # 确认 KPOINTS 合理
    echo "  KPOINTS 存在"
fi

# ---- gap=2.0 续跑 ----
echo "--- gap=2.0: 准备续跑 ---"
cd "$BASE/data/hetero_intrinsic/2-relax_2_continue"

# CONTCAR.gz → POSCAR
if [ -f "CONTCAR.gz" ]; then
    gunzip -c CONTCAR.gz > POSCAR
    echo "  CONTCAR.gz → POSCAR"
elif [ -f "CONTCAR" ]; then
    cp CONTCAR POSCAR
    echo "  CONTCAR → POSCAR"
fi

# 修改 INCAR: NSW=150, EDIFF=1E-5, ISTART=0
if [ -f "INCAR.gz" ]; then
    gunzip -c INCAR.gz > INCAR
fi
if [ -f "INCAR" ]; then
    sed -i 's/^NSW\s*=.*/NSW = 150/' INCAR
    sed -i 's/^EDIFF\s*=.*/EDIFF = 1E-5/' INCAR
    echo "  INCAR 已更新: NSW=150, EDIFF=1E-5"
fi

# ---- 提交 Slurm ----
echo ""
echo "提交作业到 B5070 Slurm..."

for job in gap1.0 gap2.0; do
    if [ "$job" = "gap1.0" ]; then
        WORKDIR="$BASE/data/gap_scan/intrinsic_gap_1.0_continue"
        JOB_NAME="gap1.0_continue"
    else
        WORKDIR="$BASE/data/hetero_intrinsic/2-relax_2_continue"
        JOB_NAME="gap2.0_continue"
    fi

    if [ ! -f "$WORKDIR/POSCAR" ]; then
        echo "  $WORKDIR/POSCAR 不存在, 跳过 $job"
        continue
    fi

    cd "$WORKDIR"
    CMD="cd $WORKDIR && module load vasp-gpu && srun vasp_std"
    sbatch --job-name="$JOB_NAME" \
           --output="$BASE/logs/${JOB_NAME}.log" \
           --mem=20G --nodes=1 --ntasks=1 --partition=partGPU \
           --time=100:00:00 --gpus-per-task=1 \
           --wrap="module load vasp-gpu; cd $WORKDIR && srun vasp_std"
    echo "  提交 $JOB_NAME → $WORKDIR"
done

echo "=== 全部提交完成 ==="
EOF
}

# ===== 查看 B5070 上作业状态 =====
check_status() {
    echo "=== B5070 作业队列 ==="
    ssh "$REMOTE" "squeue -u mcmf507" 2>/dev/null || echo "  (无法连接)"
    echo ""
    echo "=== B5070 上续跑目录结构 ==="
    ssh "$REMOTE" "ls -la $REMOTE_BASE/data/gap_scan/intrinsic_gap_1.0_continue/ && echo '---' && ls -la $REMOTE_BASE/data/hetero_intrinsic/2-relax_2_continue/" 2>/dev/null || echo "  (部分目录不存在)"
}

# ===== CLI =====
case "${1:---submit}" in
    --submit|--push)
        push_data
        submit_jobs
        ;;
    --push-only)
        push_data
        ;;
    --status)
        check_status
        ;;
    *)
        echo "用法: $0 [--submit|--push-only|--status]"
        exit 1
        ;;
esac
