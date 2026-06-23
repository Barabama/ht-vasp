#!/bin/bash
set -e

REMOTE="mcmf507@10.144.144.11"
REMOTE_BASE="/workspace/gaominliang/ht-vasp"
LOCAL_BASE="/nfs_hdd/2025/gaominliang/ht-vasp/CMCH_at_CoNi-LDH-S"

RSYNC_OPTS="-avz --progress --exclude=WAVECAR.gz --ignore-existing" #  --exclude=*_out.json

echo "=== Pulling data directories (bulk + slab results) ==="

# DATA_DIRS=(
#     CoMnH2CO5
#     # CoMnH2CO5-slab
#     CoNiOH2
#     # CoNiOH2-slab
#     CoNiOH2S-noH
#     # CoNiOH2S-noH-slab
#     # CoNiOH2S-noH-slab-flip
# )
DATA_DIRS=(
    CoNiHO
    CoNiHOS-Co3
    CoNiHOS-Ni3
    CoNiHOS-Co1Ni2
    CoNiHOS-Co2Ni1
    CoNiHOS-Co3-noH
)

for dir in "${DATA_DIRS[@]}"; do
    echo "[$(date +%H:%M:%S)] Pulling data/$dir ..."
    rsync $RSYNC_OPTS "$REMOTE:$REMOTE_BASE/data/$dir/" "$LOCAL_BASE/data/$dir/"
done

# echo ""
# echo "=== Pulling poscars ==="
# rsync $RSYNC_OPTS "$REMOTE:$REMOTE_BASE/data/poscars/" "$LOCAL_BASE/data/poscars/"

# echo ""
# echo "=== Pulling log files ==="

# LOG_FILES=(
#     CoMnH2CO5-nscf.log
#     # CoMnH2CO5-slab-nscf.log
#     CoNiOH2-nscf.log
#     # CoNiOH2-slab-nscf.log
#     CoNiOH2S-noH-nscf.log
#     # CoNiOH2S-noH-slab-nscf.log
#     # CoNiOH2S-noH-slab-flip-nscf.log
# )

# for log in "${LOG_FILES[@]}"; do
#     echo "[$(date +%H:%M:%S)] Pulling logs/$log ..."
#     scp "$REMOTE:$REMOTE_BASE/logs/$log" "$LOCAL_BASE/logs/$log"
# done

echo ""
echo "=== Done ==="
du -sh "$LOCAL_BASE/data/"*/
