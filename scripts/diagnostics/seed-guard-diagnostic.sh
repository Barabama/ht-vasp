#!/usr/bin/env bash
#SBATCH --job-name=seeddiag
#SBATCH --output=/nfs_hdd/2025/gaominliang/ht-vasp/logs/seeddiag-%j.out
#SBATCH --error=/nfs_hdd/2025/gaominliang/ht-vasp/logs/seeddiag-%j.err
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=partCPU
#SBATCH --time=00:10:00
#
# Read-only diagnostic: reproduce step 4's charge-seed guard inside a real
# Slurm allocation on the same node the pilot uses, and dump everything the
# guard looks at. No VASP, no writes to the seed directory.
set -uo pipefail
umask 027

HTVASP_ROOT=/nfs_hdd/2025/gaominliang/ht-vasp
PYTHON="${HTVASP_ROOT}/.conda/bin/python"
PREP="${HTVASP_ROOT}/scripts/legacy_dfpt_qha_pilot_bcc_ti_v.py"
V="${HTVASP_ROOT}/data/legacy-dfpt-qha-pilot/BCC-Ti-V/legacy-fixed-grid/runs/pilot-20261005T145928Z/volumes/QHA-01"
seed_root="$V/charge-seed"

module purge
module load vasp-cpu
ulimit -s unlimited
export PMG_VASP_PSP_DIR=/nfs_ssd/documents I_MPI_FABRICS=ofi
cd "$HTVASP_ROOT"

echo "hostname=$(hostname)"
echo "cwd=$(pwd)"
echo "seed_root=$seed_root"
echo "--- ls -la seed_root ---"
ls -la "$seed_root"
echo "--- stat CHGCAR ---"
stat -c 'size=%s inode=%i mtime=%y %n' "$seed_root/CHGCAR" || echo "stat FAILED"
echo "--- test -s ---"
if [[ -s "$seed_root/CHGCAR" ]]; then echo "RESULT: PRESENT"; else echo "RESULT: ABSENT"; fi
echo "--- dfpt_has_vasp_artifacts probe ---"
for f in OUTCAR vasprun.xml OSZICAR vasp.log CONTCAR DOSCAR EIGENVAL IBZKPT WAVECAR CHGCAR CHG PCDAT REPORT XDATCAR; do
    [[ -e "$seed_root/$f" ]] && echo "  artifact: $f"
done
echo "--- check-recipe rc ---"
"$PYTHON" "$PREP" check-recipe --stage charge-seed --dfpt-variant production \
    --manifest "$seed_root/input-manifest.json" > /tmp/seeddiag-check.json 2>&1
echo "rc=$?"
grep -E '"(status|changed_tags|ignored_inert_tags)"' /tmp/seeddiag-check.json | sed 's/^/  /'
echo "--- verify rc ---"
"$PYTHON" "$PREP" verify --vasprun "$seed_root/vasprun.xml" --outcar "$seed_root/OUTCAR" \
    --stage charge-seed --report /tmp/seeddiag-verify.json > /dev/null 2>&1
echo "rc=$?"
"$PYTHON" -c "
import json;d=json.load(open('/tmp/seeddiag-verify.json'));print('  status=%s failures=%s'%(d['status'],d['failures']))" 2>&1
echo "--- done ---"
