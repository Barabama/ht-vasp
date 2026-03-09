#!/bin/bash
#SBATCH --job-name=gml-qha
#SBATCH --output=job_qha.log
##SBATCH --error=job.err
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=44
#SBATCH --cpus-per-task=1
#SBATCH --partition=partCPU
#SBATCH --time=100:00:00

# Echo job info
echo "Working directory: $(pwd)"
# Set env
ulimit -s unlimited
module purge && module load vasp-cpu

source /opt/miniconda3/bin/activate
conda activate mpea-tdb
python main_qha.py

# Echo job done
echo "Job completed at: $(date)"
