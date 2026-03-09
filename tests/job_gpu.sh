#!/bin/bash
#SBATCH --job-name=gml-test
#SBATCH --output=job_gpu.log
##SBATCH --error=job.err
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-task=1
#SBATCH --partition=partGPU
#SBATCH --time=100:00:00

# Echo job info
echo "Working directory: $(pwd)"
# Set env
ulimit -s unlimited
module purge && module load vasp-gpu

WORKDIR="$(pwd)/vasprun/gpu/"
mkdir -p $WORKDIR
cp $(pwd)/POSCAR-32 $WORKDIR/POSCAR
cp $(pwd)/INCAR $WORKDIR/INCAR
cd $WORKDIR

echo -e "102\n2\n0.04\n" | vaspkit
echo -e "108\n" | vaspkit

srun stdbuf -oL vasp_std

# Echo job done
echo "Job completed at: $(date)"
