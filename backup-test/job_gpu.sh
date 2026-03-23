#!/bin/bash
#SBATCH --job-name=gml-test
#SBATCH --output=job_gpu.log
##SBATCH --error=job.err
#SBATCH --mem=10G
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-task=1
#SBATCH --partition=partGPU
#SBATCH --time=100:00:00

ROOTDIR=$(pwd)

# Echo job info
echo "Working directory: $ROOTDIR"
# Set env
ulimit -s unlimited
module purge && module load vasp-gpu

WORKDIR=$ROOTDIR/vasprun/gpu/
mkdir -p $WORKDIR
cp $ROOTDIR/POSCAR-16 $WORKDIR/POSCAR
cp $ROOTDIR/INCAR $WORKDIR/INCAR
cd $WORKDIR

echo -e "102\n1\n0.04\n" | vaspkit
echo -e "108\n" | vaspkit

srun stdbuf -oL vasp_std

# Echo job done
echo "Job completed at: $(date)"
