#!/bin/bash
#SBATCH --job-name=gml-test
#SBATCH --output=job_cpu.log
##SBATCH --error=job.err
#SBATCH --mem=20G
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --partition=partCPU
#SBATCH --time=100:00:00

ROOTDIR=$(pwd)

# Echo job info
echo "Working directory: $ROOTDIR"
# Set env
ulimit -s unlimited
module purge && module load vasp-cpu

WORKDIR=$ROOTDIR/vasprun/cpu/
mkdir -p $WORKDIR
cp $ROOTDIR/POSCAR-16 $WORKDIR/POSCAR
cp $ROOTDIR/INCAR $WORKDIR/INCAR
cd $WORKDIR

echo -e "102\n1\n0.04\n" | vaspkit
echo -e "108\n" | vaspkit

srun vasp_std

# Echo job done
echo "Job completed at: $(date)"
