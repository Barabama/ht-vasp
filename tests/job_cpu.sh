#!/bin/bash
#SBATCH --job-name=gml-test
#SBATCH --output=job_cpu.log
##SBATCH --error=job.err
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --partition=partCPU
#SBATCH --time=100:00:00

# Echo job info
echo "Working directory: $(pwd)"
# Set env
ulimit -s unlimited
module purge && module load vasp-cpu

WORKDIR="$(pwd)/vasprun/cpu/"
mkdir -p $WORKDIR
cp $(pwd)/POSCAR-32 $WORKDIR/POSCAR
cp $(pwd)/INCAR $WORKDIR/INCAR
cd $WORKDIR

echo -e "102\n2\n0.04\n" | vaspkit
echo -e "108\n" | vaspkit

srun vasp_std

# Echo job done
echo "Job completed at: $(date)"
