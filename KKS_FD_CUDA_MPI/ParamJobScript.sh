#!/bin/bash
#SBATCH -N 4
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=2
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=1450_slurm_NiAlMo
#SBATCH --partition=gpu

module load gcc/10.2.0 cuda spack gnu8 gsl

. /home/apps/spack/share/spack/setup-env.sh
spack load nvhpc

mpirun -n $SLURM_NTASKS ./microsim_kks_fd_cuda_mpi Input.in Filling.in Output

