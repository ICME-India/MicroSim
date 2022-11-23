#!/bin/sh

#SBATCH -N 2

#SBATCH --ntasks-per-node=32 

#SBATCH --time=23:00:00 

#SBATCH --gres=gpu:2

#SBATCH --job-name=GP4ku

#SBATCH --error=job.%J.err_gpu 

#SBATCH --output=job.%J.out_gpu 

#SBATCH --partition=gpumultinode 

source ../../../pravega-env-setup.sh

mpirun -np 64 ./main2d.gnu.MPI.CUDA.ex input2.in
