#!/bin/sh
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:05:00

#SBATCH --output=job.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4

mpirun -n 2 ./main2d.gnu.MPI.ex inputs

