#!/bin/sh

#SBATCH -N 4

#SBATCH --ntasks-per-node=16

#SBATCH --time=23:00:00 

#SBATCH --job-name=GP4ku

#SBATCH --error=job.%J.err_cpu 

#SBATCH --output=job.%J.out_cpu 

#SBATCH --partition=small 

source ../../../pravega-env-setup.sh

export LD_LIBRARY_PATH='/home/ext/apps/spack/opt/spack/linux-centos7-cascadelake/gcc-11.2.0/gsl-2.7-h5q52xdndmnqwjwnysbdt5c6ccbk4fv6/lib'

mpirun -np 64 ./main2d.gnu.MPI.ex input2.in

