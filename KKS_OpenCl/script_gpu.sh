#!/bin/sh
#SBATCH --account=abhik
#SBATCH --job-name=ffffa	      # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --output=%j.log   # Standard output and error log
#SBATCH --gres=gpu:1
#SBATCH --partition=debug

make 
#time ./microsim_kks_opencl i.in   Filling_isoT_rad5pts.in  KKS_OpenCl_binary_1_dT13_F2

time ./microsim_kks_opencl Input_KKS_github_tdb_und_13K.in Filling_isoT.in MPMC_F4_atr > output.log

