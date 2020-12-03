#!/bin/sh
#SBATCH --nodes=16
#SBATCH --time=1-0
#SBATCH --ntasks-per-node=8
#SBATCH --partition=med
mpirun -n 128 uc
