#!/bin/bash
# parallel job using 32 processors. and runs for 4 hours (max)
#SBATCH -N 4 # node count
#SBATCH --ntasks-per-node=8
#SBATCH -t 4:00:00
# sends mail when process begins, and
# when it ends. Make sure you define your email
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=jhlu@princeton.edu

# Load openmpi environment
module load python/2.7
module load python/2.7/scipy-mkl
module load python/2.7/numpy-mkl
module load anaconda

srun ./run_causal.sh