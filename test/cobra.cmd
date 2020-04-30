#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./ompgpu.%j.script.out
#SBATCH -e ./ompgpu.%j.script.err
# Initial working directory:
#SBATCH -D ./
#
#SBATCH -J gpu
#
#SBATCH  --partition=gpudev

# Node feature:
#SBATCH --constraint="gpu"
# Specify type and number of GPUs to use:
#   GPU type can be v100 or rtx5000
# #SBATCH --gres=gpu:v100:2         # If using both GPUs of a node
#SBATCH --gres=gpu:v100:1       # If using only 1 GPU of a shared node
#SBATCH --mem=22500             # Memory is necessary if using only 1 GPU
#
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
# #SBATCH --ntasks-per-node=40      # If using both GPUs of a node
#SBATCH --ntasks-per-node=20    # If using only 1 GPU of a shared node
#
#SBATCH --mail-type=none
#SBATCH --mail-user=<userid>@rzg.mpg.de
#
# wall clock limit:
#SBATCH --time=00:05:00
#module load cuda 

# Run the program:
#srun ./ompgpu > ompgpu.$SLURM_JOB_ID.out 2> ompgpu.$SLURM_JOB_ID.err

# profile the program:

module  load  nsight_systems
nsys profile -t cuda,nvtx srun ./ompgpu

# to analyze the resulting report1.qdrep, run "nsight-sys" and open it from there.
