#!/bin/bash
#SBATCH -o job.%j.out
#SBATCH -J wd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -A pi_hpc
#SBATCH -p gpu4Q
#SBATCH -q gpuq
#SBATCH --gres=gpu:1 

module load NVHPCSDK/23.11
module load GNU/gcc-9.3.0
export PATH=/public/software/nvhpc_2311_cuda_12.3/Linux_x86_64/23.11/comm_libs/hpcx/bin/:$PATH
srun hostname -s | sort -n >slurm.hosts

make
mpirun -n 1 -machinefile slurm.hosts ./exe
