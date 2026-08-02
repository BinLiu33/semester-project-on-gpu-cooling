#!/bin/bash -l
#SBATCH --job-name=sod_10_1
#SBATCH --account=uzh8              
#SBATCH --partition=debug
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4         
#SBATCH --cpus-per-task=72          
#SBATCH --gpus-per-task=1           
#SBATCH --uenv=prgenv-gnu
#SBATCH --view=modules
#SBATCH --output=sod_10_1.%j.out
#SBATCH --error=sod_10_1.%j.err

set -euo pipefail


module load boost cmake cray-mpich fftw gcc gsl hdf5 python cuda

cd /users/binliu/pkdgrav3/
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt

cd /users/binliu/sod_tube

PKDGRAV3=/users/binliu/pkdgrav3/build/pkdgrav3

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MPICH_GPU_SUPPORT_ENABLED=1          

srun --cpus-per-task=${SLURM_CPUS_PER_TASK} \
     "${PKDGRAV3}" ./initial/sod_10_1.par
