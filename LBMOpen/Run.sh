#!/bin/bash
#SBATCH --job-name=AllWang
#SBATCH -t 3-23:00:00
#SBATCH -o slurm.%j.out
#SBATCH --nodes 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH -e slurm.%j.err
#SBATCH --mail-type=END  
#SBATCH --mail-user=thomascarmo@outlook.com

module load gnu13/13.2.0

# Configurar o número de threads OpenMP para 20
export OMP_NUM_THREADS=20


make clean && make 

./LBMOpen
