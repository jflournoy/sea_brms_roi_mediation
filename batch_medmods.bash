#!/bin/bash
#SBATCH -J brms-med
#SBATCH --mem 5G
#SBATCH -p ncf
#SBATCH --cpus-per-task 4
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 0-12:00
#SBATCH -o %x_%A_%a.out
#SBATCH --mail-user=john_flournoy@fas.harvard.edu
#SBATCH --mail-type=ALL

module load gcc/7.1.0-fasrc01
module load R/3.5.1-fasrc01

cp ~/.R/Makevars.gcc ~/.R/Makevars

export R_LIBS_USER=/ncf/mclaughlin/users/jflournoy/R_3.5.1_GCC:$R_LIBS_USER

runme=/ncf/mclaughlin/users/jflournoy/code/

srun -c 4 `which Rscript` 