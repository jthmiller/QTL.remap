#!/bin/bash -l
#SBATCH -o /home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/array_error_out/kinship/out_%x_%A_%a.txt
#SBATCH -e /home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/array_error_out/kinship/err_%x_%A_%a.txt
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2000
#SBATCH -p high

scriptdir='/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP'

Rscript $scriptdir/compare.kinship.R --vanilla NBH $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
Rscript $scriptdir/compare.kinship.R --vanilla ELR $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
Rscript $scriptdir/compare.kinship.R --vanilla NEW $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
#Rscript $scriptdir/compare.kinship.R --vanilla $pop $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
