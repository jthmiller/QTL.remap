#!/bin/bash -l
#SBATCH -o /home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/array_error_out/ripple/out_%x_%A_%a.txt
#SBATCH -e /home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/array_error_out/ripple/err_%x_%A_%a.txt
#SBATCH -t 12:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=2000
#SBATCH -p high
#SBATCH --array=1-24

scriptdir='/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP'

pop=$1
### This requires 2 cpu with 3G memory each, but takes awhile
Rscript $scriptdir/18.8.22.MAP_qtl.3_3.R --vanilla $pop $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
