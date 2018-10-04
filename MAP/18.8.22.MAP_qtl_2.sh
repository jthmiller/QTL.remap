#!/bin/bash -l
#SBATCH -o /home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/array_error_out/pardrop/out_%x_%A_%a.txt
#SBATCH -e /home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/array_error_out/pardrop/err_%x_%A_%a.txt
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4000
#SBATCH -p med
#SBATCH --array=1-24
#SBATCH --quiet

####QTLs are on chrm '1 2 8 13 18 24'
scriptdir='/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP'

pop=$1
### This requires 12 cpus with 5G of memory (dasically an entire node, but cant run on 24 cpus with only 2.5Gs of memory)
Rscript $scriptdir/18.8.22.MAP_qtl.2_3.R --vanilla $pop $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
