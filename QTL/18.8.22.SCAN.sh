#!/bin/bash -l
#SBATCH -o /home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/array_error_out/scan/%j_out_%A_%a.txt
#SBATCH -e /home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/array_error_out/scan/%j_err_%A_%a.txt
#SBATCH --time=48:00:00
#SBATCH --job-name=QTL.Remap
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=5000
#SBATCH -p low
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jthmiller@ucdavis.edu
#SBATCH --array=1-24%5

####QTLs are on chrm '1 2 8 13 18 24'

scriptdir='/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/QTL'

pop=$1

Rscript $scriptdir/18.8.22.qtl_SCAN.R --vanilla $pop $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
