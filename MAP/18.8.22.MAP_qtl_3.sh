#!/bin/bash -l
#SBATCH -J "final_map"
#SBATCH -o /home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/array_error_out/final_map_out_%A_%a.txt
#SBATCH -e /home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/array_error_out/final_map_err_%A_%a.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=QTL.Remap
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5000
#SBATCH -p med
######SBATCH --mail-type=ALL
######SBATCH --mail-user=jthmiller@ucdavis.edu
#SBATCH --array=1-24

####QTLs are on chrm '1 2 8 13 18 24'

step=final_map
scriptdir='/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP'

pop='NBH'
### This requires 2 cpu with 3G memory each, but takes awhile
Rscript $scriptdir/18.8.22.remap_qtl.3_3.R --vanilla $pop $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK

wait

pop='ELR'

### This requires 2 cpu with 3G memory each, but takes awhile
Rscript $scriptdir/18.8.22.remap_qtl.3_3.R --vanilla $pop $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
