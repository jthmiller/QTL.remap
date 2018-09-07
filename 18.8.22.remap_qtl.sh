#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o /home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/array_error_out/array_job_out_%A_%a.txt
#SBATCH -e /home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/array_error_out/array_job_err_%A_%a.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=QTL.Remap
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=4500
#SBATCH -p med
####SBATCH -n 2
####SBATCH --mail-type=ALL
####SBATCH --mail-user=jthmiller@ucdavis.edu
#SBATCH --array=1,2,8,13,18,24%6


####QTLs are on chrm '1 2 8 13 18 24'

scriptdir='/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap'
#declare -a pops=("NBH" "ELR")

#for pop in "${pops[@]}"
#do
pop='NBH'
#pop='ELR'
Rscript $scriptdir/18.8.22.remap_all_qtl.R --vanilla $pop $SLURM_ARRAY_TASK_ID
#echo 'working on QTL chromosome number' $SLURM_ARRAY_TASK_ID 'in' $pop
#Rscript $scriptdir/test.R $pop
#done
#18.8.22.remap_all_qtl.R

## set better mem (limit mem in Rscript)
## consider splitting parallel part out
## lighten up mem requirements
## fix marker.density
