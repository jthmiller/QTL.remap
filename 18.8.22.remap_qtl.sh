#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o /home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/array_error_out/array_job_out_%A_%a.txt
#SBATCH -e /home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/array_error_out/array_job_err_%A_%a.txt
#SBATCH --time=24:00:00
#SBATCH --job-name=QTL.Remap
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=30000
#SBATCH -p med
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jthmiller@ucdavis.edu
#SBATCH --array=1,2,8,13,18,24%6

####QTLs are on chrm '1 2 8 13 18 24'

scriptdir='/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap'

pop='NBH'
#pop='ELR'

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=3000
#Rscript $scriptdir/18.8.22.remap_qtl.1_3.R --vanilla $pop $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
echo $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
wait

#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=30000
#Rscript $scriptdir/18.8.22.remap_all_qtl.2_3.R --vanilla $pop $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
echo $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
#wait

#Rscript $scriptdir/18.8.22.remap_all_qtl.2_3.R --vanilla $pop $SLURM_ARRAY_TASK_ID $SLURM_CPUS_PER_TASK
