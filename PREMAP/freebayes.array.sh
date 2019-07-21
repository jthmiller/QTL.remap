#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o /home/jmiller1/QTL_Map_Raw/QTL_scripts/out_er/array_job_out_%A_%a.txt
#SBATCH -e /home/jmiller1/QTL_Map_Raw/QTL_scripts/out_er/array_job_err_%A_%a.txt
#SBATCH --array=1-8577
#SBATCH -p med
#SBATCH --mem=30000
#SBATCH -t 48:00:00
###### number of nodes
###### number of processors
#SBATCH --cpus-per-task=6

bwagenind=/home/nreid/popgen/kfish3/killifish20130322asm.fa

my_freebayes=/home/jmiller1/bin/freebayes/bin/freebayes
my_bedtools=/home/jmiller1/bin/bedtools2/bin/bedtools
my_bamtools=/home/jmiller1/bin/bamtools-master/bin/bamtools
#module load $my_freebayes

scafnum=$(expr $SLURM_ARRAY_TASK_ID + -1)
#scaf=$(sed -n "$SLURM_ARRAY_TASK_ID p" /home/jmiller1/QTL_Map_Raw/bedfiles/unique.scaf.lg.txt)
scaf=Scaffold$scafnum
endpos=$(expr $(grep -P "$scaf\t" /home/nreid/popgen/kfish3/killifish20130322asm.fa.fai | cut -f 2) - 1)

region=$scaf:1..$endpos
echo $region

outfile=$scaf.vcf

bam_dir=/home/jmiller1/QTL_Map_Raw/align
vcf_out=/home/jmiller1/QTL_Map_Raw/vcf/freebayes.array.kfish3
bed_regions=~/QTL_Map_Raw/align/radsites.bed
bam_list=/home/jmiller1/QTL_Map_Raw/align/bamlist.txt
#merged_bams=/home/jmiller1/QTL_Map_Raw/align/SOMM0_ALL.bam
#all_scaf=/home/jmiller1/bin/code/ALLMAPS_OUT/unsplit_merge.fasta.genomefile.txt

$my_bamtools merge -list $bam_list -region $region| \
$my_bamtools filter -in stdin -mapQuality '>30' -isProperPair true | \
$my_bedtools intersect -sorted -a stdin -b $bed_regions | \
$my_freebayes -f $bwagenind --use-best-n-alleles 4 --pooled-discrete --stdin > $vcf_out/$outfile
