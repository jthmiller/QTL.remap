#!/bin/bash

plink=~/bin/plink
vcfdir=/home/jmiller1/QTL_Map_Raw/popgen/vcf
outfiles=/home/jmiller1/QTL_Map_Raw/popgen/outfiles
infiles=/home/jmiller1/QTL_Map_Raw/popgen/infiles
indpops=/home/jmiller1/QTL_Map_Raw/popgen/plinkfiles/ind.pops
pheno=$infiles/SOMM.FAM.2.txt
init_flagset='--make-bed  --allow-extra-chr --autosome-num 24 --allow-no-sex --family'
flagset='--set-missing-var-ids @:#  --allow-extra-chr --autosome-num 24 --allow-no-sex --family --chr 1-24'
geno='--geno .4'
maf='--mac 1'


module load vcftools

for X in NBH BRP NEW
do

#vcftools --gzvcf $vcfdir/SOMM.vcf.gz --keep $infiles/$X.samples --remove-filtered-all --recode --stdout | vcftools --vcf - --max-meanDP 90 --maxDP 90 --stdout --minGQ 20 $maf $maxf --remove-filtered-all --recode | gzip -c > $vcfdir/$X.vcf.gz

$plink --vcf $vcfdir/$X.vcf.gz --out $indpops/$X  $init_flagset --pheno $pheno --all-pheno --update-ids $infiles/SOMM.txt

$plink --bfile $indpops/$X --out $indpops/$X $flagset --make-bed --pheno $pheno --all-pheno --keep-cluster-names $X --make-founders

$plink --bfile $indpops/$X --out $indpops/$X $flagset --pheno $pheno --all-pheno --keep-cluster-names $X $geno $maf --recode --biallelic-only strict --snps-only just-acgt --nonfounders

done

### REFILTER ELR. It contains the parent 'BLI' code (in addition to ELR fam)

$plink --vcf $vcfdir/$X.vcf.gz --out $indpops/$X  $init_flagset --pheno $pheno --all-pheno --update-ids $infiles/SOMM.txt

$plink --bfile $indpops/$X --out $indpops/$X $flagset --make-bed --pheno $pheno --all-pheno --keep-cluster-names ELR BLI --make-founders

$plink --bfile $indpops/$X --out $indpops/$X $flagset --pheno $pheno --all-pheno --keep-cluster-names ELR BLI $geno $maf --recode --biallelic-only strict --snps-only just-acgt --nonfounders
