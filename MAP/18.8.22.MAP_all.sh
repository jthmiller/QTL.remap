#!/bin/bash

scriptdir='/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP'

# first job - no dependencies

echo "Population Chromosome Genotyping_Error_Rate" > /home/jmiller1/QTL_Map_Raw/popgen/rQTL/remap_out/genotyping_error_rate.txt

RES=$(sbatch -J "NBH.initial" --parsable $scriptdir/18.8.22.MAP_qtl_1.sh NBH)
RES2=$(sbatch -J "ELR.initial" -o $elrdir --parsable --dependency=afterany:$RES $scriptdir/18.8.22.MAP_qtl_1.sh ELR)
RES3=$(sbatch -J "NBH.pardrop" -o $nbhdir --parsable --dependency=afterany:$RES2 $scriptdir/18.8.22.MAP_qtl_2.sh NBH)
RES4=$(sbatch -J "ELR.pardrop" -o $elrdir  --parsable --dependency=afterany:$RES3 $scriptdir/18.8.22.MAP_qtl_2.sh ELR)
RES5=$(sbatch -J "NBH.final" -o $nbhdir --parsable --dependency=afterany:$RES4 $scriptdir/18.8.22.MAP_qtl_3.sh NBH)
sbatch -J "ELR.final" -o $elrdir --dependency=afterok:$RES5 $scriptdir/18.8.22.MAP_qtl_3.sh ELR
