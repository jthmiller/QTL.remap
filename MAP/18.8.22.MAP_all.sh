#!/bin/bash

scriptdir='/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP'

# first job - no dependencies

echo "Population Chromosome Genotyping_Error_Rate" > /home/jmiller1/QTL_Map_Raw/popgen/rQTL/remap_out/genotyping_error_rate.txt

RES=$(sbatch --parsable $scriptdir/18.8.22.MAP_qtl_1.sh NBH)
RES2=$(sbatch --parsable --dependency=afterany:$RES $scriptdir/18.8.22.MAP_qtl_1.sh ELR)
RES3=$(sbatch --parsable --dependency=afterany:$RES2 $scriptdir/18.8.22.MAP_qtl_2.sh NBH)
RES4=$(sbatch --parsable --dependency=afterany:$RES3 $scriptdir/18.8.22.MAP_qtl_2.sh ELR)
RES5=$(sbatch --parsable --dependency=afterany:$RES4 $scriptdir/18.8.22.MAP_qtl_3.sh NBH)
sbatch --dependency=afterok:$RES3 $scriptdir/18.8.22.MAP_qtl_3.sh ELR
