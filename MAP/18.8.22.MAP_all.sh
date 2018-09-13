echo "Population Chromosome Genotyping_Error_Rate" > /home/jmiller1/QTL_Map_Raw/popgen/rQTL/remap_out/genotyping_error_rate.txt

scriptdir='/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP'

# first job - no dependencies

RES=$(sbatch --parsable $scriptdir/18.8.22.MAP_qtl_1.sh)

RES2=$(sbatch --parsable --dependency=afterok:$RES $scriptdir/18.8.22.MAP_qtl_2.sh)

sbatch --dependency=afterok:$RES2 $scriptdir/18.8.22.MAP_qtl_3.sh
