#!/bin/bash

scriptdir='/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/'

# first job - no dependencies

echo "Population Chromosome Genotyping_Error_Rate" > /home/jmiller1/QTL_Map_Raw/popgen/rQTL/remap_out/genotyping_error_rate.txt

RES1=$(sbatch -J "NBH.initial" --parsable $scriptdir/MAP/18.8.22.MAP_qtl_1.sh NBH)
RES2=$(sbatch -J "NBH.pardrop"  --parsable --dependency=afterany:$RES1 $scriptdir/MAP/18.8.22.MAP_qtl_2.sh NBH)
RES3=$(sbatch -J "NBH.final" --parsable --dependency=afterany:$RES2 $scriptdir/MAP/18.8.22.MAP_qtl_3.sh NBH)

#E. river
### RES4=$(sbatch -J "ELR.initial" --parsable $scriptdir/MAP/18.8.22.MAP_qtl_1.sh ELR)
RES4=$(sbatch -J "ELR.initial" --parsable --dependency=afterany:$RES3 $scriptdir/MAP/18.8.22.MAP_qtl_1.sh ELR)
RES5=$(sbatch -J "ELR.pardrop"  --parsable --dependency=afterany:$RES4 $scriptdir/MAP/18.8.22.MAP_qtl_2.sh ELR)
RES6=$(sbatch -J "ELR.final"  --parsable --dependency=afterok:$RES5 $scriptdir/MAP/18.8.22.MAP_qtl_3.sh ELR)

#Newark
### RES7=$(sbatch -J "NEW.initial" --parsable $scriptdir/18.8.22.MAP_qtl_1.sh NEW)
RES7=$(sbatch -J "NEW.initial" --parsable --dependency=afterany:$RES6 $scriptdir/MAP/18.8.22.MAP_qtl_1.sh NEW)
RES8=$(sbatch -J "NEW.pardrop"  --parsable --dependency=afterany:$RES7 $scriptdir/MAP/18.8.22.MAP_qtl_2.sh NEW)
RES9=$(sbatch -J "NEW.final" --parsable --dependency=afterany:$RES8 $scriptdir/MAP/18.8.22.MAP_qtl_3.sh NEW)

## Scan Models
RES10=$(sbatch -J "NBH.scan" --parsable --dependency=afterany:$RES9 $scriptdir/QTL/18.8.22.SCAN.sh NBH)
RES11=$(sbatch -J "ELR.scan" --parsable --dependency=afterany:$RES10 $scriptdir/QTL/18.8.22.SCAN.sh ELR)
RES12=$(sbatch -J "NEW.scan" --parsable --dependency=afterany:$RES11 $scriptdir/QTL/18.8.22.SCAN.sh NEW)
