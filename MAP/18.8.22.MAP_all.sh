#!/bin/bash
scriptdir='/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/'

# first job - no dependencies
echo "Population Chromosome Genotyping_Error_Rate" > /home/jmiller1/QTL_Map_Raw/popgen/rQTL/remap_out/genotyping_error_rate.txt

NBH1=$(sbatch -J "NBH.initial" --parsable $scriptdir/MAP/18.8.22.MAP_qtl_1.sh NBH)
NBH2=$(sbatch -J "NBH.pardrop"  --parsable --dependency=afterany:$NBH1 $scriptdir/MAP/18.8.22.MAP_qtl_2.sh NBH)
NBH3=$(sbatch -J "NBH.final" --parsable --dependency=afterany:$NBH2 $scriptdir/MAP/18.8.22.MAP_qtl_3.sh NBH)
NBH4=$(sbatch -J "NBH.scan" --parsable --dependency=afterany:$NBH3 $scriptdir/QTL/18.8.22.SCAN.sh NBH)
NBH5=$(sbatch -J "NBH.scan2" --parsable --dependency=afterany:$NBH4 $scriptdir/QTL/18.8.22.SCAN2.sh NBH)

#E. river
### RES4=$(sbatch -J "ELR.initial" --parsable $scriptdir/MAP/18.8.22.MAP_qtl_1.sh ELR)
ELR1=$(sbatch -J "ELR.initial" --parsable --dependency=afterany:$NBH5 $scriptdir/MAP/18.8.22.MAP_qtl_1.sh ELR)
ELR2=$(sbatch -J "ELR.pardrop"  --parsable --dependency=afterany:$ELR1 $scriptdir/MAP/18.8.22.MAP_qtl_2.sh ELR)
ELR3=$(sbatch -J "ELR.final"  --parsable --dependency=afterany:$ELR2 $scriptdir/MAP/18.8.22.MAP_qtl_3.sh ELR)
ELR4=$(sbatch -J "ELR.scan" --parsable --dependency=afterany:$ELR3 $scriptdir/QTL/18.8.22.SCAN.sh ELR)
ELR5=$(sbatch -J "NBH.scan" --parsable --dependency=afterany:$ELR4 $scriptdir/QTL/18.8.22.SCAN2.sh NBH)

#Newark
### RES7=$(sbatch -J "NEW.initial" --parsable $scriptdir/18.8.22.MAP_qtl_1.sh NEW)
NEW1=$(sbatch -J "NEW.initial" --parsable --dependency=afterany:$ELR5 $scriptdir/MAP/18.8.22.MAP_qtl_1.sh NEW)
NEW2=$(sbatch -J "NEW.pardrop"  --parsable --dependency=afterany:$NEW1 $scriptdir/MAP/18.8.22.MAP_qtl_2.sh NEW)
NEW3=$(sbatch -J "NEW.final" --parsable --dependency=afterany:$NEW2 $scriptdir/MAP/18.8.22.MAP_qtl_3.sh NEW)
NEW4=$(sbatch -J "NEW.scan" --parsable --dependency=afterany:$NEW3 $scriptdir/QTL/18.8.22.SCAN.sh NEW)
NEW5=$(sbatch -J "NBH.scan" --parsable --dependency=afterany:$NEW4 $scriptdir/QTL/18.8.22.SCAN2.sh NBH)
