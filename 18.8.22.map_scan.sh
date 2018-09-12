#!/bin/bash

echo "Population Chromosome Genotyping_Error_Rate" > /home/jmiller1/QTL_Map_Raw/popgen/rQTL/remap_out/genotyping_error_rate.txt

# first job - no dependencies
jid1=$(sbatch 18.8.22.remap_qtl_1.sh)

jid2=$(sbatch --dependency=afterok:$jid1 118.8.22.remap_qtl_2.sh)

jid3=$(sbatch --dependency=afterok:$jid2 118.8.22.remap_qtl_3.sh)


# multiple jobs can depend on a single job
jid2=$(sbatch --dependency=afterok:$jid1 18.8.22.scan_plot_rQTL2.sh)
