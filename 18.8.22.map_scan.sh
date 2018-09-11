#!/bin/bash

# first job - no dependencies
jid1=$(sbatch 18.8.22.remap_qtl.sh)

# multiple jobs can depend on a single job
jid2=$(sbatch --dependency=afterok:$jid1 18.8.22.scan_plot_rQTL2.sh)
