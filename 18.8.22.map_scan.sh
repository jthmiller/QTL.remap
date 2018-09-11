#!/bin/bash

sbatch 18.8.22.remap_qtl.sh

wait

sbatch 18.8.22.scan_plot_rQTL2.sh
