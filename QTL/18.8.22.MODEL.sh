#!/bin/bash

scriptdir='/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/QTL'

RES1=$(sbatch -J "NBH.scan" --parsable $scriptdir/18.8.22.SCAN.sh NBH)
