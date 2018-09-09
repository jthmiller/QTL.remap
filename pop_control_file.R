#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','ELR','NEW','BP')]
X <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) ## X is equal to chrom number
slurmcore <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
## Directories
basedir <- '/home/jmiller1/QTL_Map_Raw/popgen'
plotdir <- file.path(basedir,'rQTL/plots')
indpops <- file.path(basedir,'plinkfiles/ind.pops')
popdir <- file.path(basedir,'rQTL',pop,'REMAPS')
setwd(popdir)

## Funtions for processing rQTL data
source(file.path(basedir,'rQTL/scripts/QTL_remap/removeDoubleXO.R'))
source(file.path(basedir,'rQTL/scripts/QTL_remap/QTL_map_sourcefile.R'))

## Parameters for rQTL for all datasets

## Parameters for rQTL for population specific datasets (NBH markers require at least 70% genotypes )
if (pop=='NBH'){
  inds <- c('ind15','ind89','ind88','ind14') # determined to be dropped low cov
  missing <- 0.75
  grpLod <- 6 ## Standard LG form LOD
  finLod <- 8 ## Higher final NBH LOD

}
if (pop=='ELR'){
  inds <- c('ind2') # determined to be dropped low cov
  missing <- 0.65
  grpLod <- 4 ## Standard LG form LOD
  finLod <- 6 ## Higher final NBH LOD

}
