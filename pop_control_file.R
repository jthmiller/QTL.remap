#!/bin/R
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','ELR','NEW','BP')]
X <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) ## X is equal to chrom number
slurmcore <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
## Directories
basedir <- '/home/jmiller1/QTL_Map_Raw/popgen'
plotdir <- file.path(basedir,'rQTL/plots')
indpops <- file.path(basedir,'plinkfiles/ind.pops')
popdir <- file.path(basedir,'rQTL',pop,'REMAPS')
qtldir <- file.path(basedir,'rQTL/remap_out')
errfile <- file.path(qtldir,'genotyping_error_rate.txt')
setwd(popdir)

## Funtions for processing rQTL map data
source(file.path(basedir,'rQTL/scripts/QTL_remap/removeDoubleXO.R'))
source(file.path(basedir,'rQTL/scripts/QTL_remap/QTL_map_sourcefile.R'))
source(file.path(basedir,'rQTL/scripts/QTL_remap/custom_rQTL_functions.R'))

## Libraries
packs <- c('qtl','foreach','doParallel')
lapply(packs, require, character.only = TRUE)
## Load a couple fixed rQTL functions

## Parameters for rQTL for population specific datasets (NBH markers require at least 70% genotypes )
if (pop=='NBH'){
  inds <- c('ind15','ind89','ind88','ind14','ind20') # determined to be dropped low cov
  missing <- 0.80
  grpLod <- 8 ## Standard LG form LOD
  finLod <- 12 ## Higher final NBH LOD
  grpRf <- 0.35
  finRf <- 0.25
}
if (pop=='ELR'){
  inds <- c('ind2') # determined to be dropped low cov
  missing <- 0.75
  grpLod <- 8 ## Standard LG form LOD
  finLod <- 12 ## Higher final NBH LOD
  grpRf <- 0.35
  finRf <- 0.25
}
## Try to get error exported by map
expr <- paste('tac ',errfile,' | grep -m 1 \'',pop,' ',X,'\' | awk \'{print $3}\'',sep='')
try(ers <- as.numeric(system(expr,intern=T)))

if (length(ers)==0){ print('couldnt find the error. Using 0.02')
  ers <- 0.02
} else {print(paste(ers,'genotyping error'))}
