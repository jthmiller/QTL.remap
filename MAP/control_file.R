#!/bin/R
## Each chrom/pop info
pop <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','ELR','NEW','BP')]
X <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')) ## X is equal to chrom number
slurmcore <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
popq <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','ELR','NEW','BP')]

## QTL Scans 
chrms <- c(1:24)
pops <- c('NBH','NEW')

## Only use previously mapped markers?
mapped.only=TRUE

if (mapped.only==TRUE) {
  outname <- 'NW_dropped'
} else { outname <- 'NW' }

## Directories
basedir <- '/home/jmiller1/QTL_Map_Raw/popgen'
plotdir <- file.path(basedir,'rQTL/plots')
indpops <- file.path(basedir,'plinkfiles/ind.pops')
popdir <- file.path(basedir,'rQTL',pop,'REMAPS')
qtldir <- file.path(basedir,'rQTL/remap_out')
errfile <- file.path(qtldir,'genotyping_error_rate.txt')

## Funtions for processing rQTL map data
source(file.path(basedir,'rQTL/scripts/QTL_remap/MAP/source_file.R'))

## Libraries
packs <- c('qtl','foreach','doParallel')
lapply(packs, require, character.only = TRUE)
## Load a couple fixed rQTL functions
#require(qtl2,lib.loc='/share/apps/rmodules')


dis.nbh <- c(2,13,20)
dis.elr <- c(18)
cov.nbh <- c(13,18)

## Parameters for rQTL for population specific datasets (NBH markers require at least 70% genotypes )
if (pop=='NBH'){
  inds <- c('ind15','ind89','ind88','ind14','ind20') # determined to be dropped low cov
  missing <- 0.9
  if (X %in% cov.nbh){
    missing <- 0.8
    finRf <- 0.20
  }
  grpLod <- 12 ## Standard LG form LOD
  finLod <- 14 ## Higher final NBH LOD
  grpRf <- 0.20
  finRf <- 0.10
  cutoff <- 1.0e-10
  miss <- 5
  if (X %in% dis.nbh){cutoff <- 1.0e-10}

} else if (pop=='ELR'){
  ### ELR parents are incorrectly ID'd
  inds <- c('ind2') # determined to be dropped low cov
  missing <- 0.8
  grpLod <- 8 ## Standard LG form LOD
  finLod <- 10 ## Higher final ELR LOD
  grpRf <- 0.2
  finRf <- 0.1
  cutoff <- 1.0e-4
  miss <- 10
  if (X %in% dis.elr){cutoff <- 1.0e-10}
} else if ( pop=='NEW'){
  inds <- c(NA) # determined to be dropped low cov
  missing <- 0.8
  grpLod <- 10 ## Standard LG form LOD
  finLod <- 12 ## Higher final ELR LOD
  grpRf <- 0.25
  finRf <- 0.15
  cutoff <- 1.0e-08
  miss <- 10
}

## Try to get error exported by map
expr <- paste('tac ',errfile,' | grep -m 1 \'',pop,' ',X,'\' | awk \'{print $5}\'',sep='')
try(ers <- as.numeric(system(expr,intern=T)))

if (length(ers)==0|is.null(ers)){ print('couldnt find the error. Using 0.03')
  ers <- 0.03
} else {print(paste(ers,'genotyping error'))}
