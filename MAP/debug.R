outname <- 'NW_dropped'
## Only use previously mapped markers?
mapped.only=TRUE
## Each chrom/pop info

####### DEBUG ONLY ####
pop <- 'NBH'
X <- 13
slurmcore <- 12
####### DEBUG ONLY ####

basedir <- '/home/jmiller1/QTL_Map_Raw/popgen'
plotdir <- file.path(basedir,'rQTL/plots')
indpops <- file.path(basedir,'plinkfiles/ind.pops')
popdir <- file.path(basedir,'rQTL',pop,'REMAPS')
qtldir <- file.path(basedir,'rQTL/remap_out')
errfile <- file.path(qtldir,'genotyping_error_rate.txt')
setwd(popdir)

## Funtions for processing rQTL map data
source(file.path(basedir,'rQTL/scripts/QTL_remap/MAP/source_file.R'))

## Libraries
packs <- c('qtl','foreach','doParallel')
lapply(packs, require, character.only = TRUE)
## Load a couple fixed rQTL functions

## Parameters for rQTL for population specific datasets (NBH markers require at least 70% genotypes )
if (pop=='NBH'){
  inds <- c('ind15','ind89','ind88','ind14','ind20') # determined to be dropped low cov
  missing <- 0.8
  grpLod <- 10 ## Standard LG form LOD
  finLod <- 14 ## Higher final NBH LOD
  grpRf <- 0.30
  finRf <- 0.10
}
if (pop=='ELR'){
  inds <- c('ind2') # determined to be dropped low cov
  missing <- 0.8
  grpLod <- 8 ## Standard LG form LOD
  finLod <- 12 ## Higher final ELR LOD
  grpRf <- 0.30
  finRf <- 0.10
}
## Try to get error exported by map
expr <- paste('cat ',errfile,' | grep -m 1 \'',pop,' ',X,'\' | awk \'{print $5}\'',sep='')
try(ers <- as.numeric(system(expr,intern=T)))

if (length(ers)==0){ print('couldnt find the error. Using 0.02')
  ers <- 0.05
} else {print(paste(ers,'genotyping error'))}

### Begin of code
setwd('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/')
