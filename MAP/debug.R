## Figure out if all AAxAB are lumped or out of order and plot errorlod
## ers in controlfil still giving ers. get rid of it
####### DEBUG ONLY ####
pop <- 'NEW'
X <- 1
slurmcore <- 12
setwd('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/')
####### DEBUG ONLY ####

outname <- 'NW_dropped'
## Only use previously mapped markers?
mapped.only=TRUE

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

dis.nbh <- c(2,13,20)
dis.elr <- c(18)
cov.nbh <- c(13,18)

## Parameters for rQTL for population specific datasets (NBH markers require at least 70% genotypes )
if (pop=='NBH'){
  inds <- c('ind15','ind89','ind88','ind14','ind20') # determined to be dropped low cov
  missing <- 0.9
  grpLod <- 12 ## Standard LG form LOD
  finLod <- 14 ## Higher final NBH LOD
  grpRf <- 0.20
  finRf <- 0.10
  cutoff <- 1.0e-08
  if (X %in% dis.nbh){cutoff <- 1.0e-08}
  if (X %in% cov.nbh){
    missing <- 0.8
    grpLod <- 8 ## Standard LG form LOD
    finLod <- 10 ## Higher final NBH LOD
  }
}
if (pop=='ELR'){
  inds <- c('ind2') # determined to be dropped low cov
  missing <- 0.75
  grpLod <- 8 ## Standard LG form LOD
  finLod <- 10 ## Higher final ELR LOD
  grpRf <- 0.20
  finRf <- 0.10
  cutoff <- 1.0e-08
  if (X %in% dis.elr){
    cutoff <- 1.0e-08
    grpLod <- 6 ## Standard LG form LOD
    finLod <- 8 ## Higher final ELR LOD
  }
}
## Try to get error exported by map
expr <- paste('tac ',errfile,' | grep -m 1 \'',pop,' ',X,'\' | awk \'{print $5}\'',sep='')
try(ers <- as.numeric(system(expr,intern=T)))

if (length(ers)==0|is.null(ers)){ print('couldnt find the error. Using 0.03')
  ers <- 0.03
} else {print(paste(ers,'genotyping error'))}

#Used to update phenotypes
##ind.inx <- grep('NG',cross.18$pheno$ID)
##repl <- phenos[as.character(cross.18$pheno$ID[grep('NG',cross.18$pheno$ID)]),1]
#cross.18$pheno$pheno_05[ind.inx] <- as.character(repl)
