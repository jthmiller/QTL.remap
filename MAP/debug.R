## Figure out if all AAxAB are lumped or out of order and plot errorlod
## ers in controlfil still giving ers. get rid of it
####source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/debug.R')

####### DEBUG ONLY ####
slurmcore <- 12
setwd('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/')
####### DEBUG ONLY ####
### Prompt
pop <- function()c('NBH','NEW','ELR','NEW')[menu(c('NBH','NEW','ELR','NEW'), title="Which pop")]
pop <- popq <- pop()
X <- function()c(1:24)[menu(1:24, title="Which Chromosome")]
X <- X()
## Only use previously mapped markers?
mapped.only <- function()c(TRUE,FALSE)[menu(c(TRUE,FALSE), title="Mapped markers only?")]
mapped.only <- mapped.only()
## Only use granparent confirmed markers?
confirmed <- function()c(TRUE,FALSE)[menu(c(TRUE,FALSE), title="Only use granparent confirmed markers?")]
confirmed <- confirmed()

if (mapped.only==TRUE) {
  outname <- 'NW_dropped'
  print('Excluding markers on scaffolds')
} else { outname <- 'NW'
  print('Including markers on scaffolds')
}

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
flib <- '/share/apps/rmodules'
fpacks <- c('devtools','httr','ggplot2','reshape','pheatmap','RColorBrewer')
lapply(fpacks, require, character.only = TRUE,lib.loc=flib)

mylib <- "/home/jmiller1/R/x86_64-pc-linux-gnu-library/3.5"
mpacks <- c('qtl','foreach','doParallel','qtl2','qtlTools','gplots','qgraph')
lapply(mpacks, require, character.only = TRUE,lib.loc=mylib)

### Phenotype translation
trsl.bin <- c(0,0,0,1,1,1)
names(trsl.bin) <- as.character(0:5)

chrms <- 1:24

## Parameters for rQTL for population specific datasets (NBH markers require at least 70% genotypes )
## Parameters for rQTL for population specific datasets (NBH markers require at least 70% genotypes )
if (pop=='NBH'){
  confirmed=T
  reorder<-F
  mapped.only=TRUE
  grpLod <- 12 ## Standard LG form LOD
  finLod <- 14 ## Higher final NBH LOD
  grpRf <- 0.20
  finRf <- 0.1
  cutoff <- 1.0e-8
  miss <- 8
} else if (pop=='ELR'){
  mapped.only <- TRUE
  confirmed <- FALSE
  reorder<-F
  missing <- 0.9
  grpLod <- 10 ## Standard LG form LOD
  finLod <- 12 ## Higher final ELR LOD
  grpRf <- 0.2
  finRf <- 0.1
  cutoff <- 1.0e-4
  miss <- 2 ## Higher, need more power to detect seg distortion
} else if ( pop=='NEW'){
  confirmed=T
  reorder<-F
  mapped.only=TRUE
  inds <- c(NA) # determined to be dropped low cov
  missing <- 0.8
  grpLod <- 12 ## Standard LG form LOD
  finLod <- 14 ## Higher final ELR LOD
  grpRf <- 0.20
  finRf <- 0.10
  cutoff <- 1.0e-8
  miss <- 8
}

if (mapped.only==TRUE & reorder==F) {
  outname <- 'NW_dropped_physical'
} else if (mapped.only==TRUE & reorder==T ){
   outname <- 'NW_Mapped'
} else {
   outname <- 'NW_ReMapped'
}

## Try to get error exported by map

expr <- paste('tac ',errfile,' | grep -m 1 \'',pop,' ',X,'\' | awk \'{print $5}\'',sep='')
try(ers <- as.numeric(system(expr,intern=T)))

if (length(ers)==0|is.null(ers)){ print('couldnt find the error. Using 0.03')
  ers <- 0.03
} else {print(paste(ers,'genotyping error'))}
