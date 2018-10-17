## rQTL2
slurmcore <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
popq <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','ELR','NEW','BP')]
if (popq>1){

}

## Directories
basedir <- '/home/jmiller1/QTL_Map_Raw/popgen'
plotdir <- file.path(basedir,'rQTL/plots')
indpops <- file.path(basedir,'plinkfiles/ind.pops')
qtldir <- file.path(basedir,'rQTL/remap_out')
errfile <- file.path(qtldir,'genotyping_error_rate.txt')
popdir <- file.path(basedir,'rQTL',popq,'REMAPS')
setwd(qtldir)

## Funtions for processing rQTL data
source(file.path(basedir,'rQTL/scripts/QTL_remap/QTL/source_file.R'))
source(file.path(basedir,'rQTL/scripts/QTL_remap/QTL/model_source_file.R'))

## Libraries
.libPaths(c( .libPaths(), "/home/jmiller1/R/x86_64-pc-linux-gnu-library/3.5"))
packs <- c('qtl','foreach','doParallel')
lapply(packs, require, character.only = TRUE)
require(qtl2,lib.loc='/share/apps/rmodules')

## QTL LGs to consider
#X <- c(1,2,8,13,18,24) ## When only mapping qTLs
chrms <- c(1:24)
pops <- c('NBH','ELR')

## Try to get chromosome avg error exported by map
ers <- mean(sapply(chrms,function(Z,pop='NBH'){
  expr <- paste('tac ',errfile,' | grep -m 1 \'',pop,' ',Z,'\' | awk \'{print $5}\'',sep='')
  try(ers <- as.numeric(system(expr,intern=T)))
    if (length(ers)==0){ ers <- 0.03 }
    return(round(ers,digits=4))
    }
  )
)
print(paste(ers,'genotyping error'))
