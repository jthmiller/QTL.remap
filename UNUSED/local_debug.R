## Local Debug for Full QTL model

## rQTL2
#slurmcore <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
slurmcore <- 1
#popq <- commandArgs(TRUE)[commandArgs(TRUE) %in% c('NBH','ELR','NEW','BP')]
popq <- 'NEW'
pop <- 'NEW'

## Directories
basedir <- '~/Dropbox/QTL_Paper/CODE/QTL'
#plotdir <- file.path(basedir,'rQTL/plots')
#indpops <- file.path(basedir,'plinkfiles/ind.pops')
#qtldir <- file.path(basedir,'rQTL/remap_out')
#errfile <- file.path(qtldir,'genotyping_error_rate.txt')

## Funtions for processing rQTL data
source(file.path(basedir,'source_file.R'))
source(file.path(basedir,'model_source_file.R'))

## Libraries
packs <- c('qtl','foreach','doParallel')
lapply(packs, require, character.only = TRUE)
require(qtl2,lib.loc='/share/apps/rmodules')

## QTL LGs to consider
#X <- c(1,2,8,13,18,24) ## When only mapping qTLs
chrms <- c(1:24)
pops <- c('NBH','ELR')
