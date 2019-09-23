#!/bin/R
### Map QTLs 1 of 3
#debug.cross <- T
#source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library('qtl')
pop <- 'ELR'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

file_list <- gsub('.csv','',list.files(mpath, '*downsmpl_map*'))
file_list <- list.files(mpath, '*downsmpl_map*',include.dirs = T)

i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]

filename <- paste0('ELR_gts_CHR',i,'_downsmpl_map.csv')

cross <- read.cross(file=filename,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

################################################################################

refine_maps <- function(X){

 tmp <- orderMarkers(X, window=4,verbose=FALSE,
                 use.ripple=FALSE, error.prob=0.01, sex.sp=FALSE,
                 map.function="kosambi",maxit=1000, tol=1e-3)

 tmp <- calc.errorlod(tmp, err=0.01)

 tmp_map <-  est.map(tmp, error.prob=0.01,
              map.function="kosambi",
              maxit=1000, tol=1e-4, sex.sp=FALSE,
              verbose=FALSE, n.cluster=6)

 tmp <- qtl:::replace.map(tmp,tmp_map)

 drp1 <- droponemarker(tmp, error.prob=0.01,
                    map.function="kosambi",
                    maxit=1000, tol=1e-3, sex.sp=FALSE,
                    verbose=FALSE)

 dropByDropone(tmp, drp1, endMarkerThresh = 20,
  midMarkerThresh = 20, map.function = "kosambi",
  re.est.map = T, sex.sp=F, error.prob=0.01,maxit=1000, tol=1e-4, sex.sp=FALSE,
  verbose=FALSE, n.cluster=6)

}

cross <- refine_maps(cross)

filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_gts_CHR',i,'_downsmpl_map')
write.cross(cross,chr=i,filestem=filename,format="csv")













################################################################################
### fix 13
cross <- orderMarkers(cross, window=4,verbose=FALSE,chr=13,
                 use.ripple=FALSE, error.prob=0.01, sex.sp=FALSE,
                 map.function="kosambi",maxit=1000, tol=1e-3)


cross <- calc.errorlod(cross, err=0.01)

cross_map <-  est.map(cross, chr=13, error.prob=0.01,
              map.function="kosambi",
              maxit=100, tol=1e-4, sex.sp=FALSE,
              verbose=FALSE, omit.noninformative=TRUE, n.cluster=6)

cross <- qtl:::replace.map(cross,cross_map)

crs <- droponemarker(cross, chr, error.prob=0.01,
                    map.function="kosambi",
                    maxit=1000, tol=1e-3, sex.sp=FALSE,
                    verbose=FALSE)

################################################################################

################################################################################
## SCAN
################################################################################

##cross <- switchAlleles(cross, markers = markernames(cross,chr=13))
