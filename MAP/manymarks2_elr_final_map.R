#!/bin/R
### Map QTLs 1 of 3
debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")

library('qtl')
################################################################################
## read in the QTL cross
################################################################################
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'
i <- 18
cross.r <- read.cross(
 file = file.path(mpath, paste0('ELR_chr_',i,'.csv')),
 format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
 estimate.map = FALSE
)
################################################################################

################################################################################
dups <- findDupMarkers(cross.r, exact.only = T, adjacent.only = F)
cross.final.nodups <- drop.markers(cross.r, unlist(dups))
################################################################################

cross_nodups <- orderMarkers(cross.final.nodups, window=7,verbose=FALSE,
                 use.ripple=TRUE, error.prob=0.01, sex.sp=FALSE,
                 map.function="kosambi",maxit=4000, tol=1e-4)


filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_mapped_chr_',i)
write.cross(cross_nodups,chr=i,filestem=filename,format="csv")
################################################################################
