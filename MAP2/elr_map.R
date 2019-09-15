#!/bin/R
### Map QTLs 1 of 3
debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")

library('qtl')

mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

################################################################################
## read in the QTL cross
################################################################################
fl <- file.path(mpath,'ELR_unmapped_filtered_added_markers.csv')

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
 estimate.map = FALSE
)
################################################################################

################################################################################
dups <- findDupMarkers(cross, exact.only = T, adjacent.only = F)
cross.final.nodups <- drop.markers(cross, unlist(dups))
cross.final.nodups <- subset(cross.final.nodups,chr=i)
################################################################################
cross_nodups <- orderMarkers(cross.final.nodups, window=7,verbose=FALSE,chr=i,
                 use.ripple=TRUE, error.prob=0.01, sex.sp=FALSE,
                 map.function="kosambi",maxit=2000, tol=1e-4)

filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_mapped_chr_',i)
write.cross(cross_nodups,chr=i,filestem=filename,format="csv")
################################################################################
