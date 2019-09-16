#!/bin/R
### Map QTLs 1 of 3
library('qtl')

################################################################################
## read in the QTL cross
################################################################################
i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]

print(i)

mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

fl <- file.path(mpath,'ELR_subsetted.csv')

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
 estimate.map = FALSE
)
################################################################################

nmars <- nmar(cross)

cross <- subset(cross,ind=nmissing(cross) < (nmars*.5))

################################################################################
dups <- findDupMarkers(cross, exact.only = T, adjacent.only = F)
cross <- drop.markers(cross, unlist(dups))
cross <- subset(cross,chr=i)
##cross <- est.map(cross, error.prob=0.1, map.function="kosambi",sex.sp=F,chr=i)

################################################################################
cross <- orderMarkers(cross, window=7,verbose=FALSE,chr=i,
                 use.ripple=TRUE, error.prob=0.01, sex.sp=FALSE,
                 map.function="kosambi",maxit=1, tol=1e-4)

filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_mapped_chr_',i)
write.cross(cross,chr=i,filestem=filename,format="csv")
################################################################################
