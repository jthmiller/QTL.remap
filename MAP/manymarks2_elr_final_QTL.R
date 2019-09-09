#!/bin/R
### Map QTLs 1 of 3
debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")

library('qtl')
################################################################################
## read in the QTL cross
################################################################################
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

cross <- read.cross(
 file = file.path(mpath, 'ELR.all.csv'),
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


################################################################################
for i in ELR*.csv; do
 cut -d',' -f6- $i > ${i%.csv}.nofirst
done

paste *.nofirst > ELR_marks
cut -d',' -f1-5 ELR_chr_24.csv > ELR.info

paste -d',' ELR.info ELR_marks > ELR.all.scv
awk -F',' '{print $0}' ELR.all.scv > ELR.all.csv
################################################################################


perms.bin.em <- scanone(cross.r, method = "em", model = "binary", maxit = 1000,
  n.perm = 500, pheno.col = 4, n.cluster = 1)
