#!/bin/R
### Map QTLs 1 of 3
library('qtl')

################################################################################
## read in the QTL cross
################################################################################
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]

fl <- file.path(mpath,'ELR_unmapped_filtered_added_markers.csv')

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
 estimate.map = FALSE
)
################################################################################


if (i==1){
 chr1 <- gsub("AHR2a_del","350000",markernames(cross,chr=1))
 ord <- order(as.numeric(gsub('.*:','',chr1)))
 cross <- switch.order(cross, chr = 1, ord, error.prob = 0.1, map.function = "kosambi",
    maxit = 1000, tol = 0.001, sex.sp = F)

} else if (i==2){
 chr2 <- gsub("AIP_252","37860632",markernames(cross,chr=2))
 ord <- order(as.numeric(gsub('.*:','',chr2)))
 cross <- switch.order(cross, chr = 2, ord, error.prob = 0.1, map.function = "kosambi",
     maxit = 1000, tol = 0.001, sex.sp = F)

} else {
  chr <- markernames(cross,chr=i)
  ord <- order(as.numeric(gsub('.*:','',chr)))
  cross <- switch.order(cross, chr = i, ord, error.prob = 0.1, map.function = "kosambi",
       maxit = 1000, tol = 0.001, sex.sp = F)
}


################################################################################
dups <- findDupMarkers(cross, exact.only = T, adjacent.only = F)
cross <- drop.markers(cross, unlist(dups))
cross <- subset(cross,chr=i)
cross <- est.map(cross, error.prob=0.1, map.function="kosambi",sex.sp=F,chr=i)

################################################################################
cross_nodups <- orderMarkers(cross, window=7,verbose=FALSE,chr=i,
                 use.ripple=TRUE, error.prob=0.01, sex.sp=FALSE,
                 map.function="kosambi",maxit=5, tol=1e-4)

filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_mapped_chr_',i)
write.cross(cross_nodups,chr=i,filestem=filename,format="csv")
################################################################################


for (i in 1:24){

  ord <- order(as.numeric(gsub('.*:','',chr)))
  cross <- switch.order(cross, chr = i, ord, error.prob = 0.1, map.function = "kosambi",
       maxit = 10, tol = 0.001, sex.sp = F)
}
mp <- pull.map(cross)

loca <- lapply(mp, function(X){
 setNames(gsub("AHR2a_del","350000",names(X)),names(X))
})

loca <- lapply(loca, function(X){
 gsub("AIP_252","37860632",X)
})
loca <- lapply(loca, function(X){
 gsub("AIP_261","37860640",X)
})
loca <- lapply(loca, function(X){
 setNames(as.numeric(gsub('.*:','',X)),names(X))
})

cross.ss <- replace.map(cross,loca)
subs <- sapply(pull.map(cross.ss),pickMarkerSubset,1000)
subs <- unname(unlist(subs))
drops <- markernames(cross.ss)[! markernames(cross.ss) %in% subs]
cross.ss <- drop.markers(cross.ss,drops)

mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'
write.table(markernames(cross.ss),file.path(mpath,'ER_markers_subst.table'))

fl <- file.path(mpath,'ELR_subsetted')
write.cross(cross.ss,filestem=fl,format="csv")
