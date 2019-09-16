#!/bin/R
### Map QTLs 1 of 3
debug.cross <- F
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")

library('qtl')

mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]


print(i)

mis_tol <- 3
cross.1 <- subset(cross,chr=1)
gt <- geno.table(cross.1)

png(paste0('~/public_html/ELR_miss',i,'.png'),width=1000)
plot(gsub(".*:","",rownames(gt)),gt$missing, pch=18)
points(gsub(".*:","",rownames(gt.af)),gt.af$missing, pch=18,col='red')
points(gsub(".*:","",rownames(gt3)),gt3$missing, pch=18,col='green')
dev.off()



keep_all <- rownames(gt)[gt$missing<mis_tol & -log(gt$P.value) < 5]

drop <- rownames(gt)[!rownames(gt) %in% keep_all]

cross.1.ds <- drop.markers(cross.1,drop)

gt.af <- geno.table(cross.1.ds)

cross.1.rf <- formLinkageGroups(cross, max.rf = 0.05, reorgMarkers = TRUE)
swit_18 <- markernames(cross.1.rf, chr=2)
cross.1.rf<- switchAlleles(cross.1.rf, markers = swit_18)
cross.1.rf <- formLinkageGroups(cross.1.rf, max.rf = 0.05, reorgMarkers = TRUE)
swit_18 <- markernames(cross.1.rf, chr=1)
cross.1.rf<- switchAlleles(cross.1.rf, markers = swit_18)
cross.1.rf <- formLinkageGroups(cross.1.rf, max.rf = 0.1, reorgMarkers = TRUE)
cross.1.rf <-  subset(cross.1.rf,chr=1)
################################################################################
################################################################################
fl <- file.path(mpath,'ELR_unmapped_filtered_added_markers.csv')

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
 estimate.map = FALSE
)

################################################################################
################################################################################
################################################################################
cross <- subset(cross,chr=i)
nmars <- nmar(cross)
cross <- subset(cross,ind=nmissing(cross) < (nmars*.5))

## initial order
 ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross)[[1]]))))

cross <- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
 maxit = 10, tol = 0.001, sex.sp = F)
################################################################################
################################################################################


### CROSS FOR LOCAT. BAD MARKERS #####################################
#### RND 1  #####################################
loc.these <- table(unlist(locateXO(cross)))

png(paste0('~/public_html/ELR_xo_a',i,'.png'))
hist(sort(table(unlist(locateXO(cross)))),breaks=30)
dev.off()

loc.these <- names(loc.these[as.numeric(loc.these) > 5])
drops <- sapply(as.numeric(loc.these),function(X){ find.marker(cross, chr=i, pos=X) } )
cross <- drop.markers(cross,drops)

cross <- calc.genoprob(cross)
cross <- sim.geno(cross)
cross <- calc.errorlod(cross, err=0.01)

png(paste0('~/public_html/ELR_gts_a',i,'.png'),width=2000)
plotGeno(cross)
dev.off()


################################################################################
################################################################################
#### RND 2

cross <- calc.genoprob(cross)
cross <- sim.geno(cross)
cross <- calc.errorlod(cross, err=0.01)
cross <- removeDoubleXO(cross)

loc.these <- table(unlist(locateXO(cross)))

png(paste0('~/public_html/ELR_xo_b',i,'.png'))
hist(sort(table(unlist(locateXO(cross)))))
dev.off()

loc.these <- names(loc.these[as.numeric(loc.these) > 5])
drops <- sapply(as.numeric(loc.these),function(X){ find.marker(cross, chr=i, pos=X) } )
cross <- drop.markers(cross,drops)

png(paste0('~/public_html/ELR_gts_b',i,'.png'),width=2000)
plotGeno(cross)
dev.off()

cross <- removeDoubleXO(cross)
cross <- calc.errorlod(cross, err=0.01)

png(paste0('~/public_html/ELR_gts_c',i,'.png'),width=4000,height=1000)
plotGeno(cross)
dev.off()

################################################################################
################################################################################
### THIN MARKERS ###############################################################

mp <- as.numeric(gsub(".*:",'',markernames(cross)))
names(mp) <- markernames(cross)
mp <- list('13'=mp)
cross <- replace.map(cross,mp)

gts <- geno.table(cross)
weight <- 1 - gts$missing/rowSums(gts[,c(3:5)])*10

dwnsmpl <- pickMarkerSubset(pull.map(cross)[[1]],2000, weights=weight)

drops <- markernames(cross)[! markernames(cross) %in% dwnsmpl]
cross <- drop.markers(cross,drops)

cross <- calc.genoprob(cross)
cross <- sim.geno(cross)
cross <- removeDoubleXO(cross)
cross <- calc.errorlod(cross, err=0.01)

png(paste0('~/public_html/ELR_gts_c',i,'.png'),width=1000)
plotGeno(cross)
dev.off()

################################################################################

mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'
write.table(markernames(cross.ss),file.path(mpath,'ER_markers_subst.table'))

fl <- file.path(mpath,'ELR_subsetted')
write.cross(cross.ss,filestem=fl,format="csv")

################################################################################
loglik <- err <- c(0.001, 0.005, 0.01, 0.015, 0.02)
for(i in seq(along=err)) {
 cat(i, "of", length(err), "\n")
 tempmap <- est.map(mapthis, error.prob=err[i])
 loglik[i] <- sum(sapply(tempmap, attr, "loglik"))
}
lod <- (loglik - max(loglik))/log(10)
