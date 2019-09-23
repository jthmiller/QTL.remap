#!/bin/R
### Map QTLs 1 of 3
debug.cross <- T
pop <- 'ELR'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")

library('qtl')

mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]

#cross <- read.cross.jm(file = file.path(indpops, paste0(pop, ".unphased.f2.csvr")),
#format = "csvr", geno = c(1:3), estimate.map = FALSE)



################################################################################
################################################################################
################################################################################

for (i in 3:24){

 cross <- subset(cross.all,chr=i)
 nmars <- nmar(cross)
 ## initial order
 ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross)[[1]]))))
 cross <- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
  maxit = 10, tol = 0.001, sex.sp = F)
 ################################################################################
 ################################################################################

 cross <- subset(cross,ind=!is.na(cross$pheno$Pheno))

 cross <- calc.errorlod(cross, err=0.01)
 png(paste0('~/public_html/ELR_gts_preclean',i,'.png'),height=2500,width=4500)
 plotGeno(cross)
 dev.off()

 png(paste0('~/public_html/ELR_xo_a',i,'.png'))
 hist(sort(table(unlist(locateXO(cross)))),breaks=30)
 dev.off()

 loc.xocount <- table(unlist(locateXO(cross)))

 marker <- sapply(as.numeric(names(loc.xocount)),function(X){
  find.marker(cross, chr=i, pos=X) })

 dropdf <- data.frame(loc.xocount,marker,stringsAsFactors=F)

 dropdf$tot <- sapply(dropdf$mark, function(X){ sum(table(pull.geno(cross,i)[,X]))})
 drops <- unique(dropdf[dropdf$Freq/dropdf$tot > 0.10,'marker'])

 cross <- drop.markers(cross,drops)
 cross <- calc.genoprob(cross)
 cross <- sim.geno(cross)
 cross <- calc.errorlod(cross, err=0.01)

 png(paste0('~/public_html/ELR_gts_preclean_droppedmark',i,'.png'),height=2500,width=4500)
 plotGeno(cross)
 dev.off()

 cross <- cleanGeno_jm(cross, chr=i, maxdist=100, maxmark=8, verbose=TRUE)
 cross <- calc.errorlod(cross, err=0.025)
 cross <- removeDoubleXO(cross)
 cross <- calc.errorlod(cross, err=0.025)
 cross <- cleanGeno_jm_2(cross, chr=i, maxdist=50, maxmark=4, verbose=TRUE)
 cross <- calc.errorlod(cross, err=0.025)

 png(paste0('~/public_html/ELR_clean.png'),height=2500,width=4000)
 plotGeno(cross,cex=3)
 dev.off()

 png(paste0('~/public_html/ELR_RF_clean',i,'.png'))
 plotRF(cross)
 dev.off()

 fl <- file.path(mpath,paste0(i,'ELR_unmapped_unfiltered'))
 write.cross(cross,filestem=fl,format="csv")


################################################################################
### THIN MARKERS IF NEEDED #####################################################

mp <- as.numeric(gsub(".*:",'',markernames(cross)))
names(mp) <- markernames(cross)
mp <- list(mp)
names(mp) <- i
cross <- replace.map(cross,mp)

gts <- geno.table(cross)
weight <- 1 - gts$missing/rowSums(gts[,c(3:5)])*10

dwnsmpl <- pickMarkerSubset(pull.map(cross)[[1]],2000, weights=weight)

drops <- markernames(cross)[! markernames(cross) %in% dwnsmpl]
cross.dwn <- drop.markers(cross,drops)

cross.dwn <- calc.genoprob(cross.dwn)
cross.dwn <- sim.geno(cross.dwn)
cross.dwn <- calc.errorlod(cross.dwn, err=0.01)

png(paste0('~/public_html/ELR_gts_CHR',i,'_downsmpl.png'),height=1500,width=4500)
plotGeno(cross.dwn ,cex=3)
dev.off()

#####MAP ########################################################################

cross.dwn <- subset(cross.dwn,ind=!cross$pheno$ID %in% c('ELR_10869','ELR_ER1124F','ELR_10977','ELR_10988','BLI_BI1124M'))

cross.dwn <- orderMarkers(cross.dwn, window=7,verbose=FALSE,chr=i,
                 use.ripple=TRUE, error.prob=0.025, sex.sp=FALSE,
                 map.function="kosambi",maxit=500, tol=1e-4)

cross.dwn <- calc.genoprob(cross.dwn)
cross.dwn <- sim.geno(cross.dwn)
cross.dwn <- calc.errorlod(cross.dwn, err=0.01)

##cross.dwn <- read.cross(
## file = filename,
## format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
## estimate.map = FALSE
##)
##
cross.dwn_map <-  est.map(cross.dwn,  error.prob=0.025,
            map.function="kosambi",
            maxit=10000, tol=1e-6, sex.sp=FALSE,
            verbose=FALSE, omit.noninformative=TRUE, n.cluster=6)

 cross.dwn_map <- shiftmap(cross.dwn_map, offset=0)

cross.dwn <- replace.map(cross, cross.dwn_map)

filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_gts_CHR',i,'_downsmpl_map')
write.cross(cross.dwn,chr=i,filestem=filename,format="csv")

}
################################################################################
