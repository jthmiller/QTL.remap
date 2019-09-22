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


fl <- file.path(mpath,'ELR_unmapped_filtered_added_markers.csv')

cross.all <- read.cross(
 file = fl,
 format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
 estimate.map = FALSE
)

################################################################################
################################################################################
################################################################################

 cross <- subset(cross.all,chr=i)
 nmars <- nmar(cross)
 ## initial order

if (i==1){
CHR1 <- colnames(pull.geno(cross))
CHR1[CHR1=="AHR2a_del"] <- 350000
ord <- order(as.numeric(gsub(".*:","",CHR1)))
} else if (i==2){


}

 cross <- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
 maxit = 10, tol = 0.001, sex.sp = F)

 rm(cross.all)
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

 fl <- file.path(mpath,paste0(i,'ELR_unmapped_filtered_cleaned'))
 write.cross(cross,filestem=fl,format="csv")


################################################################################
### THIN MARKERS IF NEEDED #####################################################

cross <- subset(cross,ind=!cross$pheno$ID %in% c('ELR_10869','ELR_ER1124F','ELR_10977','ELR_10988','BLI_BI1124M'))
ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross)[[1]]))))
cross <- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
maxit = 10, tol = 0.001, sex.sp = F)

mp <- as.numeric(gsub(".*:",'',markernames(cross)))
names(mp) <- markernames(cross)
mp <- list(mp)
names(mp) <- i
cross <- replace.map(cross,mp)

gts <- geno.table(cross)

cross <- calc.errorlod(cross, err=0.025)
erlod <- top.errorlod(cross)['errorlod']/10
erind <- !duplicated(top.errorlod(cross)[,'marker'])
erlod <- erlod[erind,]
names(erlod) <- top.errorlod(cross)[erind,'marker']

weight <- 1 - gts$missing/rowSums(gts[,c(3:5)])*10
weight[names(erlod)] <- weight[names(erlod)] - erlod

dwnsmpl <- pickMarkerSubset(pull.map(cross)[[1]],1000, weights=weight)

drops <- markernames(cross)[! markernames(cross) %in% dwnsmpl]

cross <- drop.markers(cross,drops)

cross <- calc.genoprob(cross)
cross <- sim.geno(cross)
cross <- calc.errorlod(cross, err=0.01)

png(paste0('~/public_html/ELR_gts_CHR',i,'_downsmpl.png'),height=1500,width=4500)
plotGeno(cross ,cex=3)
dev.off()

#####MAP ########################################################################

cross <- orderMarkers(cross, window=5,verbose=FALSE,chr=i,
                 use.ripple=TRUE, error.prob=0.01, sex.sp=FALSE,
                 map.function="kosambi",maxit=1000, tol=1e-4)

cross <- calc.genoprob(cross)
cross <- sim.geno(cross)
cross <- calc.errorlod(cross, err=0.01)

  cross_map <-  est.map(cross,  error.prob=0.01,
              map.function="kosambi",
              maxit=10000, tol=1e-4, sex.sp=FALSE,
              verbose=FALSE, omit.noninformative=TRUE, n.cluster=6)

 cross_map <- shiftmap(cross_map, offset=0)

cross <- replace.map(cross, cross_map)

filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_gts_CHR',i,'_downsmpl_map')
write.cross(cross,chr=i,filestem=filename,format="csv")

################################################################################
