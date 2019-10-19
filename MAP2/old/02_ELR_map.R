#!/bin/R
### Map QTLs 1 of 3
#debug.cross <- T
#source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library('qtl')
pop <- 'ELR'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

################################################################################
i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'
fl <- file.path('ELR_unmapped_filtered.csv')
cross <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
cross <- subset(cross,chr=i)
################################################################################

################################################################################
cross <- cleanGeno_jm(cross, chr=i, maxdist=100, maxmark=8, verbose=TRUE)
cross <- calc.errorlod(cross, err=0.025)
cross <- removeDoubleXO(cross, chr=i)
cross <- calc.errorlod(cross, err=0.025)
cross <- cleanGeno_jm_2(cross, chr=i, maxdist=50, maxmark=4, verbose=TRUE)
cross <- calc.errorlod(cross, err=0.025)
################################################################################

################################################################################
png(paste0('~/public_html/ELR_RF_',i,'.png'))
  plotRF(cross, chr=i)
dev.off()
png(paste0('~/public_html/ELR_gts_cleaned',i,'.png'),height=2500,width=4500)
  plotGeno(cross, chr=i, cex=2)
dev.off()
################################################################################

################################################################################
mp <- pull.map(cross)
vc <- as.numeric(gsub(".*:",'',names(mp[[as.character(i)]]) ))
mp <- list(vc)
names(mp) <- i

attr(mp[[as.character(i)]], "loglik") <- -1
attr(mp[[as.character(i)]],"names") <- names(pull.map(cross)[[as.character(i)]])
attr(mp[[as.character(i)]],"class") <- "A"
attr(mp, "class") <- "map"

cross <- replace.map(cross,mp)
gts <- geno.table(cross)
cross <- calc.errorlod(cross, err=0.01)

if( any( top.errorlod(cross, cutoff=0)[,4] > 4)){
 erlod <- top.errorlod(cross)['errorlod']/10
 erind <- !duplicated(top.errorlod(cross)[,'marker'])
 erlod <- erlod[erind,]
 names(erlod) <- top.errorlod(cross)[erind,'marker']
 weight <- 1 - gts$missing/rowSums(gts[,c(3:5)])*10
 weight[names(erlod)] <- weight[names(erlod)] - erlod
} else {
 weight <- 1 - gts$missing/rowSums(gts[,c(3:5)])*10
}

dwnsmpl <- pickMarkerSubset(pull.map(cross)[[as.character(i)]],1000, weights=weight)

cross <- pull.markers(cross,dwnsmpl)
################################################################################

################################################################################
cross_map <-  est.map(cross, error.prob=0.01,map.function="kosambi",maxit=10000,
  tol=1e-4, sex.sp=FALSE, verbose=FALSE)

cross_map <- shiftmap(cross_map, offset=0)

cross <- qtl:::replace.map(cross,cross_map)

cross <- orderMarkers(cross, verbose=FALSE,use.ripple=FALSE, error.prob=0.01,
  sex.sp=FALSE, map.function="kosambi", maxit=100000, tol=1e-7)

cross <- calc.errorlod(cross, err=0.01)

cross_map <-  est.map(cross, error.prob=0.01,map.function="kosambi",maxit=100000,
  tol=1e-4, sex.sp=FALSE, verbose=FALSE)

cross_map <- shiftmap(cross_map, offset=0)

cross <- qtl:::replace.map(cross,cross_map)

#drp1 <- droponemarker(cross, error.prob=0.01,map.function="kosambi",maxit=100000,
##  tol=1e-3, sex.sp=FALSE,verbose=FALSE)

##cross <- dropByDropone(cross, drp1, endMarkerThresh = 20, midMarkerThresh = 20,
##  map.function = "kosambi", re.est.map = T, error.prob=0.01, maxit=100000,
##  tol=1e-4, sex.sp=FALSE, verbose=FALSE)

png(paste0('~/public_html/ELR_RF_reordered',i,'.png'))
 plotRF(cross, chr=i)
dev.off()

filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_gts_CHR',i,'_downsmpl_reordered')
write.cross(cross,chr=i,filestem=filename,format="csv")
################################################################################
## END #########################################################################
