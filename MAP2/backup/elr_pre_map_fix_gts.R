#!/bin/R
### Map QTLs 1 of 3
i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]
debug.cross <- T
pop <- 'ELR'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")

library('qtl')

mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'


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
 ## initial order

if (i==1){
 CHR1 <- colnames(pull.geno(cross))
 CHR1[CHR1=="AHR2a_del"] <- 343745
 ord <- order(as.numeric(gsub(".*:","",CHR1)))
} else if (i==2){
 CHR2 <- colnames(pull.geno(cross))
 CHR2[CHR2=="AIP_261"] <- 29370504
 CHR2[CHR2=="AIP_252"] <- 29370500
 ord <- order(as.numeric(gsub(".*:","",CHR2)))
} else {
 ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross)[[1]]))))
}

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

png(paste0('~/public_html/ELR_cleaned_',i,'.png'),height=2500,width=4000)
plotGeno(cross,cex=3)
dev.off()

png(paste0('~/public_html/ELR_RF_clean',i,'.png'))
plotRF(cross)
dev.off()

fl <- file.path(mpath,paste0(i,'ELR_unmapped_filtered_cleaned'))
write.cross(cross,filestem=fl,format="csv")

################################################################################
### THIN MARKERS IF NEEDED #####################################################

cross.bk <- subset(cross,ind=!cross$pheno$ID %in% c('ELR_11103','ELR_10869','ELR_ER1124F','ELR_10977','ELR_10988','BLI_BI1124M'))

if (i==1){
 CHR1 <- colnames(pull.geno(cross.bk))
 CHR1[CHR1=="AHR2a_del"] <- 343745
 ord <- order(as.numeric(gsub(".*:","",CHR1)))
} else if (i==2){
 CHR2 <- colnames(pull.geno(cross.bk))
 CHR2[CHR2=="AIP_261"] <- 29370504
 CHR2[CHR2=="AIP_252"] <- 29370500
 ord <- order(as.numeric(gsub(".*:","",CHR2)))
} else {
 ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross.bk)[[1]]))))
}

cross.bk <- switch.order(cross.bk, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
maxit = 10, tol = 0.001, sex.sp = F)

mp <- pull.map(cross.bk)
vc <- as.numeric(gsub(".*:",'',names(mp[[as.character(i)]]) ))
mp <- list(vc)
names(mp) <- i

attr(mp[[as.character(i)]], "loglik") <- -1
attr(mp[[as.character(i)]],"names") <- names(pull.map(cross.bk)[[as.character(i)]])
attr(mp[[as.character(i)]],"class") <- "A"
attr(mp, "class") <- "map"

cross.bk <- replace.map(cross.bk,mp)

gts <- geno.table(cross.bk)
cross.bk <- calc.errorlod(cross.bk, err=0.01)

if( any( top.errorlod(cross.bk, cutoff=0)[,4] > 4)){
 erlod <- top.errorlod(cross.bk)['errorlod']/10
 erind <- !duplicated(top.errorlod(cross.bk)[,'marker'])
 erlod <- erlod[erind,]
 names(erlod) <- top.errorlod(cross.bk)[erind,'marker']
 weight <- 1 - gts$missing/rowSums(gts[,c(3:5)])*10
 weight[names(erlod)] <- weight[names(erlod)] - erlod
} else {
 weight <- 1 - gts$missing/rowSums(gts[,c(3:5)])*10
}

dwnsmpl <- pickMarkerSubset(pull.map(cross.bk)[[1]],1000, weights=weight)

#####MAP ########################################################################
cross <- pull.markers(cross,dwnsmpl)
cross <- subset(cross,ind=!cross$pheno$ID %in% c('ELR_11103','ELR_10869','ELR_ER1124F','BLI_BI1124M'))

cross <- calc.genoprob(cross)
cross <- sim.geno(cross)
cross <- calc.errorlod(cross, err=0.01)

png(paste0('~/public_html/ELR_gts_CHR',i,'_downsmpl.unordered.png'),height=1500,width=4500)
plotGeno(cross,cex=3)
dev.off()

cross <- orderMarkers(cross, verbose=FALSE,error.prob=0.01, sex.sp=FALSE,
                 map.function="kosambi",maxit=10000, tol=1e-3, use.ripple=FALSE)

 cross <- calc.errorlod(cross, err=0.01)

 cross_map <-  est.map(cross, error.prob=0.01,
              map.function="kosambi",
              maxit=1000, tol=1e-4, sex.sp=FALSE,
              verbose=FALSE, n.cluster=6)

 tmp <- qtl:::replace.map(cross,cross_map)

 drp1 <- droponemarker(cross, error.prob=0.01,
                    map.function="kosambi",
                    maxit=100, tol=1e-3, sex.sp=FALSE,
                    verbose=FALSE)

 cross <- dropByDropone(cross, drp1, endMarkerThresh = 20,
                        midMarkerThresh = 20, map.function = "kosambi",
                        re.est.map = T, error.prob=0.01,maxit=1, tol=1e-3, sex.sp=FALSE,
                        verbose=FALSE)

cross <- orderMarkers(cross,verbose=FALSE,chr=i,error.prob=0.01, sex.sp=FALSE,
                        map.function="kosambi",maxit=1000, tol=1e-4, use.ripple=FALSE)

cross <- calc.errorlod(cross, err=0.01)

cross_map <-  est.map(cross,  error.prob=0.01,
              map.function="kosambi",
              maxit=10000, tol=1e-3, sex.sp=FALSE,
              verbose=FALSE, omit.noninformative=TRUE)

cross <- qtl:::replace.map(cross,cross_map)

filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_gts_CHR',i,'_downsmpl_map')
write.cross(cross,chr=i,filestem=filename,format="csv")

################################################################################

ordered <- read.cross(file=filename,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
