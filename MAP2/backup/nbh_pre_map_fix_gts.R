#!/bin/R

pop <- 'NBH'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")

library('qtl')

mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]

################################################################################
## read in the QTL cross
cross <- read.cross.jm(file = file.path(indpops, paste(pop, ".unphased.f2.csvr",
  sep = "")), format = "csvr", geno = c(1:3), estimate.map = FALSE)
################################################################################

################################################################################
 cross <- subset(cross,chr=i)
 #nmars <- nmar(cross)
 ### initial order
 ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross)[[1]]))))
 cross <- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
  maxit = 10, tol = 0.001, sex.sp = F)
 #################################################################################

### Pull names from plinkfile
path <- file.path(indpops, paste(pop, ".ped", sep = ""))
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
cross$pheno$ID <- paste(popname, indname, sep = "_")
################################################################################

#### PHENO #####################################################################
cross$pheno$bin <- ifelse(cross$pheno$Pheno > 2, 1 , 0)
cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno))
################################################################################
### ENRICH FOR AAxBB
gt.cp <- geno.table(cross)
gt.cp <- rownames(gt.cp[which(gt.cp$P.value > 0.001),])
cross.par <- subset(cross,ind=c('NBH_NBH1M','NBH_NBH1F'))
cross.par <- pull.markers(cross.par,gt.cp)
gt.cross.par <- geno.table(cross.par)

## DROP DANGEROUS ABxAB cross
DROP <- rownames(gt.cross.par)[which(gt.cross.par$AB==2)]
cross <- drop.markers(cross,DROP)
cross.par <- drop.markers(cross.par,DROP)

### PHASE FROM PARENTS
NBH_NBH1M <- geno.table(subset(cross.par, ind='NBH_NBH1M'))
NBH_NBH1F <- geno.table(subset(cross.par, ind='NBH_NBH1F'))

swit <- rownames(NBH_NBH1M)[which(NBH_NBH1M$AA==1 & NBH_NBH1F$BB==1)]
cross <- switchAlleles(cross, markers = swit)

## FROM THAT SET, DROP
dp1 <- rownames(NBH_NBH1M)[which(NBH_NBH1M$AB==1 & NBH_NBH1F$AA==1)]
dp2 <- rownames(NBH_NBH1F)[which(NBH_NBH1F$AB==1 & NBH_NBH1M$BB==1)]
drops <- unique(c(dp1,dp2))
cross <- drop.markers(cross,drops)


pullgts <- pull.geno(cross)
rownames(pullgts) <- cross$pheno$ID
swit <- colnames(pullgts)[which(pullgts['NBH_NBH1M',]==1)]
cross <- switchAlleles(cross, markers = swit)
swit <- colnames(pullgts)[which(pullgts['NBH_NBH1F',]==3)]
cross <- switchAlleles(cross, markers = swit)
gtpar <- geno.table(subset(cross,ind=is.na(cross$pheno$Pheno)))
likely.par.markers <- rownames(gtpar)[which(gtpar$AA==1 & gtpar$BB==1)]
################################################################################
png(paste0('~/public_html/NBH_noseg_geno',i,'.png'),width=5000,height=2000)
plotGeno(cross ,ind=loc.xocount, cex=2)
dev.off()

cross <- subset(cross,ind=!is.na(cross$pheno$Pheno))
gts <- geno.table(cross)
not_par <- gts[!rownames(gts) %in% likely.par.markers,]
keep_seg <- rownames(not_par[-log(not_par$P.value) < 10,])
keeps <- unique(c(likely.par.markers,keep_seg))
cross <- pull.markers(cross,keeps)

 png(paste0('~/public_html/NBH_missing',i,'.png'))
 hist(gts$missing,breaks=30)
 dev.off()

 png(paste0('~/public_html/NBH_pval',i,'.png'))
 hist(-log(gts[,'P.value']),breaks=50)
 abline(v=5)
 dev.off()

crossbk <- cross
cross <- formLinkageGroups(cross, max.rf = 0.05, min.lod = 15, reorgMarkers = TRUE)
cross <- switchAlleles(cross, markers = markernames(cross,chr=2))
cross <- formLinkageGroups(cross, max.rf = 0.05, min.lod = 15, reorgMarkers = TRUE)
cross <- switchAlleles(cross, markers = markernames(cross,chr=2))
cross <- formLinkageGroups(cross, max.rf = 0.05, min.lod = 15, reorgMarkers = TRUE)
cross <- subset(cross,chr=1)

names(cross$geno) <- i

nmars <- nmar(cross)
 ## initial order
ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross)[[1]]))))
cross <- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
 maxit = 10, tol = 0.001, sex.sp = F)
################################################################################

cross <- subset(cross,ind=!is.na(cross$pheno$Pheno))

cross <- calc.errorlod(cross, err=0.01)
png(paste0('~/public_html/NBH_gts_preclean',i,'.png'),height=2500,width=4500)
plotGeno(cross)
dev.off()

png(paste0('~/public_html/NBH_xo_a',i,'.png'))
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

png(paste0('~/public_html/NBH_gts_preclean_droppedmark',i,'.png'),height=2500,width=4500)
plotGeno(cross)
dev.off()

cross <- cleanGeno_jm(cross, chr=i, maxdist=100, maxmark=8, verbose=TRUE)
cross <- calc.errorlod(cross, err=0.025)
cross <- removeDoubleXO(cross)
cross <- calc.errorlod(cross, err=0.025)
cross <- cleanGeno_jm_2(cross, chr=i, maxdist=50, maxmark=4, verbose=TRUE)
cross <- calc.errorlod(cross, err=0.025)

png(paste0('~/public_html/NBH_cleaned_',i,'.png'),height=2500,width=4000)
plotGeno(cross,cex=3)
dev.off()

png(paste0('~/public_html/NBH_RF_clean',i,'.png'))
plotRF(cross)
dev.off()

fl <- file.path(mpath,paste0(i,'NBH_unmapped_filtered_cleaned'))
write.cross(cross,filestem=fl,format="csv")


###############################################################################
### THIN MARKERS IF NEEDED #####################################################

cross.bk <- subset(cross,ind=!cross$pheno$ID %in% c())
ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross.bk)[[1]]))))
cross.bk <- switch.order(cross.bk, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
maxit = 10, tol = 0.001, sex.sp = F)

mp <- pull.map(cross.bk)
vc <- as.numeric(gsub(".*:",'',names(mp[[as.character(i)]])))
mp <- list(vc)
names(mp) <- i

##attr(mp[[as.character(i)]], "loglik") <- -1028
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

################################################################################

cross <- pull.markers(cross,dwnsmpl)
cross <- calc.genoprob(cross)
cross <- sim.geno(cross)
cross <- calc.errorlod(cross, err=0.01)

png(paste0('~/public_html/NBH_gts_CHR',i,'_downsmpl.unordered.png'),height=1500,width=4500)
plotGeno(cross,cex=3)
dev.off()

cross <- orderMarkers(cross, verbose=FALSE,error.prob=0.01, sex.sp=FALSE,
                 map.function="kosambi",maxit=10000, tol=1e-3, use.ripple=FALSE)

 cross <- calc.errorlod(cross, err=0.01)

 cross_map <-  est.map(tmp, error.prob=0.01,
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

filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/NBH_gts_CHR',i,'_downsmpl_map')
write.cross(cross,chr=i,filestem=filename,format="csv")

png(paste0('~/public_html/NBH_RF_clean_mapped_',i,'.png'))
plotRF(cross)
dev.off()

################################################################################
