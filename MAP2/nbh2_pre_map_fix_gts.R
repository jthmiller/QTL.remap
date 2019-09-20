#!/bin/R
### Map QTLs 1 of 3
debug.cross <- T
pop <- 'NBH'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")

library('qtl')

mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]

#cross <- read.cross.jm(file = file.path(indpops, paste0(pop, ".unphased.f2.csvr")),
#format = "csvr", geno = c(1:3), estimate.map = FALSE)


################################################################################
## read in the QTL cross
cross <- read.cross.jm(file = file.path(indpops, paste(pop, ".unphased.f2.csvr",
  sep = "")), format = "csvr", geno = c(1:3), estimate.map = FALSE)
################################################################################


################################################################################
################################################################################
################################################################################
 cross <- subset(cross,chr=i)
 #nmars <- nmar(cross)
 ### initial order
 #ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross)[[1]]))))
 #cross <- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
#  maxit = 10, tol = 0.001, sex.sp = F)
 #################################################################################

### Pull names from plinkfile
path <- file.path(indpops, paste(pop, ".ped", sep = ""))
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
cross$pheno$ID <- paste(popname, indname, sep = "_")
################################################################################

################################################################################
## FILTER TABLES
################################################################################
pullgts <- pull.geno(cross)
rownames(pullgts) <- cross$pheno$ID

filter_01 <- colnames(pullgts) %in% markernames(cross)
names(filter_01) <- colnames(pullgts)
marks_filt <- data.frame(filter_01, stringsAsFactors=F)

parents_01 <- !is.na(cross$pheno$Phen)
names(parents_01)  <- cross$pheno$ID
ind_filt <- data.frame(parents_01,stringsAsFactors=F)
################################################################################

#### PHENO #####################################################################
cross$pheno$bin <- ifelse(cross$pheno$Pheno > 2, 1 , 0)
cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno))
################################################################################

gt.cp <- geno.table(cross)
gt.cp <- rownames(gt.cp[which(gt.cp$P.value > 0.001),])
cross.par <- subset(cross,ind=c('NBH_NBH1M','NBH_NBH1F'))
cross.par <- pull.markers(cross.par,gt.cp)
gt.cross.par <- geno.table(cross.par)

DROP <- rownames(gt.cross.par)[which(gt.cross.par$AB==2)]

cross <- drop.markers(cross,DROP)
cross.par <- drop.markers(cross.par,DROP)

crsplot <-



ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross.par)[[1]]))))
 cross.par <- switch.order(cross.par, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
 maxit = 10, tol = 0.001, sex.sp = F)

cross.par <- calc.errorlod(cross.par, err=0.1)
png(paste0('~/public_html/NBH_par_geno',i,'.png'),width=5000)
plotGeno(cross.par, cex=2)
dev.off()


crs <- pull.markers(cross,gt.cp)
ord <- order(as.numeric(gsub(".*:","",names(pull.map(crs)[[1]]))))
 crs<- switch.order(crs, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
 maxit = 10, tol = 0.001, sex.sp = F)


 loc.xocount <- crs$pheno$ID[order(countXO(crs))]

png(paste0('~/public_html/NBH_noseg_geno',i,'.png'),width=5000,height=2000)
plotGeno(crs ,ind=loc.xocount, cex=2)
dev.off()




fix.cross <- pull.geno(cross)[,gt.cp]
one_pars <- cross$pheno$ID[order(rowSums(fix.cross==1 | fix.cross==3, na.rm=T))]

cross$pheno$ID[9]

swit <- colnames(pullgts)[which(pullgts['NBH_NBH1M',]==1)]
cross <- switchAlleles(cross, markers = swit)
swit <- colnames(pullgts)[which(pullgts['NBH_NBH1F',]==3)]
cross <- switchAlleles(cross, markers = swit)
gtpar <- geno.table(subset(cross,ind=is.na(cross$pheno$Pheno)))
likely.par.markers <- rownames(gtpar)[which(gtpar$AA==1 & gtpar$BB==1)]
################################################################################
cross <- subset(cross,ind=!is.na(cross$pheno$Pheno))
gts <- geno.table(cross)

 png(paste0('~/public_html/NBH_missing',i,'.png'))
 hist(gts$missing,breaks=30)
 dev.off()

 png(paste0('~/public_html/NBH_pval',i,'.png'))
 hist(-log(gts[,'P.value']),breaks=50)
 abline(v=5)
 dev.off()

not_par <- gts[!rownames(gts) %in% likely.par.markers,]
keep_seg <- rownames(not_par[-log(not_par$P.value) < 10,])
keeps <- unique(c(likely.par.markers,keep_seg))
cross <- pull.markers(cross,keeps)

crs <- formLinkageGroups(cross, max.rf = 0.05, min.lod = 15, reorgMarkers = TRUE)

sum(likely.par.markers %in% markernames(crs,chr=3))

crs <- switchAlleles(crs, markers = markernames(crs,chr=2))
crs <- formLinkageGroups(crs, max.rf = 0.05, min.lod = 15, reorgMarkers = TRUE)
crs <- switchAlleles(crs, markers = markernames(crs,chr=2))
crs <- formLinkageGroups(crs, max.rf = 0.05, min.lod = 15, reorgMarkers = TRUE)

 crs <- subset(crs,chr=i)
 nmars <- nmar(crs)
 ## initial order
 ord <- order(as.numeric(gsub(".*:","",names(pull.map(crs)[[1]]))))
 crs <- switch.order(crs, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
  maxit = 10, tol = 0.001, sex.sp = F)
 ################################################################################


geno.table(crs,chr=2)






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
