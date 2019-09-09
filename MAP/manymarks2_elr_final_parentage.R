#!/bin/R
### Map QTLs 1 of 3
debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
## For plotting
marker_dens <- list()

# Table of Chroms with sig QTLs
test.QTLs <- read.table(file.path(basedir, "rQTL/metadata/QTLs.txt"), sep = "\t",
  header = T)

## Get chrom number vector
test.QTLs$chrm.n <- gsub("chr", "", test.QTLs$chrom)

################################################################################
## read in the QTL cross
cross <- read.cross.jm(file = file.path(indpops, paste(pop, ".unphased.f2.csvr",
  sep = "")), format = "csvr", geno = c(1:3), estimate.map = FALSE)
################################################################################

################################################################################
### Pull names from plinkfile
path <- file.path(indpops, paste(pop, ".ped", sep = ""))
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
cross$pheno$ID <- paste(popname, indname, sep = "_")
################################################################################

#### PHENO #####################################################################
cross$pheno$bin <- ifelse(cross$pheno$Pheno > 2, 1 , 0)
cross$pheno$pheno_norm <- nqrank(cross$pheno$Pheno)
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

################################################################################
### Switch phase and keep only parent conf markers##############################
swit <- colnames(pullgts)[which(pullgts['BLI_BI1124M',]==1)]
cross <- switchAlleles(cross, markers = swit)
################################################################################

################################################################################
#### FILTER BY PARENT ALLELES ##################################################
pullgts <- pull.geno(cross)
rownames(pullgts) <- cross$pheno$ID
gts <- geno.table(cross)
pos <- as.numeric(gsub(".*:",'',rownames(gts)))
is_homs_02 <- rownames(marks_filt) %in% colnames(pullgts)[pullgts['BLI_BI1124M',]==3]
not_NA_03 <- rownames(marks_filt) %in% colnames(pullgts)[!is.na(pullgts['BLI_BI1124M',])]
marks_filt <- cbind(filter_01,is_homs_02,not_NA_03)
################################################################################

################################################################################
### TEST SAMPLE GT SIMILARITY ##################################################
drop <- rownames(marks_filt)[!rowSums(marks_filt[,1:3])==3]
cross.1 <- subset(drop.markers(cross,drop), ind = ind_filt[,1])
cpgt <- comparegeno(cross.1)
colnames(cpgt) <- cross.1$pheno$ID
rownames(cpgt) <- cross.1$pheno$ID
cpgt[cpgt==NaN] <- NA
diag(cpgt) <- NA
cpgt <- cpgt[rowSums(is.na(cpgt)) < nind(cross.1),colSums(is.na(cpgt)) < nind(cross.1)]
################################################################################

################################################################################
###### Remove the samples related by more than 60% of genotypes #####
wh <- which(cpgt > 0.6, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
mats <- cbind(rownames(wh),colnames(cpgt)[as.numeric(wh[,2])])
toss.missing <- apply(mats,1,function(X){
 X[which.max(c(nmissing(cross.1)[X[1]],nmissing(cross.1)[X[2]]))]
})
###### SAME GENOS FILTER #######################################################
same_geno_02 <- !cross$pheno$ID %in% toss.missing
ind_filt <- cbind(parents_01,same_geno_02)
################################################################################

################################################################################
#### RATIO OF HETS TO HOM. SOMETHING IS OFF IN SOME SAMPLES. Remove them #######
################################################################################
mar_index <- rownames(marks_filt)[rowSums(marks_filt[,c(1:3)])==3]
ind_index <- cross$pheno$ID %in% rownames(ind_filt)[rowSums(ind_filt[,c(1,2)])==2]

drop <- markernames(cross)[markernames(cross) %in% mar_index]
cross.2 <- subset( drop.markers(cross,drop), ind = ind_index)
pullgts.2 <- pull.geno(cross.2)

homb <- apply(pullgts.2,1,function(X){sum(X==3, na.rm=T) })
homa <- apply(pullgts.2,1,function(X){ sum(X==1, na.rm=T) })
het <- apply(pullgts.2,1,function(X){ sum(X==2, na.rm=T)})
rat <- (homa + homb)/het
names(rat) <- cross.2$pheno$ID

## HET v HOMZ RATIO FILTER #####################################################
het_ratio_3 <- rownames(ind_filt) %in% names(rat)[rat < 2.55]
ind_filt <- cbind(parents_01,same_geno_02,het_ratio_3)
################################################################################

################################################################################
### Pvalue FILTER and PLOT #####################################################
### Marker filters = 3, Inds = 3
################################################################################
drop <- rownames(marks_filt)[!rowSums(marks_filt)==3]
ind_bool <- cross$pheno$ID %in% rownames(ind_filt)[rowSums(ind_filt[,c(1:3)])==3]
cross.3 <- subset(drop.markers(cross,drop), ind = ind_bool)
pullgts.3 <- pull.geno(cross.3)
gts.3 <- geno.table(cross.3)
pos.3 <- as.numeric(gsub(".*:",'',rownames(gts.3)))

######SET MISSING AND PVALU FILTER #############################################
keep_all <- rownames(gts.3)[gts.3$missing<3 & -log(gts.3$P.value) < 2.5]
keep_5 <- rownames(gts.3)[gts.3$chr==5 & gts.3$missing<5 & -log(gts.3$P.value) < 5]
keep_11 <- rownames(gts.3)[gts.3$chr==11 & gts.3$missing<5 & -log(gts.3$P.value) < 5]
keep_17 <- rownames(gts.3)[gts.3$chr==17 & gts.3$missing<5 & -log(gts.3$P.value) < 5]
keep_2.5 <- unique(c(keep_all,keep_5,keep_11,keep_17))
################################################################################

## PVAL and MISSING FITER ######################################################
miss_pval2.5_04 <- rownames(marks_filt) %in% keep_2.5
marks_filt <- cbind(filter_01, is_homs_02, not_NA_03, miss_pval2.5_04)
################################################################################

save.image('ELR_final_markerset_unmapped.rsave')

################################################################################
## FINAL TABLES ##
################################################################################
drop <- rownames(marks_filt)[!rowSums(marks_filt)==4]
ind_bool <- cross$pheno$ID %in% rownames(ind_filt)[rowSums(ind_filt[,c(1:3)])==3]
cross.4 <- subset(drop.markers(cross,drop), ind = ind_bool)

gts.4 <- geno.table(cross.4)
pos.4 <- as.numeric(gsub(".*:",'',rownames(gts.4)))

pullgts.4 <- pull.geno(cross.4)
rownames(pullgts.4) <- cross.4$pheno$ID

post_filt_marks <- rownames(marks_filt) %in% colnames(pullgts.4)
names(post_filt_marks) <- rownames(marks_filt)
final_marks <- data.frame(post_filt_01, stringsAsFactors=T)

post_filt_ind <- cross$pheno$ID %in% cross.4$pheno$ID
names(post_filt_ind) <- rownames(cross$pheno$ID)
final_ind <- data.frame(post_filt_01, stringsAsFactors=T)
################################################################################

################################################################################
# FINAL MARKER SET ABOVE #######################################################
# ALL PHYS POS COVERED TO HERE #################################################
################################################################################
################################################################################

final_marks <- data.frame(post_filt_01, stringsAsFactors=T)
y18 <- as.numeric(gsub(".*:",'',rownames(final_marks)))
x18 <- as.numeric(gsub(":.*",'',rownames(final_marks)))

for(i in 1:24){
 inchr <- i
 reorg.lg <- formLinkageGroups(subset(cross.4, chr=inchr), max.rf = 0.1, min.lod = 15, reorgMarkers = TRUE)

 png(paste0('~/public_html/ER_RF_LG_',i,'.png'))
 plotRF(reorg.lg,chr=1:4)
 dev.off()

 nms <- markernames(reorg.lg, chr=1)
 mark <- data.frame(rownames(final_marks) %in% nms, stringsAsFactors=F)
 cur.chr <- paste0('chr',i)
 colnames(mark) <- cur.chr

 final_marks <- cbind(final_marks,mark)

 save.image('ELR_final_markerset_unmapped.rsave')

 print(i)

 yf <- as.numeric(gsub(".*:",'',rownames(final_marks)[final_marks[,cur.chr]]))
 xf <- as.numeric(gsub(":.*",'',rownames(final_marks)[final_marks[,cur.chr]]))

 png(paste0('~/public_html/ER_physLG_',i,'.png'))
 plot(x18[x18==i],y18[x18==i])
 points(xf,yf,col='green',pch=16)
 dev.off()

}
################################################################################
## FIX BY HAND ####
################################################################################
## Flip Chr 2 on 18

inchr <- 18
reorg.lg <- formLinkageGroups(subset(cross.4, chr=inchr), max.rf = 0.1, min.lod = 15, reorgMarkers = TRUE)
swit_18 <- markernames(reorg.lg, chr=2)
cross.4 <- switchAlleles(cross.4, markers = swit_18)

reorg.lg <- formLinkageGroups(subset(cross.4, chr=18), max.rf = 0.1, min.lod = 15, reorgMarkers = TRUE)
nms <- markernames(reorg.lg, chr=1)
mark <- data.frame(rownames(final_marks) %in% nms, stringsAsFactors=F)
colnames(mark) <- 'chr18'
final_marks[,'chr18'] <- mark
reorg.lg <- formLinkageGroups(subset(cross.4, chr=18), max.rf = 0.1, min.lod = 15, reorgMarkers = TRUE)
png(paste0('~/public_html/ER_RF_LG_',i,'.png'))
plotRF(reorg.lg, chr=18)
dev.off()

################################################################################
## Lower CHR 10 linkage
################################################################################
inchr <- 10
reorg.lg <- formLinkageGroups(subset(cross.4, chr=inchr), max.rf = 0.1, min.lod = 12, reorgMarkers = TRUE)
nms <- markernames(reorg.lg, chr=1)
mark <- data.frame(rownames(final_marks) %in% nms, stringsAsFactors=F)
colnames(mark) <- 'chr10'
final_marks[,'chr10'] <- mark




################################################################################
cross.final <- switchAlleles(cross, markers = swit_18)
drops <- rownames(final_marks)[!rowSums(final_marks[,c(2:25)])==1]
ind_bool <- cross$pheno$ID %in% rownames(ind_filt)[rowSums(ind_filt[,c(1:3)])==3]
cross.final <- subset(drop.markers(cross.final,drop), ind = ind_bool)
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

for (i in 1:24){

 filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_chr_',i)

 write.cross(cross.final,chr=i,filestem=filename,format="csv")

}

dups <- findDupMarkers(cross.final, exact.only = T, adjacent.only = F)
cross.final.nodups <- drop.markers(cross.final, unlist(dups))
################################################################################









cross.18 <- read.cross.jm(file = file.path(indpops, paste(pop, ".unphased.f2.csvr",
  sep = "")), format = "csvr", geno = c(1:3), estimate.map = FALSE)

### Pull names from plinkfile
path <- file.path(indpops, paste(pop, ".ped", sep = ""))
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
cross.18$pheno$ID <- paste(popname, indname, sep = "_")

indo <- cross.18$pheno$ID[!cross.18$pheno$ID %in% crs.bk$pheno$ID]
cross.18 <- subset(cross.18, ind = indo)

cross.18 <- subset(cross.18, ind = !is.na(cross.18$pheno$Phen))

gt.missing <- geno.table(cross.18)

cross.18 <- drop.markers(cross.18, rownames(gt.missing)[gt.missing$missing > 4])
cross.18 <- drop.markers(cross.18, rownames(gt.missing[gt.missing$P.value < 0.05,]))

cross.18$pheno$bin <- ifelse(cross.18$pheno$Pheno > 2, 1 , 0)

scan.bin.mr <- scanone(cross.18, method = "mr", model = "binary", pheno.col = 4)






crs.bk$pheno$pheno_norm <- nqrank(crs.bk$pheno$Pheno)




marker.warning()


gt.missing <- geno.table(cross.18)

#drops <- markernames(cross.18)[!markernames(cross.18) %in% marks.crs.bk]
#cross.18 <- drop.markers(cross.18,drops)

cross.18$pheno$bin <- ifelse(cross.18$pheno$Pheno > 2, 1 , 0)

cross.18 <- subset(cross.18, ind = !is.na(cross.18$pheno$Phen))

scan.bin.mr.18 <- scanone(cross.18, method = "mr", model = "binary", pheno.col = 4)

scan.np.em.18 <- scanone(cross.18, method = "em", model = "np", pheno.col = 3, maxit = 5000)











png('~/public_html/Geno_ER1.png')
plot(sort(rat))
dev.off()


gt.table.ratio <- apply



png('~/public_html/Geno_ER1.png')
plotGeno(cross.18.reorg1,chr=1)
dev.off()

png('~/public_html/RF_ER1.png')
plotRF(cross.18.reorg1,chr=1)
dev.off()

cross.18.reorg1 <- subset(cross.18.strict, chr=1)
cross.18.reorg1 <- formLinkageGroups(cross.18.reorg1, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)



cross.18.reorg2.rf <- formLinkageGroups(cross.18.reorg2, max.rf = 0.05, reorgMarkers = TRUE)

png('~/public_html/RF_ER2.png')
plotRF(cross.18.reorg2.rf)
dev.off()

cross.18.reorg2.rf <- switchAlleles(cross.18.reorg2.rf, markernames(cross.18.reorg2.rf, chr=4 ))



png('~/public_html/Geno_ER2.png',width=1000,height=1000)
plotGeno(cross.18.reorg2,chr=1)
dev.off()

png('~/public_html/RF_ER2.png')
plotRF(cross.18.reorg1,chr=1:2)
dev.off()



gts <- geno.table(cross.18)

cross.18 <- cross.18.reorg

##cross.18 <- subset(cross.18.reorg, chr=1)
##cross.18 <- movemarker(cross.18, markernames(cross.18), newchr=2)

crs.bk <- cross.18.reorg

for (i in lenght(1:24)){

ord <- order(as.numeric(gsub(paste(i, ":", sep = ""), "", markernames(cross.18, chr = i))))

crs.bk <<- switch.order(cross.18, chr = i, ord, error.prob = 0.05, map.function = "kosambi",
  maxit = 1000, tol = 0.001, sex.sp = F)
}

png('~/public_html/Geno_ER.png',width=1200,height=1000)
plotGeno(cross.18,chr=1)
dev.off()


cross.18 <- removeDoubleXO(cross.18, verbose = T)
print("Done removing dxo..")

png('~/public_html/Geno_ER.png',width=4000,height=2000)
plotGeno(cross.18,chr=1)
dev.off()






png('~/public_html/pval.png')
hist(log10(gt.missing$P.value),breaks=50)
dev.off()


cross.18 <- removeDoubleXO(cross.18, verbose = T)
print("Done removing dxo..")








fileConn <- file(file.path(popdir, paste(X, "BB_markes.txt", sep = "_")))
writeLines(swit, fileConn)
close(fileConn)

print("estimating map with markers at physical positions")
ord <- order(as.numeric(gsub(paste(X, ":", sep = ""), "", markernames(cross.18, chr = X))))

cross.18 <- switch.order(cross.18, chr = X, ord, error.prob = 0.01, map.function = "kosambi",
  maxit = 1000, tol = 0.001, sex.sp = F)

swit <- checkAlleles(cross.18, threshold = 6, verbose = F)
cross.18 <- switchAlleles(cross.18, swit)

cross.18 <- formLinkageGroups(cross.18, max.rf = 0.25, min.lod = 12, reorgMarkers = TRUE)
cross.18 <- subset(cross.18, chr = which.max(nmar(cross.18)))
names(cross.18$geno) <- X

print("estimating map with markers at physical positions")
ord <- order(as.numeric(gsub(paste(X, ":", sep = ""), "", markernames(cross.18, chr = X))))

cross.18 <- switch.order(cross.18, chr = X, ord, error.prob = 0.01, map.function = "kosambi",
  maxit = 1000, tol = 0.001, sex.sp = F)

marker.warning()

cross.18 <- removeDoubleXO(cross.18, verbose = T)
print("Done removing dxo..")

cross.18 <- drop.missing(cross.18, miss)

marker.warning()

POS.map.18 <- est.map(cross.18, error.prob = 0.1, map.function = "kosambi", chr = X,
  maxit = 1000)

cross.18 <- replace.map(cross.18, POS.map.18)

print(summary(pull.map(cross.18))[as.character(X), ])

print("Writing the markers to rQTL format")
write.cross(cross.18, filestem = paste(popdir, "/chr", X, "_", outname, ".manymarks.QTLmap",
  sep = ""), format = "csv", chr = X)


#####
## scratch
####

##############################
##############################
pos <- as.numeric(gsub(".*:","",rownames(gt.missing)))
names(pos) <- rownames(gt.missing)
head(sort(abs(pos -  343835)))

1:317181

1:363497

crs.bk

chr1gts <- pull.geno(crs.bk, 1)

chr1phn <- pull.pheno(crs.bk, 1)



chr1gts <- pull.geno(cross.18, 1)

chr1phn <- pull.pheno(cross.18, 1)


AHR <- cbind(chr1phn,chr1gts[,'1:317181'],chr1gts[,'1:363497'])

AHR <- cbind(chr1phn,chr1gts[,'1:317181'],chr1gts[,'1:363497'])
AHR <- AHR[order(AHR[,1]),]

table(AHR[AHR[,1]<2,3])
table(AHR[AHR[,1]>2,3])

table(AHR[AHR[,1]==0,3])
table(AHR[AHR[,1]==1,3])
table(AHR[AHR[,1]==4,3])
table(AHR[AHR[,1]==5,3])

chr1.pars <- pull.geno(cross.pars, 1)
rbind(chr1.pars[,'1:317181'], chr1.pars[,'1:363497'])

TAKE THE HOMZYGOUS GENOTYPES FOR THE ONE PARENT AND SEE IF THEY TEND TOWARD 1:2:1 compared to
het in parent.

AHR[,1] <- as.factor(AHR[,1])
AHR[,2] <- as.factor(AHR[,2])
AHR[,3] <- as.factor(AHR[,3])

png('~/public_html/ER_AHR.png')
plot(table(AHR[AHR[,1]<2,2]))
dev.off()

for(
table(AHR[,1])

print("Removing duplicates")
##dups <- findDupMarkers(cross.18, exact.only = F, adjacent.only = F)
##cross.18 <- drop.markers(cross.18, unlist(dups))
##confirm ahr2a 343745   343931 AHR2a
##mid is 343835
