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
cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno))
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
swit <- colnames(pullgts)[which(pullgts['NBH_NBH1M',]==1)]
cross <- switchAlleles(cross, markers = swit)
swit <- colnames(pullgts)[which(pullgts['NBH_NBH1F',]==3)]
cross <- switchAlleles(cross, markers = swit)
gtpar <- geno.table(subset(cross,ind=is.na(cross$pheno$Pheno)))
likely.par.markers <- rownames(gtpar)[which(gtpar$AA==1 & gtpar$BB==1)]
################################################################################
scan.bin.mr <- scanone(cross, method = "mr", model = "binary", pheno.col = 4)
scan.norm.mr <- scanone(cross, method = "mr", model = "normal", pheno.col = 5)
################################################################################
cross.qtls <-

gt14 <- geno.table(subset(cross,chr=14))
gt8 <- geno.table(subset(cross,chr=8))


8:12345508
1:5540802
14:12487649  missassembly
5:1856088


1:5540802

ns19 <- gt19[which(gt19$P.value > 3.345965e-03),]
which(rownames(ns19)=='19:2435634')
islink <- rownames(ns19[1:200,])
chr19qtls <- pull.markers(cross,islink)


'14:12556225' -> '14:12487649'
chr14 <- subset(cross,chr=14)
gt14 <- geno.table(chr14)
ns <- gt14[which(gt14$P.value >= 0.05),]
chr14 <- pull.markers(cross,rownames(ns))
chr14 <- formLinkageGroups(chr14, max.rf = 0.05, reorgMarkers = TRUE)
gt14 <- geno.table(chr14)
gt14['14:12487649',]

scan.norm.mr14 <- scanone(chr14, method = "mr", model = "normal", pheno.col = 5)


png(paste0('~/public_html/NBH_AHR2.png'))
plotPXG(cross,c('18:17816908','1:5540802'),infer=F)
dev.off()




### NEED TO CHECH ON
5:1856088
################################################################################
#### FILTER BY PARENT ALLELES ##################################################
pullgts <- pull.geno(cross)
rownames(pullgts) <- cross$pheno$ID
gts <- geno.table(cross)
pos <- as.numeric(gsub(".*:",'',rownames(gts)))

is_homs_both_02 <- rownames(marks_filt) %in% colnames(pullgts)[pullgts['NBH_NBH1M',]==3 & pullgts['NBH_NBH1F',]==1]

homz_NA_03 <- rownames(marks_filt) %in% colnames(pullgts)[pullgts['NBH_NBH1M',]==3 & is.na(pullgts['NBH_NBH1F',])]
homz_NB_03 <- rownames(marks_filt) %in% colnames(pullgts)[is.na(pullgts['NBH_NBH1M',]) & pullgts['NBH_NBH1F',]==1]

homz_either_03 <- apply(cbind(homz_NA_03,homz_NB_03),2,any)

marks_filt <- cbind(filter_01,is_homs_both_02,homz_either_03)
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
het_ratio_3 <- rownames(ind_filt) %in% names(rat)[rat < 10]
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
keep_all <- rownames(gts.3)[gts.3$missing<3 & -log(gts.3$P.value) < 3]
#keep_5 <- rownames(gts.3)[gts.3$chr==5 & gts.3$missing<5 & -log(gts.3$P.value) < 5]
#keep_11 <- rownames(gts.3)[gts.3$chr==11 & gts.3$missing<5 & -log(gts.3$P.value) < 5]
#keep_17 <- rownames(gts.3)[gts.3$chr==17 & gts.3$missing<5 & -log(gts.3$P.value) < 5]
#keep_2.5 <- unique(c(keep_all,keep_5,keep_11,keep_17))
################################################################################

## PVAL and MISSING FITER ######################################################
miss_pval2.5_04 <- rownames(marks_filt) %in% keep_all
marks_filt <- cbind(filter_01, is_homs_02, not_NA_03, miss_pval2.5_04)
################################################################################

#save.image('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_final_markerset_unmapped.rsave')
load('/home/jmiller1/QTL_Map_Raw/NBH_final_map/NBH_final_markerset_unmapped.rsave')

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
final_marks <- data.frame(post_filt_marks, stringsAsFactors=T)

post_filt_ind <- cross$pheno$ID %in% cross.4$pheno$ID
names(post_filt_ind) <- rownames(cross$pheno$ID)
final_ind <- data.frame(post_filt_marks, stringsAsFactors=T)
################################################################################

################################################################################
# FINAL MARKER SET ABOVE #######################################################
# ALL PHYS POS COVERED TO HERE #################################################
################################################################################
################################################################################

final_marks <- data.frame(post_filt_01=post_filt_marks, stringsAsFactors=T)
#y18 <- as.numeric(gsub(".*:",'',rownames(final_marks)))
#x18 <- as.numeric(gsub(":.*",'',rownames(final_marks)))

for(i in 1:24){
 inchr <- i
 reorg.lg <- formLinkageGroups(subset(cross.4, chr=inchr), max.rf = 0.1, min.lod = 15, reorgMarkers = TRUE)

 #png(paste0('~/public_html/ER_RF_LG_',i,'.png'))
 #plotRF(reorg.lg,chr=1:4)
 #dev.off()

 nms <- markernames(reorg.lg, chr=1)
 mark <- data.frame(rownames(final_marks) %in% nms, stringsAsFactors=F)
 cur.chr <- paste0('chr',i)
 colnames(mark) <- cur.chr

 final_marks <- cbind(final_marks,mark)

 save.image('/home/jmiller1/QTL_Map_Raw/NBH_final_markerset_unmapped.rsave')

 print(i)

 #yf <- as.numeric(gsub(".*:",'',rownames(final_marks)[final_marks[,cur.chr]]))
 #xf <- as.numeric(gsub(":.*",'',rownames(final_marks)[final_marks[,cur.chr]]))

 ## png(paste0('~/public_html/ER_physLG_',i,'.png'))
 ## plot(x18[x18==i],y18[x18==i])
 ## points(xf,yf,col='green',pch=16)
 ## dev.off()

}
save.image('/home/jmiller1/QTL_Map_Raw/NBH_final_markerset_unmapped.rsave')
################################################################################
## FIX BY HAND ####
################################################################################
## Flip

inchr <-
reorg.lg <- formLinkageGroups(subset(cross.4, chr=inchr), max.rf = 0.1, min.lod = 15, reorgMarkers = TRUE)
swit_18 <- markernames(reorg.lg, chr=2)
cross.4 <- switchAlleles(cross.4, markers = swit_18)

reorg.lg <- formLinkageGroups(subset(cross.4, chr=18), max.rf = 0.1, min.lod = 15, reorgMarkers = TRUE)
nms <- markernames(reorg.lg, chr=1)
mark <- data.frame(rownames(final_marks) %in% nms, stringsAsFactors=F)
colnames(mark) <- 'chr18'
final_marks[,'chr18'] <- mark
reorg.lg <- formLinkageGroups(subset(cross.4, chr=18), max.rf = 0.1, min.lod = 15, reorgMarkers = TRUE)
png(paste0('~/public_html/ER_RF_LG_redux',i,'.png'))
plotRF(reorg.lg, chr=18)
dev.off()

################################################################################
##
################################################################################
inchr <-
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
cross.final$pheno$pheno_norm <- signif(cross.final$pheno$pheno_norm,5)



## for (i in 1:24){
##
##  filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_chr_',i)
##
##  write.cross(cross.final,chr=i,filestem=filename,format="csv")
##
## }

write.table(markernames(cross.final),'/home/jmiller1/QTL_Map_Raw/ELR_final_map/goodmarks.rtable')
write.table(cross.final$pheno$ID,'/home/jmiller1/QTL_Map_Raw/ELR_final_map/goodsamps.rtable')

mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'
fl <- file.path(mpath,'ELR_unmapped_filtered')
write.cross(cross.final,filestem=fl,format="csv")

marks <- rownames(final_marks)[!rowSums(final_marks)==2]
cross.final <- drop.markers(cross,drop)
scan.bin.mr <- scanone(cross.final, method = "mr", model = "binary", pheno.col = 4)
scan.nor.mr <- scanone(cross.final, method = "mr", pheno.col = 5)

################################################################################
