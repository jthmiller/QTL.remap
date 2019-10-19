#!/bin/R
### Map QTLs 1 of 3

pop <- 'ELR'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")

################################################################################
## read in the QTL cross
cross <- read.cross.jm(file = file.path(indpops, paste0(pop, ".unphased.f2.csvr")),
format = "csvr", geno = c(1:3), estimate.map = FALSE)
################################################################################

### Pull names from plinkfile ##################################################
path <- file.path(indpops, paste(pop, ".ped", sep = ""))
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
cross$pheno$ID <- paste(popname, indname, sep = "_")
################################################################################

#### PHENO #####################################################################
cross$pheno$bin <- ifelse(cross$pheno$Pheno > 2, 1 , 0)
cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno))
################################################################################

### Switch phase and keep only parent conf markers##############################
### is 10869 the real ELR mother?
### ENRICH FOR AAxBB
##cross.bk <- cross
## DROP DANGEROUS ABxAB cross
DROP <- pull.geno(cross)[cross$pheno$ID=='BLI_BI1124M',]
DROP <- names(DROP)[which(as.numeric(DROP)==2)]
cross <- drop.markers(cross,DROP)
################################################################################

### SWITCH ALLELES THAT ARE PROB AA x BB #######################################
bfix <- pull.geno(cross)[cross$pheno$ID=='BLI_BI1124M',]
bfix_swit <- names(bfix)[which(as.numeric(bfix)==1)]
gt <- geno.table(cross)
bfix_swit <- intersect(rownames(gt[which(gt$P.value > 0.05),]) ,bfix_swit)
cross <- switchAlleles(cross, markers = bfix_swit)
################################################################################

### Get highly likely AB x AB markers ##########################################
bfix <- pull.geno(cross)[cross$pheno$ID=='BLI_BI1124M',]
bfix <- names(bfix)[which(as.numeric(bfix)==3)]
parABxAB <- intersect(rownames(gt[which(gt$P.value > 0.05),]) ,bfix)
cross.1 <- pull.markers(cross,parABxAB)
################################################################################

### TEST SAMPLE GT SIMILARITY ##################################################
cross.1 <- subset(cross.1,ind=!cross$pheno$ID%in%c('BLI_BI1124M'))
cpgt <- comparegeno(cross.1)
colnames(cpgt) <- cross.1$pheno$ID
rownames(cpgt) <- cross.1$pheno$ID
cpgt[cpgt==NaN] <- NA
diag(cpgt) <- NA
cpgt <- cpgt[rowSums(is.na(cpgt)) < nind(cross.1),colSums(is.na(cpgt)) < nind(cross.1)]
################################################################################

################################################################################
###### Remove the samples related by more than 80% of genotypes #####
wh <- which(cpgt > 0.7, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
mats <- cbind(rownames(wh),colnames(cpgt)[as.numeric(wh[,2])])
toss.missing <- apply(mats,1,function(X){
 X[which.max(c(nmissing(cross.1)[X[1]],nmissing(cross.1)[X[2]]))]
})

## USE MIS_ID'd samples for map, but not QTL
toss.missing <- c(mats[,1],"ELR_10869")
################################################################################

################################################################################
#### FILTER BY Pvalue and Missing ##############################################
cross.1 <- subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,'BLI_BI1124M','ELR_ER1124F'))
gt <- geno.table(cross.1)
bfixA <- rownames(gt[which(gt$P.value > 0.0001 & gt$missing < 5),])
################################################################################
## cross.bk <- cross
###### FILTER #######################################################
cross <- pull.markers(cross,bfixA)
cross <- subset(cross,ind=!cross$pheno$ID %in% c(toss.missing,'BLI_BI1124M','ELR_ER1124F'))
################################################################################

for(Z in 1:24){
 reorg.1 <- formLinkageGroups(subset(cross,chr=Z), max.rf = 0.2, min.lod = 20, reorgMarkers = TRUE)
 swits <- markernames(reorg.1, chr=2)
 reorg.1 <- switchAlleles(reorg.1, markers = markernames(reorg.1,chr=2))
 reorg.2 <- formLinkageGroups(reorg.1, max.rf = 0.1, min.lod = 20, reorgMarkers = TRUE)
 subs <- markernames(reorg.2, chr=1)
 drops <- markernames(reorg.1)[!markernames(reorg.1) %in% subs]
 cross <<- switchAlleles(cross, swits)
 cross <<- drop.markers(cross, drops)
}

mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'
fl <- file.path(mpath,'ELR_unmapped_filtered')
write.cross(cross,filestem=fl,format="csv")
