#!/bin/R

pop <- 'BRP'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")

################################################################################
## read in the QTL cross
cross <- read.cross.jm(file = file.path(indpops, paste0(pop, ".unphased.f2.csvr")),
format = "csvr", geno = c(1:3), estimate.map = FALSE)
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

pars <- c('BRP_BRP1M','BRP_BRP8F','BRP_BRP1F','BRP_BRP8M')

##################################################################################
##### Switch phase and keep only parent conf markers #############################
##### ENRICH FOR AAxBB ##########################################################
##
#### DROP DANGEROUS ABxAB cross ##################################################
##DROP1M <- pull.geno(cross)[cross$pheno$ID=='BRP_BRP1M',]
##DROP1F <- pull.geno(cross)[cross$pheno$ID=='BRP_BRP1F',]
##DROP8M <- pull.geno(cross)[cross$pheno$ID=='BRP_BRP8M',]
##DROP8F <- pull.geno(cross)[cross$pheno$ID=='BRP_BRP8F',]
##
##m <- names(DROP1M)[which(DROP1M == DROP8M)]
##f <- names(DROP1F)[which(DROP1F == DROP8F)]
##mf <- intersect(m,f)
##cross <- pull.markers(cross, mf)
##
##### SWITCH ALLELES THAT ARE PROB AA x BB #######################################
##bfix <- pull.geno(cross)[cross$pheno$ID=='BRP_BRP8M',]
##bfix_swit1 <- names(bfix)[which(as.numeric(bfix)==1)]
##bfix <- pull.geno(cross)[cross$pheno$ID=='BRP_BRP1F',]
##bfix_swit2 <- names(bfix)[which(as.numeric(bfix)==3)]
##bfix_swit12 <- intersect(bfix_swit1 ,bfix_swit2)
##
##cross <- switchAlleles(cross, markers = bfix_swit12)
##################################################################################
##
##
##################################################################################
##### Get highly likely AB x AB markers ##########################################
##bfix1 <- pull.geno(cross)[cross$pheno$ID=='BRP_BRP8M',]
##bfix1 <- names(bfix1)[which(as.numeric(bfix1)==3)]
##bfix2 <- pull.geno(cross)[cross$pheno$ID=='BRP_BRP1F',]
##bfix2 <- names(bfix2)[which(as.numeric(bfix2)==1)]
##parABxAB <- intersect(bfix1,bfix2)
##
##gt_nopar <- geno.table(subset(cross,ind=!cross$pheno$ID %in% pars))
##parABxAB <- intersect(rownames(gt_nopar[which(gt_nopar$P.value > 0.01),]) ,parABxAB)
##cross.1 <- pull.markers(cross,parABxAB)
##################################################################################

##### TEST SAMPLE GT SIMILARITY ##################################################
##cross.1 <- subset(cross.1,ind=!cross.1$pheno$ID%in% pars)
##cpgt <- comparegeno(cross.1)
##colnames(cpgt) <- cross.1$pheno$ID
##rownames(cpgt) <- cross.1$pheno$ID
##cpgt[cpgt==NaN] <- NA
##diag(cpgt) <- NA
##cpgt <- cpgt[rowSums(is.na(cpgt)) < nind(cross.1),colSums(is.na(cpgt)) < nind(cross.1)]
##################################################################################
png(paste0('~/public_html/BRP_relat.png'))
 hist(cpgt)
dev.off()
################################################################################
toss.missing <- c("BRP_2535","BRP_2410","BRP_2687")
################################################################################
#cross <- cros.bk
################################################################################
#### Pvalue and Missing ##############################################
gt <- geno.table(subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,pars)))
bfixA <- rownames(gt[which(gt$P.value > 0.000001 & gt$missing < 5),])
################################################################################
#cros.bk <- cross
###### FILTER #######################################################
cross <- pull.markers(cross,bfixA)
cross <- subset(cross,ind=!cross$pheno$ID %in% c(toss.missing,pars))
################################################################################

for(Z in 1:24){
 reorg.1 <- formLinkageGroups(subset(cross,chr=Z), max.rf = 0.125, min.lod = 14, reorgMarkers = TRUE)
 swits <- markernames(reorg.1, chr=2)
 reorg.1 <- switchAlleles(reorg.1, markers = markernames(reorg.1,chr=2))
 reorg.2 <- formLinkageGroups(reorg.1, max.rf = 0.125, min.lod = 14, reorgMarkers = TRUE)
 subs <- markernames(reorg.2, chr=1)
 drops <- markernames(reorg.1)[!markernames(reorg.1) %in% subs]
 cross <<- switchAlleles(cross, swits)
 cross <<- drop.markers(cross, drops)
}

mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'
fl <- file.path(mpath,'BRP_unmapped_filtered')
write.cross(cross,filestem=fl,format="csv")
