#!/bin/R
### Map QTLs 1 of 3
pop <- 'NBH'
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


cross <- subset(cross,chr=c(1,2,3,8,13,18,24))
toss.missing <- c("NBH_5525","NBH_6177","NBH_5528","NBH_5951")
cross <- subset(cross,ind=!cross$pheno$ID %in% toss.missing)

################################################################################
### Switch phase and keep only parent conf markers #############################
### ENRICH FOR AAxBB ##########################################################

## DROP DANGEROUS ABxAB cross ##################################################
DROP1 <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1M',]
DROP1 <- names(DROP1)[which(as.numeric(DROP1)==2)]
DROP2 <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1F',]
DROP2 <- names(DROP2)[which(as.numeric(DROP2)==2)]
DROP <- intersect(DROP1,DROP2)

cross <- drop.markers(cross,DROP)
################################################################################

### SWITCH ALLELES THAT ARE PROB AA x BB #######################################
bfix <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1M',]
bfix_swit1 <- names(bfix)[which(as.numeric(bfix)==1)]
bfix <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1F',]
bfix_swit2 <- names(bfix)[which(as.numeric(bfix)==3)]
bfix_swit12 <- intersect(bfix_swit1 ,bfix_swit2)

cross <- switchAlleles(cross, markers = bfix_swit12)
################################################################################


################################################################################
### Get highly likely AB x AB markers ##########################################
bfix1 <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1M',]
bfix1 <- names(bfix1)[which(as.numeric(bfix1)==3)]
bfix2 <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1F',]
bfix2 <- names(bfix2)[which(as.numeric(bfix2)==1)]
parABxAB <- intersect(bfix1,bfix2)

gt_nopar <- geno.table(subset(cross,ind=!cross$pheno$ID %in% c('NBH_NBH1M','NBH_NBH1F')))
parABxAB <- intersect(rownames(gt_nopar[which(gt_nopar$P.value > 0.01),]) ,parABxAB)
cross.1 <- pull.markers(cross,parABxAB)
cross.1 <- subset(cross.1,ind=!cross.1$pheno$ID %in% c('NBH_NBH1M','NBH_NBH1F'))
gts <- geno.table(cross.1)
bfixA <- rownames(gts[which(gts$P.value > 0.01),])
cross.1 <- pull.markers(cross.1,bfixA)

################################################################################


################################################################################
###### Remove the samples related by more than 80% of genotypes #####

## USE MIS_ID'd samples for map, but not QTL

################################################################################

################################################################################
#### Pvalue and Missing ##############################################
gt <- geno.table(subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F')))
bfixA <- rownames(gt[which(gt$missing < 1),])
################################################################################
cros.bk <- cross
###### FILTER #######################################################
cross <- pull.markers(cross,bfixA)
cross <- subset(cross,ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))

gts <- geno.table(cross)
cross <-  drop.markers(cross,rownames(gts)[which(gts$AA < 3 | gts$BB < 3)])
################################################################################

gts <- geno.table(cross)
bfixA <- rownames(gts[which(gts$P.value > 0.01),])
cross <- pull.markers(cross,bfixA)


for(Z in chrnames(cross)){
 reorg.1 <- subset(cross,chr=Z)
 reorg.1 <- markerlrt(reorg.1)
 reorg.1 <- formLinkageGroups(reorg.1, min.lod = 20, reorgMarkers = TRUE)
 reorg.2 <- subset(reorg.1,chr=1)
 subs <- markernames(reorg.2, chr=1)
 drops <- markernames(reorg.1)[!markernames(cross,chr=Z) %in% subs]
 cross <<- drop.markers(cross, drops)
}


crossmap <-tspOrder(cross = cross, hamiltonian = TRUE, method="concorde", concorde_path='/home/jmiller1/concorde_build/TSP/')

save.image('NBH_seg.dist')

 png(paste("/home/jmiller1/public_html/nbhcrossmap.png"))
  plotRF(crossmap)
 dev.off()



#### Maybe markerlrt or es.rf can help? what are the strongest between chromosome differences in recombination frequencies?



load('scantwo.scans.nbh.rsave')



 png(paste("/home/jmiller1/public_html/nbh",Z,".png"))
  plot(sort(-log10(geno.table(reorg.1)[,'P.value'])))
 dev.off()
 reorg.1 <- est.rf(reorg.1)
 reorg.1 <- formLinkageGroups(reorg.1, max.rf = 0.125, min.lod = 20, reorgMarkers = TRUE)
 swits <- markernames(reorg.1, chr=2)

 reorg.1 <- switchAlleles(reorg.1, markers = markernames(reorg.1,chr=2))
 reorg.2 <- formLinkageGroups(reorg.1, max.rf = 0.125, min.lod = 20, reorgMarkers = TRUE)

 subs <- markernames(reorg.2, chr=1)
 drops <- markernames(reorg.1)[!markernames(reorg.1) %in% subs]
 cross <<- switchAlleles(cross, swits)



mpath <- '/home/jmiller1/QTL_Map_Raw/ELRp'
fl <- file.path(mpath,'NBH_unmapped_filtered')
write.cross(cross,filestem=fl,format="csv")



addcovar=NULL, intcovar=NULL



sc2_cov <- scantwo(gg_step2, pheno.col=5, model="normal",
             method="imp",
             addcovar=gg_step2$pheno$sex, intcovar=ic , weights=NULL,
             use="complete.obs",
             incl.markers=FALSE, clean.output=T,
             clean.nmar=1, clean.distance=30,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=12)


summary(sc2_cov, what= "add", perms=sc2_normal_imp_perms, alphas=c(9.1, 7.1, 6.3, 6.3, 3.3), lodcolumn=5, pvalues=F, allpairs=T)

summary(sc2_cov, perms=sc2_normal_imp_perms, thresholds=c(9.1, 7.1, 6.3, 6.3, 3.3), lodcolumn=5)


summary(sc2_cov,what= "int", perms=sc2_normal_imp_perms, lodcolumn=5, pvalues=F, allpairs=T)
summary(sc2_cov,what= "full", perms=sc2_normal_imp_perms, lodcolumn=5, pvalues=F, allpairs=T,alphas=c(9.1, 7.1, 6.3, 6.3, 3.3))
summary(sc2_cov, what="best",perms=sc2_normal_imp_perms, lodcolumn=5, pvalues=F, allpairs=T)


scan.norm.mr <- scanone(gg_step2, method = "mr", model = "normal", pheno.col = 5)

ic <- pull.geno(gg_step2)[,'2:35718468']

c4 :c7    52   32    11.57   9.835   10.30   1.271 -0.46028
c13:c19    8   24    12.88  11.029   10.32   2.560  0.71275
