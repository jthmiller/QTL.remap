#!/bin/R
### Map QTLs 1 of 3
#debug.cross <- T
#source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
#Sys.setenv(TAR = "/bin/tar")
##install_github("kbroman/ASMap")
##install_github("jtlovell/qtlTools")
##install_github("mckaylab/TSPmap")
################################################################################
i <- 5

libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap")
suppressMessages(sapply(libs2load, require, character.only = TRUE))
pop <- 'NBH'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")

################################################################################
## read in the QTL cross
cross <- read.cross.jm(file = file.path(indpops, paste0(pop, ".unphased.f2.csvr")),
format = "csvr", geno = c(1:3), estimate.map = FALSE)
cross <- subset(cross,chr=5)
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
sex <- read.table("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/data/sex.txt",stringsAsFactors=F)
rownames(sex) <- sex$ID
sex.vec <- sex[as.character(cross$pheno$ID), 'sex']
cross$pheno$sex <- sex.vec
cross$pheno$pgm <- 0
## If the male is switched to 'B', then pgm = 0 (or the resistant strain is A)
## 0=female and 1=male
################################################################################

################################################################################
## Switch only in male ##################################################
bfixm <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1M',]
bfixf <- pull.geno(cross)[cross$pheno$ID=='NBH_NBH1F',]

bfix_switm <- names(bfixm)[which(as.numeric(bfixm)==1)]
bfix_switf <- names(bfixf)[which(as.numeric(bfixf)==3)]
bfix <- intersect(bfix_switm,bfix_switf)
cross <- switchAlleles(cross, markers = bfix)

##swit <- c(bfix_switf[!bfix_switf %in% bfix_switm],
fem_gts <- pull.geno(cross)[which(cross$pheno$ID=='NBH_NBH1F'),]
mal_gts <- pull.geno(cross)[which(cross$pheno$ID=='NBH_NBH1M'),]

bb_in_male <- names(which(mal_gts==3))
ab_in_female <- names(which(fem_gts==2))
aa_in_female <- names(which(fem_gts==1))
ab_fem_b_mal <- intersect(bb_in_male, ab_in_female)

goodmarks <- intersect(bb_in_male, aa_in_female)
################################################################################
toss.missing <- c("NBH_5525","NBH_6177")
cross <- subset(cross,ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F','NBH_5528'))
################################################################################

################################################################################
names(cross$geno) <- 'X'
class(cross$geno[['X']]) <- 'X'
################################################################################

cross.bk <- cross

cross.gm <- pull.markers(cross,goodmarks)
cross.gm$pheno$pgm <- 0
cross.gm <- formLinkageGroups(cross.gm, max.rf = 0.25, min.lod = 4, reorgMarkers = TRUE)

cross.gm <- subset(cross.gm, chr=1)
names(cross.gm$geno) <- 'X'
class(cross.gm$geno[['X']]) <- 'X'

cross.gm <- tspOrder(cross = cross.gm, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
png(paste0('~/public_html/NBH_RF_concord_gm.png'))
  plotRF(cross.gm)
dev.off()
################################################################################



################################################################################
reorg <- formLinkageGroups(cross, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
reorg <- subset(reorg , chr=1)
names(reorg$geno) <- 'X'
class(reorg$geno[['X']]) <- 'X'
reorg <- formLinkageGroups(reorg, max.rf = 0.05, min.lod = 10, reorgMarkers = TRUE)

reorg2 <- subset(reorg , chr=2)
reorg3 <- subset(reorg , chr=3)
reorg3 <- subset(reorg , chr=4)

## lod <- pull.rf(reorg, what="lod")
## png(paste0('~/public_html/NBH_RF_concord_try.png'))
##   plot(sort(lod))
## dev.off()

################################################################################
################################################################################
## cross.HOM FEMALES (AA x AB)
## cross.HET MALES (AA x BB)
cross.HET <- subset(cross, ind=cross$pheno$sex==1)
cross.HOM <- subset(cross, ind=cross$pheno$sex==0)
cross.HOM <- formLinkageGroups(cross.HOM, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
cross.HET <- formLinkageGroups(cross.HET, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
save.image('/home/jmiller1/QTL_Map_Raw/ELR_final_map/rsaves/sex_het_hom.map')

## MALES (HET) should be AB (where there is Y, chr2.. and BB elsewhere
##table(pull.geno(cross.HET,chr=3))
##table(pull.geno(cross.HET,chr=8))
hetswits <- markernames(cross.HET,chr=4)
cross.HET.2_3 <- subset(cross.HET, chr=2:7)
cross.HET.2_3  <- switchAlleles(cross.HET.2_3, markers = markernames(cross.HET.2_3,chr=4))
cross.HET.2_3 <- formLinkageGroups(cross.HET.2_3, max.rf = 0.1, min.lod = 6, reorgMarkers = TRUE)
#cross.HETS <- pull.markers(cross.HETS,markernames(cross.HET.2_3,chr=4)
#filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/cross.HET.2_3')
#write.cross(cross,filestem=filename,format="csv")
##cross.HET.2_3  <- switchAlleles(cross.HET.2_3, markers = markernames(cross.HET.2_3,chr=4))
################################################################################
################################################################################

################################################################################
################################################################################
## FEMALES (HOM) should be AB or AA (more AA, bc male non-rec X is all A)
#geno.table(cross.HOM,chr=3)
#table(pull.geno(cross.HOM,chr=4))
homswits <- markernames(cross.HOM,chr=c(5,6,8,9,13))
cross.HOM.3_4 <- subset(cross.HOM, chr=c(3:13))
cross.HOM.3_4  <- switchAlleles(cross.HOM.3_4, markers = homswits)
cross.HOM.3_4 <- formLinkageGroups(cross.HOM.3_4, max.rf = 0.1, min.lod = 6, reorgMarkers = TRUE)
##table(pull.geno(cross.HET.2_3,chr=3))
##table(pull.geno(cross.HOM.3_4,chr=7))
##geno.table(cross.HET,chr=3)
################################################################################
################################################################################

gthom <- geno.table(cross.HOM.3_4)
gthet <- geno.table(cross.HET.2_3)

homs <- rownames(gthom)[gthom$missing < 5]
hets <- rownames(gthet)[gthet$missing < 5]

cross.try <- pull.markers(cross, unique(c(homs,hets)))
cross.try <- switchAlleles(cross.try, c(homswits,hetswits))
cross.try <- formLinkageGroups(cross.try, max.rf = 0.1, min.lod = 10, reorgMarkers = TRUE)

#cross.try.bk <- cross.try
cross.try <- cross.try.bk
cross.try <- switchAlleles(cross.try, markernames(cross.try,chr=3))

cross.try1 <- subset(cross.try,chr=1)
names(cross.try1$geno)[1] <- 'X'
class(cross.try1$geno[['X']]) <- 'X'
cross.try1 <- formLinkageGroups(cross.try1, max.rf = 0.08, min.lod = 4, reorgMarkers = TRUE)

cross.try1 <- tspOrder(cross = cross.try1, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
png(paste0('~/public_html/NBH_RF_concord_try1.png'))
  plotRF(cross.try1)
dev.off()


cross.try2 <- subset(cross.try,chr=2)
names(cross.try2$geno)[1] <- 'X'
class(cross.try2$geno[['X']]) <- 'X'
cross.try2 <- formLinkageGroups(cross.try2, max.rf = 0.1, min.lod = 8, reorgMarkers = TRUE)
cross.try2 <- tspOrder(cross = cross.try2, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
png(paste0('~/public_html/NBH_RF_concord_try2.png'))
  plotRF(cross.try2)
dev.off()

cross.try3 <- subset(cross.try,chr=3)
names(cross.try3$geno)[1] <- 'X'
class(cross.try3$geno[['X']]) <- 'X'
cross.try3 <- formLinkageGroups(cross.try3, max.rf = 0.1, min.lod = 8, reorgMarkers = TRUE)
cross.try3 <- tspOrder(cross = cross.try3, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
png(paste0('~/public_html/NBH_RF_concord_try3.png'))
  plotRF(cross.try3)
dev.off()

cross.try4 <- subset(cross.try,chr=4)
names(cross.try4$geno)[1] <- 'X'
class(cross.try4$geno[['X']]) <- 'X'
cross.try4 <- formLinkageGroups(cross.try4, max.rf = 0.1, min.lod = 8, reorgMarkers = TRUE)
cross.try4 <- tspOrder(cross = cross.try4, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
png(paste0('~/public_html/NBH_RF_concord_try4.png'))
  plotRF(cross.try4)
dev.off()

################################################################################


cross.try <- tspOrder(cross = cross.try, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')



cross.try5 <- subset(cross.try, ind=cross$pheno$sex==1)
cross.try5 <- formLinkageGroups(cross.try5, max.rf = 0.08, min.lod = 4, reorgMarkers = TRUE)




cross.try5 <- subset(cross.try, ind=cross.try$pheno$sex==0)
cross.try5 <- subset(cross.try5,chr=5)
names(cross.try5$geno)[1] <- 'X'
class(cross.try5$geno[['X']]) <- 'X'
cross.try5 <- formLinkageGroups(cross.try5, max.rf = 0.05, min.lod = 5, reorgMarkers = TRUE)
names(cross.try5$geno)[1] <- 'X'
class(cross.try5$geno[['X']]) <- 'X'
cross.try5 <- tspOrder(cross = cross.try5, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
png(paste0('~/public_html/NBH_RF_concord_try5.png'))
  plotRF(cross.try5)
dev.off()





################################################################################

names(cross.try2$geno) <- 'X'
class(cross.try$geno[['X']]) <- 'X'
cross.try <- tspOrder(cross = cross.try, hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

png(paste0('~/public_html/NBH_RF_concord_try.png'))
  plotRF(cross.try)
dev.off()

##geno.table(cross.try,chr=3)

#cross.HOM <- subset(cross.HOM,chr=1)
#cross.HET <- subset(cross.HET,chr=1)
################################################################################
chr1, (AA or AB in F, AYxBY in m)
chr2, (??x?? in F, AYxBY in m)
chr3, (AA    in F, BY in m)







gt <- geno.table(cross)
png(paste0('~/public_html/NBH_nbh_sex_pval.png'))
 hist(-log10(gt$P.value),breaks=50)
dev.off()

keep <- rownames(gt)[-log10(gt$P.value) < 3]
try <- pull.markers(cross,keep)

try <- formLinkageGroups(try, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
try <- subset(try,chr=1)

ord <- order(as.numeric(gsub(".*:","",names(pull.map(try)[[as.character(1)]]))))
try <- switch.order(try, chr = 1, ord, error.prob = 0.01, map.function = "kosambi",
 maxit = 1, tol = 0.1, sex.sp = F)
names(try$geno) <- 'X'
class(try$geno[['X']]) <- 'X'
################################################################################

ret <- removeDoubleXO(try, chr='X')

try <- tspOrder(cross = try ,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

png(paste0('~/public_html/NBH_RF_concord_try.png'))
  plotRF(try)
dev.off()

################################################################################
cross.HET <- subset(cross, ind=cross$pheno$sex==1)
cross.HOM <- subset(cross, ind=cross$pheno$sex==0)
################################################################################

################################################################################
cross <- cleanGeno_jm(cross, chr=i, maxdist=100, maxmark=8, verbose=TRUE)
cross <- calc.errorlod(cross, err=0.05)
cross <- removeDoubleXO(cross, chr=i)
cross <- calc.errorlod(cross, err=0.05)
cross <- cleanGeno_jm_2(cross, chr=i, maxdist=50, maxmark=4, verbose=TRUE)
cross <- calc.errorlod(cross, err=0.05)
################################################################################

################################################################################
cross <-tspOrder(cross = cross,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
cross <- shiftmap(cross, offset=0)
cross_map <-  est.map(cross, error.prob=0.04,map.function="kosambi",maxit=100000,tol=1e-7, sex.sp=FALSE, verbose=FALSE)
cross <- qtl:::replace.map(cross,cross_map)

filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/NBH_all_mark_',i,'_tsp')
write.cross(cross,chr=i,filestem=filename,format="csv")

png(paste0('~/public_html/NBH_RF_concord',i,'_tsp.png'))
  plotRF(cross)
dev.off()
################################################################################

















mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

filename <- paste0('NBH_all_mark_HETS',i,'_tsp.csv')
reorg.HET <- read.cross.jm(file=filename,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

filename <- paste0('NBH_all_mark_HOM',i,'_tsp.csv')
reorg.HOM <- read.cross.jm(file=filename,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

### CHR ########################################################################
gt.het <- geno.table(reorg.HET)
#  10   11   12    9    8    7    6    5    4    3    2    1
#  19   19   19   20   22   64  120  140  269  338 1258 8343
##gt.het[intersect(rownames(gt.het),goodmarks),]
ht.gm <- gt.het[intersect(rownames(gt.het),goodmarks),]
ord.het <- names(sort(table(ht.gm$chr)))
get.het <- table(ht.gm$chr)[ord.het]/table(gt.het$chr)[ord.het]
view_hets <- cbind(table(ht.gm$chr)[ord.het], table(gt.het$chr)[ord.het])
##CHR2 and 5 seems to have highest proportion

### In fixed markers, chr1 is BB, there are no chr2 markers, and chr3 is AA x BB
markers.HET <- markernames(reorg.HET, chr=c(2,5))
################################################################################

################################################################################
gt.hom <- geno.table(reorg.HOM)
#   7    8    6    5    4    3    2    1
#   5    5   25   62  400 1438 1913 9823
hm.gm <- gt.hom[intersect(rownames(gt.hom),goodmarks),]
ord.hom <- names(sort(table(hm.gm$chr)))
get.hom <- table(hm.gm$chr)[ord.hom]/table(gt.hom$chr)[ord.hom]
view_homs <- cbind(table(hm.gm$chr)[ord.hom], table(gt.hom$chr)[ord.hom])

### CHR2 and 3
## In fixed markers, cross_HOM chr1 is AA, chr2 is very few BB, and chr3 is many AAxBB
markers.HOM <- markernames(reorg.HOM, chr=c(2,3))
################################################################################

################################################################################
reorg.HET <- pull.markers(cross.HET, markers.HET)
names(reorg.HET$geno) <- 'X'
class(reorg.HET$geno[['X']]) <- 'X'
reorg.HET <- formLinkageGroups(reorg.HET, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)



reorg.HOM <- pull.markers(cross.HOM, markers.HOM)
names(reorg.HOM$geno) <- 'X'
class(reorg.HOM$geno[['X']]) <- 'X'
reorg.HOM <- formLinkageGroups(reorg.HOM, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)

reorg.all <- pull.markers(cross, unique(c(markers.HET,markers.HOM)))
names(reorg.all$geno) <- 'X'
class(reorg.all$geno[['X']]) <- 'X'
################################################################################
reorg.all.flg.bk <- reorg.all.flg
reorg.all.flg <- formLinkageGroups(reorg.all, max.rf = 0.05, min.lod = 10, reorgMarkers = TRUE)

reorg.all.flg <- tspOrder(cross = reorg.all.flg ,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

reorg.all.flg.hi <- formLinkageGroups(subset(reorg.all.flg, chr=c(1:4)), max.rf = 0.15, min.lod = 8, reorgMarkers = TRUE)
reorg.all.flg.hi <- tspOrder(cross = reorg.all.flg.hi ,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

png(paste0('~/public_html/nbh_HOM_5_flg_tsp.png'))
  plotRF(reorg.all.flg.hi)
dev.off()


png(paste0('~/public_html/nbh_HOM_5_flg_tsp.png'))
  plotRF(reorg.all.flg,chr=2)
dev.off()

################################################################################

cross.all <- formLinkageGroups(cross, max.rf = 0.15, min.lod = 8, reorgMarkers = TRUE)
cross.all.bk <- cross.all
cross.all <- switchAlleles(cross.all, markers = markernames(cross.all,chr=2))
cross.all <- formLinkageGroups(cross.all, max.rf = 0.15, min.lod = 8, reorgMarkers = TRUE)

save.image('~/home/jmiller1/QTL_Map_Raw/ELR_final_map/rsaves/sex.map')

################################################################################
reorg <- tspOrder(cross = reorg.all ,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
png(paste0('~/public_html/nbh_HOM_5_COMBINED_tsp.png'))
  plotRF(reorg)
dev.off()

reorg.hm <- tspOrder(cross = reorg.HOM ,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
png(paste0('~/public_html/nbh_HOM_5_HETs_tsp.png'))
  plotRF(reorg.hm)
dev.off()

reorg.ht <- tspOrder(cross = reorg.HET ,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
png(paste0('~/public_html/nbh_HOM_5_HOMs_tsp.png'))
  plotRF(reorg.ht)
dev.off()
################################################################################
paternal grandmother (Resistant founder). If AA, then  0.



you could code females as AA and AB and males as AA and BB,
in which case this needs to be indicated through the genotypes argument in read.cross.

In my case, HOM is AA/BB and HET is

There are 200 F2 individuals typed at 94 markers, including 3 on the X chromosome. There is one quantitative phenotype,
along with an indication of sex (0=female, 1=male) and the direction of the cross (pgm = paternal grandmother,
0=A, meaning the cross was (AxB)x(AxB), and
1=B, meaning the cross was (AxB)x(BxA)).
Note that the X chromosome genotypes are coded in a special way (see read.cross).

for individuals with pgm=0, sex=0, 1=AA and 2=AB;
for individuals with pgm=0, sex=1, 1=A and 2=B (hemizygous);

for individuals with pgm=1, sex=0, 1=BB and 2=AB;
for individuals with pgm=1, sex=1, 1=A and 2=B.

## males
for individuals with pgm=0, sex=1, 1=A and 2=B (hemizygous);
for individuals with pgm=1, sex=1, 1=A and 2=B.

for individuals with pgm=0, sex=0, 1=AA and 2=AB;
for individuals with pgm=1, sex=0, 1=BB and 2=AB;

## ################################################################################
## gt <- geno.table(reorg.HOM)
## het <- rownames(gt)[which(gt$AB==1)]
## sort(rowSums(pull.geno(reorg.HOM)==2,na.rm=T))
## g <- pull.geno(reorg.HOM)
## gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
## gfreq <- t(t(gfreq) / colSums(gfreq))
## png(paste0('~/public_html/nbh_trsh.png'))
## plot(sort(gfreq[2,]))
## dev.off()
## ################################################################################

################################################################################
### Switch phase and keep only parent conf markers #############################
### ENRICH FOR AAxBB ##########################################################

################################################################################
#### Pvalue and Missing ##############################################
gt <- geno.table(subset(cross, ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F')))
bfixA <- rownames(gt[which(gt$P.value > 0.0001 & gt$missing < 5),])
################################################################################
cros.bk <- cross
###### FILTER #######################################################
cross <- pull.markers(cross,bfixA)
cross <- subset(cross,ind=!cross$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))
################################################################################


mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'
fl <- file.path(mpath,'NBH_unmapped_filtered')
write.cross(cross,filestem=fl,format="csv")

###
####
####

################################################################################

ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross)[[as.character(i)]]))))
cross <- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
 maxit = 1, tol = 0.1, sex.sp = F)

png(paste0('~/public_html/NBH_gts_preclean',i,'.png'),height=2500,width=4500)
 plotGeno(cross, chr=i, cex=2)
dev.off()

################################################################################
### SEX SPEC MAPS

if (i == 5){

cross.m <- subset(cross, ind=cross$pheno$sex==1)
cross.f <- subset(cross, ind=cross$pheno$sex==0)

reorg.m <- drop.markers(cross.m,names(which(colSums(is.na(pull.geno(cross.m))) > 1)))
reorg.m <- formLinkageGroups(reorg.m , max.rf = 0.05, min.lod = 15, reorgMarkers = TRUE)
reorg.m  <-tspOrder(cross = subset(reorg.m ,chr=1) ,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
png(paste0('~/public_html/nbh_M',i,'_tsp.png'))
  plotRF(reorg.f)
dev.off()

reorg.f <- drop.markers(cross.f,names(which(colSums(is.na(pull.geno(cross.f))) > 1)))
reorg.f <- formLinkageGroups(subset(cross.f,chr=5), max.rf = 0.5, min.lod = 15, reorgMarkers = TRUE)
reorg.f <-tspOrder(cross = subset(reorg.f,chr=1),hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')
png(paste0('~/public_html/nbh_F',i,'_tsp.png'))
  plotRF(reorg.f)
dev.off()

}
cross <- cross.m
################################################################################
cross <- cleanGeno_jm(cross, chr=i, maxdist=100, maxmark=8, verbose=TRUE)
cross <- calc.errorlod(cross, err=0.05)
cross <- removeDoubleXO(cross, chr=i)
cross <- calc.errorlod(cross, err=0.05)
cross <- cleanGeno_jm_2(cross, chr=i, maxdist=50, maxmark=4, verbose=TRUE)
cross <- calc.errorlod(cross, err=0.05)
################################################################################

cross <-tspOrder(cross = cross,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

cross <- shiftmap(cross, offset=0)

cross_map <-  est.map(cross, error.prob=0.04,map.function="kosambi",maxit=100000,tol=1e-7, sex.sp=FALSE, verbose=FALSE)

cross <- qtl:::replace.map(cross,cross_map)

filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/NBH_all_mark_',i,'_tsp')
write.cross(cross,chr=i,filestem=filename,format="csv")

png(paste0('~/public_html/NBH_RF_concord',i,'_tsp.png'))
  plotRF(cross)
dev.off()


g <- pull.geno(cross.f)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
