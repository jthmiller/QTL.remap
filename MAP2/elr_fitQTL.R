#!/bin/R
### Map QTLs 1 of 3
#debug.cross <- T
#source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library('qtl')

################################################################################
## unmapped
################################################################################
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'
fl <- file.path(mpath,'ELR_unmapped_filtered_added_markers.csv')

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
 estimate.map = FALSE
)
################################################################################

cross$pheno$bin <- ifelse(cross$pheno$Pheno > 2, 1 , 0)
cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno))

## rf.1 <- est.rf(subset(cross,chr=1), maxit=1000, tol=1e-6)
## rf.2 <- est.rf(subset(cross,chr=2), maxit=1000, tol=1e-6)
## rf.8 <- est.rf(subset(cross,chr=8), maxit=1000, tol=1e-6)
## rf.18 <- est.rf(subset(cross,chr=18), maxit=1000, tol=1e-6)
## rf.13 <- est.rf(subset(cross,chr=13), maxit=1000, tol=1e-6)
## save.image(file.path(mpath,'ER_chrFR.rsave'))

##################################################################################
#### read in the QTL cross
##################################################################################
##mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'
##fl <- file.path(mpath,'ELR_unmapped_filtered.csv')
##
##cross <- read.cross(
## file = fl,
## format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
## estimate.map = FALSE
##)
##################################################################################

## UNFILTERED
##################################################################################
#### read in the QTL cross
## cross <- read.cross.jm(file = file.path(indpops, paste(pop, ".unphased.f2.csvr",
##  sep = "")), format = "csvr", geno = c(1:3), estimate.map = FALSE)
##################################################################################
##################################################################################
##### Pull names from plinkfile
##path <- file.path(indpops, paste(pop, ".ped", sep = ""))
##popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
##indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
##cross$pheno$ID <- paste(popname, indname, sep = "_")
##################################################################################
##
###### PHENO #####################################################################
##cross$pheno$bin <- ifelse(cross$pheno$Pheno > 2, 1 , 0)
##cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno))
##################################################################################
##
## #### FITER ########################################################################
## marks <- read.table('/home/jmiller1/QTL_Map_Raw/ELR_final_map/goodmarks.rtable',stringsAsFactors=F)
## inds <- read.table('/home/jmiller1/QTL_Map_Raw/ELR_final_map/goodsamps.rtable',stringsAsFactors=F)
## ####################################################################################
## drops <- markernames(cross)[!markernames(cross) %in%  marks[,1]]
## cross <- drop.markers(cross,drops)
## cross <- subset(cross, ind= !is.na(cross$pheno$Pheno))
##
## ################################################################################
dups <- findDupMarkers(cross, exact.only = F, adjacent.only = F)
cross <- drop.markers(cross, unlist(dups))
##
cross <- subset(cross, ind=nmissing(cross)<5000)
cross <- calc.genoprob(cross,error.prob=0.05)
cross <- sim.geno(cross,error.prob=0.05)
################################################################################
perms.bin.em <- scanone(cross, method = "em", model = "binary", maxit = 100,
  n.perm = 5, pheno.col = 4, n.cluster = 2)

scan.bin.em <- scanone(cross, method = "em", model = "binary", pheno.col = 4)
scan.norm.em <- scanone(cross, method = "em", model = "normal", pheno.col = 5)
scan.bin.mr <- scanone(cross, method = "mr", model = "binary", pheno.col = 4)
scan.norm.mr <- scanone(cross, method = "mr", model = "normal", pheno.col = 5)
scan.norm.imp <- scanone(cross, method = "imp", model = "normal", pheno.col = 5)
scan.bin.imp <- scanone(cross, method = "imp", model = "binary", pheno.col = 4)


cbind(summary(scan.bin.em),a=summary(scan.norm.em)[,3],b=summary(scan.bin.mr)[,3],c=summary(scan.np.mr )[,3])
################################################################################

hyper <- sim.geno(cross.np)

qtl <- makeqtl(hyper, chr=18, pos=1587.8,  what="draws")

out.i.18 <- addqtl(hyper, qtl=qtl, formula=y~Q1*Q2, method="imp",model="binary",pheno.col=4)
out.a.18 <- addqtl(hyper, qtl=qtl, formula=y~Q1+Q2, method="imp",model="binary",pheno.col=4)

png(paste0('~/public_html/ELR_multi_imputation_interaction_lod.png'))
plotPXG(cross,c('18:20273448','15:2483654') ,1,infer=F)
dev.off()


png(paste0('~/public_html/ELR_multi_imputation_interaction_lod.png'))
plotPXG(cross,c('18:20273448','AHR2a_del') ,1,infer=F, jitter=2, pch=18)
dev.off()

png(paste0('~/public_html/ELR_multi_imputation_interaction_lod.png'))
 effectplot(cross,pheno.col=1,mname1='AHR2a_del',mname2='18:20273448')
dev.off()

png(paste0('~/public_html/ELR_multi_interaction_lod.png'))
 effectplot(cross, pheno.col=1,mname1='15:2483654', mname2='18:20273448')
dev.off()

png(paste0('~/public_html/ELR_multi_interaction_lod.png'))
 effectplot(cross, pheno.col=4,mname1='9:26045359', mname2='18:20273448')
dev.off()

geno.crosstab(cross,"18:20273448",'15:2483654')

png(paste0('~/public_html/ELR_multi_interaction_lod.png'),width=2000)
effectscan(hyper,pheno.col=5)
dev.off()


png(paste0('~/public_html/ELR_multi_interaction_lod.png'),width=2000)
plot(out.a.18)
dev.off()

out.i.18 <- addqtl(hyper, qtl=qtl, formula=y~Q1*Q2, method="imp",model="normal",pheno.col=5)
out.a.18 <- addqtl(hyper, qtl=qtl, formula=y~Q1+Q2, method="imp",model="normal",pheno.col=5)




qtl <- makeqtl(hyper, chr=1, pos=7.62e-04,  what="draws")
out.i.1 <- addqtl(hyper, qtl=qtl, formula=y~Q1*Q2, method="imp",model="binary",pheno.col=4)




















geno <- c('AA','AB','BB')
names(geno) <- c(1:3)

GxE <- function(crs,mark,phenocol,ch){
 table(
  pull.pheno(cross)[,phenocol],
  geno[pull.geno(crs,ch)[,mark]]
 )
}

GxE(crs=cross,mark='1:20299408',phenocol=4,ch=1)
GxE(crs=cross,mark="AHR2a_del",phenocol=4,ch=1)
GxE(crs=cross,mark="18:20273448",phenocol=4,ch=18)

geno.crosstab(cross,"18:20273448",'1:20299408')

pull.rf(rf, what="lod", chr=1 )["AHR2a_del",]
sort(rf.1[,"AHR2a_del"])
################################################################################
##for i in ELR*.csv; do
## cut -d',' -f6- $i > ${i%.csv}.nofirst
##done
##
##paste *.nofirst > ELR_marks
##cut -d',' -f1-5 ELR_chr_24.csv > ELR.info
##
##paste -d',' ELR.info ELR_marks > ELR.all.scv
##awk -F',' '{print $0}' ELR.all.scv > ELR.all.csv
################################################################################

dups <- findDupMarkers(cross.final, exact.only = T, adjacent.only = F)
cross.final.nodups <- drop.markers(cross.final, unlist(dups))
cross.final.nodups$pheno$pheno_norm <- signif(cross.final.nodups$pheno$pheno_norm,5)

### SCANS
scan.bin.mr <- scanone(cross.final.nodups, method = "mr", model = "binary", pheno.col = 4)
scan.bin.mr <- scanone(cross.final.nodups, method = "mr", model = "np", pheno.col = 4)



scan.norm.mr <- scanone(cross.final.nodups, method = "mr", model = "normal", pheno.col = 5)
scan.bin.imp <- scanone(cross.final.nodups, method = "imp",model = "normal", pheno.col = 5)
scan.np.em <- scanone(cross.final.nodups, method = "em", model = "np", pheno.col = 1, maxit = 5000)
scan.bin.em <- scanone(cross.final.nodups, method = "em", model = "binary", pheno.col = 4)


scan.norm.ehk <- scanone(cross.final.nodups, method = "ehk", model = "normal", maxit = 5000, pheno.col = 5)


scan.2p.imp <- scanone(cross.final.nodups, method = "em", model = "2part", pheno.col = 1)

perms.bin.em <- scanone(cross.final.nodups, method = "em", model = "binary", maxit = 5000,
  n.perm = 500, pheno.col = 4, n.cluster = 4)
summary(perms.bin.em)


save.image('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_permutations.rsave')


chr18 <- orderMarkers(cross.final.nodups,chr=18, use.ripple=F,maxit=100,sex.sp=F)
scan.np.em.18 <- scanone(chr18, method = "em", model = "np", pheno.col = 1, maxit = 5000)
pull.map(chr18, chr=18)


1:26531574

hyper <-jittermap(chr18, amount=1e-6)
hyper <- calc.genoprob(hyper, step=2.5, err=0.01)
hyper <- sim.geno(hyper, step=2.5, n.draws=25, err=0.01)

qtl <- makeqtl(hyper, chr=18, pos=608,  what="draws")
out.i.18 <- addqtl(hyper, qtl=qtl, formula=y~Q1*Q2, method="imp",model="binary",pheno.col=4)

qtl <- makeqtl(hyper, chr=1, pos=7.62e-04,  what="draws")
out.i.1 <- addqtl(hyper, qtl=qtl, formula=y~Q1*Q2, method="imp",model="binary",pheno.col=4)


int_sum <- summary(out.i)
int_sum[order(int_sum$lod),]

out.a <- addqtl(hyper, qtl=qtl, formula=y~Q1+Q2, method="imp",model="binary")
add_sum <- summary(out.i)
add_sum[order(add_sum$lod),]

png(paste0('~/public_html/ELR_multi_imputation_interaction_lod.png'))
 plot(out.i - out.a)
dev.off()


scan_18 <- scanqtl(cross, pheno.col=4, chr=18, pos=608, covar=NULL, formula=y~Q1,
            method="imp", model="binary",
            incl.markers=T, verbose=TRUE, tol=1e-4, maxit=100,
            forceXcovar=FALSE)





scan.norm.em <- scanone(cross.18, method = "em", model = "normal", maxit = 5000,
  pheno.col = 5)
### Normal scan on transformed phenotype w/Extended haley knott (better for
### selective/missing genos at non-gt'ed ind)
scan.norm.ehk <- scanone(cross.final.nodups, method = "ehk", model = "normal", maxit = 5000, pheno.col = 5)

### Normal scan on transformed phenotype fast haley knott (not robust to missing
### data. LOD inflation)
scan.norm.hk <- scanone(cross.18, method = "hk", model = "normal", maxit = 5000,
  pheno.col = 6)



scan.bin.mr <- scanone(cross.final.nodups, method = "mr", model = "binary", pheno.col = 4)
mr_sum <- summary(scan.bin.mr)
mr_sum[order(mr_sum$lod),]

geno.crosstab(cross.final.nodups,'18:20367780','1:26531574')


mr_sum[order(mr_sum$lod),]
png(paste0('~/public_html/18_1_pxg.png'))
plotPXG(cross.final.nodups, c('18:20367780','1:26531574'))
dev.off()

png(paste0('~/public_html/1pxg.png'))
 plotPXG(cross.final.nodups,pheno.col=1,'1:26531574')
dev.off()

png(paste0('~/public_html/18pxg.png'))
 plotPXG(cross.final.nodups,pheno.col=1,'18:20367780')
dev.off()

lod <- subset(scan.bin.mr, chr=13,lodcolumn=1)
xf <- as.numeric(gsub(".*:",'',rownames(lod)))

png(paste0('~/public_html/chr1_lod.png'))
 plot(xf,lod$lod)
dev.off()



### MAPPED
hyper <- sim.geno(cross.18.nodup, step=1, n.draws=256, err=0.01)
qtl <- makeqtl(hyper, chr=14, pos=106.9, what="draws")

out.i <- addqtl(hyper, qtl=qtl, formula=y~Q1*Q2, method="imp")
out.a <- addqtl(hyper, qtl=qtl, formula=y~Q1+Q2, method="imp")

png(paste0('~/public_html/one_v_two_modle_qtl.png'),width=1000)
plot(out.i - out.a)
dev.off()

mname1 <- find.marker(cross.18.nodup, 1, 1)
mname2 <- find.marker(cross.18.nodup, 18, 152.427)
mname2 <- find.marker(cross.18.nodup, 14, 106.9) # marker D13Mit147










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
