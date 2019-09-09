#!/bin/R
### Map QTLs 1 of 3
debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")

library('qtl')
################################################################################
## read in the QTL cross
################################################################################
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'
fl <- file.path(mpath,'ELR_unmapped_filtered.csv')

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
 estimate.map = FALSE
)
################################################################################

################################################################################
dups <- findDupMarkers(cross, exact.only = T, adjacent.only = F)
cross <- drop.markers(cross, unlist(dups))

################################################################################
perms.bin.em <- scanone(cross.r, method = "em", model = "binary", maxit = 1000,
  n.perm = 500, pheno.col = 4, n.cluster = 1)


################################################################################


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


scan.bin.mr.18 <- scanone(cross.18, method = "mr", model = "binary", pheno.col = 4)
scan.np.em.18 <- scanone(cross.18, method = "em", model = "np", pheno.col = 3, maxit = 5000)

perms.bin.em <- scanone(cross.r, method = "em", model = "binary", maxit = 1000,
  n.perm = 500, pheno.col = 4, n.cluster = 1)

dups <- findDupMarkers(cross.final, exact.only = T, adjacent.only = F)
cross.final.nodups <- drop.markers(cross.final, unlist(dups))
cross.final.nodups$pheno$pheno_norm <- signif(cross.final.nodups$pheno$pheno_norm,5)

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
