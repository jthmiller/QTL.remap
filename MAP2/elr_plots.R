#!/bin/R
### Map QTLs 1 of 3
#debug.cross <- T
#source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library('qtl')

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

qtl <- makeqtl(hyper, chr=18, pos=1587.8,  what="draws")
out.i.18 <- addqtl(hyper, qtl=qtl, formula=y~Q1*Q2, method="imp",model="binary",pheno.col=4)
out.a.18 <- addqtl(hyper, qtl=qtl, formula=y~Q1+Q2, method="imp",model="binary",pheno.col=4)


cbind(summary(scan.bin.em),a=summary(scan.norm.em)[,3],b=summary(scan.bin.mr)[,3],c=summary(scan.np.mr )[,3])
################################################################################




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





png(paste0('~/public_html/NBH_mr.png'))
 effectplot(cross, pheno.col=5,mname1='19:2435634')
dev.off()


png(paste0('~/public_html/NBH_mr_lod.png'))
plotPXG(cross,'19:2435634',infer=F, jitter=2, pch=18)
dev.off()






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
