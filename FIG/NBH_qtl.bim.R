#!/bin/R
debug.cross <- T
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')
library('qtlbim')
clean <- qtl::clean
genotyped.only <- T

popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NBH/REMAPS"
cross.18 <- read.cross(format = "csv", dir = popdir, file = paste(outname, ".BACKUP.QTL_chr.QTLmap.csv",
  sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", "B"))
sex <- read.table(file = file.path(dirso, "sex.txt"))
rownames(sex) <- sex$ID
cross.18$pheno$sex <- sex[as.character(cross.18$pheno$ID), 2]
cross.18$pheno$binary <- as.numeric(cross.18$pheno$pheno >= 3)

pheno <- read.csv('~/QTL_Map_Raw/popgen/rQTL/data/PhenoDist.csv')
rownames(pheno) <- paste(pheno$pop_all,pheno$Sample,sep='_')
cross.18$pheno$gt <- pheno[as.character(cross.18$pheno$ID),6]

#### IF GENOTYPED IND ONLY
if (genotyped.only == T) cross.18 <- subset(cross.18, ind = cross.18$pheno$gt ==
  'GT')
#### Pheno to GT/NG
## Remove problematic individuals (found by kinship analysis)
con <- file(file.path(popdir, "kinship.keep.ind.txt"), open = "r")
keepers <- readLines(con)
close(con)

print("Dropping kinship outliers")
cross.18 <- subset(cross.18, ind = cross.18$pheno$ID %in% keepers)
cross.18 <- subset(cross.18, ind = !is.na(cross.18$pheno$phen))

# rqtl binary scan for prior
cross.18 <- sim.geno(cross.18, error.prob = 0.025, step = 5, n.draws = 500)
scan.norm.imp <- scanone(cross.18, model = "normal", pheno.col = 1, method = "imp",
  addcovar = cross.18$pheno$sex)
perms.norm.imp <- scanone(cross.18, method = "imp", model = "normal", n.perm = 500,
  pheno.col = 1, perm.strata = as.character(cross.18$pheno$gt))
norm.qtl <- summary(scan.norm.imp, perms = perms.norm.imp, alpha = 0.05)
qtl.uns <- makeqtl(cross.18, chr = norm.qtl$chr, pos = norm.qtl$pos)

## n.cluster = slurmcore removed. Farm not working
cross.18 <- calc.genoprob(cross.18,error.prob=0.025)
scan.bin.mr <- scanone(cross.18, method = "mr", model = "binary", pheno.col = 5)
perms.norm.mr <- scanone(cross.18, method = "mr", model = "binary", n.perm = 5000,
  perm.strata = cross.18$pheno$gt, pheno.col = 5)
bin.qtl <- summary(scan.bin.mr, perms = perms.norm.mr, alpha = 0.05)
qtl.mr <- makeqtl(cross.18, chr = bin.qtl$chr, pos = bin.qtl$pos,what="prob")
## n.cluster = slurmcore removed. Farm not working


# full <- stepwiseqtl(cross.18, additive.only = T, method = 'imp', pheno.col = 1,
# scan.pairs = T)
fit <- fitqtl(cross.18, pheno.col = 1, qtl.uns, method = "imp", model = "normal",
  dropone = TRUE, get.ests = TRUE, run.checks = TRUE, tol = 1e-04, maxit = 10000)

fit.mr <- fitqtl(cross.18, pheno.col = 5, qtl.mr, method = "hk", model = "binary",
  dropone = TRUE, get.ests = TRUE, run.checks = TRUE, tol = 1e-04, maxit = 10000)

qtl <- find.marker(cross.18, qtl.uns$chr, qtl.uns$pos)
qtl.mr <- find.marker(cross.18, qtl.mr$chr, qtl.mr$pos)

### Reg and interval mapping

png("/home/jmiller1/public_html/nbh_bin_regression.png", width = 2000)
plot(scan.bin.mr, bandcol = "gray70", cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2)
abline(h = summary(perms.norm.mr)[1, 1])
dev.off()

png("/home/jmiller1/public_html/nbh_norm_imp.png", width = 2000)
plot(scan.norm.imp, bandcol = "gray70", cex.lab = 2, cex.axis = 2, cex.main = 2,
  cex.sub = 2)
abline(h = summary(perms.norm.imp)[1, 1])
dev.off()

### Effect plots

png(paste("/home/jmiller1/public_html/1_2_nbh_effectplot", pop, qtl[1], "_", qtl[2],
  "_pxg.png"))
effectplot(cross.18, mname1 = qtl[1], mname2 = qtl[2])
dev.off()

png(paste("/home/jmiller1/public_html/3_4_nbh_effectplot", pop, qtl[3], "_", qtl[4],
  "_pxg.png"))
effectplot(cross.18, mname1 = qtl[3], mname2 = qtl[4])
dev.off()

png(paste("/home/jmiller1/public_html/1_4_nbh_effectplot", pop, qtl[1], "_", qtl[4],
  "_pxg.png"))
effectplot(cross.18, mname1 = qtl[1], mname2 = qtl[4])
dev.off()

png(paste("/home/jmiller1/public_html/1_3_nbh_effectplot", pop, qtl[1], "_", qtl[3],
  "_pxg.png"))
effectplot(cross.18, mname1 = qtl[1], mname2 = qtl[3], var.flag = "pooled", main = "Genotype interaction \nChrs 2 at 27MB (AIP) and 18 (AHRb) at 20MB")
dev.off()


for (i in 1:length(qtl.mr)) {
  png(paste("/home/jmiller1/public_html/nbh_pxg", pop, qtl.mr[i], "_pxg.png"))
  plotPXG(cross.18, qtl.mr[i], pheno.col = 1, jitter = 1.5, infer = F, pch = 19, main = paste(qtl.mr[i],
    pop))
  dev.off()
}

### qtlbim
crOb <- cross.18
#crOb <- switchAlleles(crOb, markers = swits)
crOb <- qb.genoprob(crOb, step = 10)
###########
qbData.b <- qb.data(crOb, pheno.col = 5, trait = "binary",rancov = 2)
qbModel.b <- qb.model(crOb, epistasis = T, main.nqtl = 5, mean.nqtl = 6, depen = FALSE,
  max.qtl = 8,interval=rep(5,24), chr.nqtl = rep(2,nchr(crOb)))
mc.b <- qb.mcmc(crOb, qbData.b, qbModel.b, pheno.col = 5, n.iter = 300000,genoupdate=FALSE)
so.b <- qb.scanone(mc.b, epistasis = T, type.scan = "heritability", chr = 1:24)
so.b.bf <- qb.scanone(mc.b, epistasis = T, type.scan = "2logBF", chr = 1:24)
best.b <- qb.BayesFactor.jm(mc.b,items = c("pattern", "nqtl"))
two.b <- qb.scantwo(mc.b, chr = c(1,3,2,7,8,9,13,18,20,22,23))
close.b <- qb.close(mc.b)

#######
qbData <- qb.data(crOb, pheno.col = 1, trait = "ordinal", rancov = 2)
mc <- qb.mcmc(crOb, qbData, qbModel, pheno.col = 1, n.iter = 30000)
qbModel <- qb.model(crOb, epistasis = T, main.nqtl = 3, mean.nqtl = 3, depen = FALSE,
  max.qtl = 0, interval = 20)
# chr.nqtl = rep.int(3, 24))
so <- qb.scanone(mc.b, epistasis = T, type.scan = "heritability", chr = 1:24)
best <- qb.BayesFactor.jm(mc.b, items = c("pattern", "nqtl"))
two <- qb.scantwo(mc.b, scan, type.scan = "2logBF")


scan <- list(upper = "main", lower = "epistasis")
st <- qb.scantwo(mc, scan, type.scan = "nqtl", chr = c(1, 2, 6, 8, 13, 18, 20, 23,
  24), epistasis = TRUE)
two <- qb.scantwo(mc, scan, type.scan = "2logBF")
slice <- qb.sliceone(mc, type = "cellmean", chr = c(2, 8, 18, 19, 24))
so <- qb.scanone(mc, epistasis = F, type.scan = "heritability", chr = 1:24)
best <- qb.BayesFactor.jm(mc.b, items = c("pattern", "nqtl"))

#save.image('/home/jmiller1/public_html/NBH.bim.Rsave')

#### FOR FIGURE IN MS

a <- find.marker(cross.18,2,180 )
b <- find.marker(cross.18, 18, 110)
png("/home/jmiller1/public_html/NBH_2_18.png")
effectplot(crOb,pheno.col=1, mname1 = a, mname2 = b,ylim=c(0,5),main = "Genotype interaction \nChrs 2 at 27MB (AIP) and 18 (AHRb) at 20MB")
dev.off()

png("/home/jmiller1/public_html/NBH_2_pxg.png")
plotPXG(cross.18, a, pheno.col = 1, jitter = 1.5, infer = F, pch = 19, main = a)
dev.off()
png("/home/jmiller1/public_html/NBH_18_pxg.png")
plotPXG(cross.18, b, pheno.col = 1, jitter = 1.5, infer = F, pch = 19, main = b)
dev.off()


#slice on 24
temp <- qb.sliceone(mc.b, slice =3,center.type ='mode',sum.scan = "only",type.scan = "cellmean")
png("/home/jmiller1/public_html/slice.png")
plot(temp,chr=c(1,3:17,19:24))
dev.off()

slice <- qb.slicetwo(mc.b,chr = c(2,18),pos=c(180 ,80),type.scan = "LPD",width = 100,smooth=50)
png("/home/jmiller1/public_html/nbh_slice.two_2_18.png", width = 1000)
plot(slice)
dev.off()






png("/home/jmiller1/public_html/nbh.qbBayes.png", width = 600)
plot(best.b)
dev.off()

png("/home/jmiller1/public_html/nbh.scanone.so.png", width = 3000)
plot(so.b, chr = c(1:24), cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5,
  cex.sub = 2.5, xlab = NA)
dev.off()

png("/home/jmiller1/public_html/nbh.scantwo.slice.png", width = 3000)
plot(two.b, chr = c(1:24), slice = 1)
dev.off()

png("/home/jmiller1/public_html/nbh.scantwo.so.slice.png", width = 300)
plot(slice, chr = c(1:24), slice = 2)
dev.off()

png("/home/jmiller1/public_html/nbh.model.png", width = 2000)
# plot(qb.close(mc, target = so))
plot(so.b, c(1:24))
dev.off()

png("/home/jmiller1/public_html/nbh.scanone.st.png", width = 2000, height = 2000)
plot(st, c(1:24))
dev.off()

png("/home/jmiller1/public_html/nbh.coda.png", width = 2000)
plot(qb.coda(mc, variables = c("nqtl")))
dev.off()

png("/home/jmiller1/public_html/nbh_scan_diagnost.so.png", width = 2000)
plot(qb.hpdone(mc.b, chr = c(1:24), level = 0.9, scan = "sum"))
dev.off()

### qtlbim markers qtl.bm <- as.character(summary(qb.hpdone(mc))$chr)
rownames(summary(best)$pattern[2, ])
bimqtl <- summary(qb.scanone(mc, type = "heritability", chr = c(2, 8, 18, 19)))
qtl.bm <- find.marker(cross.18, rownames(bimqtl), bimqtl$pos)

# cross2 <- argmax.geno(crOb, step = 5, off.end = 5, err = 0.1) cross2 <-
# subset(cross2, ind = (!is.na(cross2$pheno$sex)))
cross2 <- convert2cross2(cross.18)

for (i in 1:length(qtl.bm)) {
  chr <- rownames(bimqtl)[i]
  mark <- qtl.bm[i]
  index <- cross2$geno[[as.character(chr)]][, mark] > 0
  png(paste("/home/jmiller1/public_html/", pop, chr, "_interation_plot.png"))
  plot_pxg(geno = cross2$geno[[as.character(chr)]][, mark][index], pheno = cross2$pheno[,
    1][index], sort = F, SEmult = 2)
  dev.off()

  png(paste("/home/jmiller1/public_html/", pop, chr, "_pxg.png"))
  plotPXG(crOb, mark, pheno.col = 1, jitter = 1, infer = F)
  dev.off()

  print(mark)
  print(table(cross2$geno[[as.character(chr)]][, mark]))
}


## Specific to NBH
png(paste("/home/jmiller1/public_html/NBH_pxg", pop, qtl.bm, "_pxg.png", sep = ""))
plotPXG(cross.18, qtl.bm, pheno.col = 1, jitter = 1.5, infer = F, pch = 19, main = paste(qtl.bm,
  pop))
dev.off()

qtl.uns <- makeqtl(cross.18, chr = rownames(bimqtl), pos = bimqtl$pos)

fit.bm <- fitqtl(cross.18, pheno.col = 1, qtl.uns, method = "imp", model = "normal",
  dropone = TRUE, get.ests = TRUE, run.checks = TRUE, tol = 1e-04, maxit = 10000)

capture.output(c(summary(fit), summary(fit.bm), summary(best), summary(so), summary(st)),
  file = "/home/jmiller1/public_html/NBH_genotyped_out.txt")


need fit for mr, imp, binary, bayes (and bayes with binary) all with and without added gts


a <- find.marker(cross.18, 1,20 )
b <- find.marker(cross.18, 13, 105)
png("/home/jmiller1/public_html/nbh_13x1.check.png")
effectplot(crOb, mname1 = a, mname2 = b,ylim=c(0,5))
dev.off()

png("/home/jmiller1/public_html/nbh_13_pxg.png")
plotPXG(cross.18, b, pheno.col = 1, jitter = 1.5, infer = F, pch = 19, main = b)
dev.off()
png("/home/jmiller1/public_html/nbh_1_pxg.png")
plotPXG(cross.18, a, pheno.col = 1, jitter = 1.5, infer = F, pch = 19, main = a)
dev.off()
