#!/bin/R
debug.cross <- T
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')
library('qtlbim')
clean <- qtl::clean
genotyped.only <- T

popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/ELR/REMAPS"
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


##flips <- which(!before.flip == after.flip)
flips <- c(1,2,3,5,8,9,12,15,17,20,21,22)
mak <- markernames(cross.18, chr =flips)
cross.18 <- switchAlleles(cross.18, markers = mak)



# rqtl binary scan for prior
cross.18 <- sim.geno(cross.18, error.prob = 0.025, step = 5, n.draws = 500)
scan.norm.imp <- scanone(cross.18, model = "normal", pheno.col = 1, method = "imp",
  addcovar = cross.18$pheno$sex)
perms.norm.imp <- scanone(cross.18, method = "imp", model = "normal", n.perm = 500,
  pheno.col = 1, perm.strata = as.character(cross.18$pheno$gt))
norm.qtl <- summary(scan.norm.imp, perms = perms.norm.imp, alpha = 0.05)
qtl.uns <- makeqtl(cross.18, chr = norm.qtl$chr, pos = norm.qtl$pos)
fit <- fitqtl(cross.18, pheno.col = 1, qtl.uns, method = "imp", model = "normal",
  dropone = TRUE, get.ests = TRUE, run.checks = TRUE, tol = 1e-04, maxit = 10000)
qtl <- find.marker(cross.18, qtl.uns$chr, qtl.uns$pos)

## n.cluster = slurmcore removed. Farm not working
cross.18 <- calc.genoprob(cross.18,error.prob=0.025)
scan.bin.mr <- scanone(cross.18, method = "mr", model = "binary", pheno.col = 5)
perms.norm.mr <- scanone(cross.18, method = "mr", model = "binary", n.perm = 5000,
  perm.strata = cross.18$pheno$gt, pheno.col = 5)
bin.qtl <- summary(scan.bin.mr, perms = perms.norm.mr, alpha = 0.05)
qtl.mr <- makeqtl(cross.18, chr = bin.qtl$chr, pos = bin.qtl$pos,what="prob")
fit.mr <- fitqtl(cross.18, pheno.col = 5, qtl.mr, method = "hk", model = "binary",
  dropone = TRUE, get.ests = TRUE, run.checks = TRUE, tol = 1e-04, maxit = 10000)
qtl.mr <- find.marker(cross.18, qtl.mr$chr, qtl.mr$pos)

## n.cluster = slurmcore removed. Farm not working
scan.bin.np <- scanone(cross.18, model = "np", pheno.col = 1)
perms.norm.np <- scanone(cross.18, model = "np", n.perm = 5000,
  perm.strata = cross.18$pheno$gt, pheno.col = 1)
nor.qtl <- summary(scan.bin.np, perms = perms.norm.np, alpha = 0.1)

qtl.np <- makeqtl(cross.18, chr = nor.qtl$chr, pos = nor.qtl$pos,what="prob")

# full <- stepwiseqtl(cross.18, additive.only = T, method = 'imp', pheno.col = 1,
# scan.pairs = T)

### Reg and interval mapping plots
png("/home/jmiller1/public_html/elr_bin_regression.png", width = 2000)
plot(scan.bin.mr, bandcol = "gray70", cex.lab = 2, cex.axis = 2, cex.main = 2, cex.sub = 2)
abline(h = summary(perms.norm.mr)[1, 1])
dev.off()

png("/home/jmiller1/public_html/elr_norm_imp.png", width = 2000)
plot(scan.norm.imp, bandcol = "gray70", cex.lab = 2, cex.axis = 2, cex.main = 2,
  cex.sub = 2)
abline(h = summary(perms.norm.imp)[1, 1])
dev.off()

### Effect plots

png(paste("/home/jmiller1/public_html/elr_effectplot", pop, qtl[1], "_", qtl[2],
  "_pxg.png", sep = ""))
effectplot(cross.18, mname1 = qtl[1], mname2 = qtl[2])
dev.off()

for (i in 1:length(qtl)) {
  png(paste("/home/jmiller1/public_html/elr_pxg", pop, qtl[i], "_pxg.png", sep = ""))
  plotPXG(cross.18, qtl[i], pheno.col = 1, jitter = 1.5, infer = F, pch = 19, main = paste(qtl[i],
    pop))
  dev.off()
}

## Specific to ELR
png(paste("/home/jmiller1/public_html/elr_pxg_mr", pop, qtl.mr[1], "_pxg.png", sep = ""))
plotPXG(cross.18, qtl.mr[1], pheno.col = 1, jitter = 1.5, infer = F, pch = 19, main = paste(qtl[i],
  pop))
dev.off()

### qtlbim
crOb <- cross.18
crOb <- switchAlleles(crOb, markers = swits)

crOb <- qb.genoprob(crOb, step = 10)
###########
qbData.b <- qb.data(crOb, pheno.col = 5, trait = "binary")
qbModel.b <- qb.model(crOb, epistasis = T, main.nqtl = 2, mean.nqtl = 4, depen = FALSE,
  max.qtl = 0)
mc.b <- qb.mcmc(crOb, qbData.b, qbModel.b, pheno.col = 5, n.iter = 30000)
so.b <- qb.scanone(mc.b, epistasis = T, type.scan = "heritability", chr = 1:24)
so.b.bf <- qb.scanone(mc.b, epistasis = T, type.scan = "2logBF", chr = 1:24)
best.b <- qb.BayesFactor.jm(mc.b,items = c("pattern", "nqtl"))
two.b <- qb.scantwo(mc.b, chr = c(1,2,7,8,9,13,18,22,23))
close.b <- qb.close(mc.b)
try <- qb.BestPattern(mc.b, epistasis = TRUE,include ='exact', category =  "nqtl", level = 5)

#######
qbData <- qb.data(crOb, pheno.col = 1, trait = "ordinal", rancov = 2)
mc <- qb.mcmc(crOb, qbData, qbModel, pheno.col = 1, n.iter = 30000)
qbModel <- qb.model(crOb, epistasis = T, main.nqtl = 3, mean.nqtl = 3, depen = FALSE,
  max.qtl = 0, interval = 20)
# chr.nqtl = rep.int(3, 24))
so <- qb.scanone(mc.b, epistasis = T, type.scan = "heritability", chr = 1:24)
best <- qb.BayesFactor.jm(mc.b, items = c("pattern", "nqtl"))
two <- qb.scantwo(mc.b, scan, type.scan = "2logBF")
#################
temp <- qb.slicetwo(mc.b, chr = c(1,18), pos = c(45,12))
qb.varcomp(qbObject, scan, aggregate = TRUE, ...)
qb.pairloci(qbObject, chr, ...)



temp <- qb.sliceone(mc.b, slice = 1, chr = 2:24, type.scan = "cellmean")
png("/home/jmiller1/public_html/slice.1.18.png", width = 3000)
plot(temp)
dev.off()


slice <- qb.slicetwo(mc.b,chr = c(13,18),pos=c( 78.1 , 130.8),width = 50)
png("/home/jmiller1/public_html/slice.1.18.png", width = 3000)
plot(slice)
dev.off()



a <- find.marker(cross.18, 1,247.41 )
b <- find.marker(cross.18, 18, 130.8)
png("/home/jmiller1/public_html/rqtl.check.png")
effectplot(crOb, mname1 = a, mname2 = b)
dev.off()




z <- 3
gt <- geno.table(crOb,chr=z)
pos <- gsub(paste(z,':',sep=''),'',rownames(gt))
png("/home/jmiller1/public_html/rqtl.check.png")
plot(pos,gt$AB,ylim=c(8,45))
points(pos,gt$AA, col='blue')
points(pos,gt$BB, col='red')
dev.off()

a <- find.marker(crOb, 1,250 )
png("/home/jmiller1/public_html/elr_pxg")
plotPXG(crOb, a, pheno.col = 1, jitter = 1.5, infer = F, pch = 19)
dev.off()




best.hpd <- qb.hpdone(mc.b,profile = "2logBF")
png("/home/jmiller1/public_html/elr.best.hpd.one.png", width = 1000)
plot(best.hpd, chr = c(1:24), cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5,
  cex.sub = 2.5, xlab = NA)
dev.off()


two.b <- qb.scantwo(mc.b, chr = c(1,2,7,8,9,18,22,23),type = "LPD")
png("/home/jmiller1/public_html/elr.two.png", width = 1000, height = 1000)
plot(two.b, cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5,
  cex.sub = 2.5, xlab = NA)
dev.off()






two <- qb.scantwo(mc.b, chr = summary(best.b)$chr,type = "LPD")
slice <- qb.slicetwo(mc.b,chr = c(2,18),pos=c(241.76, 125.77))


, summary(best.b)$pos[c(1,3)])


plot(slice, figs = c("effects", "cellmean", "effectplot"))


plot(two, chr = 6, slice = 15)

plot(two, chr = 15, slice = 6)

two.lpd <- qb.scantwo(qbHyper, chr = c(6,15),type = "LPD")
plot(two.lpd, chr = 6, slice = 15)
plot(two.lpd, chr = 15, slice = 6)

save.image("/home/jmiller1/public_html/ELR.bim.Rsave")

png("/home/jmiller1/public_html/elr.scanone.so.png", width = 1000)
plot(so, chr = c(1:24), cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5,
  cex.sub = 2.5, xlab = NA)
dev.off()

png("/home/jmiller1/public_html/elr.qbBayes.png", width = 600)
plot(best)
dev.off()

best <- qb.BayesFactor.jm(mc, items = "nqtl")
png("/home/jmiller1/public_html/elr.nqtl.png", width = 600)
plot(best)
dev.off()

png("/home/jmiller1/public_html/elr.model.png", width = 2000)
# plot(qb.close(mc, target = so))
plot(so.LPD, c(1:24))
dev.off()

png("/home/jmiller1/public_html/elr.scanone.st.png", width = 2000, height = 2000)
plot(st, c(1:24))
dev.off()

png("/home/jmiller1/public_html/elr_coda.png", width = 2000)
plot(qb.coda(mc, variables = c("nqtl")))
dev.off()

png("/home/jmiller1/public_html/elr_scan_diagnost.so.png", width = 2000)
plot(qb.hpdone(mc.b, chr = c(1:24), level = 0.9, scan = "sum"))
dev.off()

### qtlbim markers
qtl.bm <- c(1,13, 18)  ##set manual
bimqtl <- summary(qb.scanone(mc.b, type = "heritability", chr = qtl.bm))
qtl.bm <- find.marker(cross.18, rownames(bimqtl)[3], bimqtl$pos[3])

crs <- rownames(summary(best)$pattern[3, ])
crs <- unique(unlist(strsplit(crs,split=',')))
bimqtl <- summary(qb.scanone(mc, type = "heritability", chr = crs))
qtl.bm <- find.marker(cross.18, rownames(bimqtl), bimqtl$pos)

# cross2 <- argmax.geno(crOb, step = 5, off.end = 5, err = 0.1) cross2 <-
# subset(cross2, ind = (!is.na(cross2$pheno$sex)))
cross2 <- convert2cross2(crOb)

for (i in 1:length(qtl.bm)) {
  chr <- crs[i]
  mark <- qtl.bm[i]
  index <- cross2$geno[[as.character(chr)]][, mark] > 0
  png(paste("/home/jmiller1/public_html/", pop, chr, "_interation_plot.png", sep = ""))
  plot_pxg(geno = cross2$geno[[as.character(chr)]][, mark][index], pheno = cross2$pheno[,
    1][index], sort = F, SEmult = 2)
  dev.off()

  png(paste("/home/jmiller1/public_html/", pop, chr, "_pxg.png", sep = ""))
  plotPXG(crOb, mark, pheno.col = 1, jitter = 1, infer = F)
  dev.off()

  print(mark)
  print(table(cross2$geno[[as.character(chr)]][, mark]))
}

## Specific to ELR
png(paste("/home/jmiller1/public_html/elr_pxg", pop, qtl.bm, "_pxg.png", sep = ""))
plotPXG(cross.18, qtl.bm, pheno.col = 1, jitter = 1.5, infer = F, pch = 19, main = paste(qtl.bm,
  pop))
dev.off()

qtl.uns <- makeqtl(cross.18, chr = rownames(bimqtl), pos = bimqtl$pos)

fit.bm <- fitqtl(cross.18, pheno.col = 1, qtl.uns, method = "imp", model = "normal",
  dropone = TRUE, get.ests = TRUE, run.checks = TRUE, tol = 1e-04, maxit = 10000)

capture.output(c(summary(fit), summary(fit.bm), summary(best), summary(so)), file = "/home/jmiller1/public_html/ELR_out.txt")
capture.output(fit.mr, file = "/home/jmiller1/public_html/ELR_removed_mr_out.txt")
