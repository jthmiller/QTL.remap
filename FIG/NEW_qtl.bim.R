#!/bin/R
debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library("qtlbim")
clean <- qtl::clean
genotyped.only <- T

popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NEW/REMAPS"
cross.18 <- read.cross(format = "csv", dir = popdir, file = paste(outname, ".BACKUP.QTL_chr.QTLmap.csv", 
  sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", "B"))
sex <- read.table(file = file.path(dirso, "sex.txt"))
rownames(sex) <- sex$ID
cross.18$pheno$sex <- sex[as.character(cross.18$pheno$ID), 2]
cross.18$pheno$binary <- as.numeric(cross.18$pheno$pheno >= 3)

pheno <- read.csv("~/QTL_Map_Raw/popgen/rQTL/data/PhenoDist.csv")
rownames(pheno) <- paste(pheno$pop_all, pheno$Sample, sep = "_")
cross.18$pheno$gt <- pheno[as.character(cross.18$pheno$ID), 6]

# mak <- markernames(cross.18, chr = flips)
cross.18 <- switchAlleles(cross.18, markers = markernames(cross.18))


#### IF GENOTYPED IND ONLY
if (genotyped.only == T) cross.18 <- subset(cross.18, ind = cross.18$pheno$gt == 
  "GT")
#### Pheno to GT/NG Remove problematic individuals (found by kinship analysis)
con <- file(file.path(popdir, "kinship.keep.ind.txt"), open = "r")
keepers <- readLines(con)
close(con)

print("Dropping kinship outliers")
cross.18 <- subset(cross.18, ind = cross.18$pheno$ID %in% keepers)
cross.18 <- subset(cross.18, ind = !is.na(cross.18$pheno$phen))


# rqtl binary scan for prior cross.18$pheno$nqrank <-
# nqrank(cross.18$pheno$pheno)
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
cross.18 <- calc.genoprob(cross.18, error.prob = 0.025)
scan.bin.mr <- scanone(cross.18, method = "mr", model = "binary", pheno.col = 5)
perms.norm.mr <- scanone(cross.18, method = "mr", model = "binary", n.perm = 5000, 
  perm.strata = cross.18$pheno$gt, pheno.col = 5)
bin.qtl <- summary(scan.bin.mr, perms = perms.norm.mr, alpha = 0.05)
qtl.mr <- makeqtl(cross.18, chr = bin.qtl$chr, pos = bin.qtl$pos, what = "prob")
fit.mr <- fitqtl(cross.18, pheno.col = 5, qtl.mr, method = "hk", model = "binary", 
  dropone = TRUE, get.ests = TRUE, run.checks = TRUE, tol = 1e-04, maxit = 10000)
qtl.mr <- find.marker(cross.18, qtl.mr$chr, qtl.mr$pos)

### Reg and interval mapping

png("/home/jmiller1/public_html/new_bin_regression.png", width = 2000)
plot(scan.bin.mr, main = pop, bandcol = "gray70", cex.lab = 2, cex.axis = 2, cex.main = 2, 
  cex.sub = 2)
abline(h = summary(perms.norm.mr)[1, 1])
dev.off()

png("/home/jmiller1/public_html/new_norm_imp.png", width = 2000)
plot(scan.norm.imp, main = pop, bandcol = "gray70", cex.lab = 2, cex.axis = 2, cex.main = 2, 
  cex.sub = 2)
abline(h = summary(perms.norm.imp)[1, 1])
dev.off()

### Effect plots

png(paste("/home/jmiller1/public_html/new_effectplot", pop, qtl[1], "_", qtl[2], 
  "_pxg.png"))
effectplot(cross.18, mname1 = qtl[1], mname2 = qtl[2])
dev.off()

png(paste("/home/jmiller1/public_html/new_effectplot", pop, qtl[2], "_", qtl[1], 
  "_pxg.png"))
effectplot(cross.18, mname1 = qtl[2], mname2 = qtl[1])
dev.off()

for (i in 1:length(qtl)) {
  png(paste("/home/jmiller1/public_html/new_pxg", pop, qtl[i], "_pxg.png"))
  plotPXG(cross.18, qtl[i], pheno.col = 1, jitter = 1.5, infer = F, pch = 19, main = paste(qtl[i], 
    pop))
  dev.off()
}


crOb <- cross.18
# crOb <- switchAlleles(crOb, markers = swits)
crOb <- qb.genoprob(crOb, step = 10)
########### 
qbData.b <- qb.data(crOb, pheno.col = 5, trait = "binary")
qbModel.b <- qb.model(crOb, epistasis = T, main.nqtl = 2, mean.nqtl = 3, depen = FALSE, 
  max.qtl = 0)
mc.b <- qb.mcmc(crOb, qbData.b, qbModel.b, pheno.col = 5, n.iter = 30000)
so.b <- qb.scanone(mc.b, epistasis = T, type.scan = "heritability", chr = 1:24)
so.b.bf <- qb.scanone(mc.b, epistasis = T, type.scan = "2logBF", chr = 1:24)
best.b <- qb.BayesFactor.jm(mc.b, items = c("pattern", "nqtl"))
two.b <- qb.scantwo(mc.b, chr = c(1, 2, 7, 8, 9, 13, 18, 22, 23))
close.b <- qb.close(mc.b)
# try <- qb.BestPattern(mc.b, epistasis = TRUE,include ='exact', category =
# 'nqtl', level = 5)

slice <- qb.sliceone(mc.b, slice = 3, smooth = 5, type.scan = "cellmean")
png("/home/jmiller1/public_html/slice.png")
plot(temp, chr = c(12, 13, 24))
dev.off()

#### FOR FIGURE

a <- find.marker(cross.18, 2, 107)
b <- find.marker(cross.18, 18, 106)
png("/home/jmiller1/public_html/NEW_2_18.png")
effectplot(crOb, pheno.col = 1, mname1 = a, mname2 = b, ylim = c(0, 5), main = NULL)
dev.off()

png("/home/jmiller1/public_html/NEW_2_pxg.png")
plotPXG(cross.18, a, pheno.col = 1, jitter = 1.5, infer = F, pch = 19, main = a)
dev.off()
png("/home/jmiller1/public_html/NEW_18_pxg.png")
plotPXG(cross.18, b, pheno.col = 1, jitter = 1.5, infer = F, pch = 19, main = b)
dev.off()





#### norm chr.nqtl = rep.int(2, 24))
mc <- qb.mcmc(crOb, qbData, qbModel, pheno.col = 1, n.iter = 30000)
so <- qb.scanone(mc, epistasis = T, type = "2logBF")
so.LPD <- qb.scanone(mc, epistasis = T, type = "LPD")
scan <- list(upper = "main", lower = "epistasis")
st <- qb.scantwo(mc, scan, type.scan = "nqtl", chr = c(1, 2, 6, 8, 13, 18, 20, 23, 
  24), epistasis = TRUE)
two <- qb.scantwo(mc, scan, type.scan = "2logBF")
slice <- qb.sliceone(mc, type = "cellmean", chr = c(2, 8, 18, 24))
so <- qb.scanone(mc, epistasis = T, type.scan = "heritability", chr = 1:24)
best <- qb.BayesFactor.jm(mc, items = "pattern")

save.image("/home/jmiller1/public_html/NEW.bim.Rsave")

png("/home/jmiller1/public_html/new.scanone.so.png", width = 1000)
plot(so, chr = c(1:24), cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, 
  cex.sub = 2.5, xlab = NA)
dev.off()

png("/home/jmiller1/public_html/new.model.png", width = 2000)
# plot(qb.close(mc, target = so))
plot(so.LPD, c(1:24))
dev.off()

png("/home/jmiller1/public_html/new.scanone.st.png", width = 2000, height = 2000)
plot(st, c(1:24))
dev.off()

png("/home/jmiller1/public_html/new.coda.png", width = 2000)
plot(qb.coda(mc, variables = c("nqtl")))
dev.off()

png("/home/jmiller1/public_html/new_scan_diagnost.so.png", width = 2000)
plot(qb.hpdone(mc))
dev.off()

qtl.bm <- as.character(summary(qb.hpdone(mc))$chr)
bimqtl <- summary(qb.scanone(mc, type = "heritability", chr = qtl.bm))
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
  plotPXG(cross.18, mark, pheno.col = 1, jitter = 1, infer = F)
  dev.off()
  
  print(mark)
  print(table(cross2$geno[[as.character(chr)]][, mark]))
}

qtl.uns <- makeqtl(cross.18, chr = rownames(bimqtl), pos = bimqtl$pos)

fit.bm <- fitqtl(cross.18, pheno.col = 1, qtl.uns, method = "imp", model = "normal", 
  dropone = TRUE, get.ests = TRUE, run.checks = TRUE, tol = 1e-04, maxit = 10000)

capture.output(c(summary(fit), summary(fit.bm), summary(best), summary(so), summary(st)), 
  file = "/home/jmiller1/public_html/NEW_genotyped_out.txt")
