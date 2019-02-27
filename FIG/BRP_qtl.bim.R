debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library("qtlbim")
clean <- qtl::clean
genotyped.only <- T

######### MAP BRP ####

# cross.18 <- reconst(X = chrms, pop = popq, temp.dir = popdir, a = 2)
# print('Writing the merged chromosome markers to rQTL format')
# write.cross(cross.18, filestem = paste(popdir, '/', outname,
# '.BACKUP.QTL_chr.QTLmap', sep = ''), format = 'csv')

load("/home/jmiller1/public_html/BRP_remap.Rsave")
cross.18 <- brp.remap
if (genotyped.only == T) cross.18 <- subset(cross.18, ind = cross.18$pheno$gt == 
  1)
popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/BRP/REMAPS"
popq <- "BRP"
mak <- markernames(cross.18)
cross.18 <- switchAlleles(cross.18, markers = mak)
# mak <- markernames(cross.18, chr = 18) cross.18 <- switchAlleles(cross.18,
# markers = mak)


# popdir <- '/home/jmiller1/QTL_Map_Raw/popgen/rQTL/BRP/REMAPS' cross.18 <-
# read.cross(format = 'csv', dir = popdir, file = paste(outname,
# '.BACKUP.QTL_chr.QTLmap.csv', sep = ''), geno = c('AA', 'AB', 'BB'), alleles =
# c('A', 'B'))

# sex <- read.table(file = file.path(dirso, 'sex.txt')) rownames(sex) <- sex$ID
# cross.18$pheno$sex <- sex[as.character(cross.18$pheno$ID), 2]
cross.18$pheno$binary <- as.numeric(cross.18$pheno$pheno >= 3)

# rqtl binary scan for prior
cross.18$pheno$nqrank <- nqrank(cross.18$pheno$pheno)
### only for brp
cross.18 <- sim.geno(cross.18, error.prob = 0.1, step = 5, n.draws = 500)

scan.norm.imp <- scanone(cross.18, model = "normal", pheno.col = 1, method = "imp", 
  addcovar = cross.18$pheno$sex)
perms.norm.imp <- scanone(cross.18, method = "imp", model = "normal", n.perm = 500, 
  pheno.col = 1, perm.strata = as.character(cross.18$pheno$gt))
norm.qtl <- summary(scan.norm.imp, perms = perms.norm.imp, alpha = 0.05)
qtl.uns <- makeqtl(cross.18, chr = norm.qtl$chr, pos = norm.qtl$pos)

## n.cluster = slurmcore removed. Farm not working
cross.18 <- calc.genoprob(cross.18, error.prob = 0.1)
scan.bin.mr <- scanone(cross.18, method = "mr", model = "binary", pheno.col = 5)
perms.norm.mr <- scanone(cross.18, method = "mr", model = "binary", n.perm = 5000, 
  perm.strata = cross.18$pheno$gt, pheno.col = 5)
bin.qtl <- summary(scan.bin.mr, perms = perms.norm.mr, alpha = 0.05)
qtl.mr <- makeqtl(cross.18, chr = bin.qtl$chr, pos = bin.qtl$pos, what = "prob")

fit.mr <- fitqtl(cross.18, pheno.col = 5, qtl.mr, method = "hk", model = "binary", 
  dropone = TRUE, get.ests = TRUE, run.checks = TRUE, tol = 1e-04, maxit = 10000)

# full <- stepwiseqtl(cross.18, additive.only = T, method = 'imp', pheno.col = 1,
# scan.pairs = T)
qtl.mr <- find.marker(cross.18, qtl.mr$chr, qtl.mr$pos)

qtl <- find.marker(cross.18, qtl.uns$chr, qtl.uns$pos)

### Reg and interval mapping
pop <- "BRP"
png("/home/jmiller1/public_html/brp_bin_regression.png", width = 2000)
plot(scan.bin.mr, main = pop, bandcol = "gray70", cex.lab = 2, cex.axis = 2, cex.main = 2, 
  cex.sub = 2)
abline(h = summary(perms.norm.mr)[1, 1])
dev.off()

png("/home/jmiller1/public_html/brp_norm_imp.png", width = 2000)
plot(scan.norm.imp, main = pop, bandcol = "gray70", cex.lab = 2, cex.axis = 2, cex.main = 2, 
  cex.sub = 2)
abline(h = summary(perms.norm.imp)[1, 1])
dev.off()

### Effect plots

png(paste("/home/jmiller1/public_html/brp_effectplot", pop, qtl[1], "_", qtl[2], 
  "_pxg.png"))
effectplot(cross.18, mname1 = qtl[1], mname2 = qtl[2])
dev.off()

png(paste("/home/jmiller1/public_html/brp_effectplot", pop, qtl[2], "_", qtl[1], 
  "_pxg.png"))
effectplot(cross.18, mname1 = qtl[2], mname2 = qtl[1])
dev.off()

for (i in 1:length(qtl)) {
  png(paste("/home/jmiller1/public_html/brp_pxg", pop, qtl[i], "_pxg.png"))
  plotPXG(cross.18, qtl[i], pheno.col = 1, jitter = 1.5, infer = F, pch = 19, main = paste(qtl[i], 
    pop))
  dev.off()
}

crOb <- cross.18
crOb <- qb.genoprob(crOb, step = 10)
########### 
qbData.b <- qb.data(crOb, pheno.col = 5, trait = "binary")
qbModel <- qb.model(crOb, epistasis = T, main.nqtl = 2, mean.nqtl = 2, depen = FALSE, 
  max.qtl = 0)
mc.b <- qb.mcmc(crOb, qbData.b, qbModel, pheno.col = 5, n.iter = 30000)
so <- qb.scanone(mc.b, epistasis = T, type.scan = "heritability", chr = 1:24)
best <- qb.BayesFactor.jm(mc.b, items = c("pattern", "nqtl"))
two <- qb.scantwo(mc.b, scan, type.scan = "2logBF")
####### 
qbData <- qb.data(crOb, pheno.col = 1, trait = "ordinal", rancov = 2)
mc <- qb.mcmc(crOb, qbData, qbModel, pheno.col = 1, n.iter = 30000)
qbModel <- qb.model(crOb, epistasis = T, main.nqtl = 3, mean.nqtl = 3, depen = FALSE, 
  max.qtl = 0, interval = 20)
# chr.nqtl = rep.int(3, 24))
so <- qb.scanone(mc.b, epistasis = T, type.scan = "heritability", chr = 1:24)
best <- qb.BayesFactor.jm(mc.b, items = c("pattern", "nqtl"))
two <- qb.scantwo(mc.b, scan, type.scan = "2logBF")


save.image("/home/jmiller1/public_html/BRP.bim.Rsave")

png("/home/jmiller1/public_html/brp.qbBayes.png", width = 3000)
plot(best)
dev.off()

png("/home/jmiller1/public_html/brp.scanone.so.png", width = 1000)
plot(so, chr = c(1:24), cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, 
  cex.sub = 2.5, xlab = NA)
dev.off()

png("/home/jmiller1/public_html/brp.model.png", width = 2000)
# plot(qb.close(mc, target = so))
plot(so.LPD, c(1:24))
dev.off()

png("/home/jmiller1/public_html/brp.scanone.st.png", width = 2000, height = 2000)
plot(st, c(1:24))
dev.off()

png("/home/jmiller1/public_html/brp.coda.png", width = 2000)
plot(qb.coda(mc, variables = c("nqtl")))
dev.off()

png("/home/jmiller1/public_html/brp_scan_diagnost.so.png", width = 2000)
plot(qb.hpdone(mc.b))
dev.off()

# qtl.bm <- as.character(summary(qb.hpdone(mc))$chr)
bimqtl <- summary(qb.scanone(mc, type = "heritability", chr = c(2, 18)))
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

fit.bm <- fitqtl(cross.18, pheno.col = 1, qtl.uns, method = "imp", model = "normal", 
  dropone = TRUE, get.ests = TRUE, run.checks = TRUE, tol = 1e-04, maxit = 10000)

capture.output(c(summary(fit), summary(fit.bm), summary(best), summary(so)), file = "/home/jmiller1/public_html/BRP_out.txt")
save.image("/home/jmiller1/public_html/brp.rsave")
