debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library("qtlbim")
clean <- qtl::clean

if(pop=='NBH'){

  popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NBH/REMAPS"
  cross.18 <- read.cross(format = "csv", dir = popdir, file = paste(outname, ".BACKUP.QTL_chr.QTLmap.csv",
    sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", "B"))
  sex <- read.table(file = file.path(dirso, "sex.txt"))
  rownames(sex) <- sex$ID
  cross.18$pheno$sex <- sex[as.character(cross.18$pheno$ID), 2]
  cross.18$pheno$binary <- as.numeric(cross.18$pheno$pheno >= 3)
  crOb <- cross.18

  # rqtl binary scan for prior
  cross.18$pheno$nqrank <- nqrank(cross.18$pheno$pheno)
  cross.18 <- sim.geno(cross.18, error.prob = 0.025, step = 5, n.draws = 250)

  scan.norm.imp <- scanone(cross.18, model = "normal", pheno.col = 1, method = "imp", addcovar = cross.18$pheno$sex)
  scan.bin.mr <- scanone(cross.18, method = "mr", model = "binary", pheno.col = 5)

  perms.norm.imp <- scanone(cross.18, method = "imp", model = "normal", n.perm = 50,
     pheno.col = 1,perm.strata = as.character(cross.18$pheno$gt))

  perms.norm.mr <- scanone(cross.18, method = "mr", model = "binary", n.perm = 500,
    perm.strata = cross.18$pheno$gt, pheno.col = 5)

    ## n.cluster = slurmcore removed. Farm not working

  norm.qtl <- summary(scan.norm.imp, perms = perms.norm.imp, alpha = 0.05)
  bin.qtl <- summary(scan.bin.mr, perms = perms.norm.mr, alpha = 0.05)

  qtl.uns <- makeqtl(cross.18, chr = norm.qtl$chr, pos = norm.qtl$pos)
  qtl.mr <- makeqtl(cross.18, chr = bin.qtl$chr, pos = bin.qtl$pos)

  #full <- stepwiseqtl(cross.18, additive.only = T, method = "imp", pheno.col = 1, scan.pairs = T)
  fit <- fitqtl(cross.18, pheno.col=1, qtl.uns,method="imp",
    model="normal", dropone=TRUE, get.ests=TRUE,
    run.checks=TRUE, tol=1e-4, maxit=10000)

  cross.nomis <- subset(cross.18,ind=cross.18$pheno$gt==1)
  fit.nomis <- fitqtl(cross.nomis, pheno.col=1, qtl.uns,method="imp",
    model="normal", dropone=TRUE, get.ests=TRUE,
    run.checks=TRUE, tol=1e-4, maxit=10000)

  qtl <- find.marker(cross.18, qtl.uns$chr,  qtl.uns$pos)
  qtl.mr <- find.marker(cross.18, qtl.mr$chr,  qtl.mr$pos)

  ### Reg and interval mapping

  png("/home/jmiller1/public_html/nbh_bin_regression.png", width=2000)
  plot(scan.bin.mr,
     bandcol="gray70",cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  abline(h=summary(perms.norm.mr)[1,1])
  dev.off()

  png("/home/jmiller1/public_html/nbh_norm_imp.png", width=2000)
  plot(scan.norm.imp,
     bandcol="gray70",cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2)
  abline(h=summary(perms.norm.imp)[1,1])
  dev.off()

  ### Effect plots

  png(paste("/home/jmiller1/public_html/1_2_nbh_effectplot", pop, qtl[1],"_",qtl[2], "_pxg.png"))
    effectplot(cross.18, mname1 = qtl[1], mname2 = qtl[2])
  dev.off()

  png(paste("/home/jmiller1/public_html/3_4_nbh_effectplot", pop, qtl[3],"_",qtl[4], "_pxg.png"))
    effectplot(cross.18, mname1 = qtl[3], mname2 = qtl[4])
  dev.off()

  png(paste("/home/jmiller1/public_html/1_4_nbh_effectplot", pop, qtl[1],"_",qtl[4], "_pxg.png"))
    effectplot(cross.18, mname1 = qtl[1], mname2 = qtl[4])
  dev.off()

  png(paste("/home/jmiller1/public_html/1_4_nbh_effectplot", pop, qtl[1],"_",qtl[2], "_pxg.png"))
    effectplot(cross.18, mname1 = qtl[1], mname2 = qtl[3])
  dev.off()


for (i in 1:length(qtl)){
  png(paste("/home/jmiller1/public_html/nbh_pxg", pop, qtl[i],"_pxg.png"))
    plotPXG(cross.18, qtl[i], pheno.col = 1, jitter = 1.5, infer = F, pch = 19,main=paste(qtl[i],pop))
  dev.off()
  }

  ### qtlbim
  crOb <- cross.18
  crOb <- qb.genoprob(crOb, step = 10)
  qbData <- qb.data(crOb, pheno.col = 1, trait = "ordinal", rancov = 2)
  qbModel <- qb.model(crOb, epistasis = T, main.nqtl = 4, mean.nqtl = 4, depen = FALSE)
  mc <- qb.mcmc(crOb, qbData, qbModel, pheno.col = 1, n.iter = 30000)
  so <- qb.scanone(mc, epistasis = T, type = "2logBF")
  so.LPD <- qb.scanone(mc, epistasis = T, type = "LPD")
  scan <- list(upper="main",lower='epistasis')
  st <- qb.scantwo(mc,scan, type.scan="nqtl",chr = c(1, 2,6, 8,13, 18,20, 23,24),epistasis = TRUE)

  png("/home/jmiller1/public_html/nbh.scanone.so.png", width = 1000)
  plot(so, chr = c(1:24), cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5,
    cex.sub = 2.5, xlab = NA)
  dev.off()

  png("/home/jmiller1/public_html/nbh.model.png", width = 2000)
  #plot(qb.close(mc, target = so))
  plot(so.LPD, c(1:24))
  dev.off()

  png("/home/jmiller1/public_html/nbh.scanone.st.png", width = 2000, height = 2000)
  plot(st, c(1:24))
  dev.off()

  png("/home/jmiller1/public_html/nbh.coda.png", width = 2000)
  plot(qb.coda(mc, variables = c('nqtl')))
  dev.off()

  png("/home/jmiller1/public_html/nbh_scan_diagnost.so.png", width = 2000)
  plot(qb.hpdone(mc))
  dev.off()

  ### qtlbim markers
  qtl.bm <- as.character(summary(qb.hpdone(mc))$chr)
  bimqtl <- summary(qb.scanone(mc,type = "heritability",chr=qtl.bm))
  qtl.bm <- find.marker(cross.18, rownames(bimqtl),  bimqtl$pos)


  #cross2 <- argmax.geno(crOb, step = 5, off.end = 5, err = 0.1)
  #cross2 <- subset(cross2, ind = (!is.na(cross2$pheno$sex)))
  cross2 <- convert2cross2(cross.18)

for (i in 1:length(qtl.bm)){
  chr <- rownames(bimqtl)[i]
  mark <- qtl.bm[i]
  index <- cross2$geno[[as.character(chr)]][,mark]>0
  png(paste("/home/jmiller1/public_html/",pop,chr,"_interation_plot.png"))
    plot_pxg(geno=cross2$geno[[as.character(chr)]][,mark][index],pheno=cross2$pheno[,1][index],sort=F,SEmult=2)
  dev.off()

  png(paste("/home/jmiller1/public_html/",pop,chr,"_pxg.png"))
    plotPXG(crOb , mark, pheno.col=1, jitter=1, infer=F)
  dev.off()

  print(mark)
  print(table(cross2$geno[[as.character(chr)]][,mark]))
}


## Specific to NBH
png(paste("/home/jmiller1/public_html/NBH_pxg", pop, qtl.bm,"_pxg.png",sep=''))
  plotPXG(cross.18, qtl.bm, pheno.col = 1, jitter = 1.5, infer = F, pch = 19,main=paste(qtl.bm,pop))
dev.off()

qtl.uns <- makeqtl(cross.18, chr = rownames(bimqtl), pos = bimqtl$pos)

fit.bm <- fitqtl(cross.18, pheno.col=1, qtl.uns,method="imp",
  model="normal", dropone=TRUE, get.ests=TRUE,
  run.checks=TRUE, tol=1e-4, maxit=10000)

capture.output(c(summary(fit),summary(fit.bm),summary(fit.nomis)), file = "/home/jmiller1/public_html/NBH_out.txt")
