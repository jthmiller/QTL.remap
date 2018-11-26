#!/bin/R
debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library("qtlbim")

if (pop == "ELR") {
  
  popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/ELR/REMAPS"
  cross.18 <- read.cross(format = "csv", dir = popdir, file = paste(outname, ".BACKUP.QTL_chr.QTLmap.csv", 
    sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", "B"))
  sex <- read.table(file = file.path(dirso, "sex.txt"))
  rownames(sex) <- sex$ID
  cross.18$pheno$sex <- sex[as.character(cross.18$pheno$ID), 2]
  cross.18$pheno$binary <- as.numeric(cross.18$pheno$pheno >= 3)
  crOb <- cross.18
  crOb <- qb.genoprob(crOb, step = 10)
  qbData <- qb.data(crOb, pheno.col = 1, trait = "ordinal", rancov = 2)
  qbModel <- qb.model(crOb, epistasis = T, main.nqtl = 2, mean.nqtl = 2, depen = FALSE)
  mc <- qb.mcmc(crOb, qbData, qbModel, pheno.col = 1, n.iter = 30000)
  so <- qb.scanone(mc, epistasis = T, type = "2logBF")
  so.LPD <- qb.scanone(mc, epistasis = T, type = "LPD")
  scan <- list(upper = "main", lower = "epistasis")
  st <- qb.scantwo(mc, scan, type.scan = "nqtl", chr = c(1, 2, 6, 8, 13, 18, 20, 
    23, 24), epistasis = TRUE)
  
  png("/home/jmiller1/public_html/elr.scanone.so.png", width = 1000)
  plot(so, chr = c(1:24), cex = 1.5, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, 
    cex.sub = 2.5, xlab = NA)
  dev.off()
  
  png("/home/jmiller1/public_html/elr.model.png", width = 2000)
  plot(qb.close(mc, target = so))
  # plot(so.LPD, c(1:24))
  dev.off()
  
  png("/home/jmiller1/public_html/elr.scanone.st.png", width = 2000, height = 2000)
  plot(st, c(1:24))
  dev.off()
  
  # locus specific
  
  png("/home/jmiller1/public_html/elr_13_interation_plot.png")
  plotPXG(cross.18, "13:12791943", pheno.col = 1, jitter = 1.5, infer = F, pch = 19)
  dev.off()
  
  ### Plot chr1, chr6, chr2, chr13, chr20
  
  one <- find.marker(crOb, chr = 13, pos = 84.959)
  two <- find.marker(crOb, chr = 8, pos = 149.64)
  
  png("/home/jmiller1/public_html/elr_interation_plot.png")
  effectplot(crOb, mname1 = one, mname2 = two)
  dev.off()
  
  png("/home/jmiller1/public_html/elr_interation_plot.png")
  plot(qb.pairloci(mc, c(13, 18)))
  dev.off()
  
  chrp <- c(13, 18, 23)
  pos <- summary(so)[as.character(chrp), 2]
  qtl <- find.marker(crOb, chr = chrp, pos = pos)
  names(chrp) <- qtl
  
  fake.f2 <- argmax.geno(crOb, step = 5, off.end = 5, err = 0.1)
  cross2 <- subset(fake.f2, ind = (!is.na(fake.f2$pheno$sex)))
  cross2 <- convert2cross2(cross2)
  
  for (i in 1:length(chrp)) {
    chr <- chrp[i]
    mark <- names(chrp)[i]
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
  
}



cross2 <- convert2cross2(cross2)



### rQTL
cross.18 <- sim.geno(cross.18, error.prob = 0.1, step = 5, n.draws = 250)
scon <- scanone(cross.18, model = "normal", pheno.col = 1, method = "imp", addcovar = cross.18$pheno$sex)

fake.f2 <- argmax.geno(crOb, step = 5, off.end = 5, err = 0.1)
scon <- scanone(fake.f2, model = "binary", pheno.col = 5, method = "imp", addcovar = fake.f2$pheno$sex)

sims <- pull.draws(cross.18)
