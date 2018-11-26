debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library("qtlbim")

if (pop == "NEW") {
  
  popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NEW/REMAPS"
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
  
  png("/home/jmiller1/public_html/new_coda.png", width = 2000)
  # plot(qb.coda(mc, variables = c('nqtl'))) plot(qb.epistasis(mc))
  plot(qb.diag(mc))
  dev.off()
  
  png("/home/jmiller1/public_html/new.scanone.so.png", width = 2000)
  plot(so, chr = c(1:24))
  dev.off()
  
  png("/home/jmiller1/public_html/new.close.png", width = 2000)
  plot(qb.close(mc, target = so))
  # plot(so.LPD, c(1:24))
  dev.off()
  
  png("/home/jmiller1/public_html/new.scanone.st.png", width = 2000, height = 2000)
  plot(st, c(1:24))
  dev.off()
  
  png("/home/jmiller1/public_html/new.close.png", width = 2000, height = 2000)
  plot(qb.close(mc, target = so))
  dev.off()
  
  png("/home/jmiller1/public_html/new.hpdone.png", width = 2000, height = 2000)
  plot(qb.hpdone(mc, effects = "epistasis"))
  dev.off()
  
  chrp <- c(2, 18)
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
