#!/bin/R
debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library("qtlbim")
clean <- qtl::clean
cross.18 <- reconst(X = chrms, pop = popq, temp.dir = popdir, a = 2)

# print('Writing the merged chromosome markers to rQTL format')
write.cross(cross.18, filestem = paste(popdir, "/", outname, ".BACKUP.QTL_chr.QTLmap", 
  sep = ""), format = "csv")

popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NEW/REMAPS"
popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NBH/REMAPS"
popdir <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/ELR/REMAPS"

cross.18 <- read.cross(format = "csv", dir = popdir, file = paste(outname, ".BACKUP.QTL_chr.QTLmap.csv", 
  sep = ""), geno = c("AA", "AB", "BB"), alleles = c("A", "B"))

## NBH Single and epistasis scans
crOb <- cross.18
crOb <- qb.genoprob(crOb, step = 1)
qbData.nbh <- qb.data(crOb, pheno.col = 1, trait = "ordinal")
qbModel.nbh <- qb.model(crOb, epistasis = T, main.nqtl = 2, interval = 1200, chr.nqtl = rep(1, 
  nchr(crOb)), mean.nqtl = 2, depen = FALSE)
qbModel.nbh <- qb.model(crOb, epistasis = T, main.nqtl = 4, mean.nqtl = 5, depen = FALSE)

mc <- qb.mcmc(crOb, qbData.nbh, qbModel.nbh, pheno.col = 1, n.iter = 3000)
so <- qb.scanone(mc, epistasis = T, type = "2logBF")
so.LPD <- qb.scanone(mc, epistasis = T, type = "LPD")
st <- qb.scantwo(mc, chr = c(1, 2, 8, 18, 24))
st.LPD <- qb.scantwo(mc, chr = c(1, 2, 8, 18, 24), type = "LPD")

png("/home/jmiller1/public_html/scan_diagnost.so.png", width = 2000)
plot(qb.coda(mc, variables = c("nqtl")))
plot(qb.hpdone(mc))
plot(qb.epistasis(mc))
dev.off()


png("/home/jmiller1/public_html/scanone.so.png", width = 2000)
plot(so, chr = c(1:24))
dev.off()

png("/home/jmiller1/public_html/scanone.so.LPD.png", width = 2000)
plot(so.LPD, c(1:24))
dev.off()

png("/home/jmiller1/public_html/scanone.st.png", width = 2000, height = 2000)
plot(st, c(1:24))
dev.off()
