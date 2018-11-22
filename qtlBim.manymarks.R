#!/bin/R
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library("qtlbim")

cross.18 <- reconst(X = chrms, pop = popq, temp.dir = popdir, a = 2)

NBH.x <- 
## NBH Single and epistasis scans
crOb <- cross.18
crOb <- qtl::clean(crOb)
crOb <- qb.genoprob(crOb, step = 1)
qbData.nbh <- qb.data(crOb, pheno.col = 1, trait = "binary")
qbModel.nbh <- qb.model(crOb, epistasis = T, main.nqtl = 5, interval = 50, chr.nqtl = rep(2, 
  nchr(crOb)), mean.nqtl = 4, depen = FALSE)
mc <- qb.mcmc(crOb, qbData.nbh, qbModel.nbh, pheno.col = 2)
so <- qb.scanone(nbh.mc, epistasis = T, type = "2logBF")
so.LPD <- qb.scanone(nbh.mc, epistasis = T, type = "LPD")
st <- qb.scantwo(nbh.mc, chr = c(1, 2, 8, 18, 24))
st.LPD <- qb.scantwo(nbh.mc, chr = c(1, 2, 8, 18, 24), type = "LPD")

png("/home/jmiller1/public_html/scanone.nbh.so.png", width = 2000)
plot(nbh.so, chr = 1:24)
dev.off()

png("/home/jmiller1/public_html/scanone.nbh.so.LPD.png", width = 2000)
plot(nbh.so.LPD, chr = 1:24)
dev.off()

png("/home/jmiller1/public_html/scanone.nbh.st.png", width = 2000, height = 2000)
plot(nbh.st, chr = 1:24)
dev.off()
