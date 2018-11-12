#!/bin/R

debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library("ggridges")

qtl.uns
full
out.fq <- fitqtl(ELR$cross.18, qtl = ELR$qtl.uns, formula = y ~ Q1 + Q2 + Q3 * Q4)
out.fq <- fitqtl(NBH$cross.18, qtl = NBH$qtl.uns)

mark <- find.marker(ELR$cross.18, chr = 13, pos = 73)
plot.genos(ELR$cross.18, chr = 13)


cross2 <- convert2cross2(ELR$cross.18)


plt <- cross2$geno$"13"[, mark]
plt[plt == 0] <- NA
plt[plt == 1] <- "AA"
plt[plt == 2] <- "AB"
plt[plt == 3] <- "BB"
plt <- factor(plt, levels = c("AA", "AB", "BB"))


drop <- c("ELR_10869", "ELR_10967", "ELR_11587")

er <- subset(ELR$cross.18, ind = (!ELR$cross.18$pheno$ID %in% drop))

table(paste(plt, cross2$pheno[, 2]))

png("/home/jmiller1/public_html/elr.pheno.qtl.png", width = 600)
plot_pxg(geno = plt, pheno = cross2$pheno[, 2], sort = F, SEmult = 2)
dev.off()

png("/home/jmiller1/public_html/elr.pheno.qtl.png", width = 1600, height = 1100)
plotGeno(ELR$cross.18, chr = 13)
dev.off()


png("/home/jmiller1/public_html/elr.pheno.qtl.png", width = 1000, height = 1000)
rf13 <- pull.rf(er, chr = 11)
rf13.lod <- pull.rf(er, chr = 11, what = "lod")
plot(as.numeric(rf13), as.numeric(rf13.lod))
dev.off()
