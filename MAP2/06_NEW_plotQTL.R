new#!/bin/R
### Map QTLs 1 of 3
#debug.cross <- T
#source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library('qtl')
pop <- 'NEW'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

load(file.path(mpath,'scans.nbh.rsave'))

bins <- data.frame(
 em=summary(scan.bin.em),
 imp=summary(scan.bin.imp)[,'lod'],
 mr=summary(scan.bin.mr)[,'lod'],
 np=summary(scan.np.em.b)[,'lod'])

binpo <- data.frame(
 empo=rownames(summary(scan.bin.em)),
 impo=rownames(summary(scan.bin.imp)),
 mrpo=rownames(summary(scan.bin.mr)),
 nppo=rownames(summary(scan.np.em.b)))

norms <- data.frame(
 em=summary(scan.norm.em),
 imp=summary(scan.norm.imp)[,'lod'],
 mr=summary(scan.norm.mr)[,'lod'],
 np=summary(scan.np.em.n)[,'lod'],
 ehk=summary(scan.norm.ehk)[,'lod'])

normpo <- data.frame(
 empo=rownames(summary(scan.norm.em)),
 impo=rownames(summary(scan.norm.imp)),
 mrpo=rownames(summary(scan.norm.mr)),
 nppo=rownames(summary(scan.np.em.n)))


png("/home/jmiller1/public_html/NBH_map.png")
plot(pull.map(cross))
dev.off()

png("/home/jmiller1/public_html/NBH_scan.norm.mr.png")
plot(scan.norm.mr)
dev.off()

png("/home/jmiller1/public_html/NBH_scan.sex.png")
plot(scan.bin.sex)
dev.off()


qtl.int <- cbind(
  a2=summary(out.a.2)[,c('pos','lod')],
  i2=summary(out.i.2)[,c('pos','lod')],
  a8=summary(out.a.8)[,c('pos','lod')],
  i8=summary(out.i.8)[,c('pos','lod')],
  a18=summary(out.a.18)[,c('pos','lod')],
  i18=summary(out.i.18)[,c('pos','lod')])




a <- find.marker(cross, 1, 8)
b <- find.marker(cross, 18, 27)
png("/home/jmiller1/public_html/NBH_1_18.png")
effectplot(cross, pheno.col = 4, mname1 = b, mname2 = a, ylim = c(0, 1), main = "Genotype interaction \nChrs 1 at 8cm (AHR) and 18 (AHRb) at 27cm")
dev.off()


a <- find.marker(cross, 1, 8)
b <- find.marker(cross, 18, 27)
png("/home/jmiller1/public_html/NBH_2_18.png")
effectplot(cross, pheno.col = 4, mname1 = b, mname2 = a, ylim = c(0, 1), main = "Genotype interaction \nChrs 1 at 27MB (AIP) and 18 (AHRb) at 20MB")
dev.off()


a <- find.marker(cross, 2, 30)
b <- find.marker(cross, 18, 18)
png("/home/jmiller1/public_html/NBH_2_18.png")
effectplot(cross, pheno.col = 5, mname1 = b, mname2 = a, ylim = c(0, 5), main = "Genotype interaction \nChrs 2 at 27MB (AIP) and 18 (AHRb) at 20MB")
dev.off()


a <- find.marker(cross, 2, 80)
b <- find.marker(cross, 13, 38)
png("/home/jmiller1/public_html/NBH_2_13.png")
effectplot(cross, pheno.col = 4, mname1 = b, mname2 = a, ylim = c(0, 1), main = "Genotype interaction \nChrs 2 at 27MB (AIP) and 13 (arnt) at 20MB")
dev.off()

a <- find.marker(cross, 2, 80)
b <- find.marker(cross, 13, 38)
png("/home/jmiller1/public_html/NBH_2_13_norm.png")
effectplot(cross, pheno.col = 5, mname1 = b, mname2 = a,ylim = c(0, 5), main = "Genotype interaction \nChrs 2 at 27MB (AIP) and 13 (arnt) at 20MB")
dev.off()



qtl <- makeqtl(cross, chr=c(2,13,18), pos=c(30,33,28),  what="draws")
fitted <- fitqtl(cross,qtl=qtl, formula=y~Q1+Q3+Q1*Q2)



a <- find.marker(cross, 2, 80)
b <- find.marker(cross, 18, 26)
png("/home/jmiller1/public_html/NBH_2_pxg.png")
plotPXG(cross, a, pheno.col = 1, jitter = 1.5, infer = F, pch = 19, main = a)
dev.off()
png("/home/jmiller1/public_html/NBH_18_pxg.png")
plotPXG(cross, b, pheno.col = 1, jitter = 1.5, infer = F, pch = 19, main = b)
dev.off()


png("/home/jmiller1/public_html/NBH_coef.png")
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred", "green3")
plot(c2eff, map["2"], columns=1:3, col=col)
dev.off()

g <- maxmarg(pr, map, chr=2, pos=80.24051)
png("/home/jmiller1/public_html/NBH_bin_pxg.png")
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_pxg(g, cross2$pheno[,'bin'], ylab="binary phenotype",jitter = 0.5)
dev.off()

g <- maxmarg(pr, map, chr=2, pos=80.24051)
png("/home/jmiller1/public_html/NBH_normal_pxg.png")
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_pxg(g, cross2$pheno[,'Pheno'], ylab="normal phenotype",jitter = 0.5)
dev.off()





plot(out_gwas$lod, out_gwas$snpinfo, altcol="green4", gap=0)
