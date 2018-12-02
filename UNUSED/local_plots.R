


cross.NBH <- read.cross(format = "csv", file = "NBH.csv", geno = c("AA", "AB", "BB"),
  alleles = c("A", "B"))

cross.ELR <- read.cross(format = "csv", file = "ELR.csv", geno = c("AA", "AB", "BB"),
  alleles = c("A", "B"))

cross.NEW <- read.cross(format = "csv", file = "NEW.csv", geno = c("AA", "AB", "BB"),
  alleles = c("A", "B"))



cross.NBH$pheno$pheno_norm <- nqrank(cross.NBH$pheno$pheno_05)
scan.norm.imp <- scanone(cross.NBH, method = "imp", model = "normal", pheno.col = 5)



cross.ELR$pheno$pheno_norm <- nqrank(cross.ELR$pheno$pheno_05)
cross.ELR <- sim.geno(cross.ELR, n.draws = 200, step = 1, off.end = 0, error.prob = 0.01,
  map.function = "kosambi", stepwidth = "fixed")

scan.norm.imp.elr <- scanone(cross.ELR, method = "imp", model = "normal", pheno.col = 5)
elr.q1.pos <- summary(scan.norm.imp.elr)[c(2,18),2]
elr.q1.chr <- summary(scan.norm.imp.elr)[c(13,18),1]
elr.2.qtl <- makeqtl(cross.ELR, chr = elr.q1.chr, pos = elr.q1.pos)
elr.2.fit <- fitqtl(cross.ELR, qtl = elr.2.qtl, formula = y ~ Q1 + Q2)
marks <- find.marker(cross.ELR,elr.2.qtl$chr,elr.2.qtl$pos)

plotPXG(cross.ELR,marks[1], pheno.col=5, jitter=1, infer=TRUE,pch=19)
plotPXG(cross.ELR,marks[2], pheno.col=5, jitter=1, infer=TRUE,pch=19)

effectplot(cross.ELR, pheno.col=5, mname1=marks[1],mname2=marks[2])

               geno2, main, ylim, xlab, ylab, col, add.legend=TRUE,
               legend.lab, draw=TRUE, var.flag=c("pooled","group"))
