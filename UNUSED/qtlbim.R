## qtlbim analysis

# install_github('cran/qtlbim')
library("qtlbim")

load("~/NEW.NBH.MAP.QTLmap.Rsave")

## NBH Single and epistasis scans
crOb <- NBH$cross.18
#crOb <- cross.nbh
#crOb <- NBH.x
crOb <- qb.genoprob(crOb, step=1)
qbData.nbh <- qb.data(crOb, pheno.col = 1, trait = "ordinal")
qbModel.nbh <- qb.model(crOb, epistasis = T, main.nqtl = 5, interval = 50, chr.nqtl = rep(2,
  nchr(crOb)), mean.nqtl = 4, depen = FALSE)
nbh.mc <- qb.mcmc(crOb, qbData.nbh, qbModel.nbh, pheno.col = 2)
nbh.so <- qb.scanone(nbh.mc, epistasis = T, type = "2logBF")
nbh.so.LPD <- qb.scanone(nbh.mc, epistasis = T, type="LPD")
nbh.st <- qb.scantwo(nbh.mc, chr = c(1,2,8,18,24))
nbh.st.LPD <- qb.scantwo(nbh.mc, chr = c(1,2,8,18,24),type="LPD")

png("/home/jmiller1/public_html/scanone.nbh.so.png", width = 2000)
plot(nbh.so, chr = 1:24)
dev.off()

png("/home/jmiller1/public_html/scanone.nbh.so.LPD.png", width = 2000)
plot(nbh.so.LPD, chr = 1:24)
dev.off()

png("/home/jmiller1/public_html/scanone.nbh.st.png", width = 2000,height=2000)
plot(nbh.st, chr = 1:24)
dev.off()



#######################################
bf <- qb.bf(nbh.mc, item = "nqtl")
qb.bf(nbh.mc, item = "pairs")
qb.BayesFactor(qbObject, items = c("nqtl","pattern","chrom","pairs"),
       cutoff.pattern, cutoff.pairs = 1, nmax = 15, epistasis = TRUE, ...)
bf <- qb.bf(nbh.mc, item = "pattern")
#######################################

crOb <- NEW$cross.18
#crOb <- cross.new
#crOb <- NEW.x
crOb <- qb.genoprob(crOb, step=1)
qbData.new <- qb.data(crOb, pheno.col = 1, trait = "binary")
qbModel.new <- qb.model(crOb, epistasis = T, main.nqtl = 3, interval = 50, chr.nqtl = rep(2,
  nchr(cross)), mean.nqtl = 3, depen = FALSE)
new.mc <- qb.mcmc(crOb, qbData.new, qbModel.new, pheno.col = 2)
new.so <- qb.scanone(new.mc, epistasis = T, type = "2logBF")
new.so.LPD <- qb.scanone(new.mc, epistasis = T, type="LPD")
new.st <- qb.scantwo(new.mc, chr = c(1,2,8,18,24))
new.st.LPD <- qb.scantwo(new.mc, chr = c(1,2,8,18,24),type="LPD")

png("/home/jmiller1/public_html/scanone.new.so.png", width = 2000)
plot(new.so, chr = 1:24)
dev.off()


png("/home/jmiller1/public_html/scanone.new.so.LPD.png", width = 2000)
plot(new.so.LPD, chr = 1:24)
dev.off()

png("/home/jmiller1/public_html/scanone.new.st.png", width = 2000,height=2000)
plot(new.st, chr = 1:24)
dev.off()





effectplot(multitrait, mname1="PVV4", mname2="GH.117C")




#######################################
best <- qb.best(new.mc)
best <- qb.best(new.mc)
new.bf <- qb.bf(new.mc, item =  c("nqtl","chrom"))
new.bf <- qb.bf(new.mc, item = "pattern")
qbModel.new <- qb.model(NEW.x, epistasis = T,chr.nqtl = rep(3,nchr(cross)), mean.nqtl = 10, depen = FALSE)
new.mc <- qb.mcmc(NEW.x, qbData.new, qbModel.new, pheno.col = 2)
qb.BestPattern(new.mc, epistasis = TRUE,category = "nqtl")

   cutoff, score.type =
  c("sq.atten","attenuation","variance","recombination","distance"),
  include = c("nested","all","exact"),
  center = c("median","mean"), level = 5, ...)
#######################################


#####################
crOb <- cross.ELR
#crOb <- cross.elr
#crOb <- ELR.x
crOb <- qb.genoprob(crOb, step=1)
qbData.elr <- qb.data(crOb, pheno.col = 2, trait = "ordinal")
qbModel.elr <- qb.model(crOb, epistasis = T, main.nqtl = 2, interval = 100, chr.nqtl = rep(2,
  nchr(cross)), depen = FALSE)

elr.mc <- qb.mcmc(crOb, qbData.elr, qbModel.elr, pheno.col = 2, genoupdate = TRUE)
elr.so <- qb.scanone(elr.mc, epistasis = T, type = "2logBF")
elr.st <- qb.scantwo(elr.mc, chr = c(1,2,8,13,18,24))
elr.so.LPD <- qb.scanone(elr.mc, epistasis = T, type="LPD")
elr.st.LPD <- qb.scantwo(elr.mc, chr = c(1,2,8,13,18,24),type="LPD")

png("/home/jmiller1/public_html/scanone.elr.so.png", width = 2000)
plot(elr.so, chr = 1:24)
dev.off()
png("/home/jmiller1/public_html/scanone.elr.so.LPD.png", width = 2000)
plot(elr.so.LPD, chr = 1:24)
dev.off()




NEW$cross.18 <- switchAlleles(NEW$cross.18,markernames(NEW$cross.18,chr=18))
cross.all <- c(NBH$cross.18, NEW$cross.18)


crOb <- qb.genoprob(cross.all, step=1)
qbData <- qb.data(crOb, pheno.col = 2, trait = "ordinal",fixcov=5)
qbModel <- qb.model(crOb, epistasis = T, main.nqtl = 5, interval = 75, chr.nqtl = rep(1,
  nchr(crOb)), depen = FALSE, max.nqtl=10)
qbObject.all <- qb.mcmc(crOb, qbData, qbModel, pheno.col = 2, genoupdate = TRUE)

so <- qb.scanone(qbObject, epistasis = T, type = "2logBF")
st <- qb.scantwo(qbObject, chr = c(1,2,8,13,18,24))
po <- qb.scanone(qbObject, type = "posterior")

so.LPD <- qb.scanone(qbObject, epistasis = T, type="LPD")
st.LPD <- qb.scantwo(qbObject, chr = c(1,2,8,13,18,24),type="LPD")

png("/home/jmiller1/public_html/scanone.so.png", width = 2000)
plot(so, chr = 1:24)
dev.off()
png("/home/jmiller1/public_html/scanone.so.LPD.png", width = 2000)
plot(so.LPD, chr = 1:24)
dev.off()












qb.scanone(qbObject, epistasis = TRUE, scan, type.scan, covar, adjust.covar, chr,
  sum.scan = "yes", min.iter = 1, aggregate = TRUE, smooth = 3, weight = c("sqrt",
    "count", "none"), split.chr, center.type = c("mode", "mean", "scan"), half = FALSE,
  verbose = FALSE)
