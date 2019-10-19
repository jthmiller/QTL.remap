#!/bin/R
### Map QTLs 1 of 3
#debug.cross <- T
#source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library('qtl')
pop <- 'ELR'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

file_list <- list.files(mpath, 'ELR_all_mark_?[0-9]?[0-9]_tsp.csv')

chr <- gsub("ELR_all_mark_",'',file_list)
chr <- as.numeric(gsub("_tsp.csv",'',chr))

elr <- lapply(file_list,function(X){ read.cross(file=X,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)})

gnos <- lapply(elr,function(X){
  data.frame(X[[1]][[1]][['data']],stringsAsFactors=F)
})
gnos <- do.call(cbind,gnos)
gnos <- cbind(elr[[1]]$pheno,gnos)
gnos$ID <- as.character(gnos$ID)

m_names <- unlist(sapply(elr,function(X){
  markernames(X)
}))

colnames(gnos) <- c('Pheno','sex','ID','bin','pheno_norm',m_names)
rownames(gnos) <- elr[[1]]$pheno$ID

map <- c(colnames(elr[[1]]$pheno),unname(unlist(sapply(elr,pull.map))))
zd <- as.numeric(gsub(":.*","",m_names))

zd[is.na(zd)] <- c(1,2,2)
chr <- c(colnames(elr[[1]]$pheno),zd)
info <- c(colnames(elr[[1]]$pheno),m_names)
headers <- rbind(info,chr,map)
colnames(headers) <- headers[1,]
headers[2:3,1:5] <- ''

headers.u <- unname(data.frame(headers,row.names=NULL,stringsAsFactors=FALSE))
gnos.u <- unname(data.frame(lapply(gnos, as.character),row.names=NULL,stringsAsFactors=FALSE))
colnames(headers.u) <- colnames(gnos.u) <- headers.u[1,]
to_write <- rbind(headers.u,gnos.u)
write.table(to_write,file.path(mpath,'elr.mapped.tsp.csv'),sep=',',row.names=F,quote=F,col.names = F)

################################################################################
## scan
################################################################################

fl <- file.path(mpath,'elr.mapped.tsp.csv')

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

cross <- sim.geno(cross)
cross <- calc.genoprob(cross,step=1,error.prob=0.01,off.end=5)

## binary
scan.bin.em <- scanone(cross, method = "em", model = "binary", pheno.col = 4)
scan.bin.mr <- scanone(cross, method = "mr", model = "binary", pheno.col = 4)

## normal
scan.norm.em <- scanone(cross, method = "em", model = "normal", pheno.col = 1)
scan.norm.mr <- scanone(cross, method = "mr", model = "normal", pheno.col = 1)
scan.norm.imp <- scanone(cross, method = "imp", model = "normal", pheno.col = 1)
scan.norm.ehk <- scanone(cross, method = "ehk", model = "normal", maxit = 5000, pheno.col = 1)

## normal transform
scan.normT.em <- scanone(cross, method = "em", model = "normal", pheno.col = 5)
scan.normT.mr <- scanone(cross, method = "mr", model = "normal", pheno.col = 5)
scan.normT.imp <- scanone(cross, method = "imp", model = "normal", pheno.col = 5)
scan.normT.ehk <- scanone(cross, method = "ehk", model = "normal", maxit = 5000, pheno.col = 5)

## non-parametric
scan.np.em.b <- scanone(cross, method = "em", model = "np", pheno.col = 4, maxit = 5000)
scan.np.em.n <- scanone(cross, method = "em", model = "np", pheno.col = 5, maxit = 5000)

## PERMS 5% LOD 4.07
perms.norm.imp.perms <- scanone(cross, method = "imp", model = "normal", maxit = 1000,
  n.perm = 10000, pheno.col = 5, n.cluster = 10)

perms.bin.em.perms <- scanone(cross, method = "em", model = "binary", maxit = 1000,
  n.perm = 10000, pheno.col = 4, n.cluster = 10)

## step-wise
full.norm.add_only <- stepwiseqtl(cross, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=8)
full.bin.add_only <- stepwiseqtl(cross, additive.only = T, model='binary', method = "imp", pheno.col = 4, scan.pairs = T, max.qtl=8)
save.image(file.path(mpath,'qtls_scans.elr.rsave'))

## manual
qtl <- makeqtl(cross, chr=c(13,18), pos=c(62,9),  what="draws")
fitted.twoQTL <- fitqtl(cross,qtl=qtl, formula=y~Q1+Q2, pheno.col=4)
out.ap13 <- addqtl(cross, qtl=qtl, formula = y~Q1 + Q1*Q3 + Q2,  model='normal', method = "imp", pheno.col = 5)
out.ap23 <- addqtl(cross, qtl=qtl, formula = y~Q1 + Q2 + Q2*Q3,  model='normal', method = "imp", pheno.col = 5)

save.image(file.path(mpath,'manual_qtl.scans.elr.rsave'))

## LOD thresholds (1002 permutations)
##      lod
## 5%  4.26
## 10% 3.84
################################################################################
## SCANTWO ON SUBSET

gg <- sim.geno(cross,step=4)
gg_step2 <- reduce2grid(gg)

sc2_normal_imp <- scantwo(gg_step2, pheno.col=5, model="normal",
             method="imp",
             addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs",
             incl.markers=FALSE, clean.output=FALSE,
             clean.nmar=1, clean.distance=0,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=6)

sc2_normal_imp_perms <- scantwo(gg_step2, pheno.col=5, model="normal",
             method="imp", addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs", n.perm=1000,
             incl.markers=FALSE, clean.output=FALSE,
             clean.nmar=1, clean.distance=0,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=6)

sc2_bin_em <- scantwo(gg_step2, pheno.col=4, model="binary",
             method="em",addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs",
             incl.markers=FALSE, clean.output=FALSE,
             clean.nmar=1, clean.distance=0,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=6)

sc2_bin_em_perms <- scantwo(gg_step2, pheno.col=4, model="binary",
             method="em",addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs", n.perm=1000,
             incl.markers=FALSE, clean.output=FALSE,
             clean.nmar=1, clean.distance=0,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=4)

## sm <- summary(object, thresholds,what=c("best", "full", "add", "int"),
##             perms=sc2_normal_imp_perms, alphas, lodcolumn=1,
##             pvalues=FALSE, allpairs=TRUE)
##
save.image(file.path(mpath,'scantwo.scans.elr.rsave'))

################################################################################
full.norm <- stepwiseqtl(gg_step2, additive.only = F, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=6)
save.image(file.path(mpath,'pairwise_scans.elr.rsave'))

full.bin <- stepwiseqtl(gg_step2, additive.only = F, model='binary', method = "imp", pheno.col = 4, scan.pairs = T, max.qtl=6)
save.image(file.path(mpath,'pairwise_scans.elr.rsave'))
################################################################################

fit.full.norm <- fitqtl(gg_step2,qtl=full.norm, pheno.col=4)

################################################################################
# MQM
crossaug <- mqmaugment(gg_step2,strategy='drop')  # Augmentation

result <- mqmscan(crossaug)    # Scan
    # show LOD interval of the QTL on chr 3

mq.add <- mqmscan(crossaug, cofactors=NULL, pheno.col = 5,
      model=c("additive"), forceML=FALSE,
      cofactor.significance=0.05, em.iter=1000,
      window.size=25.0, step.size=5.0,
      logtransform = FALSE, estimate.map = FALSE,
      plot=FALSE, verbose=FALSE, outputmarkers=TRUE,
      multicore=TRUE, batchsize=10, n.clusters=6,
      test.normality=FALSE,off.end=0)

mq.dom <- mqmscan(crossaug, cofactors=NULL, pheno.col = 5,
      model="dominance", forceML=FALSE,
      cofactor.significance=0.05, em.iter=1000,
      window.size=25.0, step.size=5.0,
      logtransform = FALSE, estimate.map = FALSE,
      plot=FALSE, verbose=FALSE, outputmarkers=TRUE,
      multicore=TRUE, batchsize=10, n.clusters=6,
      test.normality=FALSE,off.end=0)

save.image(file.path(mpath,'mqm.scans.elr.rsave'))
################################################################################
