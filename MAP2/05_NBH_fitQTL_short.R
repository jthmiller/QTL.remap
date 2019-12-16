#!/bin/R
library('qtl')
pop <- 'NBH'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

################################################################################
## SCANTWO ON SUBSET

fl <- file.path(mpath,'nbh.mapped.tsp.csv')

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

cross <- sim.geno(cross)
cross <- calc.genoprob(cross,step=1,error.prob=0.01,off.end=5)

gg <- sim.geno(cross,step=4)
gg_step2 <- reduce2grid(gg)

## PERMS 5% LOD 4.07
sc2_normal_imp <- scantwo(gg_step2, pheno.col=5, model="normal",
             method="imp",
             addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs",
             incl.markers=FALSE, clean.output=FALSE,
             clean.nmar=1, clean.distance=0,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=12)

sc2_normal_imp_perms <- scantwo(gg_step2, pheno.col=5, model="normal",
             method="imp", addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs", n.perm=1000,
             incl.markers=FALSE, clean.output=FALSE,
             clean.nmar=1, clean.distance=0,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=12)

save.image(file.path(mpath,'scantwo.scans.nbh.rsave'))

sc2_bin_em <- scantwo(gg_step2, pheno.col=4, model="binary",
             method="em",addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs",
             incl.markers=FALSE, clean.output=FALSE,
             clean.nmar=1, clean.distance=0,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=12)

sc2_bin_em_perms <- scantwo(gg_step2, pheno.col=4, model="binary",
             method="em",addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs", n.perm=1000,
             incl.markers=FALSE, clean.output=FALSE,
             clean.nmar=1, clean.distance=0,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=12)

## sm <- summary(object, thresholds,what=c("best", "full", "add", "int"),
##             perms=sc2_normal_imp_perms, alphas, lodcolumn=1,
##             pvalues=FALSE, allpairs=TRUE)

save.image(file.path(mpath,'scantwo.scans.nbh.rsave'))
################################################################################

################################################################################
add.norm <- stepwiseqtl(gg_step2, additive.only = T, model='normal', method = "imp", pheno.col = 5, scan.pairs = T, max.qtl=6)
full.norm <- stepwiseqtl(gg_step2, incl.markers=F, scan.pairs=T, keeptrace=T, additive.only = F, model='normal', method = "imp", pheno.col = 5, max.qtl=6)

save.image(file.path(mpath,'stepwise_grid_scans.nbh.rsave'))
################################################################################

##################################################################################
### MQM ##########################################################################
##mq <- subset(cross, ind=!cross$pheno$ID %in% 'NBH_5528')
##gt <- geno.table(mq)
##bfixA <- rownames(gt[which(gt$P.value > 0.001 & gt$missing < 5),])
##mq <- pull.markers(mq,bfixA)
##
##marks <- sapply(1:24,function(i){
##    n.missing <- nmissing(subset(mq, chr=i), what="mar")
##    # weight by -log(prop'n missing), but don't let 0 missing go to +Inf
##    wts <- -log( (n.missing+1) / (nind(mq)+1) )
##    # subset of markers on chr 4 spaced >= 5 cM, with weights = -log(prop'n missing)
##    pickMarkerSubset(pull.map(mq)[[i]],1, wts)
##})
##
##pl <- unlist(marks)
##mq <- pull.markers(mq,pl)
##
##crossaug <- mqmaugment(mq,strategy='default',verbose=T)  # Augmentation
##
##result <- mqmscan(crossaug)    # Scan
##    # show LOD interval of the QTL on chr 3
##
##mq.add <- mqmscan(crossaug, cofactors=NULL, pheno.col = 5,
##      model=c("additive"), forceML=FALSE,
##      cofactor.significance=0.02, em.iter=1000,
##      window.size=25.0, step.size=5.0,
##      logtransform = FALSE, estimate.map = FALSE,
##      plot=FALSE, verbose=FALSE, outputmarkers=TRUE,
##      multicore=TRUE, batchsize=10, n.clusters=6, test.normality=FALSE,off.end=0
##      )
##
##mq.dom <- mqmscan(crossaug, cofactors=NULL, pheno.col = 5,
##      model="dominance", forceML=FALSE,
##      cofactor.significance=0.02, em.iter=1000,
##      window.size=25.0, step.size=5.0,
##      logtransform = FALSE, estimate.map = FALSE,
##      plot=FALSE, verbose=FALSE, outputmarkers=TRUE,
##      multicore=TRUE, batchsize=10, n.clusters=6, test.normality=FALSE,off.end=0
##      )
##
##save.image(file.path(mpath,'mqm_scans.nbh.rsave'))
##################################################################################
##
