#!/bin/R
library('qtl')
pop <- 'ELR'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

################################################################################

################################################################################
## SCANTWO ON SUBSET (this map made in tspmap.R
fl <- file.path(mpath,'elr.mapped.tsp.csv')

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE)

cross <- sim.geno(cross)
cross <- calc.genoprob(cross,step=1,error.prob=0.01,off.end=5)


gg <- sim.geno(cross,step=4)
gg_step2 <- reduce2grid(gg)

## Done
sc2_normal_imp <- scantwo(gg_step2, pheno.col=5, model="normal",
             method="imp",
             addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs",
             incl.markers=FALSE, clean.output=FALSE,
             clean.nmar=1, clean.distance=0,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=6)

## Done
sc2_normal_imp_perms <- scantwo(gg_step2, pheno.col=5, model="normal",
             method="imp", addcovar=NULL, intcovar=NULL, weights=NULL,
             use="complete.obs", n.perm=1000,
             incl.markers=FALSE, clean.output=FALSE,
             clean.nmar=1, clean.distance=0,
             maxit=1000, tol=1e-4,
             verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
             assumeCondIndep=FALSE, batchsize=250, n.cluster=6)

## sm <- summary(object, thresholds,what=c("best", "full", "add", "int"),
##             perms=sc2_normal_imp_perms, alphas, lodcolumn=1,
##             pvalues=FALSE, allpairs=TRUE)
