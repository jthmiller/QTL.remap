
load('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NEW/REMAPS/QTLmap.Rsave', envir=NEW)

NBH <- dropSimilarMarkers(NBH$cross.18, rf.threshold = 0.02)
NBH <- repRipple(NBH, error.prob=0.01, map.function="kosambi",window = 6)

## Single
NBH.prem.imp <- scanone(NBH, method="imp",model='normal',n.perm=2000,perm.strata=NBH$pheno$stata,pheno.col=6, n.cluster=24)
NBH.norm.imp <- scanone(NBH, method="imp",model='normal',pheno.col=6)

### Multi-QTL models

NBH.qtl <- summary(NBH.norm.imp, perms=NBH.prem.imp, alpha=0.05)
NBH.uns <- makeqtl(NBH, chr=NBH.qtl$chr, pos=NBH.qtl$pos)
NBH.full <- stepwiseqtl(NBH, additive.only=T, qtl=NBH.uns, method="imp", pheno.col=6, scan.pairs=T)

NBH.scan2.imp <- scantwo(NBH, pheno.col=6, model="normal",
  method="imp",addcovar=NULL, intcovar=NULL, weights=NULL,
  use="complete.obs",incl.markers=FALSE, clean.output=FALSE,
  clean.nmar=1, clean.distance=0,maxit=4000, tol=1e-4,
  verbose=TRUE, perm.Xsp=FALSE, perm.strata=NBH$pheno$stata,
  assumeCondIndep=FALSE, batchsize=250, n.cluster=slurmcore)

registerDoParallel(slurmcore)
perms2.NBH.imp <- foreach(i=50,
  .combine=c,.packages = "qtl") %dopar% {

 operm <- scantwo(NBH, pheno.col=6, model="normal",
  method="imp",addcovar=NULL, intcovar=NULL, weights=NULL,
  use="complete.obs",incl.markers=FALSE, clean.output=FALSE,
  clean.nmar=1, clean.distance=0,maxit=4000, tol=1e-4,
  verbose=TRUE, n.perm=1 , perm.Xsp=FALSE, perm.strata=NBH$pheno$stata,
  assumeCondIndep=FALSE, batchsize=250, n.cluster=1)

}
NBH.pen <- calc.penalties(perms2.NBH.imp)

save.image(paste('~/QTLmap.Rsave',sep=''))
