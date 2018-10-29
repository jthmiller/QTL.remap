source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')

load(paste('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/',pop,'/REMAPS/QTLmap.Rsave',sep=''))

### Permute the null distribution of LODs
operm <- scantwo(cross.18, pheno.col=6, model="normal",
  method="imp",addcovar=NULL, intcovar=NULL, weights=NULL,
  use="all.obs",incl.markers=FALSE, clean.output=FALSE,
  clean.nmar=1, clean.distance=0,maxit=2000, tol=1e-4,
  verbose=TRUE, n.perm=12 , perm.Xsp=FALSE, perm.strata=cross.18$pheno$stata,
  assumeCondIndep=FALSE, batchsize=250, n.cluster=slurmcore)

pen <- calc.penalties(operm)

### Stepwise model
DS.stepwise.model <- stepwiseqtl(cross.18, pheno.col=6, qtl=qtl.scan1, max.qtl=15,
  covar=NULL,method="imp", model="normal", incl.markers=TRUE, refine.locations=TRUE,
  additive.only=FALSE, scan.pairs=FALSE, penalties=pen,keeplodprofile=TRUE, keeptrace=TRUE,
  verbose=TRUE,tol=1e-4, maxit=1000, require.fullrank=FALSE)

### Test interactions between pairs of loci
DS.scan2.imp <- scantwo(cross.18, pheno.col=6, model="normal",
  method="imp",addcovar=NULL, intcovar=NULL, weights=NULL,
  use="all.obs",incl.markers=FALSE, clean.output=FALSE,
  clean.nmar=1, clean.distance=0,maxit=4000, tol=1e-4,
  verbose=TRUE, perm.Xsp=FALSE, perm.strata=cross.18$pheno$stata,
  assumeCondIndep=FALSE, batchsize=250, n.cluster=slurmcore)

## Summary, citing significance levels and so estimating thresholds from the permutation results
summary(DS.scan2.imp, perms=operm, alpha=rep(0.05, 5))

## Similar, but also getting genome-scan-adjusted p-values
summary(DS.scan2.imp, perms=operm, alpha=c(0.05, 0.05, 0, 0.05, 0.05),pvalues=TRUE)

save.image(paste('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/',pop,'/REMAPS/QTLmap.Rsave',sep=''))
