load('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NBH/REMAPS/QTLmap.Rsave')


scan2.norm.imp <- scantwo(cross.18, pheno.col=6, model="normal",
  method="imp",addcovar=NULL, intcovar=NULL, weights=NULL,
  use="complete.obs",incl.markers=FALSE, clean.output=FALSE,
  clean.nmar=1, clean.distance=0,maxit=4000, tol=1e-4,
  verbose=TRUE, perm.Xsp=FALSE, perm.strata=NULL,
  assumeCondIndep=FALSE, batchsize=250, n.cluster=24)


registerDoParallel(slurmcore)
perms2.norm.imp <- foreach(i=50,
  .combine=c,.packages = "qtl") %dopar% {

 operm <- scantwo(cross.18, pheno.col=6, model="normal",
  method="imp",addcovar=NULL, intcovar=NULL, weights=NULL,
  use="complete.obs",incl.markers=FALSE, clean.output=FALSE,
  clean.nmar=1, clean.distance=0,maxit=4000, tol=1e-4,
  verbose=TRUE, n.perm=1 , perm.Xsp=FALSE, perm.strata=cross.18$pheno$stata,
  assumeCondIndep=FALSE, batchsize=250, n.cluster=1)

}
pen.em <- calc.penalties(operm.hk)
save.image(paste('QTLmap.Rsave',sep=''))
