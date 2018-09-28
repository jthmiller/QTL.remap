#!/bin/bash
setwd('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/QTL/')
source('control_file.R')

cross.18 <- reconst(X=chrms,pop=popq,temp.dir=popdir)
#cross.18 <- reconst(X,pop=popq, dir=qtldir)
print('Writing the merged chromosome markers to rQTL format')
write.cross(cross.18,filestem=paste(qtldir,'BACKUP.QTL_chr.QTLmap',sep=''),format="csv",chr=X)


## Binary phenotype scan

cross.18$pheno$stata <- gsub('[0-9]','',cross.18$pheno$ID)
cross.18 <- calc.genoprob(cross.18, step=1,error.prob=ers,map.function='kosambi')
perms.em <- scanone(cross.18, method="em",model='binary',maxit=4000,
  n.perm=100,n.cluster=10,perm.strata=cross.18$pheno$stata)
scan.em <- scanone(cross.18, method="em",model='binary',maxit=4000,n.cluster=slurmcore)


scan.ehk <- scanone(cross.18, method="ehk",model='binary',maxit=4000,n.cluster=1,perm.strata=cross.18$pheno$stata)




lod.05 <- summary(perms.em)['5%',]
lod.chr <- summary(scan.em)
Qs <- lod.chr$chr[which(lod.chr$lod > lod.05)]
Ps <- lod.chr$pos[which(lod.chr$lod > lod.05)]

qtl.em <- makeqtl(cross.18, chr=Qs, pos=Ps, what="prob")

fitqtl(cross.18, pheno.col=1, qtl.em, covar=NULL, formula=y~Q1, method=("imp"),
          model="binary", dropone=TRUE, get.ests=FALSE,
          run.checks=TRUE, tol=1e-4, maxit=1000, forceXcovar=FALSE)

scan.mr <- scanone(cross.18, method="mr",model='binary',n.cluster=slurmcore)
perms.mr <- scanone(cross.18, method="mr",model='binary',maxit=4000,
  n.perm=1200,n.cluster=slurmcore,perm.strata=cross.18$pheno$stata)

imp.aug <- mqmaugment(cross.18, minprob=1.0,strategy="default")
scan.mqm <- mqmscan(imp.aug,n.cluster=slurmcore)





## Take the hist two QTLs from 2 chromsomes and make QTL
qtl.uns <- makeqtl(cross, chr=rownames(uns)[1:2], pos=uns$pos[1:2],what='prob')
qtl.rf.uns <- refineqtl(cross, qtl=qtl.uns, formula=y~Q1+Q2,method="hk",verbose=FALSE)
out.fit.uns <- fitqtl(cross, qtl=qtl.rf.uns, formula=y~Q1+Q2,method="hk",get.ests=T)

qtl<- makeqtl(cross, chr=rownames(str)[1:2], pos=str$pos[1:2],what='prob')
qtl.rf <- refineqtl(cross, qtl=qtl, formula=y~Q1+Q2,method="hk",verbose=FALSE)
out.fit <- fitqtl(cross, qtl=qtl.rf, formula=y~Q1+Q2,method="hk",get.ests=T)


### After fitting single addative QTLs, search model space for addition addative
stepout1.uns <- stepwiseqtl(cross,qtl=qtl.rf.uns,pheno.col=1,model="binary",
  additive.only=TRUE, max.qtl=7, verbose=FALSE)

stepout1 <- stepwiseqtl(cross,qtl=qtl.rf,pheno.col=1,model="binary",
  additive.only=TRUE, max.qtl=7, verbose=FALSE)

save.image(paste('QTLmap.Rsave',sep=''))

try(
cross.no2 <- subset(cross,chr=c(1,3:24))
registerDoParallel(slurmcore)
operm.hk <- foreach(i=50,
  .combine=c,.packages = "qtl") %dopar% {
    operm <- scantwo(cross.no2, method="hk", n.perm=1,perm.strata=strata)
}
pen.em <- calc.penalties(operm.hk)
save.image(paste('QTLmap.Rsave',sep=''))
)
try(
### qtl scan2
registerDoParallel(slurmcore)
operm.em <- foreach(i=50,
  .combine=c,.packages = "qtl") %dopar% {
    operm <- scantwo(cross.no2, method="em", n.perm=1,perm.strata=strata)
}
pen.em <- calc.penalties(operm.em)
)
try(
registerDoParallel(slurmcore)
operm.uns <- foreach(i=50,
  .combine=c,.packages = "qtl") %dopar% {
    operm <- scantwo(cross.no2, method="em")
}
pen.em <- calc.penalties(operm.uns)
)
print('Done scanning. Saving...')
save.image(paste('QTLmap.Rsave',sep=''))
## probs_map <- interp_genoprob(pr, map)
