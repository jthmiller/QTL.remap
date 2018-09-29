#!/bin/bash
#setwd('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/QTL/')
#source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')

#cross.18 <- reconst(X=chrms,pop=popq,temp.dir=popdir)
#print('Writing the merged chromosome markers to rQTL format')
#write.cross(cross.18,filestem=paste(qtldir,'BACKUP.QTL_chr.QTLmap',sep=''),format="csv",chr=X)

cross.18 <- read.cross.jm(file=file.path(popdir,'tempout'),format='csv',
  geno=c(1:3),estimate.map=FALSE)

cross.18$pheno$stata <- gsub('[0-9]','',cross.18$pheno$ID)

#### Nonparametric scan
cross.18 <- calc.genoprob(cross.18, step=1,error.prob=ers,map.function='kosambi')
scan.np.em <- scanone(cross.18, model="np", pheno.col=2 ,method='em')
perms.np.em <- scanone(cross.18, model="np", pheno.col=2, n.perm=1000, method='em',perm.strata=cross.18$pheno$stata)
summary(scan.np.em, perms=perms.np.em, alpha=0.05, pvalues=TRUE)


### no good way to fit a model with nonparamtric/binary data with selective genotyping
### There are ways to scan for single QTLs- but all multi models use hk (not good for Sel Gt) and imp (not coded for binary or np models)


scan.np.hk <- scanone(cross.18, model="np",pheno.col=2, method='hk')

plot(scan.np, ylab="LOD score", alternate.chrid=TRUE)

pval.np <- summary(operm.np, 0.05)
summary(out.np, perms=operm.np, alpha=0.05, pvalues=TRUE)

#### Two-part model
out.2p <- scanone(listeria, model="2part", upper=TRUE, pheno.col="logsurv")

# Plot the two-part model results
plot(out.2p, lodcolumn=1:3, ylab="LOD score", alternate.chrid=TRUE)

# Permutation test for the two-part model
operm.2p <- scanone(listeria, model="2part", upper=TRUE,pheno.col="logsurv", n.perm=1000, perm.Xsp=TRUE)

# LOD thresholds by 2-part model
summary(operm.2p, alpha=0.05)

# Summary of the significant loci
summary(out.2p, perms=operm.2p, alpha=0.05, pvalues=TRUE)

# Alternate summary
summary(out.2p, perms=operm.2p, alpha=0.05, pvalues=TRUE,
        format="allpeaks")




## Binary phenotype scan

cross.18$pheno$stata <- gsub('[0-9]','',cross.18$pheno$ID)
cross.18 <- calc.genoprob(cross.18, step=1,error.prob=ers,map.function='kosambi')
perms.em <- scanone(cross.18, method="em",model='binary',maxit=4000,
  n.perm=100,n.cluster=10,perm.strata=cross.18$pheno$stata)
scan.em <- scanone(cross.18, method="em",model='binary',maxit=4000,n.cluster=slurmcore)


scan.ehk <- scanone(cross.18, method="ehk",model='binary',maxit=4000,n.cluster=1,perm.strata=cross.18$pheno$stata)



lod.05 <- summary(perms.em, alpha=0.05)
binary.QTLs <- summary(out.bin, perms=operm.bin, alpha=0.05, pvalues=TRUE)


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
