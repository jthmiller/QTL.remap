#!/bin/bash
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')

cross.18 <- reconst(X=chrms,pop=popq,temp.dir=popdir)
print('Writing the merged chromosome markers to rQTL format')
write.cross(cross.18,filestem=paste(qtldir,'BACKUP.QTL_chr.QTLmap',sep=''),format="csv",chr=X)

### Ids for strata permutations
cross.18$pheno$stata <- gsub('[0-9]','',cross.18$pheno$ID)

### Error rate and genoprobs
ers <- 0.002
cross.18 <- calc.genoprob(cross.18, step=1,error.prob=ers,map.function='kosambi')

### Cross object with only genotyped individuals for model comp
indy <- Z$pheno$ID[grep('ind',Z$pheno$ID)]
cross.gt <- subset(Z,ind=as.character(indy))

#### Transform the phenotype for model flexibility. Imputations are better for selective genotyping than HK
Z <- transformPheno(cross.18, pheno.col=2, transf=nqrank)
Z <- sim.geno(Z,error.prob=0.002)





#### Nonparametric scan
scan.np.em <- scanone(cross.18, method='em', model="np", pheno.col=2, maxit=5000)
### 2 part model
scan.2p.em <- scanone(cross.18, method='em', model="2part", pheno.col=1, maxit=5000)
### Binary model with EM
scan.bin.em <- scanone(cross.18, method="em",model='binary', pheno.col=1, maxit=5000)
### Binary model with marker regression
scan.bin.mr <- scanone(cross.18, method="mr",model='binary',pheno.col=1)
### Binary model with haley knott (heavy inflation of LOD)
scan.bin.hk <- scanone(cross.18, method="hk",model='binary',pheno.col=1)
### Normal scan on transformed phenotype
scan.norm.em <- scanone(Z, method="em",model='normal',maxit=500, pheno.col=2)
### Normal scan on transformed phenotype w/Extended haley knott (better for selective/missing genos at non-gt'ed ind)
scan.norm.ehk <- scanone(Z, method="ehk",model='normal',maxit=500, pheno.col=2)
### Normal scan on transformed phenotype fast haley knott (not robust to missing data. LOD inflation)
scan.norm.hk <- scanone(Z, method="hk",model='normal',maxit=500, pheno.col=2)
### Imputation on transformed phenotype
scan.norm.imp <- scanone(Z, method="imp",model='normal',maxit=500, pheno.col=2)
### Imputation on transformed phenotype
scan.norm.imp.gt <- scanone(cross.gt, method="imp",model='normal',maxit=500, pheno.col=2)
### Imputation on un-transformed scores
scan.norm.imp.05 <- scanone(cross.18, method="imp",model='normal',maxit=500, pheno.col=2)
### EM scan without ungenotyped
scan.norm.em.gt <- scanone(cross.gt, method="em",model='normal',maxit=500, pheno.col=2)



## Find LOD thresholds to get a null dist. of
### remove markers that are exactly the same to speed up (same results)
dups <- findDupMarkers(cross.18, exact.only=FALSE, adjacent.only=FALSE)
cross.18 <- drop.markers(cross.18,unlist(dups))

perms.np.em <- scanone(cross.18, model="np", pheno.col=2, n.perm=2000, method='em',perm.strata=cross.18$pheno$stata, n.cluster=24)
perms.2p.em <- scanone(cross.18, model="2part",n.perm=2000,perm.strata=cross.18$pheno$stata,pheno.col=2, n.cluster=24)
perms.bin.em <- scanone(cross.18, method="em",model='binary',maxit=2000,n.perm=500,perm.strata=cross.18$pheno$stata,pheno.col=1, n.cluster=24)
perms.bin.mr <- scanone(cross.18, method="mr",model='binary', n.perm=2000, perm.strata=cross.18$pheno$stata, n.cluster=24)
perms.bin.em <- scanone(cross.18, method="hk",model='binary',maxit=2000,n.perm=500,perm.strata=cross.18$pheno$stata,pheno.col=1, n.cluster=24)
perms.norm.em <- scanone(Z, method="em",model='normal',maxit=500,n.perm=2000,perm.strata=cross.18$pheno$stata,pheno.col=2, n.cluster=24)
perms.norm.ehk <- scanone(Z, method="ehk",model='normal',maxit=500,n.perm=2000,perm.strata=cross.18$pheno$stata,pheno.col=2, n.cluster=24)
perms.norm.hk <- scanone(Z, method="hk",model='normal',maxit=500,n.perm=2000,perm.strata=cross.18$pheno$stata,pheno.col=2, n.cluster=24)
perms.norm.imp <- scanone(Z, method="imp",model='normal',n.perm=2000,perm.strata=cross.18$pheno$stata,pheno.col=2, n.cluster=24)
perms.norm.imp.2 <- scanone(Z, method="imp",model='normal',chr=-2,n.perm=2000,perm.strata=cross.18$pheno$stata,pheno.col=2,n.cluster=24)

### Multi-QTL models
th <- summary(perms.norm.imp)[1,]
norm.qtl <- summary(scan.norm.imp, perms=perms.norm.imp, alpha=0.05)
qtl.uns <- makeqtl(Z, chr=norm.qtl$chr, pos=norm.qtl$pos)
full <- stepwiseqtl(Z, additive.only=T, method="imp", pheno.col=2, scan.pairs=T)






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
