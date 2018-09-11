#!/bin/bash
## Permutations take a long time
## Map QTLs 4
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/qtl_control_file.R')
slurmcore <- detectCores()
## rQTL2 to viz single QTL models
## rQTL for 2 qtls for:
### permutations
### multi QTL models
out <- file.path(qtldir,'out')
cross.18 <- reconst(X,pop='NBH',out=out)
print('Writing the merged chromosome markers to rQTL format')
write.cross(cross.18,filestem=paste(qtldir,'BACKUP.QTL_chr.QTLmap',sep=''),format="csv",chr=X)

### rQTL2
ers <- 0.02
cross2 <- convert2cross2(cross)
map <- insert_pseudomarkers(cross2$gmap, step=1)
pr <- calc_genoprob(cross2, map, err=ers, cores=slurmcore)
pr <- clean_genoprob(pr)
apr <- genoprob_to_alleleprob(pr)
## Scan for QTLs in each pop
perms <- scan1perm(pr,cross2$pheno, model="binary", cores=slurmcore,n_perm=2000,perm_strata=cross2$pheno[,1])
cutoff <- summary(perms)['0.05',]
#perms.unstrat <- scan1perm(pr,cross2$pheno, model="binary", cores=0,n_perm=2000)
#cutoff.us <- summary(perms.unstrat)['0.05',]
out_bin <- scan1(pr,cross2$pheno[,1], model="binary", cores=slurmcore)
## uses a lod10 likelyhood function that results in a posterior dist of the QTL location
bayes_int(out_bin, map, lodcolumn=1, prob=0.95, threshold=cutoff)
single <- find_peaks(out_bin, map, threshold=cutoff, peakdrop=2, prob=0.95)

out_coef <- scan1coef(pr[,2],cross2$pheno[,1],model = 'binary', maxit = 1000,
  contrasts=cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)))

## Fit QTLs from rQTL2 with rQTL and genoprobs (Haley knott)
cross <- calc.genoprob(cross.18, step=1,error.prob=0.02,map.function='kosambi')
qtl <- makeqtl(cross, chr=single$chr, pos=single$pos,what='prob')
qtl.rf <- refineqtl(cross, qtl=qtl, formula=y~Q1+Q2+Q3,method="hk",verbose=FALSE)
out.fit <- fitqtl(cross, qtl=qtl.rf, formula=y~Q1+Q2+Q3,method="hk",get.ests=T)

### After fitting single addative QTLs, search model space for addition addative
stepout1 <- stepwiseqtl(cross,qtl=qtl.rf,pheno.col=1,model="binary",
  additive.only=TRUE, max.qtl=6, verbose=FALSE)

### qtl1 with em
strata <- cross$pheno$pheno==1
names(strata) <- cross$pheno$ID
out_scan.em <- scanone(cross, pheno.col=1, model="binary", method="em",perm.strata=strata)
out_perm.em <- scanone(cross, pheno.col=1, model="binary", method="em",perm.strata=strata,n.perm=100)

save.image(paste('QTLmap.Rsave',sep=''))

cross.no2 <- subset(cross,chr=c(1,3:24))
registerDoParallel(slurmcore)
operm.hk <- foreach(i=500,
  .combine=c,.packages = "qtl") %dopar% {
    operm <- scantwo(cross.no2, method="hk", n.perm=1,
                       perm.strata=strata)
}
save.image(paste('QTLmap.Rsave',sep=''))

### qtl scan2
registerDoParallel(slurmcore)
operm.em <- foreach(i=500,
  .combine=c,.packages = "qtl") %dopar% {
    operm <- scantwo(cross.no2, method="em", n.perm=1,
                       perm.strata=strata)
}
pen <- calc.penalties(operm)

##contrasts > to get mean, additive effect, and dominance effect.
## contrasts matrix (arg provided)

print('Done scanning. Saving...')
save.image(paste('QTLmap.Rsave',sep=''))
## probs_map <- interp_genoprob(pr, map)
