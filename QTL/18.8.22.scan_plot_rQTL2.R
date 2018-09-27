#!/bin/bash
setwd('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/QTL/')
source('control_file.R')

cross.18 <- reconst(X=chrms,pop=popq,temp.dir=popdir)
#cross.18 <- reconst(X,pop=popq, dir=qtldir)
print('Writing the merged chromosome markers to rQTL format')
write.cross(cross.18,filestem=paste(qtldir,'BACKUP.QTL_chr.QTLmap',sep=''),format="csv",chr=X)

### rQTL2
cross2 <- convert2cross2(cross.18)
map <- insert_pseudomarkers(cross2$gmap, step=1)
pr <- calc_genoprob(cross2, map, err=ers, cores=slurmcore)
pr <- clean_genoprob(pr)
apr <- genoprob_to_alleleprob(pr)
## Scan for QTLs in each pop
pr.sub <- calc_genoprob(subset(cross2,chr=1:24[-2]), map, err=ers, cores=slurmcore)
perms <- scan1perm(pr.sub,cross2$pheno, model="binary", cores=slurmcore,n_perm=2000,perm_strata=cross2$pheno[,1])
cutoff <- summary(perms)['0.05',]
perms.unstrat <- scan1perm(pr,cross2$pheno, model="binary", cores=0,n_perm=2000)
### Figure out why permutations is usually returning the max LOD from each Chr
cutoff.us <- summary(perms.unstrat)['0.05',]
out_bin <- scan1(pr,cross2$pheno[,1], model="binary", cores=slurmcore)
out_coef <- scan1coef(pr[,2],cross2$pheno[,1],model = 'binary', maxit = 1000,
  contrasts=cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(0, 1, 0)))
## uses a lod10 likelyhood function that results in a posterior dist of the QTL location
## bayes_int(out_bin, map, lodcolumn=1, prob=0.95, threshold=cutoff)
single <- find_peaks(out_bin, map, threshold=cutoff, peakdrop=3)
single.us <- find_peaks(out_bin, map, threshold=cutoff.us, peakdrop=3)
str <- as.character(unique(single$chr))
uns <- as.character(unique(single.us$chr))
uns <- lapply(split(single.us, single.us$chr),function(r){r[which.max(r$lod),4:5]})
str <- lapply(split(single, single$chr),function(r){r[which.max(r$lod),4:5]})
uns <- as.data.frame(cbind(lod=unlist(lapply(uns, "[[", 2)),pos=unlist(lapply(uns, "[[", 1))))
str <- as.data.frame(cbind(lod=unlist(lapply(str, "[[", 2)),pos=unlist(lapply(str, "[[", 1))))
uns <- uns[order(as.numeric(uns$lod),decreasing=T),]
str <- str[order(as.numeric(str$lod),decreasing=T),]

## Fit QTLs from rQTL2 with rQTL and genoprobs (Haley knott)
cross <- calc.genoprob(cross.18, step=1,error.prob=ers,map.function='kosambi')
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
