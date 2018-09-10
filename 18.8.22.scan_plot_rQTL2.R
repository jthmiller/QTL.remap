#!/bin/bash
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/qtl_control_file.R')

slurmcore <- detectCores()
## rQTL2 to viz single QTL models
## rQTL for 2 qtls for:
### permutations
### multi QTL models
out <- file.path(qtldir,'out')
cross <- reconst(X,pop='NBH',out=out)
print('Writing the markers to rQTL format')
write.cross(cross.18,filestem=paste(plotdir,'BACKUP.QTL_chr.QTLmap',sep=''),format="csv",chr=X)



### rQTL2
ers <- 0.02
cross2 <- convert2cross2(cross.18)
map <- insert_pseudomarkers(cross2$gmap, step=1)
pr <- calc_genoprob(cross2, map, err=ers, cores=0)
pr <- clean_genoprob(pr)
apr <- genoprob_to_alleleprob(pr)

## Get genotype probabilities for inserted psuedomarkers
## Two populations can be compared at CM position rather than markers
## probs_map <- interp_genoprob(pr, map)

## Scan for QTLs
perms <- scan1perm(pr,cross2$pheno, model="binary", cores=0,n_perm=2000,perm_strata=cross2$pheno)
cutoff <- summary(perms)['0.05',]
perms.unstrat <- scan1perm(pr,cross2$pheno, model="binary", cores=0,n_perm=2000)
cutoff.us <- summary(perms.unstrat)['0.05',]
out_bin <- scan1(pr,cross2$pheno, model="binary", cores=0)
out_coef <- scan1coef(pr,cross2$pheno,model = 'binary')
max_pos <- rownames(max(out_bin, map['18']))
fit <- fit1(pr[['18']][,,max_pos], cross2$pheno, model="binary")

### Single QTLs
bayes_int(out_bin, map, lodcolumn=1, prob=0.95, threshold=cutoff)
a <- find_peaks(out_bin, map, threshold=cutoff, peakdrop=3, prob=0.95)

## Make QTL in rQTL1 with rQTL2 locarions
GP <- calc.genoprob(cross.18, step=2.5)
GP <- sim.geno(GP,n.draws=1000, step=2, err=0.02)
### This uses single qtls from rQTL2
qtl <- makeqtl(GP, chr=a$chr, pos=a$pos, what="prob")
out.fq <- fitqtl(GP, qtl=qtl, method="hk")
rqtl <- refineqtl(GP, qtl=qtl, formula=y~Q1+Q2+Q3, verbose=FALSE)
stepout1 <- stepwiseqtl(GP,qtl=rqtl, additive.only=TRUE, max.qtl=6, verbose=FALSE)




find_peaks(out_bin, map, threshold=cutoff, peakdrop=3, drop=1.5)




sapply(X,function(X){
  bayes_int(out_bin, map, lodcolumn=1, chr=X, prob=0.95)
  }
)

interp_genoprob "to get two sets onto the same map for comparison purposes"


### Plot max lod on each chr #####################################
png(paste(plotdir,'/maxlod.qtl.png',sep='' ),width =1000)
par(mar=c(5.1, 4.1, 1.1, 1.1))
plot(out_bin, map, lodcolumn=1, col="slateblue", ylim=c(0, 10))
dev.off()

### Perform permutation to determine significance ####
operm <- scan1perm(pr, bin_pheno, n_perm=10000, cores=0)
summary(operm, alpha=c(0.2, 0.05))

print(find_peaks(out_bin, map, threshold=3.56, peakdrop=1.8, drop=1.5))





print('Scanning for a single QTL')
GP <- calc.genoprob(cross.18, step=2.5)
GP <- sim.geno(GP,n.draws=1000, step=2, err=0.02)
scanQTL <- scanone(GP, pheno.col=1, model="binary", method="hk")

print('Done scanning. Saving...')
save.image(paste('chr',X,'.QTLmap.Rsave',sep=''))
print(paste('done with chrom',X,'in pop',pop))
