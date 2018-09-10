#!/bin/bash
#load(paste(popdir,'/chr',X,'.QTLmap.Rsave',sep=''))

cross.18 <- reconst(X,pop='NBH',out=out)
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
#out_blup <- scan1blup(pr,cross2$pheno)
max_pos <- rownames(max(out_bin, map['18']))
fit <- fit1(pr[['18']][,,max_pos], cross2$pheno, model="binary")

find_peaks(out_bin, map, threshold=cutoff, peakdrop=1.8, drop=1.5)
find_peaks(out_bin, map, threshold=cutoff, peakdrop=1, prob=0.95)

bayes_int(out_bin, map, lodcolumn=1, chr=2, prob=0.95)


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
## For NBH
#      Pheno
# 0.2   3.56
# 0.05  4.33
print(find_peaks(out_bin, map, threshold=3.56, peakdrop=1.8, drop=1.5))




## Scan for QTL with rQTL

print('Scanning for a single QTL')
GP <- calc.genoprob(cross.18, step=2.5)
GP <- sim.geno(GP,n.draws=1000, step=2, err=0.02)
scanQTL <- scanone(GP, pheno.col=1, model="binary", method="hk")

print('Done scanning. Saving...')
save.image(paste('chr',X,'.QTLmap.Rsave',sep=''))
print(paste('done with chrom',X,'in pop',pop))
