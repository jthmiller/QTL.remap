#!/bin/bash
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/pop_control_file.R')
require(qtl2,lib.loc='/share/apps/rmodules')

pheno.all <- read.table('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/metadata/ALL_phenotype_Dist.txt',header=T)
pheno.all$pheno_all[which(pheno.all$pheno_all<2)] <- 0
pheno.all$pheno_all[which(pheno.all$pheno_all>1)] <- 1
pheno.all$sex <- 0




pheno.all$ID <- paste('NG',1:length(pheno.all$sex),sep='_')
table(pheno.all$pheno_all)-table(cross.18$pheno$Pheno)


chr18.NBH <- read.cross(format='csv',dir=file.path(basedir,'rQTL'),
  file=paste(pop,'chr',X,'.QTLmap.csv',sep=''),
  geno=c('AA','AB','BB'),alleles=c("A","B"))

load(paste('chr',X,'.QTLmap.Rsave',sep=''))

## QTL2 plotting

cross2 <- convert2cross2(chr18.NBH)
map <- insert_pseudomarkers(cross2$gmap, step=1)
pr <- calc_genoprob(cross2, map, err=ers, cores=0)
pr <- clean_genoprob(pr)
apr <- genoprob_to_alleleprob(pr)

## Get genotype probabilities for inserted psuedomarkers
## Two populations can be compared at CM position rather than markers
probs_map <- interp_genoprob(pr, map)



## Scan for QTLs
perms <- scan1perm(pr,cross2$pheno, model="binary", cores=0,n_perm=1000,perm_strata=cross2$pheno)
cutoff <- summary(perms)['0.05',]

out_bin <- scan1(pr,cross2$pheno, model="binary", cores=0)
out_coef <- scan1coef(pr,cross2$pheno,model = 'binary')
out_blup <- scan1blup(pr,cross2$pheno)

find_peaks(out_bin, map, threshold=cutoff, peakdrop=1.8, drop=1.5)
bayes_int(out_bin, map, lodcolumn=1, chr=18, prob=0.95)
find_peaks(out_bin, map, threshold=cutoff, peakdrop=1, prob=0.95)

max_pos <- rownames(max(out_bin, map['18']))
fit <- fit1(pr[['18']][,,max_pos], cross2$pheno, model="binary")




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
