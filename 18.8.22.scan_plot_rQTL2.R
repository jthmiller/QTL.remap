#!/bin/bash
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/pop_control_file.R')
require(qtl2,lib.loc='/share/apps/rmodules')

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
out_bin <- scan1(pr,cross2$pheno, model="binary", cores=0)
out_coef <- scan1coef(pr,cross2$pheno,model = 'binary')

find_peaks(out_bin, map, threshold=4, peakdrop=1.8, drop=1.5)
bayes_int(out_bin, map, lodcolumn=1, chr=18, prob=0.95)

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
