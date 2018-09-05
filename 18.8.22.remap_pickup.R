#!/bin/R

source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/pop_control_file.R')

load(paste('chr',X,'.QTLmap.Rsave',sep=''))

print('Dropping 5% of markers that inflate the map. Takes a long time...')
print(paste('# of markers =',nmar(cross.18)))
## Drop one marker, p is proportion  of worst markers to drop
cross.18 <- dropone.par(cross.18,p=0.05,chr=X,maxit=2,
  sex.sp = F,verbose=F,parallel=T)

print('saving...')
save.image(paste('chr',X,'.QTLmap.Rsave',sep=''))

print('Re-setimating map from filtered data')
print(paste('# of markers =',nmar(cross.18)))
POS.map.18 <- est.map(cross.18,error.prob=0.002,map.function="kosambi",n.cluster=6, chr=X)
cross.18 <- replace.map(cross.18, POS.map.18)

print('Re-write the markers to rQTL formate')
write.cross(file=paste('chr',X,'.QTLmap',sep=''),format="csv",filestem=plotdir,chr=X)

print('saving...')
save.image(paste('chr',X,'.QTLmap.Rsave',sep=''))

## Scan for QTL

print('scanning for a single QTL')
GP <- calc.genoprob(cross.18, step=2.5)

GP <- sim.geno(GP,n.draws=1000, step=2, err=0.02)

scanQTL <- scanone(GP, pheno.col=1, model="binary", method="hk")

print('saving...')
save.image(paste('chr',X,'.QTLmap.Rsave',sep=''))
print(paste('done with chrom',X,'in pop',P))
