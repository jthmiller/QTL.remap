#!/bin/R
### Map QTLs 3 of 3
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/pop_control_file.R')

load(paste(popdir,'/chr',X,'.QTLmap.Rsave',sep=''))

cross.18 <- read.cross(format='csv',dir=popdir,
   file=paste('chr',X,'.QTLmap.csv',sep=''),
   geno=c('AA','AB','BB'),alleles=c("A","B"))

marker.warning()

print('Re-setimating map from filtered data on')

cross.18 <- orderMarkers(cross.18,chr=X,window=5,use.ripple=T,
  error.prob=ers, map.function='kosambi',sex.sp=F,maxit=1000,tol=1e-3)

POS.map.18 <- est.map(cross.18,error.prob=ers,map.function="kosambi",n.cluster=12, chr=X)

cross.18 <- replace.map(cross.18, POS.map.18)

print('Done mapping..')
print(summary(pull.map(cross.18))[X,])

print('Re-writing the markers to rQTL format')
write.cross(cross.18,filestem=paste(plotdir,'chr',X,'.QTLmap',sep=''),format="csv",chr=X)

print('saving...')
rm(cross.18)
save.image(paste('chr',X,'.QTLmap.Rsave',sep=''))
