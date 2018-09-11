#!/bin/R
### Map QTLs 2 of 3
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/pop_control_file.R')

cross.18 <- read.cross(format='csv',dir=popdir,
   file=paste('chr',X,'.QTLmap.csv',sep=''),
   geno=c('AA','AB','BB'),alleles=c("A","B"))

marker.warning()

print('Dropping 2.5% of markers that inflate the map. Takes a long time...')
## Drop one marker, p is proportion  of worst markers to drop
cross.18 <- dropone.par(cross=cross.18,p=0.025,chr=X,maxit=1,
  error.prob=0.02,sex.sp = F,verbose=F,parallel=T)

marker.warning()

print(summary(pull.map(cross.18))[as.character(X),])

print('determine error rate')
ers <- er.rate(cross.18)
print(paste(ers,' error rate'))

## fix for those that do not have below thr error
print('dropping markers by error lod')
cross.18 <- drop.errlod(cross.18,lod=4,ers=0.02)

print('Writing the markers to rQTL format')
write.cross(cross.18,filestem=paste(popdir,'/chr',X,'.QTLmap',sep=''),format="csv",chr=X)

print('saving...')
rm(cross.18)
save.image(paste(popdir,'/chr',X,'.QTLmap.Rsave',sep=''))
