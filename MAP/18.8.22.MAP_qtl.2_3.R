#!/bin/R
### Map QTLs 2 of 3
setwd('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/')
source('control_file.R')

cross.18 <- read.cross(format='csv',dir=popdir,
   file=paste('chr',X,'.QTLmap.csv',sep=''),
   geno=c('AA','AB','BB'),alleles=c("A","B"))

marker.warning()

print('Dropping 2.5% of markers that inflate the map. Takes a long time...')
## Drop one marker, p is proportion  of worst markers to drop
cross.18d <- dropone.par(cross=cross.18,prop=0.025,chr=X, maxit=2 ,map.function = 'kosambi',
  length.imp = 1, LOD.imp = 0,error.prob=0.03,sex.sp = F,verbose=T,parallel=T,cores=slurmcore)

marker.warning()

print(summary(pull.map(cross.18))[as.character(X),])

## fix for those that do not have below thr error
print('dropping markers by error lod')
cross.18 <- drop.errlod(cross.18,lod=4,ers=ers)

## Error rate
print('determine error rate for last round of mapping')
ers <- er.rate(cross.18)
print(paste(pop,'error rate for chromosome',X,'is',ers,))
system(paste('echo',pop,X,ers,'>> /home/jmiller1/QTL_Map_Raw/popgen/rQTL/remap_out/genotyping_error_rate.txt'))

print('saving...')
save.image(paste(popdir,'/chr',X,'.QTLmap.Rsave',sep=''))

clean <- ls()[!ls() %in% c('X','ers','popdir','cross.18','marker.warning','marker.density')]
rm(list=clean)

print('Re-estimating the map')
POS.map.18 <- est.map(cross.18,error.prob=ers,map.function="kosambi",n.cluster=12, chr=X,maxit=1000)
cross.18 <- replace.map(cross.18, POS.map.18)

print('Writing the markers to rQTL format')
write.cross(cross.18,filestem=paste(popdir,'/chr',X,'.QTLmap',sep=''),format="csv",chr=X)
