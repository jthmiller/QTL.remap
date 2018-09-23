#!/bin/R
### Map QTLs 2 of 3
setwd('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/')
source('control_file.R')

cross.18 <- read.cross(format='csv',dir=popdir,
   file=paste('chr',X,'_',outname,'.QTLmap.csv',sep=''),
   geno=c('AA','AB','BB'),alleles=c("A","B"))

zero.map <- shiftmap(pull.map(cross.18))
cross.18 <- replacemap(cross.18, zero.map)

## keep genotypes for QTL markers
tokeep <- readLines(paste(popdir,'/chr',X,'_',outname,'.keepmarkers.csv',sep=''))
if(!tokeep==''){
  gi <- pull.geno(cross.18)[,tokeep]
} else {
  tokeep <- NULL
}

marker.warning()
print(summary(pull.map(cross.18))[as.character(X),])
print('Dropping 2.5% of markers that inflate the map. Takes a long time...')

## Drop one marker, p is proportion of worst markers to drop
system.time(
  cross.18 <- dropone.par(cross=cross.18, prop=0.025,chr=X, maxit=2000, map.function = 'kosambi',
    length.imp = 1, LOD.imp = 0, error.prob=0.05, sex.sp = F, verbose=F, parallel=T, cores=slurmcore)
)

marker.warning()
return.dropped.markers()
marker.warning()

print('Re-order markers')
system.time(
cross.18 <- orderMarkers(cross.18,chr=X,window=5,use.ripple=T,
  error.prob=ers, map.function='kosambi',sex.sp=F,maxit=2000,tol=1e-3)
)

print('Re-estimating the map')
system.time(
POS.map.18 <- est.map(cross.18,error.prob=0.05,map.function="kosambi", chr=X,maxit=2000)
)
cross.18 <- replace.map(cross.18, POS.map.18)

## Error rate
print('determine error rate for last round of mapping')
ers <- er.rate(cross=cross.18,cpus=slurmcore,maxit=1000)

print(paste(pop,'error rate for chromosome',X,'is',ers))
system(paste('echo',pop,X,'_',outname,ers,'>> /home/jmiller1/QTL_Map_Raw/popgen/rQTL/remap_out/genotyping_error_rate.txt'))

print('saving...')
save.image(paste(popdir,'/chr',X,'_',outname,'.QTLmap.Rsave',sep=''))

### Write Map information
system(paste('echo ',pop,X,outname,'>> /home/jmiller1/QTL_Map_Raw/popgen/rQTL/remap_out/map.txt',sep='\t'))
line <- unlist(summary(pull.map(cross.18))[as.character(X),])
write('mar  length  avesp  max',file="/home/jmiller1/QTL_Map_Raw/popgen/rQTL/remap_out/map.txt",append=TRUE)
write(line,file="/home/jmiller1/QTL_Map_Raw/popgen/rQTL/remap_out/map.txt",append=TRUE)

### Write cross to file
print('Writing the markers to rQTL format')
write.cross(cross.18,filestem=paste(popdir,'/chr',X,'_',outname,'.QTLmap',sep=''),format="csv",chr=X)
