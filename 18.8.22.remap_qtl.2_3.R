#!/bin/R
### Map QTLs 2 of 3
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/pop_control_file.R')

cross.18 <- read.cross(format='csv',dir=file.path(basedir,'rQTL'),
   file=paste(pop,'chr',X,'.QTLmap.csv',sep=''),
   geno=c('AA','AB','BB'),alleles=c("A","B"))

marker.warning()

print('Dropping 2.5% of markers that inflate the map. Takes a long time...')
## Drop one marker, p is proportion  of worst markers to drop
cross.18 <- dropone.par(cross.18,p=0.025,chr=X,maxit=3,
  sex.sp = F,verbose=F,parallel=T)

marker.warning()

print('determine error rate')
ers <- er.rate(cross.18)
print(paste(ers,' error rate'))

print('Writing the markers to rQTL format')
write.cross(cross.18,filestem=paste(plotdir,'chr',X,'.QTLmap',sep=''),format="csv",chr=X)

save.image(paste('chr',X,'.QTLmap_2.Rsave',sep=''))
