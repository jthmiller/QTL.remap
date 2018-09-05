#!/bin/R
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/pop_control_file.R')
load(paste('chr',X,'.QTLmap.Rsave',sep=''))

print('Re-write the markers to rQTL formate')
write.cross(cross.18,filestem=paste(plotdir,'chr',X,'.QTLmap',sep=''),format="csv",chr=X)

print('saving...')
save.image(paste('chr',X,'.QTLmap.Rsave',sep=''))

## Scan for QTL

print('scanning for a single QTL')
GP <- calc.genoprob(cross.18, step=2.5)

GP <- sim.geno(GP,n.draws=1000, step=2, err=0.02)

scanQTL <- scanone(GP, pheno.col=1, model="binary", method="hk")

print('saving...')
save.image(paste('chr',X,'.QTLmap.Rsave',sep=''))
print(paste('done with chrom',X,'in pop',pop))
