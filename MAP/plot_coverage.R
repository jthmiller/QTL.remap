#!/bin/R

setwd('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/')
source('control_file.R')

for (X in 1:2){
  load(paste(popdir,'/chr',X,'_',outname,'.QTLmap.Rsave',sep=''))

  ## Plot fitration step
  png(file.path(popdir,paste(X,'_pval.png',sep='')))
  hist.geno(gt.missing$P.value)
  abline(v=log10(cutoff))
  dev.off()

  png(file.path(popdir,paste(X,'_pos.png',sep='')))
  par(mfrow=c(4,1),mar=c(1,2,3,1))
  plot.geno(gt,gen.main=paste('All Mapped Markers on Chr',X))
  plot.geno(gt.missing,gen.main='Markers with < 5 missing GT, Grandparent Confirmed markers in Green',par.gt=par.confirm.marks)
  abline(h=log10(cutoff),col='red')
  plot.geno(gt.pval,gen.main=paste('Filter distortion > ',cutoff))
  plot.geno(gt.fin,gen.main='Final Markers')
  dev.off()
}
