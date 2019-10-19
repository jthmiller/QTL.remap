#!/bin/R
### Map QTLs 1 of 3
#debug.cross <- T
#source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")

#screen -m srun -t 12:00:00 -p high --mem=60000 --pty R

library('qtl')
pop <- 'NBH'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

################################################################################
fl <- file.path(mpath,'nbh.mapped.1_24.csv')

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

mrks <- sapply(1:24,function(X){
  n.missing <- nmissing(subset(cross, chr=X), what="mar")
  wts <- -log( (n.missing+1) / (nind(cross)+1) )
  pickMarkerSubset(pull.map(cross)[[as.character(X)]], 0.75, wts)
})
# subset of markers on chr 4 spaced >= 5 cM, with weights = -log(prop'n missing)

dwnsmpl <- pull.markers(cross,unlist(mrks))

gg <- sim.geno(cross,step=1)
gg_step2 <- reduce2grid(gg)

cgb <- calc.genoprob(cross, step=1)
cgb_step2 <- reduce2grid(cgb)

################################################################################
# MQM ##########################################################################
## a possible rule of thumb may be to set minprob to the percentage of data missing

crossaug <- mqmaugment(dwnsmpl, minprob=0.75, strategy="impute")  # Augmentation
result <- mqmscan(crossaug)    ## scan show LOD interval of the QTL on chr 3

#multitoset <- find.markerindex(maug, "GH.117C")
setcofactors <- mqmsetcofactors(maug, cofactors=multitoset)
mqm_co1 <- mqmscan(maug, setcofactors)

mq.add <- mqmscan(crossaug, cofactors=NULL, pheno.col = 5,
      model="additive", forceML=FALSE,
      cofactor.significance=0.02, em.iter=1000,
      window.size=25.0, step.size=5.0,
      logtransform = FALSE, estimate.map = FALSE,
      plot=FALSE, verbose=FALSE, outputmarkers=TRUE,
      multicore=TRUE, batchsize=10, n.clusters=6,
      test.normality=FALSE,off.end=0)

mq.dom <- mqmscan(crossaug, cofactors=NULL, pheno.col = 5,
      model="dominance", forceML=FALSE,
      cofactor.significance=0.02, em.iter=1000,
      window.size=25.0, step.size=5.0,
      logtransform = FALSE, estimate.map = FALSE,
      plot=FALSE, verbose=FALSE, outputmarkers=TRUE,
      multicore=TRUE, batchsize=10, n.clusters=6,
      test.normality=FALSE,off.end=0)

save.image(file.path(mpath,'scans.nbh.rsave'))
################################################################################
