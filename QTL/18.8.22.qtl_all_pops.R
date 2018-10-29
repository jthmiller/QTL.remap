#!/bin/bash
#source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')

NBH <- new.env()
load('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NBH/REMAPS/QTLmap.Rsave', envir=NBH)
ELR <- new.env()
#ELR <- load('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/ELR/REMAPS/QTLmap.Rsave', envir=ELR)
NEW <- new.env()
load('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NEW/REMAPS/QTLmap.Rsave', envir=NEW)

cbind(
  NBH=summaryMap(NBH$cross.18)$length,
  NEW=summaryMap(NBH$cross.18)$length,
  ELR=summaryMap(NBH$cross.18)$length,

)

sapply(1:24,function(X){
    cmdist <-
}
)

sim.geno(cross, n.draws=16, step=0, off.end=0, error.prob=0.0001,
map.function=c("haldane","kosambi","c-f","morgan"),
stepwidth=c("fixed", "variable", "max"))



dir <- '/home/jmiller1/QTL_Map_Raw/popgen/rQTL'

cross.18 <- read.cross(format='csv',dir=dir,BC.gen=0, F.gen=2,
   file='remap_outBACKUP.QTL_chr.QTLmap.csv',
   geno=c('AA','AB','BB'),alleles=c("A","B"))





NBH$cross.18$pheno$ID <- paste('NBH',NBH$cross.18$pheno$ID,sep='')
NEW$cross.18$pheno$ID <- paste('NEW',NEW$cross.18$pheno$ID,sep='')

trym <- combineMap(NBH$cross.18,NEW$cross.18, keep.all = FALSE, merge.by = "marker" )
tryg <- combineMap(NBH$cross.18,NEW$cross.18, keep.all = TRUE, merge.by = "genotypes")




test3 <- mstmap.cross(cross.18, chr= 1:24, id = "ID", bychr = TRUE, suffix = "numeric",
     anchor = FALSE, dist.fun = "kosambi", objective.fun = "COUNT",
     p.value = 1e-06, noMap.dist = 15, noMap.size = 0, miss.thresh = .75,
     mvest.bc = FALSE, detectBadData = FALSE, return.imputed = FALSE,
     trace = FALSE)

test4 <- mstmap(test3, bychr = TRUE, id = "ID", dist.fun = "kosambi", p.value = 2, trace = FALSE,miss.thresh = .75)



BC.gen=2, F.gen=3


try <- reduce2grid(NBH$cross.18)


hyper <- calc.genoprob(hyper, step=2)

try <- argmax.geno(NBH$cross.18, step=0, off.end=0, error.prob=0.002,map.function="kosambi",stepwidth="fixed")
#try <- pull.argmaxgeno(try)
try <- reduce2grid(try)



hypersub <- reduce2grid(hyper)
## Not run: out <- scanone(hypersub)
plot(out, incl.markers=FALSE)
## End(Not run)







interp_genoprob(probs, map, cores = 1)

NBH.cr <- convert2cross2(NBH$cross.18)
NEW.cr <- convert2cross2(NBH$cross.18)

map.h <- insert_pseudomarkers(NBH.cr$gmap, step=1)
map.w <- insert_pseudomarkers(NEW.cr$gmap, step=1)

pr.h <- calc_genoprob(NBH.cr, map.h, err=0.002, cores=slurmcore)
pr.w <- calc_genoprob(NEW.cr, map.w, err=0.002, cores=slurmcore)

probs_map.h <- interp_genoprob(pr.h, map.h, cores = slurmcore)
probs_map.w <- interp_genoprob(pr.w, map.h, cores = slurmcore)




out_bin <- scan1(probs_map.h,NBH.cr$pheno, model="normal", cores=slurmcore)

find_peaks(out_bin, map.h, threshold=4, drop=1.5)



probs_map.w <- interp_genoprob(pr.w, map.w, cores = slurmcore)








pr <- clean_genoprob(pr)
apr <- genoprob_to_alleleprob(pr)



iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2"))

probs <- calc_genoprob(iron, iron$gmap, error_prob=0.002)

# you generally wouldn't want to do this, but this is an illustration
map <- insert_pseudomarkers(iron$gmap, step=1)
probs_map <- interp_genoprob(probs, map)







### These markers worked in NBH,NEW. Maybe ELR?
marks <- unique(c(markernames(NBH$cross.18),markernames(NBH$cross.18)))

cross.18 <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
                format='csvr', geno=c(1:3),estimate.map=FALSE)

cross.18 <- drop.markers(cross.18,markernames(cross.18)[!markernames(cross.18) %in% marks])

gt <- geno.table(cross.18)
cutoff <- 1.0e-4
cross.18 <- drop.markers(cross.18,rownames(gt[gt$P.value<cutoff,]))
