#!/bin/bash
#source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/debug.R')

NBH <- new.env()
load('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NBH/REMAPS/QTLmap.Rsave', envir=NBH)
ELR <- new.env()
ELR <- load('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/ELR/REMAPS/QTLmap.Rsave', envir=ELR)
NEW <- new.env()
load('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NEW/REMAPS/QTLmap.Rsave', envir=NEW)

cm.len <- cbind(NBH=summaryMap(NBH$cross.18)$length,
  NEW=summaryMap(NEW$cross.18)$length
#  ELR=summaryMap(ELR$cross.18)$length
)
cm.len.max <- apply(cm.len,1,max)
cm.len <- cm.len.max-cm.len[,1:2]
cm.len <- cm.len[1:24,]

NBH$cross.18$ID <- paste('NBH',NBH$cross.18$ID,sep='_')
NEW$cross.18$ID <- paste('NEW',NEW$cross.18$ID,sep='_')
NBH.x <- qtl::clean(NBH$cross.18)
NEW.x <- qtl::clean(NEW$cross.18)
cross <- c(NBH$cross.18,NEW$cross.18)

cross <- sim.geno(cross, n.draws=250, step=1, off.end=0,
  error.prob=0.01,map.function="kosambi",stepwidth="fixed")

cross.nbh <- subset(cross,ind=comp$pheno$cross1==1)
cross.new <- subset(cross,ind=comp$pheno$cross2==1)

scan.norm.imp.NBH <- scanone(cross.nbh, method="imp",model='normal',pheno.col=6)
scan.norm.imp.NEW <- scanone(cross.new, method="imp",model='normal',pheno.col=6)

#scan.norm.imp.NBH <- scanone(comp, method="imp",model='normal',pheno.col=6,ind.noqtl=comp$pheno$cross1==1)
#scan.norm.imp.NEW <- scanone(comp, method="imp",model='normal',pheno.col=6,ind.noqtl=comp$pheno$cross2==1)
#scan.norm.imp.cov <- scanone(comp, method="imp",model='normal',pheno.col=6,addcovar=comp$pheno$cross1)

NBH.perm <- NBH$perms.norm.imp
NEW.perm <- NEW$perms.norm.imp

th.NBH <- summary(NBH.perm)[1,]
th.NEW <- summary(NEW.perm)[1,]

norm.qtl.NBH <- summary(scan.norm.imp.NBH, perms=NBH.perm, alpha=0.05)
norm.qtl.NEW <- summary(scan.norm.imp.NEW, perms=NEW.perm, alpha=0.05)
qtl.uns.NBH <- makeqtl(cross, chr=norm.qtl.NBH$chr, pos=norm.qtl.NBH$pos,what="draws")
qtl.uns.NEW <- makeqtl(cross, chr=norm.qtl.NEW$chr, pos=norm.qtl.NEW$pos,what="draws")
chrs.model <- unique(qtl.uns.NEW$chr, qtl.uns.NEW$chr)


scan.all <- c(scan.norm.imp.NBH,scan.norm.imp.NEW)

png('~/NBH_NEW.png')
plot(scan.norm.imp.NBH, scan.norm.imp.NEW, chr=chrs.model, lodcolumn=1,
  incl.markers=FALSE,lty=1, col=c("violetred","black"), lwd=2,
  add=FALSE, gap=25, mtick = "line",bandcol="gray70",
  show.marker.names=FALSE, alternate.chrid=T,
  type="l", cex=1, pch=1, bg="grey",
  bgrect=NULL)
dev.off()

png('~/NBH_NEW_ALL.png',width = 1000)
plot(scan.norm.imp.NBH, scan.norm.imp.NEW, chr=1:24, lodcolumn=1,
  incl.markers=FALSE,lty=1, col=c("violetred","black"), lwd=2,
  add=FALSE, gap=25, mtick = "line",bandcol="gray70",
  show.marker.names=FALSE, alternate.chrid=T,
  type="l", cex=1, pch=1, bg="grey",
  bgrect=NULL)
dev.off()


png('~/NBH_RF_2.png')
plotRF(cross.nbh, chr = 2)
dev.off()


## Positionmaps instead
for (X in 1:24){
  ord <- order(as.numeric(gsub(paste(X,":",sep=''),'',markernames(NBH$cross.18,chr=X))))
  NBH$cross.18 <- switch.order(NBH$cross.18, chr=X, ord, error.prob=0.01,
    map.function="kosambi",maxit=100, tol=1e-6, sex.sp=F)
  POS.map.18 <- est.map(NBH$cross.18,error.prob=0.01,map.function="kosambi", chr=X,maxit=10)
  NBH$cross.18 <- replace.map(NBH$cross.18, POS.map.18)
}

for (X in 1:24){
  ord <- order(as.numeric(gsub(paste(X,":",sep=''),'',markernames(NEW$cross.18,chr=X))))
  NEW$cross.18 <- switch.order(NEW$cross.18, chr=X, ord, error.prob=0.01,
    map.function="kosambi",maxit=100, tol=1e-6, sex.sp=F)
  POS.map.18 <- est.map(NEW$cross.18,error.prob=0.01,map.function="kosambi", chr=X,maxit=10)
  NEW$cross.18 <- replace.map(NEW$cross.18, POS.map.18)
}


nbh <- subset(cross.nbh,chr=2)
#Z <- repRipple(nbh, chr=2, error.prob=0.01, map.function="kosambi",window = 3)
Z <- repRipple.jm(nbh, chr=chrnames(nbh)[1], map.function="kosambi",window = 3,re.est.map=FALSE)

POS.map.18 <- est.map(Z,error.prob=ers,map.function="kosambi", chr=X,maxit=1)



save.image(file.path('~/NEW.NBH.QTLmap.Rsave'))
#perms.norm.imp <- scanone(comp, method="imp",model='normal',n.perm=5,perm.strata=cross.18$pheno$stata,pheno.col=6, n.cluster=slurmcore)






scan.cat <- c(NBH$scan.norm.imp,NEW$scan.norm.imp)


sapply(1:24,function(X){
    cmdist <-
}
)




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
