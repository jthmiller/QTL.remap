## Confirming the Grandparents in elr
packs <- c('qtl','foreach','doParallel')
lapply(packs, require, character.only = TRUE)
setwd('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/QTL/')
slurmcore <- 12
chrms <- c(1:24)
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')

# indpops
pop <- 'NBH'
cutoff <- 1.0e-5
indsa.nbh <- c('ind87','ind88','ind15')



if(pop=='NBH'){
  miss <- 3
  cross <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
          format='csvr', geno=c(1:3),estimate.map=FALSE)
  ## Take out pars
  cross.pars <- subset(cross,ind=is.na(cross$pheno$Pheno))
  cross.pars$pheno$ID <-  c('par1','par2')
  cross <- subset(cross,ind=!is.na(cross$pheno$Pheno))
  cross$pheno$ID <- paste("ind", 1:nind(cross), sep="")
  cross <- subset(cross,ind=!cross$pheno$ID %in% indsa.nbh)
  gt.cross.par <- geno.table(cross)
  cross <- drop.markers(cross,rownames(gt.cross.par[gt.cross.par$missing>miss,]))
  cross <- drop.markers(cross,rownames(gt.cross.par[gt.cross.par$P.value<cutoff,]))
  drop <- markernames(cross.pars)[!markernames(cross.pars) %in% markernames(cross)]
  cross.pars  <- drop.markers(cross.pars,drop)
  cross <- c(cross.pars,cross)
  cross <- c(cross.pars,cross)
  rela <- comparegeno(cross)
  colnames(rela) <- cross$pheno$ID
  rownames(rela) <- cross$pheno$ID
  rela[rela==NaN] <- NA
  diag(rela) <- NA

}
if(pop=='ELR'){

  cross <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
                format='csvr', geno=c(1:3),estimate.map=FALSE)

  cross.pars <- subset(cross,ind=is.na(cross$pheno$Pheno))
  cross.pars$pheno$ID <- c('par1','par2')

  cross <- subset(cross,ind=!is.na(cross$pheno$Pheno))
  cross$pheno$ID <- paste("ind", 1:nind(cross), sep="")
  #cross <- subset(cross,ind=!cross$pheno$ID %in% indsa.ELR)
  gt.cross.par <- geno.table(cross)
  miss <- 8
  cutoff <- 1.0e-06
  cross <- drop.markers(cross,rownames(gt.cross.par[gt.cross.par$missing>miss,]))
  cross <- drop.markers(cross,rownames(gt.cross.par[gt.cross.par$P.value<cutoff,]))

  drop <- markernames(cross.pars)[!markernames(cross.pars) %in% markernames(cross)]
  cross.pars  <- drop.markers(cross.pars,drop)

  cross.b <- c(cross.pars,cross)

  rela <- comparegeno(cross)
  colnames(rela) <- cross$pheno$ID
  rownames(rela) <- cross$pheno$ID
  rela[rela==NaN] <- NA
  diag(rela) <- NA

  #Drop dup inds
  wh <- which(rela > 0.75, arr=TRUE)
  wh <- wh[wh[,1] < wh[,2],]
  wh2 <- wh[,2]
  names(wh2) <- rownames(rela)[wh2]
  wh <- wh[,1]
  numiss <- nmissing(cross)
  in.drop <- sapply(1:length(wh),function(X){
    gdn <- which.max(c(numiss[wh[X]],numiss[wh2[X]]))
    return(names(c(wh[X],wh2[X]))[gdn])
    }
  )
  ###cross <- subset(cross, ind=!cross$pheno$ID %in% in.drop)

  save.image('~/ELR.kinship.Rsave')


#load('NBH.kinship.Rsave')
load('ELR.kinship.Rsave')
library(qtl)

rela <- comparegeno(cross.b)
colnames(rela) <- cross.b$pheno$ID
rownames(rela) <- cross.b$pheno$ID
rela[rela==NaN] <- NA
diag(rela) <- NA

hist(rela)
##resize
heatmap(rela,symm=T)


gt.cross.par <- geno.table(cross)
miss <- 5
cutoff <- 1.0e-02
cross <- drop.markers(cross,rownames(gt.cross.par[gt.cross.par$missing>miss,]))
cross <- drop.markers(cross,rownames(gt.cross.par[gt.cross.par$P.value<cutoff,]))

rela <- comparegeno(cross)
colnames(rela) <- cross$pheno$ID
rownames(rela) <- cross$pheno$ID
rela[rela==NaN] <- NA
diag(rela) <- NA

dups <- findDupMarkers(cross, exact.only=FALSE, adjacent.only=FALSE)
### remove markers that are exactly the same.
cross <- drop.markers(cross, unlist(dups))
write.cross(cross,filestem=paste(popdir,'/ERL_prefilter.QTLmap',sep=''),format="csv")
