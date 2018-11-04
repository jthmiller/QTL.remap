#!/bin/R
### Map QTLs 1 of 3
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')
## For plotting
marker_dens <- list()

# Table of Chroms with sig QTLs
test.QTLs <- read.table(file.path(basedir,'rQTL/metadata/QTLs.txt'),
              sep='\t',header=T)

## Get chrom number vector
test.QTLs$chrm.n <- gsub('chr','',test.QTLs$chrom)

print(paste(pop,X,sep=' '))
############

## read in the QTL cross
cross.18 <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
                format='csvr', geno=c(1:3),estimate.map=FALSE)

### Pull names from plinkfile
path <- file.path(indpops,paste(pop,'.ped',sep=''))
popname <- system(paste('cut -f1 -d\' \'',path),intern = TRUE)
indname <- system(paste('cut -f2 -d\' \'',path),intern = TRUE)
cross.18$pheno$ID <- paste(popname,indname,sep='_')

## Subset and drop parents
cross.pars <- subset(cross.18,ind=is.na(cross.18$pheno$Phen))

## Remove problematic individuals (found by kinship analysis)
con <- file(file.path(popdir,'kinship.keep.ind.txt'), open='r')
keepers <- readLines(con)
close(con)

print('Dropping kinship outliers')
cross.18 <- subset(cross.18,ind=cross.18$pheno$ID %in% keepers)
cross.18 <- subset(cross.18,ind=!is.na(cross.18$pheno$Phen))

print('Dropping all chromosomes except the one to map')
## Map each QTL chro independently
if (mapped.only==T){
  allbut <- c(1:24)[-X]
  subset.qtl <- chrnames(cross.18)[!chrnames(cross.18) %in% allbut]
  cross.18 <- subset(cross.18, chr=subset.qtl)
  cross.pars <- subset(cross.pars, chr=subset.qtl)
  marker.warning()
}

### grandparent confirmed markers
if(confirmed==F){
  print('Not using GP genos for filtering')
  dups <- findDupMarkers(cross.18, exact.only=FALSE, adjacent.only=FALSE)
  cross.18 <- drop.markers(cross.18, unlist(dups))
} else if (confirmed==T & pop=='ELR'){
  print('Filtering with one G.parent')
  cross.bli <- subset(cross.pars, ind='ELR_ER1124F')
  gt.bli <- geno.table(cross.bli)
  gt.bli <- rownames(gt.bli)[which(gt.bli$AA==1 | gt.bli$BB==1)]
  marks.drop <- markernames(cross.18)[!markernames(cross.18) %in% gt.bli]
  cross.18 <- drop.markers(cross.18,marks.drop)
  par.pos <- as.numeric(gsub(paste(X,':',sep=''),'',markernames(cross.18)))
} else {
  print('Filtering with both G.parents')
#  gt.pars <- geno.table(cross.pars)
#  gt.pars.set2 <- gt.pars[which(gt.pars$missing==1),]
#  gt.pars.set2 <- rownames(gt.pars.set2)[which(gt.pars.set2$AA==1 | gt.pars.set2$BB==1)]
#  gt.pars <- rownames(gt.pars)[which(gt.pars$AA==1 & gt.pars$BB==1)]
#  gt.pars <- c(gt.pars,gt.pars.set2)

  if(ns=='South'){
    con <- file(file.path(dirso,'fixed.south.txt'))
    gt.pars <- readLines(con)
    close(con)
  } else {
    con <- file(file.path(dirso,'fixed.north.txt'))
    gt.pars <- readLines(con)
  close(con)
  }
  marks.drop <- markernames(cross.18)[!markernames(cross.18) %in% gt.pars]
  cross.18 <- drop.markers(cross.18,marks.drop)
  par.pos <- as.numeric(gsub(paste(X,':',sep=''),'',markernames(cross.18)))

  print('Dropped markers that do not follow expectations from GP genotypes')

}

cross.18 <- formLinkageGroups(cross.18, min.lod=12, reorgMarkers=TRUE)
cross.18 <- subset(cross.18,chr=c('1','2'))
cross.18 <- switchAlleles(cross.18,markernames(cross.18,chr=chrnames(cross.18)[1]))
cross.18 <- formLinkageGroups(cross.18, min.lod=12, reorgMarkers=TRUE)
gt <- geno.table(cross.18)
#### Filter Conservative
print(paste('Dropping markers with',miss,'or more genotypes missing'))
cross.18 <- drop.missing(cross.18,miss1)
gt.missing <- geno.table(cross.18)
cross.18 <- formLinkageGroups(cross.18, max.rf=0.15,min.lod=12, reorgMarkers=TRUE)
cross.18 <- subset(cross.18,chr=c('1'))
gt.missing <- geno.table(cross.18)
cross.18 <- drop.markers(cross.18,rownames(gt.missing[gt.missing$P.value<cutoff,]))
gt.pval <- geno.table(cross.18)

print('estimating map with markers at physical positions')
ord <- order(as.numeric(gsub(paste(X,":",sep=''),'',markernames(cross.18,chr=X))))
cross.18 <- switch.order(cross.18, chr=X,ord , error.prob=0.01,
  map.function="kosambi",maxit=1000, tol=1e-3, sex.sp=F)

POS.map.18 <- est.map(cross.18,error.prob=0.10,map.function="kosambi", chr=X,maxit=1000)
cross.18 <- replace.map(cross.18, POS.map.18)

print(summary(pull.map(cross.18))[as.character(X),])

print('Writing the markers to rQTL format')
write.cross(cross.18,filestem=paste(popdir,'/chr',X,'_',outname,'.manymarks.QTLmap',sep=''),format="csv",chr=X)

cross.18 <- read.cross(format='csv',dir=popdir,
   file=paste('chr',X,'_',outname,'.manymarks.QTLmap.csv',sep=''),
   geno=c('AA','AB','BB'),alleles=c("A","B"))
