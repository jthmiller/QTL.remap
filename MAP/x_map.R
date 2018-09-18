#!/bin/R
### Map Chromosome 5 (potentiall sex)
setwd('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/')
source('control_file.R')

## For plotting
marker_dens <- list()

# Table of Chroms with sig QTLs
test.QTLs <- read.table(file.path(basedir,'rQTL/metadata/QTLs.txt'),
              sep='\t',header=T)

## Get chrom number vector
test.QTLs$chrm.n <- gsub('chr','',test.QTLs$chrom)

print(pop)
print(X)
## read in the QTL cross
cross.18 <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
                format='csvr', geno=c(1:3),estimate.map=FALSE)

## Pheno (Dev Score 0,1) -> 0 and (Dev Score 3,4,5) -> 1
cross.18 <- fix.pheno(cross.18)

## Remove problematic individuals
subset.ind <- cross.18$pheno$ID[!cross.18$pheno$ID %in% inds]
cross.18 <- subset(cross.18, ind=subset.ind)

## Map each QTL chro independently
allbut <- c(1:24)[-X]
subset.qtl <- chrnames(cross.18)[!chrnames(cross.18) %in% allbut]
cross.18 <- subset(cross.18, chr=subset.qtl)

## Starting marker number
marker.warning()

## Specific to 'outname'
if(mapped.only==T){
cross.18 <- drop.markers(cross.18,markernames(cross.18)[grep('NW',markernames(cross.18))])
}

### Conservative for sex specific coverage
##print('Dropping marker with less than 60 genotypes')
##cross.18 <- drop.missing(cross.18,35)

## Conservative
print('Dropping markers with segregation distortion < 0.0005')
cross.18 <- distort(cross.18,0.0005)

cross.18.all <- formLinkageGroups(cross.18, max.rf=grpRf,min.lod=grpLod, reorgMarkers=TRUE)

keep <- sapply(1:nchr(cross.18.all),function(i){
      return(sum(X==gsub('\\:.*','',markernames(cross.18.all,chr=i))) > 1)
      }
    )
keep <- names(cross.18.all$geno)[keep]
cross.18 <- subset(cross.18.all, chr=keep)

print('forming initial linkage groups')
cross.18 <- formLinkageGroups(cross.18, max.rf=finRf, min.lod=finLod, reorgMarkers=TRUE)

## iteritively fix phase by inverting phase of the growing LG and re-eval linkage
chrom.b4 <- nchr(cross.18)
if (chrom.b4 > 1){
  print('fixing phase')
  chrom.after <- 0
  while (chrom.b4 > chrom.after){
    chrom.b4 <- nchr(cross.18)
    cross.18 <- switchAlleles(cross.18,markernames(cross.18,chr=1))
    cross.18 <- formLinkageGroups(cross.18, max.rf=finRf, min.lod=grpLod, reorgMarkers=TRUE)
    chrom.after <- nchr(cross.18)
    }
  print('done fixing phase... switching phase added no new markers to the LG')
} else {
  print('did not need to fix phase')
}

cross.18a <- formLinkageGroups(cross.18, max.rf=finRf, min.lod=15, reorgMarkers=TRUE)
cross.18a <- subset(cross.18, chr=1)
cross.18$pheno$sex <- as.numeric(nmissing(cross.18a)[cross.18a$pheno$ID]>130)
y <- pull.pheno(cross.18, 2)

cross.18m <- subset(cross.18,ind=(y==1))
cross.18f <- subset(cross.18,ind=(y==0))

cross.18f <- drop.missing(cross.18f,5)
cross.18f <- formLinkageGroups(cross.18f, max.rf=finRf, min.lod=finLod, reorgMarkers=TRUE)
cross.18a <- orderMarkers(cross.18a,chr=X,window=5,use.ripple=T,
  error.prob=0.1, map.function='kosambi',sex.sp=F,maxit=100,tol=1e-2)
  print('Dropping marker with less than 60 genotypes')




### Try mapping each LG in 5 to see if cov is different between top/bottom half

#### indiv linkage groups (use to find diffs in cov)
### Split high missing
cross.18a <- subset(cross.18, chr=1)
cross.18a <- drop.missing(cross.18a,10)
cross.18a <- formLinkageGroups(cross.18a, max.rf=finRf, min.lod=finLod, reorgMarkers=TRUE)
cross.18a <- orderMarkers(cross.18a,chr=X,window=5,use.ripple=T,
  error.prob=0.1, map.function='kosambi',sex.sp=F,maxit=100,tol=1e-2)
  print('Dropping marker with less than 60 genotypes')


cross.18b <- subset(cross.18, chr=c(2,3))
cross.18b <- switchAlleles(cross.18b,markernames(cross.18b,chr=2))
cross.18b <- formLinkageGroups(cross.18b, max.rf=finRf, min.lod=finLod, reorgMarkers=TRUE)
cross.18b <- orderMarkers(cross.18b,chr=1,window=5,use.ripple=T,
  error.prob=0.1, map.function='kosambi',sex.sp=F,maxit=1000,tol=1e-2)
################

#### Low missing in whole group
### map low missing (if something is still not in linkage here, should distinguish)
cross.18c <- drop.missing(cross.18,82)
cross.18c <- formLinkageGroups(cross.18c, max.rf=finRf, min.lod=finLod, reorgMarkers=TRUE)
cross.18c <- subset(cross.18c, chr=1)
cross.18c <- orderMarkers(cross.18c,chr=1,window=5,use.ripple=T,
  error.prob=0.1, map.function='kosambi',sex.sp=F,maxit=1000,tol=1e-2)

cross.18d <- subset(cross.18, chr=1)
cross.18d <- formLinkageGroups(cross.18d, max.rf=finRf, min.lod=finLod, reorgMarkers=TRUE)

cross.18d <- subset(cross.18d,ind=names(nmissing(cross.18)[which(nmissing(cross.18)<100)]))
cross.18d <- drop.missing(cross.18d,10)
cross.18d <- subset(cross.18d, chr=1)
cross.18d <- orderMarkers(cross.18d,chr=1,window=5,use.ripple=T,
  error.prob=0.1, map.function='kosambi',sex.sp=F,maxit=10000,tol=1e-2)

### any markers in lower half that have lots of data?
cross.18z <- subset(cross.18,ind=names(nmissing(cross.18)[which(nmissing(cross.18)>100)]))
cross.18z <- formLinkageGroups(cross.18z, max.rf=finRf, min.lod=finLod, reorgMarkers=TRUE)
cross.18z <- drop.missing(cross.18z,10)
cross.18z <- subset(cross.18z, chr=1)
cross.18z <- formLinkageGroups(cross.18z, max.rf=finRf, min.lod=finLod, reorgMarkers=TRUE)
cross.18z <- orderMarkers(cross.18z,chr=1,window=5,use.ripple=T,
    error.prob=0.1, map.function='kosambi',sex.sp=F,maxit=10000,tol=1e-2)


save.image('~/debug.Rsave')






### Keep all with linkage to marker mapped to the LG
keep <- sapply(1:nchr(cross.18),function(i){
      return(sum(X==gsub('\\:.*','',markernames(cross.18,chr=i))) > 1)
      }
    )
keep <- names(cross.18$geno)[keep]
cross.18 <- subset(cross.18, chr=keep)
cross.18 <- formLinkageGroups(cross.18, max.rf=finRf, min.lod=grpLod, reorgMarkers=TRUE)
LGtable <- formLinkageGroups(cross.18, max.rf=finRf, min.lod=grpLod)

## form linkage groups on phase-fixed data
cross.18 <- subset(cross.18, chr=which.max(table(LGtable$LG)))
