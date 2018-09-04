#!/bin/R
### Map QTLs seperate
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/pop_control_file.R')

# Table of Chroms with sig QTLs
test.QTLs <- read.table(file.path(basedir,'rQTL/metadata/QTLs.txt'),
              sep='\t',header=T)

## Get chrom number vector
test.QTLs$chrm.n <- gsub('chr','',test.QTLs$chrom)

## inactivate karl's summary.cross(). It takes way too long

## read in the QTL cross
cross.18 <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
                format='csvr', geno=c(1:3),estimate.map=FALSE)

## Change phenotype to 0/1
cross.18 <- fix.pheno(cross.18) ## Pheno (Dev Score 0,1) -> 0 and (Dev Score 3,4,5) -> 1

## Remove problematic individuals
subset.ind <- cross.18$pheno$ID[!cross.18$pheno$ID %in% inds]
cross.18 <- subset(cross.18, ind=subset.ind)

## Map each QTL chro independently
allbut <- c(1:24)[-X]
subset.qtl <- chrnames(cross.18)[!chrnames(cross.18) %in% allbut]
cross.18 <- subset(cross.18, chr=subset.qtl)

## Conservative
cross.18 <- drop.missing(cross.18,50)

## Conservative
cross.18 <- distort(cross.18,0.005)

## Keep markers close to QTLs
qtl.index <-  which(test.QTLs$chrm.n == X)

## These markers should not be filtered (close to a QTL)
tokeep <- unlist(sapply(qtl.index,function(Z){
    markerList <- list()
    markerList[[Z]] <- keepQTL(Z,i=cross.18)
    return(markerList)
    }
  )
)
print('forming all linkage groups')
cross.18.all <- formLinkageGroups(cross.18, max.rf=0.35, min.lod=10, reorgMarkers=TRUE)

keep <- sapply(1:nchr(cross.18.all),function(i){
      l <- sum(nmar(cross.18.all))*.01
      return(sum(X==gsub('\\:.*','',markernames(cross.18.all,chr=i))) > l)
      }
    )
keep <- names(cross.18.all$geno)[keep]

cross.18.Z <- subset(cross.18.all, chr=keep)
cross.18.Z <- formLinkageGroups(cross.18.Z, max.rf=0.5, min.lod=4, reorgMarkers=TRUE)

## fix phase
chrom.b4 <- nchr(cross.18.Z)
if (chrom.b4 > 1){
  print('fixing phase')
  chrom.after <- 0
  while (chrom.b4 > chrom.after){
    chrom.b4 <- nchr(cross.18.Z)
    cross.18.Z <- switchAlleles(cross.18.Z,markernames(cross.18.Z,chr=1))
    cross.18.Z <- formLinkageGroups(cross.18.Z, max.rf=0.5, min.lod=grpLod, reorgMarkers=TRUE)
    chrom.after <- nchr(cross.18.Z)
    }
  print('done fixing phase')
} else {
  print('did not need to fix phase')
}

## form linkage groups on phase-fixed data
LGtable <- formLinkageGroups(cross.18.Z, max.rf=0.35, min.lod=finLod)
cross.18 <- subset(cross.18.Z, chr=which.max(table(LGtable$LG)))

## rename to the correct LG
names(cross.18$geno) <- X

print('saving...')
save.image(paste('chr',X,'.QTLmap.Rsave',sep=''))

print('second distortion filter')
cross.18 <- distort.18(cross.18,0.05)

cross.18 <- drop.missing.18(cross.18,missing=missing)

print('order filtered markers. Takes awhile...')
cross.18 <- orderMarkers(cross.18,chr=X,window=7,use.ripple=T,error.prob=0.002, map.function='kosambi',sex.sp=F)

print('removing double cross-overs')
cross.18 <- removeDoubleXO(cross.18, verbose=F)

print('Dropping 5% of markers that inflate the map. Takes a long time...')
## Drop one marker, p is proportion  of worst markers to drop
cross.18 <- dropone.par(cross.18,p=0.05,chr=X,maxit=2,
  sex.sp = F,verbose=F,parallel=T)

print('saving...')
save.image(paste('chr',X,'.QTLmap.Rsave',sep=''))

print('Re-setimating map from filtered data')
POS.map.18 <- est.map(cross.18,error.prob=0.002,map.function="kosambi",n.cluster=6, chr=X)
cross.18 <- replace.map(cross.18, POS.map.18)

print('Re-write the markers to rQTL formate')
write.cross(file=paste('chr',X,'.QTLmap',sep=''),format="csv",filestem=plotdir,chr=X)

print('saving...')
save.image(paste('chr',X,'.QTLmap.Rsave',sep=''))

## Scan for QTL

print('scanning for a single QTL')
GP <- calc.genoprob(cross.18, step=2.5)

GP <- sim.geno(GP,n.draws=1000, step=2, err=0.02)

scanQTL <- scanone(GP, pheno.col=1, model="binary", method="hk")

print('saving...')
save.image(paste('chr',X,'.QTLmap.Rsave',sep=''))
print(paste('done with chrom',X,'in pop',P))
