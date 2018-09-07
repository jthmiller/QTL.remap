#!/bin/R
### Map QTLs seperate
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/pop_control_file.R')

## For plotting
marker_dens <- list()

# Table of Chroms with sig QTLs
test.QTLs <- read.table(file.path(basedir,'rQTL/metadata/QTLs.txt'),
              sep='\t',header=T)

## Get chrom number vector
test.QTLs$chrm.n <- gsub('chr','',test.QTLs$chrom)

print(pop)

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

## Starting marker number
marker.warning()

## Conservative
print('Dropping marker with less than 50 genotypes')
cross.18 <- drop.missing(cross.18,50)

marker.warning()

## Conservative
print('Dropping markers with segregation distortion < 0.000005')
cross.18 <- distort(cross.18,0.005)

marker.warning()

## Keep markers close to QTLs. Should not be filtered
qtl.index <-  which(test.QTLs$chrm.n == X)
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

cross.18 <- subset(cross.18.all, chr=keep)
cross.18 <- formLinkageGroups(cross.18, max.rf=0.5, min.lod=4, reorgMarkers=TRUE)
rm(cross.18.all) ##keep memory light

## fix phase
chrom.b4 <- nchr(cross.18)
if (chrom.b4 > 1){
  print('fixing phase')
  chrom.after <- 0
  while (chrom.b4 > chrom.after){
    chrom.b4 <- nchr(cross.18)
    cross.18 <- switchAlleles(cross.18,markernames(cross.18,chr=1))
    cross.18 <- formLinkageGroups(cross.18, max.rf=0.5, min.lod=grpLod, reorgMarkers=TRUE)
    chrom.after <- nchr(cross.18)
    }
  print('done fixing phase')
} else {
  print('did not need to fix phase')
}

## form linkage groups on phase-fixed data
LGtable <- formLinkageGroups(cross.18, max.rf=0.35, min.lod=finLod)
cross.18 <- subset(cross.18, chr=which.max(table(LGtable$LG)))

marker.warning()

## rename to the correct LG
names(cross.18$geno) <- X

print('second distortion filter')
cross.18 <- distort.18(cross.18,0.05)

marker.warning()

print('second distortion filter')
cross.18 <- drop.missing.18(cross.18,missing=missing)

marker.warning()

print('order filtered markers. Takes awhile...')

cross.18 <- orderMarkers(cross.18,chr=X,window=5,use.ripple=T,
  error.prob=0.002, map.function='kosambi',sex.sp=F,maxit=2000,tol=1e-3)

## rename to the correct LG
names(cross.18$geno) <- X

print('removing double cross-overs')
cross.18 <- removeDoubleXO(cross.18, verbose=F)

marker.warning()

print('Dropping 2.5% of markers that inflate the map. Takes a long time...')
## Drop one marker, p is proportion  of worst markers to drop
cross.18 <- dropone.par(cross.18,p=0.025,chr=X,maxit=3,
  sex.sp = F,verbose=F,parallel=T)

marker.warning()

print('determine error rate')
ers <- er.rate(cross.18)
print(paste(ers,' error rate'))

print('Re-setimating map from filtered data on')
cross.18 <- orderMarkers(cross.18,chr=X,window=5,use.ripple=T,
  error.prob=ers, map.function='kosambi',sex.sp=F,maxit=2000,tol=1e-3)
POS.map.18 <- est.map(cross.18,error.prob=ers,map.function="kosambi",n.cluster=12, chr=X)
cross.18 <- replace.map(cross.18, POS.map.18)

print('Done mapping..')
print(summary(pull.map(cross.18))[X,])

print('Re-writing the markers to rQTL format')
write.cross(cross.18,filestem=paste(plotdir,'chr',X,'.QTLmap',sep=''),format="csv",chr=X)

print('saving...')
save.image(paste('chr',X,'.QTLmap.Rsave',sep=''))
