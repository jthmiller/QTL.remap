#!/bin/R
### Map QTLs 1 of 3
setwd('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/')
source('control_file.R')

## For plotting
marker_dens <- list()

# Table of Chroms with sig QTLs
test.QTLs <- read.table(file.path(basedir,'rQTL/metadata/QTLs.txt'),
              sep='\t',header=T)

## Get chrom number vector
test.QTLs$chrm.n <- gsub('chr','',test.QTLs$chrom)

print(paste(pop,X,sep=' '))


### parents ###
cross.pars <- read.cross.jm(file=file.path(indpops,paste(pop,'.parents.csvr',sep='')),
                format='csvr', geno=c(1:3),estimate.map=FALSE)


allbut <- c(1:24)[-X]
subset.qtl <- chrnames(cross.pars)[!chrnames(cross.pars) %in% allbut]
cross.pars <- subset(cross.pars, chr=subset.qtl)

if(mapped.only==T){
cross.pars <- drop.markers(cross.pars,markernames(cross.pars)[grep('NW',markernames(cross.pars))])
}

gt.par <- geno.table(cross.pars)

par.confirm.marks <- rownames(gt.par[which(gt.par$AA==1 & gt.par$BB==1),])

############

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

### Table before missing filter
gt <- geno.table(cross.18)
pos <- as.numeric(gsub(paste(X,':',sep=''),'',rownames(gt)))
pval <- log10(gt$P.value)

#### Filter Conservative
print('Dropping markers with more than 5 genotypes missing')
cross.18 <- drop.missing(cross.18,miss)
gt.missing <- geno.table(cross.18)
marker.warning()

### grandparent confirmed markers
marks.confirm <- !rownames(gt.missing) %in% par.confirm.marks
cross.18 <- drop.markers(cross.18,rownames(gt.missing)[marks.confirm])
par.pos <- as.numeric(gsub(paste(X,':',sep=''),'',par.confirm.marks))
gt.cross.par <- geno.table(cross.18)

### invariants
cross.18 <- drop.markers(cross.18,rownames(gt.cross.par[gt.cross.par$P.value<cutoff,]))
gt.pval <- geno.table(cross.18)


#a <- intersect(which(gt$neglog10P > 2.5),which(gt$AA < 0.15 & gt$BB > 0.15))
#b <- intersect(which(gt$neglog10P > 2.5),which(gt$BB < 0.15 & gt$AA > 0.15))
#c <- intersect(which(gt$neglog10P > 2.5),which(gt$BB < 0.15 & gt$AA < 0.15))
#cross.18 <- drop.markers(cross.18,rownames(gt[unique(c(a,b,c)),]))
#gt.dis <- geno.table(cross.18)
#print(paste('Dropping markers with segregation distortion < ',cutoff))
#cross.18 <- distort(cross.18,cutoff)

marker.warning()

print('Finding markers that are near known QTLs')
qtl.index <-  which(test.QTLs$chrm.n == X)
tokeep <- unlist(sapply(qtl.index,function(Z){
    markerList <- list()
    markerList[[Z]] <- keepQTL(Z,i=cross.18)
    return(markerList)
    }
  )
)

write(file=paste(popdir,'/chr',X,'_',outname,'.keepmarkers.csv',sep=''),
  tokeep,sep = ",")

if (mapped.only==T){print('keeping only previously mapped markers')
} else {print('Keeping markers that show even low linkage to any previously known mapped marker')}
print(paste('Using an initial lod of 6 to keep markers '))

cross.18.all <- formLinkageGroups(cross.18, max.rf=grpRf,min.lod=grpLod, reorgMarkers=TRUE)

keep <- sapply(1:nchr(cross.18.all),function(i){
    a <- sum(X==gsub('\\:.*','',markernames(cross.18,chr=i)))
    b <- sum(tokeep %in% markernames(cross.18,chr=i))
    return(a+b > 1)
  }
)
keep <- names(cross.18.all$geno)[keep]
cross.18 <- subset(cross.18.all, chr=keep)
rm(cross.18.all) ##keep memory light

print(paste(length(markernames(cross.18)),'markers in',nchr(cross.18),'linkage groups going into phase fix'))

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

## Pull genotypes for QTL markers
gi <- pull.geno(cross.18)[,tokeep]

### Keep all with linkage to marker mapped to the LG
keep <- sapply(1:nchr(cross.18),function(i){
    a <- sum(X==gsub('\\:.*','',markernames(cross.18,chr=i)))
    b <- sum(tokeep %in% markernames(cross.18,chr=i))
    return(a+b > 1)
  }
)
keep <- names(cross.18$geno)[keep]
cross.18 <- subset(cross.18, chr=keep)

cross.18 <- formLinkageGroups(cross.18, max.rf=finRf, min.lod=grpLod, reorgMarkers=TRUE)
LGtable <- formLinkageGroups(cross.18, max.rf=finRf, min.lod=grpLod)

keep <- sapply(1:nchr(cross.18),function(i){
    b <- sum(tokeep %in% markernames(cross.18,chr=i))
    return(b > 1)
  }
)
keep <- which.max(table(LGtable$LG))
## form linkage groups on phase-fixed data
cross.18 <- subset(cross.18, chr=keep)

## rename to the correct LG
names(cross.18$geno) <- X

marker.warning()

print('second missing filter')
cross.18 <- drop.missing.18(cross.18,missing=missing)

marker.warning()

print('initial order filtered markers with 0.1 errorprob. Takes awhile...')  ###

cross.18 <- formLinkageGroups(cross.18, max.rf=finRf, min.lod=finLod,reorgMarkers=TRUE)
LGtable <- formLinkageGroups(cross.18, max.rf=finRf, min.lod=finLod)
cross.18 <- subset(cross.18, chr=which.max(table(LGtable$LG)))
names(cross.18$geno) <- X

### add QTL markers to cross that may have been removed
##return.dropped.markers()

cross.18 <- orderMarkers(cross.18,chr=X,window=5,use.ripple=T,
  error.prob=0.5, map.function='kosambi',sex.sp=F,maxit=1000,tol=1e-2)

## rename to the correct LG
names(cross.18$geno) <- X

cross.18 <- sim.geno(cross.18,step=5,error.prob=0.1,map.function='kosambi')
cross.18 <- calc.genoprob(cross.18,step=5,error.prob=0.1,map.function='kosambi')

print('removing double cross-over genotypes')
cross.18 <- removeDoubleXO(cross.18, verbose=T)
print('Done removing dxo..')

print('Estimating the initial map with high errorprob')
POS.map.18 <- est.map(cross.18,error.prob=0.1,map.function="kosambi", chr=X,maxit=100)
cross.18 <- replace.map(cross.18, POS.map.18)
names(cross.18$geno) <- X
print(summary(pull.map(cross.18))[as.character(X),])

gt.fin <- geno.table(cross.18)

print('Writing the markers to rQTL format')
write.cross(cross.18,filestem=paste(popdir,'/chr',X,'_',outname,'.QTLmap',sep=''),format="csv",chr=X)

print('saving...')
rm(cross.18)
save.image(paste(popdir,'/chr',X,'_',outname,'.QTLmap.Rsave',sep=''))
