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

##########scratch #############


cross.pars.NBH <- cross.pars
cross.pars.NEW <- cross.pars
cross.pars.ELR <- cross.pars

all.pars <- c(cross.pars.NBH,cross.pars.NEW,cross.pars.ELR)
all.pars.gt <- geno.table(all.pars)

sum()


par1 <- subset(cross.pars,ind=cross.pars$pheno$ID[1])
par2 <- subset(cross.pars,ind=cross.pars$pheno$ID[2])
par1.gt <- geno.table(par1)
par2.gt <- geno.table(par2)
par12.gt <- geno.table(cross.pars)
ofs <- geno.table(cross.18)




gp1.het.gp2.hom <- rownames(par1.gt )[which(par1.gt$AB==1 & par12.gt$missing<1)]

ofs[gp1.het.gp2.hom,]
par1.offs <- drop.markers(cross.18,(!rownames(cross.18)%in%gp1.het.gp2.hom))
cross.18 <- formLinkageGroups(par1.offs, min.lod=finLod, reorgMarkers=TRUE)

drop <- rownames(par12.gt)[which(par12.gt$AA==2 | par12.gt$BB==2 | par12.gt$missing==2)]
cross.18 <- drop.markers(cross.18,drop)
##########################

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
  gt.pars <- geno.table(cross.pars)
  gt.pars.set2 <- gt.pars[which(gt.pars$missing==1),]
  gt.pars.set2 <- rownames(gt.pars.set2)[which(gt.pars.set2$AA==1 | gt.pars.set2$BB==1)]
  gt.pars <- rownames(gt.pars)[which(gt.pars$AA==1 & gt.pars$BB==1)]
  print('Dropped markers that do not follow expectations from GP genotypes')
  gt.pars <- c(gt.pars,gt.pars.set2)
  marks.drop <- markernames(cross.18)[!markernames(cross.18) %in% gt.pars]
  cross.18 <- drop.markers(cross.18,marks.drop)
  par.pos <- as.numeric(gsub(paste(X,':',sep=''),'',markernames(cross.18)))
}

marker.warning()

print('Removing duplicates')
dups <- findDupMarkers(cross.18, exact.only=F, adjacent.only=T)
cross.18 <- drop.markers(cross.18, unlist(dups))

marker.warning()

### Table before missing filter
gt <- geno.table(cross.18)
pos <- as.numeric(gsub(paste(X,':',sep=''),'',rownames(gt)))
pval <- log10(gt$P.value)

#### Filter Conservative
print(paste('Dropping markers with',miss,'or more genotypes missing'))
cross.18 <- drop.missing(cross.18,miss)
gt.missing <- geno.table(cross.18)

marker.warning()

### invariants
gt.cross.par <- NA
cross.18 <- drop.markers(cross.18,rownames(gt.missing[gt.missing$P.value<cutoff,]))
gt.pval <- geno.table(cross.18)

marker.warning()

#write(file=paste(popdir,'/chr',X,'_',outname,'.keepmarkers.csv',sep=''),
#  tokeep,sep = ",")

print('forming initial linkage groups to fix phase...')
cross.18 <- formLinkageGroups(cross.18, max.rf=grpRf, min.lod=grpLod, reorgMarkers=TRUE)

## iteritively fix phase by inverting phase of the growing LG and re-eval linkage
chrom.b4 <- nchr(cross.18)
if (chrom.b4 > 1){
  print('fixing phase')
  chrom.after <- 0
  while (chrom.b4 > chrom.after){
    chrom.b4 <- nchr(cross.18)
    cross.18 <- switchAlleles(cross.18,markernames(cross.18,chr=chrnames(cross.18)[1]))
    cross.18 <- formLinkageGroups(cross.18, max.rf=finRf, min.lod=finLod, reorgMarkers=TRUE)
    chrom.after <- nchr(cross.18)
    }
  print('done fixing phase... switching phase added no new markers to the LG')
} else {
  print('did not need to fix phase')
}

## Pull genotypes for QTL markers
#gi <- pull.geno(cross.18)[,tokeep]

### Continue to keep that are linked to markers that have been previously mapped
#keep <- sapply(1:nchr(cross.18),function(i){
#    a <- sum(X==gsub('\\:.*','',markernames(cross.18,chr=i)))
#    b <- sum(tokeep %in% markernames(cross.18,chr=i))
#    return(a+b > 2)
#  }
#)
#keep <- names(cross.18$geno)[keep]
#cross.18 <- subset(cross.18, chr=keep)

cross.18 <- formLinkageGroups(cross.18, max.rf=grpRf, min.lod=grpLod, reorgMarkers=TRUE)
LGtable <- formLinkageGroups(cross.18, max.rf=grpRf, min.lod=grpLod)
keep <- which.max(table(LGtable$LG))
## form linkage groups on phase-fixed data
cross.18 <- subset(cross.18, chr=keep)
## rename to the correct LG
names(cross.18$geno) <- X

marker.warning()

if (reorder.marks==T){
  print('initial order filtered markers with 0.01 errorprob. Takes awhile...')
  cross.18 <- orderMarkers(cross.18,chr=X,window=5,use.ripple=T,
    error.prob=0.05, map.function='kosambi',sex.sp=F,maxit=1000,tol=1e-2)
  print('Estimating the initial map with high errorprob')
  POS.map.18 <- est.map(cross.18,error.prob=0.01,map.function="kosambi", chr=X,maxit=100)
  cross.18 <- replace.map(cross.18, POS.map.18)
} else {
  print('estimating map with markers at physical positions')
  ord <- order(as.numeric(gsub(paste(X,":",sep=''),'',markernames(cross.18,chr=X))))
  cross.18 <- switch.order(cross.18, chr=X,ord , error.prob=0.01,
    map.function="kosambi",maxit=1000, tol=1e-3, sex.sp=F)
}

print(summary(pull.map(cross.18))[as.character(X),])

print('removing double cross-over genotypes')
## rename to the correct LG
names(cross.18$geno) <- X
cross.18 <- removeDoubleXO(cross.18, verbose=T)
print('Done removing dxo..')

marker.warning()

print('Removing duplicates induced by DXO')
dups <- findDupMarkers(cross.18, exact.only=F, adjacent.only=T)
cross.18 <- drop.markers(cross.18, unlist(dups))
print('Done..')

marker.warning()

print('Dropping markers with too many genotypes missing after removing DXO')
cross.18 <- drop.missing(cross.18,miss)
gt.missing <- geno.table(cross.18)
print('Done dropping markers')

marker.warning()

print('Final map estimation')

POS.map.18 <- est.map(cross.18,error.prob=0.01,map.function="kosambi", chr=X,maxit=1000)
cross.18 <- replace.map(cross.18, POS.map.18)

print(summary(pull.map(cross.18))[as.character(X),])

gt.fin <- geno.table(cross.18)

print('Writing the markers to rQTL format')
write.cross(cross.18,filestem=paste(popdir,'/chr',X,'_',outname,'.QTLmap',sep=''),format="csv",chr=X)

print('saving...')
rm(cross.18)
save.image(paste(popdir,'/chr',X,'_',outname,'.QTLmap.Rsave',sep=''))

print('plotting')
png(file.path(popdir,paste(X,'_pval.png',sep='')))
hist.geno(gt.missing$P.value)
abline(v=log10(cutoff))
dev.off()

png(file.path(popdir,paste(X,'_pos.png',sep='')))
par(mfrow=c(3,1),mar=c(1,2,3,1))
plot.geno(gt,gen.main=paste('All Mapped Markers on Chr',X))
#plot.geno.2(gt.missing,wtf=T, gen.main='Markers with < 5 missing GT, Grandparent Confirmed markers in Green')
abline(h=log10(cutoff),col='red')
plot.geno(gt.pval,gen.main=paste('Filter distortion > ',cutoff))
plot.geno(gt.fin,gen.main='Final Markers')
dev.off()
