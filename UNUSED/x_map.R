#!/bin/R
### Map Chromosome 5 (potentiall sex)
###setwd('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/')
###source('control_file.R')

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

cross.18.all <- drop.missing(cross.18.all,40)
cross.18.all <- formLinkageGroups(cross.18.all, min.lod=16, reorgMarkers=TRUE)
cross.18.all <- subset(cross.18.all, chr=c(1:10))
cross.18.all <- formLinkageGroups(cross.18.all, min.lod=16,max.rf=0.05, reorgMarkers=TRUE)
cross.18.all$pheno$cov <- nmissing(subset(cross.18.all,chr=1))



cross.18.all.AAxAB <- subset(cross.18.all, chr=c(1:5,7,9:11))
cross.18.all.AAxAB$pheno$hets1 <- rowSums(cross.18.all.AAxAB$geno$'1'[[1]]=='2',na.rm=T)
cross.18.all.AAxAB$pheno$hets2 <- rowSums(cross.18.all.AAxAB$geno$'2'[[1]]=='2',na.rm=T)
cross.18.all.AAxAB$pheno$hets3 <- rowSums(cross.18.all.AAxAB$geno$'3'[[1]]=='2',na.rm=T)
cross.18.all.AAxAB$pheno$hom1 <- rowSums(cross.18.all.AAxAB$geno$'1'[[1]]=='1',na.rm=T)
cross.18.all.AAxAB$pheno$hom2 <- rowSums(cross.18.all.AAxAB$geno$'1'[[1]]=='3',na.rm=T)
cross.18.all.ABxAB <- orderMarkers(cross.18.all.ABxAB,window=5,chr=3,use.ripple=T,
  error.prob=0.2, map.function='kosambi',sex.sp=F,maxit=5,tol=1e-2)



cross.18.all.ABxAB <- subset(cross.18.all, chr=c(6,8))
cross.18.all.ABxAB <- switchAlleles(cross.18.all.ABxAB,markernames(cross.18.all.ABxAB,chr=8))
cross.18.all.ABxAB <- formLinkageGroups(cross.18.all.ABxAB, min.lod=finLod,max.rf=finRf, reorgMarkers=TRUE)
cross.18.all.ABxAB <- orderMarkers(cross.18.all.ABxAB,window=5,use.ripple=T,
  error.prob=0.2, map.function='kosambi',sex.sp=F,maxit=100,tol=1e-2)
cross.18.all.ABxAB$pheno$homP <- rowSums(cross.18.all.ABxAB$geno$'1'[[1]]=='2',na.rm=T)/rowSums(!is.na(cross.18.all.ABxAB$geno$'1'[[1]]))
cross.18.all.ABxAB$pheno$hets1 <- rowSums(cross.18.all.ABxAB$geno$'1'[[1]]=='2',na.rm=T)/rowSums(!is.na(cross.18.all.ABxAB$geno$'1'[[1]]))
cross.18.all.ABxAB$pheno$hom1 <- rowSums(cross.18.all.ABxAB$geno$'1'[[1]]=='1',na.rm=T)
cross.18.all.ABxAB$pheno$hom2 <- rowSums(cross.18.all.ABxAB$geno$'1'[[1]]=='3',na.rm=T)
cross.18.all.ABxAB$pheno$homP1 <- rowSums(cross.18.all.ABxAB$geno$'1'[[1]]=='1',na.rm=T)/rowSums(!is.na(cross.18.all.ABxAB$geno$'1'[[1]]))
cross.18.all.ABxAB$pheno$homP2 <- rowSums(cross.18.all.ABxAB$geno$'1'[[1]]=='3',na.rm=T)/rowSums(!is.na(cross.18.all.ABxAB$geno$'1'[[1]]))
cross.18.all.ABxAB$pheno$homP3 <- rowSums(cross.18.all.ABxAB$geno$'1'[[1]]=='1' |
  cross.18.all.ABxAB$geno$'1'[[1]]=='3',na.rm=T) /rowSums(!is.na(cross.18.all.ABxAB$geno$'1'[[1]]))

cross.18.all.ABxAB$pheno$XO <- countXO(cross.18.all.ABxAB)

cross.18.all.ABxAB2 <- drop.missing(cross.18.all.ABxAB,66)


write.cross(cross.18.all.AB,filestem='~/debug.csv',format="csv")
save.image('~/debug.Rsave')





g <- pull.geno(cross.18.all.ABxAB2))
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,3), las=1)
for (i in 1:3){
 plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i], ylim=c(0,1))
}



pr <- convert2cross2(cross.18.all.ABxAB)
map <- insert_pseudomarkers(pr$gmap, step=1)
pr <- calc_genoprob(pr, map, error_prob=0.02, cores=1)
kinship <- calc_kinship(pr)
grid <- calc_grid(map, step=1)
pr_grid <- probs_to_grid(pr, grid)
kinship_grid <- calc_kinship(pr_grid)



## Conservative
print('Dropping markers with segregation distortion < 0.0005')
cross.18 <- distort(cross.18,0.0005)

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

#### Set sex phenotype
cross.18$pheno$sex <- as.numeric(nmissing(cross.18a)[cross.18a$pheno$ID]>850)
y <- pull.pheno(cross.18, 2)

cross.18f <- subset(cross.18,ind=(y==0))
cross.18f <- formLinkageGroups(cross.18f, max.rf=finRf, min.lod=finLod, reorgMarkers=TRUE)
cross.18a <- subset(cross.18, chr=1)

cross.18f <- drop.missing(cross.18f,5)
#names(cross.18f$geno) <- X
cross.18f <- orderMarkers(cross.18f,chr=X,window=5,use.ripple=T,
  error.prob=0.1, map.function='kosambi',sex.sp=F,maxit=10000,tol=1e-2)
#######################



cross.18m <- subset(cross.18,ind=(y==1))
cross.18m <- formLinkageGroups(cross.18m, max.rf=finRf, min.lod=finLod, reorgMarkers=TRUE)
cross.18m <- drop.missing(cross.18m,10)
names(cross.18m$geno) <- X
cross.18m <- orderMarkers(cross.18m,chr=X,window=5,use.ripple=T,
    error.prob=0.1, map.function='kosambi',sex.sp=F,maxit=10000,tol=1e-2)

## rename to the correct LG
names(cross.18$geno) <- X

cross.18 <- sim.geno(cross.18,step=5,error.prob=0.1,map.function='kosambi')
cross.18 <- calc.genoprob(cross.18,step=5,error.prob=0.1,map.function='kosambi')

print('removing double cross-overs')
cross.18 <- removeDoubleXO(cross.18, verbose=T)
print('Done removing dxo..')

print('Estimating the initial map with high errorprob')
POS.map.18 <- est.map(cross.18,error.prob=0.1,map.function="kosambi", chr=X,maxit=100)
cross.18 <- replace.map(cross.18, POS.map.18)

print(summary(pull.map(cross.18))[as.character(X),])

print('Writing the markers to rQTL format')
write.cross(cross.18,filestem=paste(popdir,'/chr',X,'_',outname,'.QTLmap',sep=''),format="csv",chr=X)

print('saving...')
rm(cross.18)
save.image(paste(popdir,'/chr',X,'_',outname,'.QTLmap.Rsave',sep=''))
