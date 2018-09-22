#####

##############
for (X in 1:24){
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

## Specific to 'outname'
if(mapped.only==T){
cross.18 <- drop.markers(cross.18,markernames(cross.18)[grep('NW',markernames(cross.18))])
}

gt.b4 <- geno.table(cross.18)
pos <- as.numeric(gsub(paste(X,':',sep=''),'',rownames(gt.b4)))
pval <- log10(gt.b4$P.value)

## Conservative
print('Dropping markers with more than 40 genotypes missing')
cross.18 <- drop.missing(cross.18,40)

gt.af <- geno.table(cross.18)

## based on distribution of P-val for seg dist from ABxAB cross
cutoff <- 0.0001
print(paste('Dropping markers with segregation distortion < ',cutoff))
#cross.18 <- distort(cross.18,0.01)
cross.18 <- distort(cross.18,cutoff)

gt.dis <- geno.table(cross.18)

png(paste(X,'pval.png',sep=''))
hist.geno(gt.af$P.value)
abline(v=log10(cutoff))
dev.off()

png(paste(X,'.png',sep=''))
par(mfrow=c(4,1),mar=c(1,2,3,1))
plot.geno(gt.b4,gen.main='All Markers')
plot.geno(gt.af,gen.main='Filter loci missing > 40')
abline(h=log10(cutoff))
plot.geno(gt.dis,gen.main=paste('Filter distortion > ',cutoff))
plot.geno(gt.nomis,gen.main='Filter loci missing > 5')
dev.off()

}
