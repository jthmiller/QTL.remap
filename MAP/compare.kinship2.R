source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')

load(paste('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/',pop,'/REMAPS/QTLmap.Rsave',sep=''))

##### Re-do kinship analysis with markers that mapped
##### Drop idividuals that are to closely/distantly related and remap everything

cross <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
        format='csvr', geno=c(1:3),estimate.map=FALSE)

cross.18 <- subset(cross.18,ind=cross.18$pheno$stata=='ind')

path <- file.path(indpops,paste(pop,'.ped',sep=''))
popname <- system(paste('cut -f1 -d\' \'',path),intern = TRUE)
indname <- system(paste('cut -f2 -d\' \'',path),intern = TRUE)
cross$pheno$ID <- paste(popname,indname,sep='_')

drop <- markernames(cross)[!markernames(cross) %in% markernames(cross.18)]
cross  <- drop.markers(cross,drop)

### Taken from Karl Broman rQTL discussion
### https://groups.google.com/forum/#!searchin/rqtl-disc/Marker$20names$20don$27t$20match%7Csort:date/rqtl-disc/wcZZh0lfTiI/PLcpTV42yjQJ
# locations of x2 markers in cross x1
pos <- find.markerpos(cross.18, markernames(cross))
# for those found, move to position in x1
for(i in which(!is.na(pos[,1]))) cross <- movemarker(cross, rownames(pos)[i], pos[i,1], pos[i,2])
# perhaps drop the markers that weren't in x1
cross <- drop.markers(cross, rownames(pos)[is.na(pos[,1])])
# probably need to sort the chromosomes
cross$geno <- cross$geno[names(cross.18$geno)]
# use map to replace parents map
map <- pull.map(cross.18)
cross <- replace.map(cross,map)

rela <- comparegeno(cross)
colnames(rela) <- cross$pheno$ID
rownames(rela) <- cross$pheno$ID
rela[rela==NaN] <- NA
diag(rela) <- NA

## remove any without data
rela <- rela[rowSums(is.na(rela)) < nind(cross),colSums(is.na(rela)) < nind(cross)]

main = paste(pop,'kinship of library (proportion of shared genotypes, 0-1)')

png(file.path(popdir,paste(pop,'_kinship_histogram_all.png',sep='')))
hist(rela,main=main,sub=sum(nmar(cross)))
dev.off()

png(file.path(popdir,paste(pop,'_kinship_heatmap_all.png',sep='')))
heatmap(rela,symm=T,main=main,sub="red 0 to 1.0 yellow/white" )
dev.off()


cross.pars <- subset(cross,ind=is.na(cross$pheno$Pheno))
cross.18 <- c(cross.18,cross.pars)

rela <- comparegeno(cross.18)
colnames(rela) <- cross.18$pheno$ID
rownames(rela) <- cross.18$pheno$ID
rela[rela==NaN] <- NA
diag(rela) <- NA

rela <- rela[rowSums(is.na(rela)) < nind(cross),colSums(is.na(rela)) < nind(cross)]

main = paste(pop,'kinship in final (proportion of shared genotypes, 0-1)')

png(file.path(popdir,paste(pop,'_kinship_histogram_final.png',sep='')))
hist(rela,main=main,sub=sum(nmar(cross.18)))
dev.off()

png(file.path(popdir,paste(pop,'_kinship_heatmap_final.png',sep='')))
heatmap(rela,symm=T,main=main,sub="red 0 to 1.0 yellow/white" )
dev.off()

save.image(paste('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/',pop,'/REMAPS/kinship_QTLmap.Rsave',sep=''))








save.image(paste('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/',pop,'/REMAPS/kinship_QTLmap.Rsave',sep=''))

cross.18 <- orderMarkers(cross.Z,window=5,use.ripple=T,
  error.prob=0.05, map.function='kosambi',sex.sp=F,maxit=2000,tol=1e-2)


save.image(paste('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/',pop,'/REMAPS/kinship_QTLmap.Rsave',sep=''))



#drop <- markernames(cross.18)[!markernames(cross.18) %in% markernames(cross)]
#cross.18 <- drop.markers(cross.18,drop)
#cross.Z <- removeDoubleXO(cross.Z, verbose=T)
#cross.18 <- orderMarkers(cross.Z,window=5,use.ripple=T,
#  error.prob=0.05, map.function='kosambi',sex.sp=F,maxit=2000,tol=1e-2)
#cross.18 <- repRipple(cross.18, error.prob=ers, map.function="kosambi",window = 6)
#POS.map.18 <- est.map(cross.18,error.prob=ers,map.function="kosambi",maxit=1000)
#cross.18 <- replace.map(cross.18, POS.map.18)
#cross.18 <- jittermap(cross.18, amount=1e-6)
#cross.18 <- sim.geno(cross.18,error.prob=ers)
#cross.pars <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
#        format='csvr', geno=c(1:3),estimate.map=FALSE)
#cross.pars <- subset(cross,ind=is.na(cross$pheno$Pheno))
#cross.pars$pheno$ID <-  c('par1','par2')
#drop <- markernames(cross.pars)[!markernames(cross.pars) %in% markernames(cross.18)]
#cross.pars  <- drop.markers(cross.pars,drop)
#cross.pars <- c(cross.pars,cross.18)

## Only use genotyped individuals to select markers
cross <- subset(cross.18,ind=cross.18$pheno$stata=='ind')
cross <- dropSimilarMarkers(cross,chr=X,rf.threshold = 0.002,byChr = TRUE,re.est.map = FALSE)
cross <- repRipple(cross, error.prob=ers, map.function="kosambi",window = 6)
cross <- removeDoubleXO(cross, verbose=T)
cross.map <- est.map(cross,error.prob=ers,map.function="kosambi",maxit=1000)
cross <- replace.map(cross, cross.map)
#cross <- est.rf(cross)
cross.2 <- repRipple(cross, error.prob=ers, map.function="kosambi",window = 6)
cross.18 <- subset(cross.18,ind=cross.18$pheno$stata=='ind')
