#!/bin/R
### Map QTLs 3 of 3
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')

#load(paste(popdir,'/chr',X,'_',outname,'.QTLmap.Rsave',sep=''))

cross.18 <- read.cross(format='csv',dir=popdir,
   file=paste('chr',X,'_',outname,'_2.QTLmap.csv',sep=''),
   geno=c('AA','AB','BB'),alleles=c("A","B"))

marker.warning()

print('dropping markers by error lod')

print('2nd time removing double cross-overs')
  cross.18 <- removeDoubleXO(cross.18, verbose=T)
print('Done removing dxo..')

dups <- findDupMarkers(cross.18, exact.only=FALSE, adjacent.only=FALSE)
### remove markers that are exactly the same.
cross.18 <- drop.markers(cross.18, unlist(dups))

if (reorder.marks==T){
print('Re-setimating map from filtered data on')
cross.18 <- orderMarkers(cross.18,chr=X,window=5,use.ripple=T,
  error.prob=ers, map.function='kosambi',sex.sp=F,maxit=10000,tol=1e-3)
}
#else {
#  print('Ripple at physical positions')
#  cross.18 <- repRipple.jm(cross, chr = chrnames(cross)[1], window = 7,method = "countxo",
#  verbose = T,map.function = "kosambi", sex.sp=F, clean1st = T, ripVerb = TRUE)
    #ripVanWink <- ripple(cross.18, X, window=6, method="likelihood",
    #   error.prob=ers, map.function="kosambi",maxit=100, tol=1e-6,
    #   sex.sp=F, verbose=TRUE, n.cluster=slurmcore)
    #cross.18 <- switch.order(cross.18, X, ripVanWink[1,])
#}

print('plotting LOD matrix')
png(file.path(popdir,paste(X,'_RF_FINAL.png',sep='')))
plotRF(cross.18,chr=chrnames(cross.18)[1],what='both',mark.diagonal=T,col.scheme="redblue")
dev.off()

print('Re-estimating the final map with many iterations...')
POS.map.18 <- est.map(cross.18,error.prob=ers,map.function="kosambi", chr=X,maxit=10000)
cross.18 <- replace.map(cross.18, POS.map.18)
print('Done mapping..')

print(summary(pull.map(cross.18))[as.character(X),])

print('Re-estimating error rate for QTL mapping')
ers <- er.rate(cross=cross.18,cpus=slurmcore,maxit=1000)
print(paste(ers,' error rate'))

### Add phenotype
cross.18$pheno$Pheno_05 <- cross.18$pheno$Pheno
cross.18$pheno$Pheno <- trsl.bin[as.character(cross.18$pheno$Pheno_05)]

print('Re-writing the markers to rQTL format')
write.cross(cross.18,filestem=paste(popdir,'/chr',X,'_',outname,'_3.QTLmap',sep=''),format="csv",chr=X)

print('Adding un-genotyped individuals for stratified analysis')
pheno.all <- phen <- read.table('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/metadata/ALL_phenotype_Dist.txt',header=T)
phen$Pheno_05 <- phen$pheno_all
index <- which(phen$pop_all==pop)

count.pheno <- sapply(0:5, function(pt){
      pt <- as.character(pt)
      total <- sum(phen$Pheno_05[index]==pt)
      incross <- sum(cross.18$pheno$Pheno_05==pt)
      return(as.numeric(total-incross))

  }
)
names(count.pheno) <- as.character(0:5)
count.pheno <- count.pheno[!is.na(count.pheno)]
count.pheno <- rep.int(names(count.pheno), times=as.numeric(count.pheno))

no_genos <- data.frame(Pheno=as.numeric(trsl.bin[as.character(count.pheno)]),sex=0,
              ID=paste('NG',1:length(count.pheno),sep=''),Pheno_05=count.pheno,
              markers= matrix('-',nrow=length(count.pheno),ncol=as.numeric(nmar(cross.18))))

phenos <- data.frame(count.pheno,Pheno=as.numeric(trsl.bin[as.character(count.pheno)]))
rownames(phenos) <- paste('NG',1:length(count.pheno),sep='')

write.table(no_genos,file=file.path(popdir,'no_genos.csv'),
  col.names=F,row.names=F,quote=F,sep=',')

system(paste('cat ',popdir,'/chr',X,'_',outname,'_3.QTLmap.csv ',popdir,'/no_genos.csv > ',popdir,'/temp.',X,sep=''))

print('saving... done with mapping ind chromosomes')

save.image(paste(popdir,'/chr',X,'_',outname,'.QTLmap.Rsave',sep=''))
