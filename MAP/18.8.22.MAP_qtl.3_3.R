#!/bin/R
### Map QTLs 3 of 3
setwd('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/')
source('control_file.R')

#load(paste(popdir,'/chr',X,'_',outname,'.QTLmap.Rsave',sep=''))

cross.18 <- read.cross(format='csv',dir=popdir,
   file=paste('chr',X,'_',outname,'_2.QTLmap.csv',sep=''),
   geno=c('AA','AB','BB'),alleles=c("A","B"))

marker.warning()

print('dropping markers by error lod')

print('2nd time removing double cross-overs once more')
  cross.18 <- removeDoubleXO(cross.18, verbose=T)
print('Done removing dxo..')

print('Re-setimating map from filtered data on')
cross.18 <- orderMarkers(cross.18,chr=X,window=5,use.ripple=T,
  error.prob=ers, map.function='kosambi',sex.sp=F,maxit=5000,tol=1e-3)

print('Re-estimating the final map with many iterations...')
POS.map.18 <- est.map(cross.18,error.prob=ers,map.function="kosambi", chr=X,maxit=5000)
cross.18 <- replace.map(cross.18, POS.map.18)
print('Done mapping..')

print(summary(pull.map(cross.18))[as.character(X),])

print('Re-estimating error rate for QTL mapping')
ers <- er.rate(cross=cross.18,cpus=slurmcore,maxit=1000)
print(paste(ers,' error rate'))

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

trsl.bin <- c(0,0,0,1,1,1)
names(trsl.bin) <- as.character(0:5)

no_genos <- data.frame(Pheno=as.numeric(trsl.bin[as.character(count.pheno)]),sex=0,
              ID=paste('NG',1:length(count.pheno),sep=''),Pheno_05=count.pheno,
              markers= matrix('-',nrow=length(count.pheno),ncol=as.numeric(nmar(cross.18))))

phenos <- data.frame(count.pheno,Pheno=as.numeric(trsl.bin[as.character(count.pheno)]))
rownames(phenos) <- paste('NG',1:length(count.pheno),sep='')

write.table(no_genos,file=file.path(popdir,'no_genos.csv'),
  col.names=F,row.names=F,quote=F,sep=',')

system(paste('cat ',popdir,'/chr',X,'_',outname,'_3.QTLmap.csv ',popdir,'/no_genos.csv > ',popdir,'/temp.',X,sep=''))

print('saving... done with mapping ind chromosomes')
rm(cross.18)
save.image(paste(popdir,'/chr',X,'_',outname,'.QTLmap.Rsave',sep=''))
