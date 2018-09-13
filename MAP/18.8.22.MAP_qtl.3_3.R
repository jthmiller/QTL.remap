  #!/bin/R
### Map QTLs 3 of 3
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/pop_control_file.R')

load(paste(popdir,'/chr',X,'.QTLmap.Rsave',sep=''))

cross.18 <- read.cross(format='csv',dir=popdir,
   file=paste('chr',X,'.QTLmap.csv',sep=''),
   geno=c('AA','AB','BB'),alleles=c("A","B"))

marker.warning()

print('Re-setimating map from filtered data on')

cross.18 <- orderMarkers(cross.18,chr=X,window=5,use.ripple=T,
  error.prob=ers, map.function='kosambi',sex.sp=F,maxit=1000,tol=1e-3)

print('2nd time removing double cross-overs once more')
  cross.18 <- removeDoubleXO(cross.18, verbose=T)
  print('Done removing dxo..')

print('Re-estimating the final map with many iterations...')
POS.map.18 <- est.map(cross.18,error.prob=ers,map.function="kosambi",n.cluster=12, chr=X,maxit=10000)
cross.18 <- replace.map(cross.18, POS.map.18)

print('Done mapping..')
print(summary(pull.map(cross.18))[as.character(X),])

print('Re-writing the markers to rQTL format')
write.cross(cross.18,filestem=paste(popdir,'/chr',X,'.QTLmap',sep=''),format="csv",chr=X)

print('Re-estimating final error rate for QTL mapping')
ers <- er.rate(cross.18)
print(paste(ers,' error rate'))

print('Re-adding un-genotyped individuals for stratified analysis')
pheno.all <- phen <- read.table('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/metadata/ALL_phenotype_Dist.txt',header=T)
phen$pheno_all[which(phen$pheno_all<2)] <- 0
phen$pheno_all[which(phen$pheno_all>1)] <- 1
index <- which(phen$pop_all==pop)
zeros <- as.numeric((table(phen$pheno_all[index])-table(cross.18$pheno$Pheno))['0'])
ones <- as.numeric((table(phen$pheno_all[index])-table(cross.18$pheno$Pheno))['1'])
no_genos <- data.frame(pheno=c(rep(0,times=zeros),rep(1,times=ones)),
              sex=rep(0,times=zeros+ones),ind=paste('NG',1:as.numeric(zeros+ones),sep=''),
              markers= matrix('-',nrow=zeros+ones,ncol=as.numeric(nmar(cross.18))))

write.table(no_genos,file=file.path(popdir,'no_genos.csv'),
  col.names=F,row.names=F,quote=F,sep=',')

system(paste('cat chr',Z,'.QTLmap.csv no_genos.csv > temp.',Z,sep=''))

print('saving... done with mapping ind chromosomes')
rm(cross.18)
save.image(paste('chr',X,'.QTLmap.Rsave',sep=''))
