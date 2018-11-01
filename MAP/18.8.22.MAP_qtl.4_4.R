#!/bin/R
### Map QTLs 4 of 3
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')
load(paste(popdir,'/chr',X,'_',outname,'.QTLmap.Rsave',sep=''))

cross.18 <- read.cross(format='csv',dir=popdir,
   file=paste('/chr',X,'_',outname,'_3.QTLmap.csv',sep=''),
   geno=c('AA','AB','BB'),alleles=c("A","B"))

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

write.table(no_genos,file=file.path(popdir,paste(X,'no_genos.csv',sep='_')),
  col.names=F,row.names=F,quote=F,sep=',')

write.cross(cross.18,filestem=paste(popdir,'/chr',X,'_',outname,'_3.QTLmap',sep=''),format="csv",chr=X)

system(paste('cat ',popdir,'/chr',X,'_',outname,'_3.QTLmap.csv ',popdir,'/',paste(X,'no_genos.csv',sep='_'),' > ',popdir,'/temp.',X,sep=''))

print('saving... done with mapping ind chromosomes')

save.image(paste(popdir,'/chr',X,'_',outname,'.QTLmap.Rsave',sep=''))
