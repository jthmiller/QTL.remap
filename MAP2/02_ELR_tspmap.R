#!/bin/R
### Map QTLs 1 of 3
#debug.cross <- T
#source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
#Sys.setenv(TAR = "/bin/tar")
##install_github("kbroman/ASMap")
##install_github("jtlovell/qtlTools")
##install_github("mckaylab/TSPmap")

i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]

libs2load<-c('devtools','qtl',"ASMap","qtlTools","TSP","TSPmap")
suppressMessages(sapply(libs2load, require, character.only = TRUE))

pop <- 'ELR'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

fl <- file.path('ELR_unmapped_filtered.csv')
cross <- read.cross(file=fl,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)
cross <- subset(cross,chr=i)

################################################################################
## SEX

sex <- read.table("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/data/sex.txt",stringsAsFactors=F)
rownames(sex) <- sex$ID
sex.vec <- sex[as.character(cross$pheno$ID), 'sex']
cross$pheno$sex <- sex.vec

################################################################################

ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross)[[as.character(i)]]))))
cross <- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
 maxit = 1, tol = 0.1, sex.sp = F)

png(paste0('~/public_html/ELR_gts_preclean',i,'.png'),height=2500,width=4500)
 plotGeno(cross, chr=i, cex=2)
dev.off()

################################################################################

if (any(i %in% c(1,2))){

 fla <-file.path(mpath, 'ER_ahr_aip_whoi_gt.csv')
 cross.ahr <- read.cross(file = fla,format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)

 cross.ahr$pheno$bin <- ifelse(cross.ahr$pheno$Pheno > 2, 1 , 0)
 cross.ahr$pheno$pheno_norm <- round(nqrank(cross.ahr$pheno$Pheno))

 ahr2 <- pull.geno(cross.ahr)[,"AHR2a_del"]
 aip252 <- pull.geno(cross.ahr)[,"AIP_252"]
 aip261 <- pull.geno(cross.ahr)[,"AIP_261"]

 add_gts <- data.frame(ahr2,aip252,aip261,stringsAsFactors=F)
 rownames(add_gts) <- cross.ahr$pheno$ID
 add_gts <- add_gts[as.character(cross$pheno$ID),]

 if (i==1){ cross <- addmarker(cross,add_gts[,'ahr2'],'ahr2a',chr=1,pos=1) }
 if (i==2){
  cross <- addmarker(cross,add_gts[,'aip261'],'aip261',chr=2,pos=1)
  cross <- addmarker(cross,add_gts[,'aip252'],'aip252',chr=2,pos=1.5)
 }
}

################################################################################
cross <- cleanGeno_jm(cross, chr=i, maxdist=100, maxmark=8, verbose=TRUE)
cross <- calc.errorlod(cross, err=0.03)
cross <- removeDoubleXO(cross, chr=i)
cross <- calc.errorlod(cross, err=0.03)
cross <- cleanGeno_jm_2(cross, chr=i, maxdist=50, maxmark=4, verbose=TRUE)
cross <- calc.errorlod(cross, err=0.03)
################################################################################

################################################################################
gt <- geno.table(cross)
bfixA <- rownames(gt[which(gt$P.value > 0.0001 & gt$missing < 6),])
cross <- pull.markers(cross,bfixA)
################################################################################

cross <-tspOrder(cross = cross,hamiltonian = TRUE, method="concorde",concorde_path='/home/jmiller1/concorde_build/TSP/')

cross <- shiftmap(cross, offset=0)

cross_map <-  est.map(cross, error.prob=0.03,map.function="kosambi",maxit=100000,tol=1e-7, sex.sp=FALSE, verbose=FALSE)

cross <- qtl:::replace.map(cross,cross_map)

filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_all_mark_',i,'_tsp')
write.cross(cross,chr=i,filestem=filename,format="csv")

png(paste0('~/public_html/ELR_RF_concord',i,'_tsp.png'))
  plotRF(cross)
dev.off()
