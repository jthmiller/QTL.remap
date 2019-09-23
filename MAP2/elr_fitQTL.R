#!/bin/R
### Map QTLs 1 of 3
#debug.cross <- T
#source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library('qtl')
pop <- 'ELR'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")

################################################################################
## unmapped
################################################################################
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

file_list <- gsub('.csv','',list.files(mpath, '*downsmpl_map*'))
file_list <- list.files(mpath, '*downsmpl_map*',include.dirs = T)

chr <- gsub("ELR_gts_CHR",'',file_list)
chr <- as.numeric(gsub("_downsmpl_map",'',chr))

elr <- lapply(file_list[-1],function(X){ read.cross(file=X,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)})

gnos <- lapply(elr,function(X){
 data.frame(X[[1]][[1]][['data']],stringsAsFactors=F)
})
gnos <- do.call(cbind,gnos)
gnos <- cbind(elr[[1]]$pheno,gnos)
gnos$ID <- as.character(gnos$ID)

ty <- as.character(gnos)

m_names <- unlist(sapply(elr,function(X){
 markernames(X)
}))

colnames(gnos) <- c('Pheno','sex','ID','bin','pheno_norm',m_names)
rownames(gnos) <- elr[[1]]$pheno$ID

map <- c(colnames(elr[[1]]$pheno),unname(unlist(sapply(elr,pull.map))))
chr <- c(colnames(elr[[1]]$pheno),gsub(":.*","",m_names))
info <- c(colnames(elr[[1]]$pheno),m_names)
headers <- rbind(info,chr,map)
colnames(headers) <- headers[1,]
headers[2:3,1:5] <- ''

headers.u <- unname(data.frame(headers,row.names=NULL,stringsAsFactors=FALSE))
gnos.u <- unname(data.frame(lapply(gnos.u, as.character),row.names=NULL,stringsAsFactors=FALSE))
colnames(headers.u) <- colnames(gnos.u) <- headers.u[1,]
to_write <- rbind(headers.u,gnos.u)


write.table(to_write,file.path(mpath,'elr.mapped.1_24.csv'),sep=',',row.names=F,quote=F,col.names = F)


fl <- file.path(mpath,'elr.mapped.1_24.csv')

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)



################################################################################

cross$pheno$bin <- ifelse(cross$pheno$Pheno > 2, 1 , 0)
cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno))

## rf.1 <- est.rf(subset(cross,chr=1), maxit=1000, tol=1e-6)
## rf.2 <- est.rf(subset(cross,chr=2), maxit=1000, tol=1e-6)
## rf.8 <- est.rf(subset(cross,chr=8), maxit=1000, tol=1e-6)
## rf.18 <- est.rf(subset(cross,chr=18), maxit=1000, tol=1e-6)
## rf.13 <- est.rf(subset(cross,chr=13), maxit=1000, tol=1e-6)
## save.image(file.path(mpath,'ER_chrFR.rsave'))

##################################################################################
#### read in the QTL cross
##################################################################################
##mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'
##fl <- file.path(mpath,'ELR_unmapped_filtered.csv')
##
##cross <- read.cross(
## file = fl,
## format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
## estimate.map = FALSE
##)
##################################################################################

## UNFILTERED
##################################################################################
#### read in the QTL cross
## cross <- read.cross.jm(file = file.path(indpops, paste(pop, ".unphased.f2.csvr",
##  sep = "")), format = "csvr", geno = c(1:3), estimate.map = FALSE)
##################################################################################
##################################################################################
##### Pull names from plinkfile
##path <- file.path(indpops, paste(pop, ".ped", sep = ""))
##popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
##indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
##cross$pheno$ID <- paste(popname, indname, sep = "_")
##################################################################################
##
###### PHENO #####################################################################
##cross$pheno$bin <- ifelse(cross$pheno$Pheno > 2, 1 , 0)
##cross$pheno$pheno_norm <- round(nqrank(cross$pheno$Pheno))
##################################################################################
##
## #### FITER ########################################################################
## marks <- read.table('/home/jmiller1/QTL_Map_Raw/ELR_final_map/goodmarks.rtable',stringsAsFactors=F)
## inds <- read.table('/home/jmiller1/QTL_Map_Raw/ELR_final_map/goodsamps.rtable',stringsAsFactors=F)
## ####################################################################################
## drops <- markernames(cross)[!markernames(cross) %in%  marks[,1]]
## cross <- drop.markers(cross,drops)
## cross <- subset(cross, ind= !is.na(cross$pheno$Pheno))
##
## ################################################################################
dups <- findDupMarkers(cross, exact.only = F, adjacent.only = F)
cross <- drop.markers(cross, unlist(dups))
##
cross <- subset(cross, ind=nmissing(cross)<5000)
cross <- calc.genoprob(cross,error.prob=0.05)
cross <- sim.geno(cross,error.prob=0.05)
################################################################################
perms.bin.em <- scanone(cross, method = "em", model = "binary", maxit = 100,
  n.perm = 5, pheno.col = 4, n.cluster = 2)

scan.bin.em <- scanone(cross, method = "em", model = "binary", pheno.col = 4)
scan.norm.em <- scanone(cross, method = "em", model = "normal", pheno.col = 5)
scan.bin.mr <- scanone(cross, method = "mr", model = "binary", pheno.col = 4)
scan.norm.mr <- scanone(cross, method = "mr", model = "normal", pheno.col = 5)
scan.norm.imp <- scanone(cross, method = "imp", model = "normal", pheno.col = 5)
scan.bin.imp <- scanone(cross, method = "imp", model = "binary", pheno.col = 4)


cbind(summary(scan.bin.em),a=summary(scan.norm.em)[,3],b=summary(scan.bin.mr)[,3],c=summary(scan.np.mr )[,3])
################################################################################

qtl <- makeqtl(cross, chr=18, pos=1587.8,  what="draws")

out.i.18 <- addqtl(cross, qtl=qtl, formula=y~Q1*Q2, method="imp",model="binary",pheno.col=4)
out.a.18 <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2, method="imp",model="binary",pheno.col=4)

out.i.18 <- addqtl(cross, qtl=qtl, formula=y~Q1*Q2, method="imp",model="normal",pheno.col=5)
out.a.18 <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2, method="imp",model="normal",pheno.col=5)

################################################################################
################################################################################

scan.norm.em <- scanone(cross.18, method = "em", model = "normal", maxit = 5000,
  pheno.col = 5)
### Normal scan on transformed phenotype w/Extended haley knott (better for
### selective/missing genos at non-gt'ed ind)
scan.norm.ehk <- scanone(cross.final.nodups, method = "ehk", model = "normal", maxit = 5000, pheno.col = 5)

### Normal scan on transformed phenotype fast haley knott (not robust to missing
### data. LOD inflation)
scan.norm.hk <- scanone(cross.18, method = "hk", model = "normal", maxit = 5000,
  pheno.col = 6)

################################################################################
################################################################################























##############################
##############################
pos <- as.numeric(gsub(".*:","",rownames(gt.missing)))
names(pos) <- rownames(gt.missing)
head(sort(abs(pos -  343835)))

1:317181

1:363497

crs.bk

chr1gts <- pull.geno(crs.bk, 1)

chr1phn <- pull.pheno(crs.bk, 1)



chr1gts <- pull.geno(cross.18, 1)

chr1phn <- pull.pheno(cross.18, 1)


AHR <- cbind(chr1phn,chr1gts[,'1:317181'],chr1gts[,'1:363497'])

AHR <- cbind(chr1phn,chr1gts[,'1:317181'],chr1gts[,'1:363497'])
AHR <- AHR[order(AHR[,1]),]

table(AHR[AHR[,1]<2,3])
table(AHR[AHR[,1]>2,3])

table(AHR[AHR[,1]==0,3])
table(AHR[AHR[,1]==1,3])
table(AHR[AHR[,1]==4,3])
table(AHR[AHR[,1]==5,3])

chr1.pars <- pull.geno(cross.pars, 1)
rbind(chr1.pars[,'1:317181'], chr1.pars[,'1:363497'])

TAKE THE HOMZYGOUS GENOTYPES FOR THE ONE PARENT AND SEE IF THEY TEND TOWARD 1:2:1 compared to
het in parent.

AHR[,1] <- as.factor(AHR[,1])
AHR[,2] <- as.factor(AHR[,2])
AHR[,3] <- as.factor(AHR[,3])

png('~/public_html/ER_AHR.png')
plot(table(AHR[AHR[,1]<2,2]))
dev.off()

for(
table(AHR[,1])

print("Removing duplicates")
##dups <- findDupMarkers(cross.18, exact.only = F, adjacent.only = F)
##cross.18 <- drop.markers(cross.18, unlist(dups))
##confirm ahr2a 343745   343931 AHR2a
##mid is 343835
