#!/bin/R
### Map QTLs 1 of 3
debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")

library('qtl')

mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

################################################################################
## ADD AHR GENOTYPES ##
################################################################################
fl <- file.path(mpath,'ELR_unmapped_filtered.csv')
cross.df <- read.csv(fl,header=FALSE,stringsAsFactors=F)
marks <- cross.df[1,6:length(cross.df[1,])]
gts <- cross.df[4:length(cross.df[,1]),6:length(cross.df[1,])]
rownames(gts) <- cross.df[c(4:length(cross.df[,1])),'V3']
nmars <- length(gts[1,])
phen <- cross.df[,1:5]
rownames(phen) <- c('info','chr','map',c(cross.df[4:length(cross.df[,1]),'V3']))


fla <-file.path(mpath, 'ER_ahr_aip_whoi_gt.csv')
cross.df.ahr <- read.csv(fla,header=FALSE,stringsAsFactors=F)
ahr.marks <- cross.df.ahr[1,4:length(cross.df.ahr[1,])]
ahr.gts <- cross.df.ahr[4:length(cross.df.ahr[,1]),4:length(cross.df.ahr[1,])]
rownames(ahr.gts) <- cross.df.ahr[c(4:length(cross.df.ahr[,1])),'V3']

phen.ah <- cross.df.ahr[,1:3]
rownames(phen.ah) <- c('info','chr','map',phen.ah[4:length(phen.ah[,1]),'V3'])
phen.ah

gts.2 <- cbind(gts, ahr.gts[rownames(gts),])

new <- rownames(ahr.gts)[!rownames(ahr.gts) %in% rownames(gts)]
a <- matrix('-',ncol=nmars,nrow=length(new))
a <- cbind(a,ahr.gts[new,])
a <- cbind(phen.ah[rownames(a),],NA,NA,a)

colnames(a) <- c('Pheno','sex','ID','bin','pheno_norm',final.marks)
gts.2 <- cbind(phen[rownames(gts.2),],gts.2)
colnames(gts.2) <- c('Pheno','sex','ID','bin','pheno_norm',final.marks)
final.gts <- rbind(gts.2,a)

row1 <- colnames(gts.2)
row2 <- c(cross.df[2,],c(1,2,2))
names(row2) <- colnames(gts.2)
row3 <- c(cross.df[3,],c(0,0,0))
names(row3) <- colnames(gts.2)

final.gts <- rbind(row1,row2,row3,final.gts)

fl <- file.path(mpath,'ELR_unmapped_filtered_added_markers.csv')
write.table(final.gts, fl,col.names=F,row.names=F,quote=F,sep=',')

################################################################################

chr1 <- gsub("AHR2a_del","350000",markernames(cross,chr=1))
ord <- order(as.numeric(gsub('.*:','',chr1)))
cross <- switch.order(cross, chr = 1, ord, error.prob = 0.1, map.function = "kosambi",
    maxit = 1000, tol = 0.001, sex.sp = F)
cross <- est.map(cross, error.prob=0.1, map.function="kosambi",sex.sp=F,chr=1)

chr2 <- gsub("AIP_252","37860632",markernames(cross,chr=2))
ord <- order(as.numeric(gsub('.*:','',chr2)))
cross <- switch.order(cross, chr = 2, ord, error.prob = 0.11, map.function = "kosambi",
    maxit = 1000, tol = 0.001, sex.sp = F)

cross <- est.map(cross, error.prob=0.1, map.function="kosambi",sex.sp=F,chr=18)


################################################################################
for (i in 3:24){

 ord <- order(as.numeric(gsub('.*:','', markernames(cross.reord,chr=i))))

 cross.reord <- switch.order(cross.reord, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
     maxit = 1000, tol = 0.001, sex.sp = F)
}

save.image(file.path(mpath,'ER_chrFR.rsave'))
################################################################################
################################################################################
for (i in 1:24){

 cross <- est.map(cross, error.prob=0.1, map.function="kosambi",sex.sp=F,chr=i)
 print(i)
 save.image(file.path(mpath,'ER_chrFR.rsave'))

}
################################################################################
