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

marks_nms <- cross.df[1,6:length(cross.df[1,])]
gts <- cross.df[4:length(cross.df[,1]),6:length(cross.df[1,])]
## 88 X 19856
rownames(gts) <- cross.df[c(4:length(cross.df[,1])),'V3']
colnames(gts) <- as.character(cross.df[1,6:length(cross.df[1,])])
## gts 88x19856
nmars <- length(gts[1,])

phenotpyes <- cross.df[,1:5]

rownames(phenotpyes ) <- c('info','chr','map',c(cross.df[4:length(cross.df[,1]),'V3']))
### phen 89x5

fla <-file.path(mpath, 'ER_ahr_aip_whoi_gt.csv')
cross.df.ahr <- read.csv(fla,header=FALSE,stringsAsFactors=F)
ahr_mark_nms <- cross.df.ahr[1,4:length(cross.df.ahr[1,])]

## AHR genotypes
ahr_gts <- cross.df.ahr[4:length(cross.df.ahr[,1]),4:length(cross.df.ahr[1,])]
rownames(ahr_gts) <- cross.df.ahr[c(4:length(cross.df.ahr[,1])),'V3']
colnames(ahr_gts) <- cross.df.ahr[1,c(4:6)]
### ahr.gts 87x3


##pheno data with blank rows included
phen.ah <- cross.df.ahr[,1:3]
rownames(phen.ah) <- c('info','chr','map',phen.ah[4:length(phen.ah[,1]),'V3'])
## phen.ah 90x3

gts.2 <- cbind(gts, ahr_gts[rownames(gts),])
## 88x19859  CONTAINS ATSM and AHR marks

new <- rownames(ahr_gts)[!rownames(ahr_gts) %in% rownames(gts.2)]

a <- matrix('-',ncol=nmars,nrow=length(new))
a <- cbind(a,ahr_gts[new,])
## 2x19859  CONTAINS ATSM and AHR marks
a <- cbind(phen.ah[rownames(a),],NA,NA,a)
colnames(a) <- c('Pheno','sex','ID','bin','pheno_norm',colnames(gts.2))
# 2x19864 (a incluedes ahr genos)

gts.2 <- cbind(phen.ah[rownames(gts.2),],gts.2)
colnames(gts.2)[1:5] <- c('Pheno','sex','ID','bin','pheno_norm')

final.gts <- rbind(gts.2[colnames(gts.2)],a)

row1 <- colnames(gts.2)
row2 <- c(cross.df[2,],c(1,2,2))
names(row2) <- colnames(gts.2)
row3 <- c(cross.df[3,],c(0,0,0))
names(row3) <- colnames(gts.2)

final.gts <- rbind(row1,row2,row3,final.gts)

fl <- file.path(mpath,'ELR_unmapped_filtered_added_markers.csv')
write.table(final.gts, fl,col.names=F,row.names=F,quote=F,sep=',')

################################################################################
################################################################################
fl <- file.path(mpath,'ELR_unmapped_filtered_added_markers.csv')

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
 estimate.map = FALSE
)

################################################################################
################################################################################
################################################################################
cross.all <- cross

for (i in 1:24){
cross <- subset(cross,chr=i)
nmars <- nmar(cross)
cross <- subset(cross,ind=nmissing(cross) < (nmars*.5))

## initial order
 ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross)[[1]]))))

cross <- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
 maxit = 10, tol = 0.001, sex.sp = F)

### CROSS FOR LOCAT. BAD MARKERS #####################################
loc.these <- table(unlist(locateXO(cross)))
loc.these <- names(loc.these[as.numeric(loc.these) > 10])
drops <- sapply(as.numeric(loc.these),function(X){ find.marker(cross, chr=i, pos=X) } )

cross <- drop.markers(cross,drops)
cross <- calc.genoprob(cross)
cross <- sim.geno(cross)
cross <- calc.errorlod(cross, err=0.01)

png(paste0('~/public_html/ELR_gts_2.png'),width=2000)
plotGeno(cross.nodup)
dev.off()

#### REMOVE DUPS
cross <- removeDoubleXO(cross)
cross <- calc.genoprob(cross)
cross <- sim.geno(cross)
cross <- calc.errorlod(cross, err=0.01)

png(paste0('~/public_html/ELR_gts_1.png'),width=2000)
plotGeno(cross.nodup)
dev.off()

loc.these <- table(unlist(locateXO(cross)))
loc.these <- names(loc.these[as.numeric(loc.these) > 4])
drops <- sapply(as.numeric(loc.these),function(X){ find.marker(cross, chr=i, pos=X) } )

cross <- drop.markers(cross,drops)
cross <- removeDoubleXO(cross)
cross <- calc.genoprob(cross)
cross <- sim.geno(cross)
cross <- calc.errorlod(cross, err=0.01)

png(paste0('~/public_html/ELR_gts_1.png'),width=4000,height=1000)
plotGeno(cross)
dev.off()







######
cross <- calc.genoprob(cross)
cross <- sim.geno(cross)
cross <- calc.errorlod(cross, err=0.01)

png(paste0('~/public_html/ELR_gts.png'),width=2000)
plotGeno(cross.nodup)
dev.off()

png(paste0('~/public_html/ELR_xo_1.png'))
hist(sort(table(unlist(locateXO(cross)))))
dev.off()

cross <- calc.errorlod(cross, err=0.01)
tel <- top.errorlod(cross)
drops <- unique(tel[tel$errorlod>10,3])



cross <- drop.markers(cross,drops)

cross <- clean(cross)

cross <- calc.errorlod(cross, err=0.01)
head(top.errorlod(cross),50)

sort(table(unlist(locateXO(cross))))


png(paste0('~/public_html/ELR_gts.png'),width=2000)
plotGeno(cross)
dev.off()



cross <- removeDoubleXO(cross)





mp <- pull.map(cross)


nmp <- names(mp[['1']])
rp <- gsub(".*:",'',nmp)
mp[["1"]] <- as.numeric(rp)
names(mp[["1"]]) <- nmp


cross <- replace.map(cross,mp)


loca <- lapply(mp, function(X){
 setNames(gsub("AHR2a_del","350000",names(X)),names(X))
})
loca <- lapply(loca, function(X){
 gsub("AIP_252","37860632",X)
})
loca <- lapply(loca, function(X){
 gsub("AIP_261","37860640",X)
})
loca <- lapply(loca, function(X){
 setNames(as.numeric(gsub('.*:','',X)),names(X))
})





ty <- markernames(cross)
ty <- as.numeric(gsub(".*:",'',ty))
names(ty) <- markernames(cross)

mp <- pull.map(cross)
#cross.ss <- cross
cross <- replace.map(cross,ty)
subs <- sapply(pull.map(cross.ss),pickMarkerSubset,2000)
subs <- unname(unlist(subs))
drops <- markernames(cross.ss)[! markernames(cross.ss) %in% subs]
cross.ss <- drop.markers(cross.ss,drops)

################################################################################
for (i in 1:24){

 ord <- order(as.numeric(pull.map(cross.ss)[[i]]))

 cross.ss <- switch.order(cross.ss, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
     maxit = 1000, tol = 0.001, sex.sp = F)
}
################################################################################

mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'
write.table(markernames(cross.ss),file.path(mpath,'ER_markers_subst.table'))

fl <- file.path(mpath,'ELR_subsetted')
write.cross(cross.ss,filestem=fl,format="csv")

################################################################################
