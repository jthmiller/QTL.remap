#!/bin/R

################################################################################
## ### Hist plot
png('~/public_html/cg.nbh.png')
hist(cpgt[lower.tri(cpgt)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cpgt[lower.tri(cpgt)])
dev.off()
################################################################################

for(i in 1:24){
tchr <- i
png(paste0('~/public_html/ER_phys_pos_dist_',i,'.png'),width=2000)
plot(pos[gts$chr==tchr],-log10(gts[gts$chr==tchr,'P.value']),pch=16)
points(pos.3[gts.3$chr==tchr & gts.3$missing<3],-log10(gts.3[gts.3$chr==tchr & gts.3$missing<3,'P.value']),pch=16, col='green')
points(pos.3[gts.3$chr==tchr & gts.3$AA<3],-log10(gts.3[gts.3$chr==tchr & gts.3$AA<3,'P.value']),pch=16, col='blue')
points(pos.3[gts.3$chr==tchr & gts.3$BB<3],-log10(gts.3[gts.3$chr==tchr & gts.3$BB<3,'P.value']),pch=16, col='red')
abline(h=2.5, col='red')
dev.off()
}

################################################################################
# FINAL MARKER SET ABOVE ###
################################################################################

png(paste0('~/public_html/ER_phys.png'),height=1000,width=1000)
plot(x0-0.75,y0)
points(xf2-0.5,yf2,col='purple',pch=16)
points(xf3-0.25,yf3,col='blue',pch=16)
points(xf4,yf4,col='green',pch=16)
dev.off()

## ALL PHYS POS COVERED TO HERE ################################################
################################################################################
png(paste0('~/public_html/RF_LOD_ER',inchr,'.png'))
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")
dev.off()

cross.18.reorg.lg <- formLinkageGroups(cross.18.reorg, max.rf = 0.1, min.lod = 15, reorgMarkers = TRUE)
cross.18.reorg.lg <- subset(cross.18.reorg.lg, chr=1)

pullgts <- pull.geno(cross.18.reorg.lg)
yf <- as.numeric(gsub(".*:",'',colnames(pullgts)))
xf <- as.numeric(gsub(":.*",'',colnames(pullgts)))

png(paste0('~/public_html/ER_phys.png'))
plot(xp18,yp18)
points(xf,yf,col='green',pch=16)
dev.off()

png('~/public_html/RF_ER1.png')
plotRF(cross.18.reorg.lg,chr=1:4)
dev.off()
print(inchr)

#switch8 <- markernames(cross.18.reorg.lg,chr=2)
##cross.18.reorg.lg <- switchAlleles(cross.18.reorg.lg, switch)
chr8 <- markernames(cross.18.reorg.lg,chr=1)
################################################################################










######################################################################################################
### QUICKSCAN MARKER REGRESSION
cross.18$pheno$bin <- ifelse(cross.18$pheno$Pheno > 2, 1 , 0)
scan.bin.mr <- scanone(cross.18, method = "mr", model = "binary", pheno.col = 4)
mr <- summary(scan.bin.mr)
mr[order(mr$lod),]
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################



matc <- sapply(as.numeric(wh[,2]),function(X){
 k <- arrayInd(X, dim(cpgt))
 colnames(cpgt)[k[,2]]
})

cbind(rownames(wh),matc)

s2 <- colnames(cpgt)[k[,2]]


s1 <- rownames(cpgt)[k[,1]]
### Drop same ind ran twice
d1 <- names(which.max(nmissing(cross.18)[c('ELR_10978','ELR_10987')]))
d2 <- names(which.max(nmissing(cross.18)[c('ELR_10982','ELR_10972')]))
d3 <- names(which.max(nmissing(cross.18)[c('ELR_10986','ELR_10977')]))
d4 <- names(which.max(nmissing(cross.18)[c('ELR_10973','ELR_10983')]))
drop.ind <- !cross.18$pheno$ID %in% c(d1,d2,d3,d4)

cross.18.d <- subset(cross.18, ind=drop.ind)
gts <- geno.table(cross.18.d)
pos <- as.numeric(gsub(".*:",'',rownames(gts)))

png(paste0('~/public_html/ER_pval8.png'),width=1000)
plot(pos[gts$chr==8], -log10(gts[gts$chr==8,'P.value']),pch=16)
dev.off()















pullgts <- pull.geno(cross.18)

hom <- apply(pullgts,1,function(X){ sum(X==1, na.rm=T) + sum(X==3, na.rm=T) })

homb <- apply(pullgts,1,function(X){sum(X==3, na.rm=T) })
names(homb) <- cross.18$pheno$ID
homa <- apply(pullgts,1,function(X){ sum(X==1, na.rm=T) })
names(homa) <- cross.18$pheno$ID
het <- apply(pullgts,1,function(X){ sum(X==2, na.rm=T)})
names(het) <- cross.18$pheno$ID
rat <- hom/het
names(rat) <- cross.18$pheno$ID

ord <- names(sort(rat))

png(paste0('~/public_html/ER_hom_het_rat.png'))
plot(het[ord])
points(homa[ord],col='red')
points(homb[ord],col='blue')
dev.off()

cross.18.b <- subset(cross.18, ind=rat < 2.55)

pullgts <- pull.geno(cross.18.b)
yf <- as.numeric(gsub(".*:",'',colnames(pullgts)))
xf <- as.numeric(gsub(":.*",'',colnames(pullgts)))

pullgts <- pull.geno(cross.18)
yp18 <- as.numeric(gsub(".*:",'',colnames(pullgts)))
xp18 <- as.numeric(gsub(":.*",'',colnames(pullgts)))

png(paste0('~/public_html/ER_phys.png'))
plot(xp18,yp18)
points(xf,yf,col='red')
dev.off()

gts <- geno.table(cross.18.b)
pos <- as.numeric(gsub(".*:",'',rownames(gts))

png(paste0('~/public_html/ER_phys.png'),width=1000)
plot(pos[gts$chr==8],gts[gts$chr==8,'missing'],pch=16)
dev.off()

png(paste0('~/public_html/ER_phys.png'),width=1000)
plotMissing(cross.18.b,chr=2)
dev.off()

cross.18.b <- subset(cross.18.b, ind = !is.na(cross.18.b$pheno$Phen))
gts <- geno.table(cross.18.b)
pos <- as.numeric(gsub(".*:",'',rownames(gts))

png(paste0('~/public_html/ER_pval8.png'),width=1000)
plot(pos[gts$chr==8], -log10(gts[gts$chr==8,'P.value']),pch=16)
dev.off()

cutoff <- 1e-02
cross.18 <- drop.markers(cross.18.b, rownames(gts[gts$P.value < cutoff,]))
gts <- geno.table(cross.18.c)
pos <- as.numeric(gsub(".*:",'',rownames(gts)))

png(paste0('~/public_html/ER_pval8.png'),width=1000)
plot(pos[gts$chr==8], -log10(gts[gts$chr==8,'P.value']),pch=16)
dev.off()


#### FIXING PHASE OF EACH LG ####################################################
inchr <- 8
cross.18.reorg <- subset(cross.18.d, chr=inchr)

rf <- pull.rf(cross.18.reorg)
lod <- pull.rf(cross.18.reorg, what="lod")

png('~/public_html/RF_LOD_ER8.png')
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")
dev.off()

cross.18.reorg.lg <- formLinkageGroups(cross.18.reorg, max.rf = 0.15, min.lod = 12, reorgMarkers = TRUE)

png('~/public_html/RF_ER1.png')
plotRF(cross.18.reorg.lg,chr=1:4)
dev.off()
print(inchr)

#switch8 <- markernames(cross.18.reorg.lg,chr=2)
##cross.18.reorg.lg <- switchAlleles(cross.18.reorg.lg, switch)
chr8 <- markernames(cross.18.reorg.lg,chr=1)
################################################################################

cross.18.d$pheno$bin <- ifelse(cross.18.d$pheno$Pheno > 2, 1 , 0)

scan.bin.mr <- scanone(cross.18.d, method = "mr", model = "binary", pheno.col = 4)
mr <- summary(scan.bin.mr)
mr[order(mr$lod),]


#### FIXING PHASE OF EACH LG ####################################################
inchr <- 1
cross.18.reorg <- subset(cross.18.d, chr=inchr)

rf <- pull.rf(cross.18.reorg)
lod <- pull.rf(cross.18.reorg, what="lod")

png('~/public_html/RF_LOD_ER1.png')
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")
dev.off()

cross.18.reorg.lg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 12, reorgMarkers = TRUE)

png('~/public_html/RF_ER1.png')
plotRF(cross.18.reorg.lg, chr=1:5)
dev.off()
print(inchr)

#switch8 <- markernames(cross.18.reorg.lg,chr=2)
##cross.18.reorg.lg <- switchAlleles(cross.18.reorg.lg, switch)
chr1 <- markernames(cross.18.reorg.lg,chr=1)














## Subset and drop parents
cross.pars <- subset(cross.18, ind = is.na(cross.18$pheno$Phen))
cross.pars.BI <- subset(cross.18, ind = 'BLI_BI1124M')

## Remove problematic individuals (found by kinship analysis)
#con <- file(file.path(popdir, "kinship.keep.ind.txt"), open = "r")
#keepers <- readLines(con)
#close(con)

#print("Dropping kinship outliers")
#cross.18 <- subset(cross.18, ind = cross.18$pheno$ID %in% keepers)

cross.18 <- subset(cross.18, ind = !is.na(cross.18$pheno$Phen))
gt.8 <- geno.table(cross.18,chr=8)

##############################
## Get parent fixed gts
cross.pars.BI <- subset(cross.pars, ind = 'BLI_BI1124M')
gt.cross.BI <- geno.table(cross.pars.BI)
drop.BI <- rownames(gt.cross.BI)[which(gt.cross.BI$AA==1 | gt.cross.BI$BB==1)]

## Drop everything not fixed in BI
gt.missing <- geno.table(cross.18)
drop.BI <- rownames(gt.missing)[!rownames(gt.missing) %in% drop.BI]
cross.18 <- drop.markers(cross.18, drop.BI)

## Drop from parent cross
cross.pars.BI <- drop.markers(cross.pars.BI, drop.BI)
gt.cross.BI <- geno.table(cross.pars.BI)
gt.missing <- geno.table(cross.18)

## Switch AA BI genotypes to Bs (BI phase is now B)
gt.cross.pars <- geno.table(cross.pars.BI)[markernames(cross.18), ]
swit <- rownames(gt.cross.pars[gt.cross.pars$AA == 1, ])

cross.pars.BI <- switchAlleles(cross.pars.BI, markers = swit)
cross.18 <- switchAlleles(cross.18, markers = swit)

gt.cross.BI <- geno.table(cross.pars.BI)
gt.missing <- geno.table(cross.18)

apply(gt.cross.BI[gt.cross.BI$chr==8,c(3:5)],2,table)

##############################
table(gt.cross.BI$chr


par.gts <- pull.geno(cross.18)









pullgts <- pull.geno(cross.18)
yp18 <- as.numeric(gsub(".*:",'',colnames(pullgts)))
xp18 <- as.numeric(gsub(":.*",'',colnames(pullgts)))

pullgts <- pull.geno(cross.pars.BI)
yp <- as.numeric(gsub(".*:",'',colnames(pullgts)))
xp <- as.numeric(gsub(":.*",'',colnames(pullgts)))

pullgts <- pull.geno(crs.bk)
yf <- as.numeric(gsub(".*:",'',colnames(pullgts)))
xf <- as.numeric(gsub(":.*",'',colnames(pullgts)))


png(paste0('~/public_html/ER_phys.png'))
plot(xp18,yp18)
points(xf,yf,col='red')
points(xp,yp,col='green')
dev.off()


g

png(paste0('~/public_html/ER_phys.png'))
plot(xp,yp)
dev.off()\]

save.image('~/ER_toss_many.rsave')

## how related are samples to the BI founder
cpgt <- comparegeno(cross.18)
colnames(cpgt) <- cross.18$pheno$ID
rownames(cpgt) <- cross.18$pheno$ID
cpgt[cpgt==NaN] <- NA
diag(cpgt) <- NA
cpgt <- cpgt[rowSums(is.na(cpgt)) < nind(cross.18),colSums(is.na(cpgt)) < nind(cross.18)]

drops <- names(which(cpgt['BLI_BI1124M',] > .6 | cpgt['BLI_BI1124M',] < .4))

cross.18 <- subset(cross.18, ind = cross.18$pheno$ID[!cross.18$pheno$ID %in% drops])
############################################################
## missing
##gt.missing <- geno.table(cross.18)

png('~/public_html/ER_missing.png')
hist(gt.missing$missing, breaks=20)
dev.off()

cross.18 <- drop.markers(cross.18, rownames(gt.missing)[gt.missing$missing > 2])
gt.missing <- geno.table(cross.18)
############################################################

cutoff <- 1e-03

cross.18 <- drop.markers(cross.18, rownames(gt.missing[gt.missing$P.value < cutoff,
  ]))



dups <- findDupMarkers(cross.18, exact.only = F, adjacent.only = F)
cross.18.nodup <- drop.markers(cross.18, unlist(dups))






rela <- cpgt
name <- paste(pop, "_kinship_mapped_markers.pdf", sep = "")
main <- paste(pop, "kinship of good markers (proportion of shared genotypes, 0-1)")

feet <- function(X,Y,Z){
  pdf(file.path('~/public_html/',Y),width = 26, height = 21)
  heatmap.2(X,symm=T,main=NA,notecol="black",labRow=NULL,na.rm=T,cexRow=2,labCol=NULL,key.title=NA,key.xlab=NA,key.ylab=NA,
  cellnote = round(rela,digits=2),notecex=1.2,srtCol=45,cexCol=2,dendrogram="row", margins = c(15,10),trace="none",keysize=.5,col=my_palette,breaks=col_breaks)
  title(Z,cex.main=3, line = -1)
  dev.off()
}
feet(rela[stricter_keep,stricter_keep], name, main)


png('~/public_html/ER_parent.png')
plot(sort(cpgt['BLI_BI1124M',]))
dev.off()


strict_toss <- rownames(cpgt)[!rownames(cpgt) %in% strict_keep]

rela[strict_keep,strict_keep]
rela[strict_keep,strict_toss]
rela[strict_toss,strict_toss]

head(gt.missing,100)

chr1 <- pull.geno(cross.18,chr=1)
which(chr1[,'1:1230855']==2)


stricter_keep <- names(which(rat < 1.2))










cross.18 <- subset(cross.18, ind=strict_keep)

png('~/public_html/ER_missing.png')
hist(gt.missing$missing, breaks=20)
dev.off()

cross.18 <- drop.markers(cross.18, rownames(gt.missing)[gt.missing$missing > 2])
gt.missing <- geno.table(cross.18)

png('~/public_html/pval.png')
hist(log10(gt.missing$P.value),breaks=100)
dev.off()

cross.18.bk <- cross.18

gt.cross.pars <- geno.table(cross.pars.BI)[markernames(cross.18), ]
swit <- rownames(gt.cross.pars[gt.cross.pars$AA == 1, ])

cross.pars.BI <- switchAlleles(cross.pars.BI, markers = swit)
cross.18 <- switchAlleles(cross.18, markers = swit)

gt.swit <- geno.table(cross.18)

### PARENT IS CONSIDERED BB

chr1 <- pull.geno(cross.18,chr=1)
which(chr1[,'1:1230855']==2)
###############
test
drop_I <- cross.18$pheno$ID[which(chr1[,'1:1230855']==2)]
cross.no <- subset(cross.18, ind=strict_keep[!strict_keep %in% drop_I])
gt.dp <- geno.table(cross.no,chr=1)

###############
cutoff <- 1e-01

cross.18 <- drop.markers(cross.18, rownames(gt.missing[gt.missing$P.value < cutoff,
  ]))

marker.warning()

gt.cross.pars <- geno.table(cross.pars.BI)[markernames(cross.18), ]
swit <- rownames(gt.cross.pars[gt.cross.pars$AA == 1, ])

cross.pars.BI <- switchAlleles(cross.pars.BI, markers = swit)
cross.18 <- switchAlleles(cross.18, markers = swit)

gt.cross.pars <- geno.table(cross.pars.BI)[markernames(cross.18), ]
gt.missing <- geno.table(cross.18)

cross.18.bk <- cross.18

##########################################

## cross.18.reorg1 <- subset(cross.18, chr=1)
## cross.18.reorg1 <- formLinkageGroups(cross.18.reorg1, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
## ##ELR_11592   ELR_10981   ELR_10980   ELR_10971   ELR_10990   ELR_10869
## sort(countXO(cross.18.reorg1))
##
## cross.18.reorg2 <- subset(cross.18, chr=2)
## cross.18.reorg2 <- formLinkageGroups(cross.18.reorg2, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
## sort(countXO(cross.18.reorg2))
## ## ELR_ER1124F   ELR_11592   ELR_11115   ELR_10967   ELR_11593   ELR_10869
##
## ## drop ELR_10869
## pullgts <- pull.geno(cross.18)
## rownames(pullgts) <- cross.18$pheno$ID
## gt.table <- apply(pullgts,1,table)
##
## homs <- gt.table[1,] + gt.table[3,]
## rat <- homs/gt.table[2,]
##
## strict_keep <- names(which(rat < 1.5))
## cross.18.strict <- subset(cross.18, ind=strict_keep)

inchr <- 1
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER1.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=2)
cross.18.strict <- switchAlleles(cross.18.strict, switch)

inchr <- 2
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER2.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=2)
cross.18.strict <- switchAlleles(cross.18.strict, switch)

inchr <- 3
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER3.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=2)
cross.18.strict <- switchAlleles(cross.18.strict, switch)

inchr <- 4
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER4.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=2)
cross.18.strict <- switchAlleles(cross.18.strict, switch)

inchr <- 5
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER5.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=2)
cross.18.strict <- switchAlleles(cross.18.strict, switch)

##6
inchr <- 6
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER6.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=2)
cross.18.strict <- switchAlleles(cross.18.strict, switch)

##7
inchr <- 7
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER7.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=2)
cross.18.strict <- switchAlleles(cross.18.strict, switch)
print(paste(inchr, 'switched'))

##8
inchr <- 8
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER8.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=c(2,3))
cross.18.strict <- switchAlleles(cross.18.strict, switch)
print(paste(inchr, 'switched'))


inchr <- 9
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER9.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=2)
cross.18.strict <- switchAlleles(cross.18.strict, switch)
print(paste(inchr, 'switched'))

inchr <- 10
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER10.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=c(2,3))
cross.18.strict <- switchAlleles(cross.18.strict, switch)
print(paste(inchr, 'switched'))

inchr <- 11
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER11.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=c(2,3))
cross.18.strict <- switchAlleles(cross.18.strict, switch)
print(paste(inchr, 'switched'))

inchr <- 12
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER12.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=c(2))
cross.18.strict <- switchAlleles(cross.18.strict, switch)
print(paste(inchr, 'switched'))


inchr <- 13
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER13.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=c(2))
cross.18.strict <- switchAlleles(cross.18.strict, switch)
print(paste(inchr, 'switched'))

inchr <- 14
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER14.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=c(2))
cross.18.strict <- switchAlleles(cross.18.strict, switch)
print(paste(inchr, 'switched'))

inchr <- 15
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER15.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=c(2))
cross.18.strict <- switchAlleles(cross.18.strict, switch)
print(paste(inchr, 'switched'))

inchr <- 16
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER16.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=c(2,3))
cross.18.strict <- switchAlleles(cross.18.strict, switch)
print(paste(inchr, 'switched'))

cross.18.strict.bk <- cross.18.strict

inchr <- 17
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER17.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=c(2,4))
cross.18.strict <- switchAlleles(cross.18.strict, switch)
print(paste(inchr, 'switched'))

inchr <- 18
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER18.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=c(2))
cross.18.strict <- switchAlleles(cross.18.strict, switch)
print(paste(inchr, 'switched'))

inchr <- 19
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER19.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=c(2))
cross.18.strict <- switchAlleles(cross.18.strict, switch)
print(paste(inchr, 'switched'))

inchr <- 20
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER20.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=c(2))
cross.18.strict <- switchAlleles(cross.18.strict, switch)
print(paste(inchr, 'switched'))

inchr <- 21
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER21.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=c(2,4))
cross.18.strict <- switchAlleles(cross.18.strict, switch)
print(paste(inchr, 'switched'))

inchr <- 22
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER22.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=c(2))
cross.18.strict <- switchAlleles(cross.18.strict, switch)
print(paste(inchr, 'switched'))

cross.18.strict.bk.2 <- cross.18.strict


inchr <- 23
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER23.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=c(2))
cross.18.strict <- switchAlleles(cross.18.strict, switch)
print(paste(inchr, 'switched'))

inchr <- 24
cross.18.reorg <- subset(cross.18.strict, chr=inchr)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)
png('~/public_html/RF_ER24.png')
plotRF(cross.18.reorg,chr=1:4)
dev.off()
print(inchr)
switch <- markernames(cross.18.reorg,chr=c(2))
cross.18.strict <- switchAlleles(cross.18.strict, switch)
print(paste(inchr, 'switched'))

cross.18.strict.bk.2 <- cross.18.strict

save.image('~/ER_fix.rsave')

cross.18.strict.drop <- cross.18.strict

for(i in 1:24){
cross.18.reorg <- subset(cross.18.strict, chr=i)
cross.18.reorg <- formLinkageGroups(cross.18.reorg, max.rf = 0.03, min.lod = 10, reorgMarkers = TRUE)
todrop <- markernames(cross.18.reorg,chr=c(2:nchr(cross.18.reorg)))
cross.18.strict.drop <<- drop.markers(cross.18.strict.drop,todrop)
}

save.image('~/ER_fix.rsave')

png('~/public_html/RF_ER18_F.png')
plotRF(cross.18.strict.drop,chr=c(1))
dev.off()


save.image('~/ER_fix.rsave')

cross.18.strict.drop.order <- cross.18.strict.drop

for(X in 1:24){
  cross.18.strict.drop.order <<- orderMarkers(cross.18.strict.drop.order, chr = X, window = 5, use.ripple = T, error.prob = 0.05,
    map.function = "kosambi", sex.sp = F, maxit = 3000, tol = 0.001)
save.image('~/ER_loopsave.rsave')
}


load('~/ER_fix.rsave')

crs.bk <- cross.18.strict.drop

for (i in 1:24){

ord <- order(as.numeric(gsub(paste(i, ":", sep = ""), "", markernames(crs.bk, chr = i))))

crs.bk <<- switch.order(crs.bk, chr = i, ord, error.prob = 0.05, map.function = "kosambi",
  maxit = 1000, tol = 0.001, sex.sp = F)

save.image('~/ER_loopsave_physical_position.rsave')

print(i)
}


load('~/ER_loopsave_physical_position.rsave')

for (i in 1:24){

png(paste0('~/public_html/ER_phys_',i,'_F.png'))
plotRF(crs.bk,chr=i)
dev.off()

print(i)
}



pullgts <- pull.geno(crs.bk)
yp <- as.numeric(gsub(".*:",'',colnames(pullgts)))
xp <- as.numeric(gsub(":.*",'',colnames(pullgts)))


png(paste0('~/public_html/ER_phys.png'))
plot(xp,yp)
dev.off()



##dups <- findDupMarkers(cross.18, exact.only = F, adjacent.only = F)
##cross.18 <- drop.markers(cross.18, unlist(dups))



pullgts <- pull.geno(crs.bk)
rownames(pullgts) <- crs.bk$pheno$ID
gt.table <- apply(pullgts,1,table)

homs <- gt.table[1,] + gt.table[3,]
rat <- homs/gt.table[2,]





crs.bk <- removeDoubleXO(crs.bk, verbose = T)
dups <- findDupMarkers(crs.bk, exact.only = F, adjacent.only = F)
crs.bk <- drop.markers(crs.bk, unlist(dups))

crs.bk$pheno$bin <- ifelse(crs.bk$pheno$Pheno > 2, 1 , 0)
cross.18.nodup$pheno$bin <- ifelse(cross.18.nodup$pheno$Pheno > 2, 1 , 0)
crs.bk$pheno$pheno_norm <- nqrank(crs.bk$pheno$Pheno)

scan.norm.em <- scanone(crs.bk, method = "em", model = "normal", maxit = 5000,
  pheno.col = 5)
em <- summary(scan.norm.em)
em[order(em$lod),]

scan.bin.em <- scanone(crs.bk, method = "em", model = "binary", pheno.col = 4,
  maxit = 5000)
bin.em <- summary(scan.bin.em )
bin.em[order(bin.em$lod),]

scan.np.em <- scanone(crs.bk, method = "em", model = "np", pheno.col = 4, maxit = 5000)
np.em <- summary(scan.np.em)
np.em[order(np.em$lod),]

scan.bin.mr <- scanone(crs.bk, method = "mr", model = "binary", pheno.col = 4)
mr <- summary(scan.bin.mr)
mr[order(mr$lod),]

### interactions?
> data(hyper)
hyper <- sim.geno(cross.18.nodup, step=1, n.draws=256, err=0.01)
qtl <- makeqtl(hyper, chr=14, pos=106.9, what="draws")

out.i <- addqtl(hyper, qtl=qtl, formula=y~Q1*Q2, method="imp")
out.a <- addqtl(hyper, qtl=qtl, formula=y~Q1+Q2, method="imp")

png(paste0('~/public_html/one_v_two_modle_qtl.png'),width=1000)
plot(out.i - out.a)
dev.off()

mname1 <- find.marker(cross.18.nodup, 1, 1)
mname2 <- find.marker(cross.18.nodup, 18, 152.427)
mname2 <- find.marker(cross.18.nodup, 14, 106.9) # marker D13Mit147

cross.18.nodup <- switchAlleles(cross.18.nodup, markers = mname2)

png(p
ste0('~/public_html/pxg_18_14_qtl.png'))
plotPXG(cross.18.nodup, c('18:18614940', mname1),pheno.col=4)
dev.off()

png(paste0('~/public_html/pxg_18_14_qtl.png'))
plotPXG(cross.18.nodup, '1:39331705',pheno.col=4)
dev.off()

14:29488677"

plotPXG(listeria, c(mname2, mname))


out.a <- addqtl(hyper, qtl=qtl, formula=y~Q1, method="imp")
png(paste0('~/public_html/one_modle_qtl.png'),width=1000)
plot(out.a)
dev.off()



gt.missing <- geno.table(crs.bk,chr=1)

perms.bin.em <- scanone(cross.18, method = "em", model = "binary", maxit = 2000,
  n.perm = 500, pheno.col = 4, n.cluster = 5)

perms.np.em <- scanone(cross.18, model = "np", pheno.col = 2, n.perm = 2000, method = "em",
  perm.strata = cross.18$pheno$stata, n.cluster = slurmcore)






th <- summary(perms.norm.imp)[1, ]
norm.qtl <- summary(scan.norm.imp, perms = perms.norm.imp, alpha = 0.05)
qtl.uns <- makeqtl(cross.18, chr = norm.qtl$chr, pos = norm.qtl$pos)
full <- stepwiseqtl(cross.18, additive.only = T, method = "imp", pheno.col = 6, scan.pairs = T)







cross.18 <- read.cross.jm(file = file.path(indpops, paste(pop, ".unphased.f2.csvr",
  sep = "")), format = "csvr", geno = c(1:3), estimate.map = FALSE)

### Pull names from plinkfile
path <- file.path(indpops, paste(pop, ".ped", sep = ""))
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
cross.18$pheno$ID <- paste(popname, indname, sep = "_")

indo <- cross.18$pheno$ID[!cross.18$pheno$ID %in% crs.bk$pheno$ID]
cross.18 <- subset(cross.18, ind = indo)

cross.18 <- subset(cross.18, ind = !is.na(cross.18$pheno$Phen))

gt.missing <- geno.table(cross.18)

cross.18 <- drop.markers(cross.18, rownames(gt.missing)[gt.missing$missing > 4])
cross.18 <- drop.markers(cross.18, rownames(gt.missing[gt.missing$P.value < 0.05,]))

cross.18$pheno$bin <- ifelse(cross.18$pheno$Pheno > 2, 1 , 0)

scan.bin.mr <- scanone(cross.18, method = "mr", model = "binary", pheno.col = 4)






crs.bk$pheno$pheno_norm <- nqrank(crs.bk$pheno$Pheno)




marker.warning()


gt.missing <- geno.table(cross.18)

#drops <- markernames(cross.18)[!markernames(cross.18) %in% marks.crs.bk]
#cross.18 <- drop.markers(cross.18,drops)

cross.18$pheno$bin <- ifelse(cross.18$pheno$Pheno > 2, 1 , 0)

cross.18 <- subset(cross.18, ind = !is.na(cross.18$pheno$Phen))

scan.bin.mr.18 <- scanone(cross.18, method = "mr", model = "binary", pheno.col = 4)

scan.np.em.18 <- scanone(cross.18, method = "em", model = "np", pheno.col = 3, maxit = 5000)











png('~/public_html/Geno_ER1.png')
plot(sort(rat))
dev.off()


gt.table.ratio <- apply



png('~/public_html/Geno_ER1.png')
plotGeno(cross.18.reorg1,chr=1)
dev.off()

png('~/public_html/RF_ER1.png')
plotRF(cross.18.reorg1,chr=1)
dev.off()

cross.18.reorg1 <- subset(cross.18.strict, chr=1)
cross.18.reorg1 <- formLinkageGroups(cross.18.reorg1, max.rf = 0.05, min.lod = 8, reorgMarkers = TRUE)



cross.18.reorg2.rf <- formLinkageGroups(cross.18.reorg2, max.rf = 0.05, reorgMarkers = TRUE)

png('~/public_html/RF_ER2.png')
plotRF(cross.18.reorg2.rf)
dev.off()

cross.18.reorg2.rf <- switchAlleles(cross.18.reorg2.rf, markernames(cross.18.reorg2.rf, chr=4 ))



png('~/public_html/Geno_ER2.png',width=1000,height=1000)
plotGeno(cross.18.reorg2,chr=1)
dev.off()

png('~/public_html/RF_ER2.png')
plotRF(cross.18.reorg1,chr=1:2)
dev.off()



gts <- geno.table(cross.18)

cross.18 <- cross.18.reorg

##cross.18 <- subset(cross.18.reorg, chr=1)
##cross.18 <- movemarker(cross.18, markernames(cross.18), newchr=2)

crs.bk <- cross.18.reorg

for (i in lenght(1:24)){

ord <- order(as.numeric(gsub(paste(i, ":", sep = ""), "", markernames(cross.18, chr = i))))

crs.bk <<- switch.order(cross.18, chr = i, ord, error.prob = 0.05, map.function = "kosambi",
  maxit = 1000, tol = 0.001, sex.sp = F)
}

png('~/public_html/Geno_ER.png',width=1200,height=1000)
plotGeno(cross.18,chr=1)
dev.off()


cross.18 <- removeDoubleXO(cross.18, verbose = T)
print("Done removing dxo..")

png('~/public_html/Geno_ER.png',width=4000,height=2000)
plotGeno(cross.18,chr=1)
dev.off()






png('~/public_html/pval.png')
hist(log10(gt.missing$P.value),breaks=50)
dev.off()


cross.18 <- removeDoubleXO(cross.18, verbose = T)
print("Done removing dxo..")








fileConn <- file(file.path(popdir, paste(X, "BB_markes.txt", sep = "_")))
writeLines(swit, fileConn)
close(fileConn)

print("estimating map with markers at physical positions")
ord <- order(as.numeric(gsub(paste(X, ":", sep = ""), "", markernames(cross.18, chr = X))))

cross.18 <- switch.order(cross.18, chr = X, ord, error.prob = 0.01, map.function = "kosambi",
  maxit = 1000, tol = 0.001, sex.sp = F)

swit <- checkAlleles(cross.18, threshold = 6, verbose = F)
cross.18 <- switchAlleles(cross.18, swit)

cross.18 <- formLinkageGroups(cross.18, max.rf = 0.25, min.lod = 12, reorgMarkers = TRUE)
cross.18 <- subset(cross.18, chr = which.max(nmar(cross.18)))
names(cross.18$geno) <- X

print("estimating map with markers at physical positions")
ord <- order(as.numeric(gsub(paste(X, ":", sep = ""), "", markernames(cross.18, chr = X))))

cross.18 <- switch.order(cross.18, chr = X, ord, error.prob = 0.01, map.function = "kosambi",
  maxit = 1000, tol = 0.001, sex.sp = F)

marker.warning()

cross.18 <- removeDoubleXO(cross.18, verbose = T)
print("Done removing dxo..")

cross.18 <- drop.missing(cross.18, miss)

marker.warning()

POS.map.18 <- est.map(cross.18, error.prob = 0.1, map.function = "kosambi", chr = X,
  maxit = 1000)

cross.18 <- replace.map(cross.18, POS.map.18)

print(summary(pull.map(cross.18))[as.character(X), ])

print("Writing the markers to rQTL format")
write.cross(cross.18, filestem = paste(popdir, "/chr", X, "_", outname, ".manymarks.QTLmap",
  sep = ""), format = "csv", chr = X)


#####
## scratch
####

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
