#!/bin/R
### Map QTLs 1 of 3
debug.cross <- T
pop <- 'ELR'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")

library('qtl')

mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]

#cross <- read.cross.jm(file = file.path(indpops, paste0(pop, ".unphased.f2.csvr")),
#format = "csvr", geno = c(1:3), estimate.map = FALSE)


fl <- file.path(mpath,'ELR_unmapped_filtered_added_markers.csv')

cross.all <- read.cross(
 file = fl,
 format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
 estimate.map = FALSE
)

################################################################################
################################################################################
################################################################################

for (i in 3:24){

 cross <- subset(cross.all,chr=i)
 nmars <- nmar(cross)
 ## initial order
 ord <- order(as.numeric(gsub(".*:","",names(pull.map(cross)[[1]]))))
 cross <- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
  maxit = 10, tol = 0.001, sex.sp = F)
 ################################################################################
 ################################################################################

 cross <- subset(cross,ind=!is.na(cross$pheno$Pheno))

 cross <- calc.errorlod(cross, err=0.01)
 png(paste0('~/public_html/ELR_gts_preclean',i,'.png'),height=2500,width=4500)
 plotGeno(cross)
 dev.off()

 png(paste0('~/public_html/ELR_xo_a',i,'.png'))
 hist(sort(table(unlist(locateXO(cross)))),breaks=30)
 dev.off()

 loc.xocount <- table(unlist(locateXO(cross)))

 marker <- sapply(as.numeric(names(loc.xocount)),function(X){
  find.marker(cross, chr=i, pos=X) })

 dropdf <- data.frame(loc.xocount,marker,stringsAsFactors=F)

 dropdf$tot <- sapply(dropdf$mark, function(X){ sum(table(pull.geno(cross,i)[,X]))})
 drops <- unique(dropdf[dropdf$Freq/dropdf$tot > 0.10,'marker'])

 cross <- drop.markers(cross,drops)
 cross <- calc.genoprob(cross)
 cross <- sim.geno(cross)
 cross <- calc.errorlod(cross, err=0.01)

 png(paste0('~/public_html/ELR_gts_preclean_droppedmark',i,'.png'),height=2500,width=4500)
 plotGeno(cross)
 dev.off()

 cross <- cleanGeno_jm(cross, chr=i, maxdist=100, maxmark=8, verbose=TRUE)
 cross <- calc.errorlod(cross, err=0.025)
 cross <- removeDoubleXO(cross)
 cross <- calc.errorlod(cross, err=0.025)
 cross <- cleanGeno_jm_2(cross, chr=i, maxdist=50, maxmark=4, verbose=TRUE)
 cross <- calc.errorlod(cross, err=0.025)

 png(paste0('~/public_html/ELR_clean.png'),height=2500,width=4000)
 plotGeno(cross,cex=3)
 dev.off()

 png(paste0('~/public_html/ELR_RF_clean',i,'.png'))
 plotRF(cross)
 dev.off()

 fl <- file.path(mpath,paste0(i,'ELR_unmapped_unfiltered'))
 write.cross(cross,filestem=fl,format="csv")


################################################################################
### THIN MARKERS IF NEEDED #####################################################

mp <- as.numeric(gsub(".*:",'',markernames(cross)))
names(mp) <- markernames(cross)
mp <- list(mp)
names(mp) <- i
cross <- replace.map(cross,mp)

gts <- geno.table(cross)
weight <- 1 - gts$missing/rowSums(gts[,c(3:5)])*10

dwnsmpl <- pickMarkerSubset(pull.map(cross)[[1]],2000, weights=weight)

drops <- markernames(cross)[! markernames(cross) %in% dwnsmpl]
cross.dwn <- drop.markers(cross,drops)

cross.dwn <- calc.genoprob(cross.dwn)
cross.dwn <- sim.geno(cross.dwn)
cross.dwn <- calc.errorlod(cross.dwn, err=0.01)

png(paste0('~/public_html/ELR_gts_CHR',i,'_downsmpl.png'),height=1500,width=4500)
plotGeno(cross.dwn ,cex=3)
dev.off()

#####MAP ########################################################################

cross.dwn <- subset(cross.dwn,ind=!cross$pheno$ID %in% c('ELR_10869','ELR_ER1124F','ELR_10977','ELR_10988','BLI_BI1124M'))

cross.dwn <- orderMarkers(cross.dwn, window=7,verbose=FALSE,chr=i,
                 use.ripple=TRUE, error.prob=0.025, sex.sp=FALSE,
                 map.function="kosambi",maxit=500, tol=1e-4)

cross.dwn <- calc.genoprob(cross.dwn)
cross.dwn <- sim.geno(cross.dwn)
cross.dwn <- calc.errorlod(cross.dwn, err=0.01)

##cross.dwn <- read.cross(
## file = filename,
## format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
## estimate.map = FALSE
##)
##
cross.dwn_map <-  est.map(cross.dwn,  error.prob=0.025,
            map.function="kosambi",
            maxit=10000, tol=1e-6, sex.sp=FALSE,
            verbose=FALSE, omit.noninformative=TRUE, n.cluster=6)

 cross.dwn_map <- shiftmap(cross.dwn_map, offset=0)

cross.dwn <- replace.map(cross, cross.dwn_map)

filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_gts_CHR',i,'_downsmpl_map')
write.cross(cross.dwn,chr=i,filestem=filename,format="csv")

}
################################################################################









################################################################################
################################################################################
### THIN MARKERS IF NEEDED #####################################################

mp <- as.numeric(gsub(".*:",'',markernames(cross)))
names(mp) <- markernames(cross)
mp <- list(get(i)=mp)
cross <- replace.map(cross,mp)

gts <- geno.table(cross)
weight <- 1 - gts$missing/rowSums(gts[,c(3:5)])*10

dwnsmpl <- pickMarkerSubset(pull.map(cross)[[1]],2000, weights=weight)

drops <- markernames(cross)[! markernames(cross) %in% dwnsmpl]
cross <- drop.markers(cross,drops)

cross <- calc.genoprob(cross)
cross <- sim.geno(cross)
cross <- removeDoubleXO(cross)
cross <- calc.errorlod(cross, err=0.01)

png(paste0('~/public_html/ELR_gts_c',i,'.png'),width=1000)
plotGeno(cross)
dev.off()

################################################################################

mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'
write.table(markernames(cross.ss),file.path(mpath,'ER_markers_subst.table'))

fl <- file.path(mpath,'ELR_subsetted')
write.cross(cross.ss,filestem=fl,format="csv")

################################################################################
loglik <- err <- c(0.001, 0.005, 0.01, 0.015, 0.02)
for(i in seq(along=err)) {
 cat(i, "of", length(err), "\n")
 tempmap <- est.map(mapthis, error.prob=err[i])
 loglik[i] <- sum(sapply(tempmap, attr, "loglik"))
}
lod <- (loglik - max(loglik))/log(10)

#### Further improve chr1

cross <- subset(cross.all,chr=1)

lots_missing <- names(sort(nmissing(cross))[sort(nmissing(cross)) > 25])
lots_xo <- names(sort(countXO(cross))[sort(countXO(cross)) > 100])
both <- intersect(lots_missing, lots_xo)

table(pull.geno(cross)[which(cross$pheno$ID=="ELR_10991"),])
table(pull.geno(cross)[which(cross$pheno$ID=="ELR_10974"),])


"ELR_10991", set all 3 to 2
"ELR_10871", set all 3 to 2
"ELR_10989" set all to 2
"ELR_10974" 1 and 3 markers after 16000 should be all 2
"ELR_10953" small region of transition from 1 to 2 to 3. change the 1/2 mix to 2
"ELR_10998" before 10000, change all to hets (or drop)
"ELR_10882" just before 10000, is 3, rest should be 2
"ELR_10969" after stretch of 3, change to 2
"ELR_10924" change all 3 to 2
"ELR_10981" change all to 2
"ELR_10990" change all to 2
"ELR_10980" change all to 2
"ELR_10967" change all to 2 until reigion of 3 starts to end
"ELR_10869" change all to 2 (really bad samp)
"ELR_10971" change all to 2 (really bad samp)
"ELR_11593" after stretch of 3, change to 2 for the rest
"ELR_11592" change all except stretch of 1


all_to_2 <- c("ELR_10991","ELR_10871","ELR_10989","ELR_10981","ELR_10990","ELR_10980","ELR_10967","ELR_10869","ELR_10971")
### ALL HET
for(id in all_to_2){
 indv <- cross$pheno$ID==id
 tozero <- which(cross$geno[[1]]$data[indv,]!=2)
 cross$geno[[1]]$data[indv,tozero] <- 2
}

### 1 to 2
id <- "ELR_10974"
indv <- cross$pheno$ID==id
###pull.geno(cross)[indv,]
start <- which(names(pull.geno(cross)[indv,])=='1:31892692')
end <- length(names(pull.geno(cross)[indv,]))
cross$geno[[1]]$data[indv,c(start,end)] <- 2

### 1 to 2
id <- "ELR_10953"
indv <- cross$pheno$ID==id
pull.geno(cross)[indv,]
start <- which(names(pull.geno(cross)[indv,])=='1:14891665')
end <- which(names(pull.geno(cross)[indv,])=='1:20393972')
chng <- names(pull.geno(cross)[indv,c(start:end)]==1)
cross$geno[[1]]$data[indv,c(start,end)] <- 2


all_1_2 <- c("ELR_10974","ELR_10998","ELR_10882","ELR_10869","ELR_11593","ELR_11593","ELR_10969")
for(id in all_1_2){
 indv <- cross$pheno$ID==id
 tozero <- which(cross$geno[[1]]$data[indv,]==1)
 cross$geno[[1]]$data[indv,tozero] <- 2
}

all_3_2 <- c("ELR_10998","ELR_10882","ELR_10869","ELR_11593","ELR_10924","ELR_11592")
for(id in all_1_2){
 indv <- cross$pheno$ID==id
 tozero <- which(cross$geno[[1]]$data[indv,]==3)
 cross$geno[[1]]$data[indv,tozero] <- 2
}

### 3 to 2
id <- "ELR_10974"
indv <- cross$pheno$ID==id
###pull.geno(cross)[indv,]
start <- which(names(pull.geno(cross)[indv,])=='1:31564972')
end <- length(names(pull.geno(cross)[indv,]))
cross$geno[[1]]$data[indv,c(start:end)] <- 2


######## PLOT #####################
cross <- calc.genoprob(cross)
png(paste0('~/public_html/ELR_gts_fix_these_genos',indo,'.png'),height=1000,width=4000)
plotGeno(subset(cross,ind=both),cex=3)
abline(v=2.862)
dev.off()
########################










indo <- "ELR_10871"

table(pull.geno(cross)[which(cross$pheno$ID==indo),])

png(paste0('~/public_html/ELR_gts_c',indo,'.png'),width=3000)
plotGeno(subset(cross,ind=indo))
dev.off()


fake.f2 <- argmax.geno(cross, step=2, off.end=5, err=0.05)

png(paste0('~/public_html/ELR_gts_c',indo,'.png'),width=3000)
plotGeno(subset(cross,ind=both))
dev.off()


fake.f2 <- reduce

png(paste0('~/public_html/ELR_gts_fix_these_genos',indo,'.png'),height=1000,width=4000)
plotGeno(subset(cross,ind=both),cex=3)
dev.off()


fake.f2 <- fill.geno(fake.f2, method=c("argmax"))

par <- subset(cross.all,chr=1,ind="BLI_BI1124M")

png(paste0('~/public_html/ELR_gts_c',indo,'.png'),height=1000,width=4000)
plotGeno(cross.all,chr=1,cex=3)
dev.off()
