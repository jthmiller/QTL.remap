#!/bin/R
### Map QTLs 1 of 3
debug.cross <- F
pop <- 'ELR'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")

library('qtl')

mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

i <- commandArgs(TRUE)[commandArgs(TRUE) %in% c(1:24)]

################################################################################
################################################################################

cross <- read.cross.jm(file = file.path(indpops, paste0(pop, ".unphased.f2.csvr")),
format = "csvr", geno = c(1:3), estimate.map = FALSE)


fl <- file.path(mpath,'ELR_unmapped_filtered_added_markers.csv')

cross.all <- read.cross(
 file = fl,
 format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
 estimate.map = FALSE
)

################################################################################
################################################################################
################################################################################
i <- 1

cross <- subset(cross.all,chr=i)
nmars <- nmar(cross)
CHR1 <- colnames(pull.geno(cross))
CHR1[CHR1=="AHR2a_del"] <- 350000
ord <- order(as.numeric(gsub(".*:","",CHR1)))

cross <- switch.order(cross, chr = i, ord, error.prob = 0.01, map.function = "kosambi",
maxit = 10, tol = 0.001, sex.sp = F)

##cross <- subset(cross,ind=!is.na(cross$pheno$Pheno))

### CROSS FOR LOCAT. BAD MARKERS #####################################
#### RND 1  #####################################
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
cross <- calc.errorlod(cross, err=0.10)

png(paste0('~/public_html/ELR_gts_preclean_droppedmark',i,'.png'),height=2500,width=4500)
plotGeno(cross)
dev.off()
##cross.bk <- cross
cross <- cleanGeno_jm(cross, chr=i, maxdist=100, maxmark=8, verbose=TRUE)
cross <- calc.errorlod(cross, err=0.05)
cross <- removeDoubleXO(cross)
cross <- calc.errorlod(cross, err=0.025)
cross <- cleanGeno_jm_2(cross, chr=i, maxdist=50, maxmark=4, verbose=TRUE)
cross <- calc.errorlod(cross, err=0.01)

png(paste0('~/public_html/ELR_clean.png'),height=2500,width=4000)
plotGeno(cross,cex=3)
dev.off()

png(paste0('~/public_html/ELR_RF_clean',i,'.png'))
plotRF(cross)
dev.off()

fl <- file.path(mpath,paste0(i,'ELR_unmapped_unfiltered'))
write.cross(cross,filestem=fl,format="csv")


cross.bk <- cross

z <- pull.geno(cross)


################################################################################
### THIN MARKERS IF NEEDED #####################################################

mp <- as.numeric(gsub(".*:",'',markernames(cross)))
names(mp) <- markernames(cross)
mp <- list('1'=mp)
cross <- replace.map(cross,mp)

gts <- geno.table(cross)
weight <- 1 - gts$missing/rowSums(gts[,c(3:5)])*10

dwnsmpl <- pickMarkerSubset(pull.map(cross)[[1]],2000, weights=weight)

drops <- markernames(cross)[! markernames(cross) %in% dwnsmpl]
cross.dwn <- drop.markers(cross,drops)

cross.dwn <- calc.genoprob(cross.dwn)
cross.dwn <- sim.geno(cross.dwn)
cross.dwn <- calc.errorlod(cross.dwn, err=0.01)

png(paste0('~/public_html/ELR_gts_CHR1_downsmpl.png'),height=1500,width=4500)
plotGeno(cross.dwn ,cex=3)
dev.off()

#####MAP ########################################################################

cross.dwn <- subset(cross.dwn,ind=!cross.dwn$pheno$ID %in% c('ELR_10869','ELR_ER1124F','BLI_BI1124M'))

cross.dwn <- orderMarkers(cross.dwn, window=7,verbose=FALSE,chr=i,
                 use.ripple=TRUE, error.prob=0.025, sex.sp=FALSE,
                 map.function="kosambi",maxit=500, tol=1e-4)


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

filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_gts_CHR1_downsmpl_map')
write.cross(cross.dwn,chr=i,filestem=filename,format="csv")





try.imp <- fill.geno(cross.dwn,method=c("imp"),error.prob=0.01)
argmax <- fill.geno(cross.dwn,method=c("argmax"),error.prob=0.01)
png(paste0('~/public_html/ELR_gts_CHR1_downsmpl.png'),height=1500,width=4500)
plotGeno(argmax  ,cex=3)
dev.off()









####MANUAL FIX #####
pgt_sort <- lapply(as.character(cross$pheno$ID[1:88]),function(id){
 indv <- cross$pheno$ID==id
 gts <- z[indv,1:50]
 phs <- names(which.max(table(gts)))
 ngt <- max(table(gts))/sum(table(gts))
 xo <- countXO(cross)[id]
 cbind(phs,xo,ngt)
})
pgt_sort <- data.frame(do.call(rbind,pgt_sort),stringsAsFactors=F)
pgtord <- rownames(pgt_sort)[order(as.numeric(pgt_sort$xo))]


png(paste0('~/public_html/ELR_gts_CHR1_prefix.png'),height=1500,width=4500)
plotGeno(cross,ind=pgtord ,cex=2)
dev.off()

png(paste0('~/public_html/ELR_RF_a',i,'.png'),width=2000)
plotRF(cross)
dev.off()

msiord <- sort(nmissing(cross))

png(paste0('~/public_html/ELR_gts_mis.png'),height=1500,width=4500)
plotGeno(cross,ind=names(msiord),cex=2)
dev.off()

################################################################################
lots_missing <- names(sort(nmissing(cross))[sort(nmissing(cross)) > 25])
lots_xo <- names(sort(countXO(cross))[sort(countXO(cross)) > 100])
both <- intersect(lots_missing, lots_xo)
################################################################################
################################################################################
"ELR_10991", set all 3 to 2 y(all 3 to 2)
"ELR_10871", set all 3 to 2 y(NR)
"ELR_10989" set all to 2 y(NR)
"ELR_10974" 1 and 3 markers after 16000 should be all 2 y(all 1 to 2)
"ELR_10953" small region of transition from 1 to 2 to 3. change the 1/2 mix to 2
"ELR_10998" before 10000, change all to hets (or drop) y(all 1 to 2)
"ELR_10882" just before 10000, is 3, rest should be 2  y(all 1 to 2)
"ELR_10969" after stretch of 3, change to 2 y(all 1 to 2)
"ELR_10924" change all 3 to 2 y(all 3 to 2)
"ELR_10981" change all to 2 y(NR)
"ELR_10990" change all to 2 y(NR)
"ELR_10980" change all to 2 y(NR)
"ELR_10967" change all to 2 until reigion of 3 starts to end y(all 1 to 2)
"ELR_10869" really bad samp n (?)
"ELR_10971" change all to 2 y
"ELR_11593" after stretch of 3, change to 2 for the rest (all 1 to 2)
"ELR_11592" change all except stretch of 1 (all 3 to 2)
"ELR_10878" y(all 3 to 2)
"ELR_10984" y(all 1 to 2)
not sure: "ELR_10869"

all_1_2 <- c("ELR_10871","ELR_10989","ELR_10974","ELR_10998","ELR_10882","ELR_10969","ELR_10981","ELR_10990","ELR_10980","ELR_10967","ELR_10971","ELR_11593","ELR_10984","ELR_11587")
for(id in all_1_2){
 indv <- cross$pheno$ID==id
 tozero <- which(cross$geno[[1]]$data[indv,]==1)
 cross$geno[[1]]$data[indv,tozero] <- 2
}
all_3_2 <- c("ELR_10991","ELR_10871","ELR_10989","ELR_10924","ELR_10981","ELR_10990","ELR_10980","ELR_10971","ELR_11592","ELR_10878")
for(id in all_3_2 ){
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

### 1 to 2
id <- "ELR_10953"
indv <- cross$pheno$ID==id
pull.geno(cross)[indv,]
start <- which(names(pull.geno(cross)[indv,])=='1:14891665')
end <- which(names(pull.geno(cross)[indv,])=='1:20393972')
chng <- names(pull.geno(cross)[indv,c(start:end)]==1)
cross$geno[[1]]$data[indv,c(start:end)] <- 2

id <- "ELR_11592"
indv <- cross$pheno$ID==id
pull.geno(cross)[indv,]
start <- which(names(pull.geno(cross)[indv,])=='1:31917407')
end <- which(names(pull.geno(cross)[indv,])=='1:39331864')
chng <- names(pull.geno(cross)[indv,c(start:end)]==1)
cross$geno[[1]]$data[indv,c(start:end)] <- 2

id <- "ELR_11593"
indv <- cross$pheno$ID==id
pull.geno(cross)[indv,]
start <- which(names(pull.geno(cross)[indv,])=='1:12634024')
end <- which(names(pull.geno(cross)[indv,])=='1:39331864')
chng <- names(pull.geno(cross)[indv,c(start:end)]==3)
cross$geno[[1]]$data[indv,c(start:end)] <- 2


id <- "ELR_10924"
indv <- cross$pheno$ID==id
pull.geno(cross)[indv,]
start <- which(names(pull.geno(cross)[indv,])=='1:210602')
end <- which(names(pull.geno(cross)[indv,])=='1:34613326')
chng <- names(pull.geno(cross)[indv,c(start:end)]==1)
cross$geno[[1]]$data[indv,c(start:end)] <- 2

id <- "ELR_10969"
indv <- cross$pheno$ID==id
pull.geno(cross)[indv,]
start <- which(names(pull.geno(cross)[indv,])=='1:20675753')
end <- which(names(pull.geno(cross)[indv,])=='1:39331864')
chng <- names(pull.geno(cross)[indv,c(start:end)]==3)
cross$geno[[1]]$data[indv,c(start:end)] <- 2

id <- "ELR_10882"
indv <- cross$pheno$ID==id
pull.geno(cross)[indv,]
start <- which(names(pull.geno(cross)[indv,])=='1:22458748')
end <- which(names(pull.geno(cross)[indv,])=='1:39331864')
chng <- names(pull.geno(cross)[indv,c(start:end)]==3)
cross$geno[[1]]$data[indv,c(start:end)] <- 2

id <- "ELR_11587"
indv <- cross$pheno$ID==id
pull.geno(cross)[indv,]
start <- which(names(pull.geno(cross)[indv,])=='1:210602')
end <- which(names(pull.geno(cross)[indv,])=='1:3846643')
chng <- names(pull.geno(cross)[indv,c(start:end)]==3)
cross$geno[[1]]$data[indv,c(start:end)] <- 2

id <- "ELR_11115"
indv <- cross$pheno$ID==id
pull.geno(cross)[indv,]
start <- which(names(pull.geno(cross)[indv,])=='1:210602')
end <- which(names(pull.geno(cross)[indv,])=='1:14268870')
chng <- names(pull.geno(cross)[indv,c(start:end)]!=2)
cross$geno[[1]]$data[indv,c(start:end)] <- 2

id <- "ELR_10967"
indv <- cross$pheno$ID==id
pull.geno(cross)[indv,]
start <- which(names(pull.geno(cross)[indv,])=='1:239958')
end <- which(names(pull.geno(cross)[indv,])=='1:12605060')
chng <- names(pull.geno(cross)[indv,c(start:end)]==3)
cross$geno[[1]]$data[indv,c(start:end)] <- 2

cross <- subset(cross,ind=!cross$pheno$ID %in% c('ELR_10869','ELR_ER1124F'))

######## PLOT #####################
cross <- calc.genoprob(cross)
png(paste0('~/public_html/ELR_gts_fix_these_genos.png'),height=1500,width=4500)
plotGeno(cross,ind=names(sort(countXO(cross))),cex=2)
abline(v=2.862)
dev.off()
########################

png(paste0('~/public_html/ELR_gts_a',i,'.png'),height=1500,width=4500)
plotGeno(cross,ind=pgtord ,cex=2)
dev.off()

################################################################################
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'
fl <- file.path(mpath,'ELR_CHR1_HAND_IMPUTED_HAPLOTYPES')
write.cross(cross,filestem=fl,format="csv")

################################################################################
### THIN MARKERS IF NEEDED #####################################################

mp <- as.numeric(gsub(".*:",'',markernames(cross)))
names(mp) <- markernames(cross)
mp <- list('1'=mp)
cross <- replace.map(cross,mp)

gts <- geno.table(cross)
weight <- 1 - gts$missing/rowSums(gts[,c(3:5)])*10

dwnsmpl <- pickMarkerSubset(pull.map(cross)[[1]],2000, weights=weight)

drops <- markernames(cross)[! markernames(cross) %in% dwnsmpl]
cross.dwn <- drop.markers(cross,drops)

cross.dwn <- calc.genoprob(cross.dwn)
cross.dwn <- sim.geno(cross.dwn)
cross.dwn <- calc.errorlod(cross.dwn, err=0.01)

png(paste0('~/public_html/ELR_gts_CHR1_downsmpl.png'),height=1500,width=4500)
plotGeno(no_dbl ,cex=3)
dev.off()

'ELR_10977,ELR_10988'

cross.dwn <- subset(cross.dwn,ind=!cross$pheno$ID %in% c('ELR_10869','ELR_ER1124F','ELR_10977','ELR_10988','BLI_BI1124M'))

try.imp <- fill.geno(cross.dwn,method=c("imp"),error.prob=0.05)
no_dbl <- fill.geno(cross.dwn,method=c("no_dbl_XO"),error.prob=0.05)
argmax <- fill.geno(cross.dwn,method=c("argmax"),error.prob=0.05)
maxm <- fill.geno(cross.dwn,method=c("maxmarginal"),error.prob=0.05)

no_dbl <- removeDoubleXO(no_dbl)

no_dbl <- orderMarkers(no_dbl, window=7,verbose=FALSE,chr=i,
                 use.ripple=TRUE, error.prob=0.025, sex.sp=FALSE,
                 map.function="kosambi",maxit=500, tol=1e-4)



no_dbl <- read.cross(
 file = filename,
 format = "csv", genotypes=c("AA","AB","BB"), alleles=c("A","B"),
 estimate.map = FALSE
)

no_dbl_map <-  est.map(no_dbl,  error.prob=0.025,
            map.function="kosambi",
            maxit=1000, tol=1e-6, sex.sp=TRUE,
            verbose=FALSE, omit.noninformative=TRUE, n.cluster=6)

 no_dbl_map <- shiftmap(no_dbl_map, offset=0)

no_dbl <- replace.map(cross, no_dbl_map)

filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_gts_CHR1_downsmpl_map')
write.cross(no_dbl,chr=i,filestem=filename,format="csv")


cleanGeno(no_dbl,  maxmark = 1)


cleanGeno(cross, chr, maxdist=2.5, maxmark=2, verbose=TRUE)
################################################################################
################################################################################
loglik <- err <- c(0.001, 0.005, 0.01, 0.015, 0.02)
for(i in seq(along=err)) {
 cat(i, "of", length(err), "\n")
 tempmap <- est.map(mapthis, error.prob=err[i])
 loglik[i] <- sum(sapply(tempmap, attr, "loglik"))
}
lod <- (loglik - max(loglik))/log(10)
