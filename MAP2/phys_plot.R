


library('qtl')
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'




new <- file.path(mpath,'new.mapped.tsp.csv')
nbh <- file.path(mpath,'nbh.mapped.tsp.csv')
elr <- file.path(mpath,'elr.mapped.tsp.csv')
brp <- file.path(mpath,'brp.mapped.tsp.csv')

new <- read.cross(file = new, format = "csv", genotypes=c("1","2","3"),estimate.map = FALSE)
nbh <- read.cross(file = nbh, format = "csv", genotypes=c("1","2","3"),estimate.map = FALSE)
elr <- read.cross(file = elr, format = "csv", genotypes=c("1","2","3"),estimate.map = FALSE)
brp <- read.cross(file = brp, format = "csv", genotypes=c("1","2","3"),estimate.map = FALSE)

gts.new <- geno.table(new)
gts.nbh <- geno.table(nbh)
gts.elr <- geno.table(elr)
gts.brp <- geno.table(brp)

nbh.pos <- gsub(".*:",'',rownames(gts.nbh))
elr.pos <- gsub(".*:",'',rownames(gts.elr))
brp.pos <- gsub(".*:",'',rownames(gts.brp))
new.pos <- gsub(".*:",'',rownames(gts.new))

nbh.chr <- gsub(":.*",'',rownames(gts.nbh))
elr.chr <- gsub(":.*",'',rownames(gts.elr))
brp.chr <- gsub(":.*",'',rownames(gts.brp))
new.chr <- gsub(":.*",'',rownames(gts.new))

####### Resistant

nbh.res <- subset(nbh,ind=nbh$pheno$bin==0)
gts.nbh.res <- geno.table(nbh.res)
nbh.pos.res <- gsub(".*:",'',rownames(gts.nbh.res))
nbh.chr.res <- gsub(":.*",'',rownames(gts.nbh.res))
#######
nbh.sen <- subset(nbh,ind=nbh$pheno$bin==1)
gts.nbh.sen <- geno.table(nbh.sen)
nbh.pos.sen <- gsub(".*:",'',rownames(gts.nbh.sen))
nbh.chr.sen <- gsub(":.*",'',rownames(gts.nbh.sen))



elr.res <- subset(elr,ind=elr$pheno$bin==0)
gts.elr.res <- geno.table(elr.res)
elr.pos.res <- gsub(".*:",'',rownames(gts.elr.res))
elr.chr.res <- gsub(":.*",'',rownames(gts.elr.res))
#######
elr.sen <- subset(elr,ind=elr$pheno$bin==1)
gts.elr.sen <- geno.table(elr.sen)
elr.pos.sen <- gsub(".*:",'',rownames(gts.elr.sen))
elr.chr.sen <- gsub(":.*",'',rownames(gts.elr.sen))



new.sen <- subset(new,ind=new$pheno$bin==1)
gts.new.sen <- geno.table(new.sen)
new.pos.sen <- gsub(".*:",'',rownames(gts.new.sen))
new.chr.sen <- gsub(":.*",'',rownames(gts.new.sen))

new.res <- subset(new,ind=new$pheno$bin==0)
gts.new.res <- geno.table(new.res)
new.pos.res <- gsub(".*:",'',rownames(gts.new.res))
new.chr.res <- gsub(":.*",'',rownames(gts.new.res))





png(paste0('~/public_html/NBH_gt_1.png'),height=2500,width=4500)
 plot(1:length(nbh.pos[ind]), nbh.pos[ind], cex=2)
 points(1:length(nbh.pos[ind]), as.numeric(gts.nbh$AA[ind]), col='blue')
 points(1:length(nbh.pos[ind]), as.numeric(gts.nbh$AA[ind]),col='red')
dev.off()







png(paste0('~/public_html/NBH_gt_seg1.png'),height=900,width=900)
par(mfrow=c(5,5))
for (i in 1:24){
 ind <- which(nbh.chr==i)
 plot(1:length(nbh.pos[ind]), as.numeric(gts.nbh$AA[ind]), col='blue',pch=16, ylim=c(0,60))
 points(1:length(nbh.pos[ind]), as.numeric(gts.nbh$BB[ind]),col='red',pch=16)
 points(1:length(nbh.pos[ind]), as.numeric(gts.nbh$AB[ind]),col='green',pch=16)
 text(33,45,i)
}
dev.off()


png(paste0('~/public_html/NEW_gt_seg1.png'),height=900,width=900)
par(mfrow=c(5,5))
for (i in 1:24){
 ind <- which(new.chr==i)
 plot(1:length(new.pos[ind]), as.numeric(gts.new$AA[ind]), col='blue',pch=16, ylim=c(0,60))
 points(1:length(new.pos[ind]), as.numeric(gts.new$BB[ind]),col='red',pch=16)
 points(1:length(new.pos[ind]), as.numeric(gts.new$AB[ind]),col='green',pch=16)
 text(33,45,i)
}
dev.off()




png(paste0('~/public_html/ELR_gt_seg1.png'),height=900,width=900)
par(mfrow=c(5,5))
for (i in 1:24){
 ind <- which(elr.chr==i)
 plot(1:length(elr.pos[ind]), as.numeric(gts.elr$AA[ind]), col='blue',pch=16, ylim=c(0,60))
 points(1:length(elr.pos[ind]), as.numeric(gts.elr$BB[ind]),col='red',pch=16)
 points(1:length(elr.pos[ind]), as.numeric(gts.elr$AB[ind]),col='green',pch=16)
 text(33,45,i)
}
dev.off()



png(paste0('~/public_html/BRP_gt_seg1.png'),height=900,width=900)
par(mfrow=c(5,5))
for (i in 1:24){
 ind <- which(brp.chr==i)
 plot(1:length(brp.pos[ind]), as.numeric(gts.brp$AA[ind]), col='blue',pch=16, ylim=c(0,60))
 points(1:length(brp.pos[ind]), as.numeric(gts.brp$BB[ind]),col='red',pch=16)
 points(1:length(brp.pos[ind]), as.numeric(gts.brp$AB[ind]),col='green',pch=16)
 text(33,45,i)
}
dev.off()


#######phys pos

png(paste0('~/public_html/NBH_gt_seg_phys_pos.png'),height=900,width=900)
par(mfrow=c(5,5))
for (i in 1:24){
 ind <- which(nbh.chr==i)

 plot(nbh.pos.res[ind], as.numeric(gts.nbh.res$AA[ind]), type="n", pch=16, ylim=c(0,100))

 abline(h=23,col='black')
 abline(h=24+50,col='black')
 abline(h=11.5,col='black')
 abline(h=12+50,col='black')

 points(nbh.pos.res[ind], as.numeric(gts.nbh.res$AA[ind]), col='blue',pch=16)
 points(nbh.pos.sen[ind], as.numeric(gts.nbh.sen$AA[ind])+50, col='blue',pch=16)
 points(nbh.pos.res[ind], as.numeric(gts.nbh.res$BB[ind]),col='red',pch=16)
 points(nbh.pos.sen[ind], as.numeric(gts.nbh.sen$BB[ind])+50,col='red',pch=16)
 points(nbh.pos.res[ind], as.numeric(gts.nbh.res$AB[ind]),col='green',pch=16)
 points(nbh.pos.sen[ind], as.numeric(gts.nbh.sen$AB[ind])+50,col='green',pch=16)

 text(150,40,i)
}
dev.off()


###############################################################
png(paste0('~/public_html/NBH_gt_seg_pval_pos.png'),height=900,width=900)
par(mfrow=c(5,5))
for (i in 1:24){
 ind <- which(nbh.chr==i)

 plot(nbh.pos.res[ind], -log(as.numeric(gts.nbh.res$P.value[ind])), col='blue', pch=16, ylim=c(0,20))
 points(nbh.pos.sen[ind], -log(as.numeric(gts.nbh.sen$P.value[ind])), col='red', pch=16)
 points(nbh.pos[ind], -log(as.numeric(gts.nbh$P.value[ind])), col='green', pch=16)

 text(2550,15,i)
}
dev.off()




scan.bin.mr <- scanone(nbh, method = "mr", model = "binary", pheno.col = 4)



head(scan.bin.mr$lod)





png(paste0('~/public_html/NBH_lodxSeg.png'),height=900,width=900)
par(mfrow=c(5,5))
for (i in 1:24){
 ind <- which(nbh.chr==i)

 plot( -log(as.numeric(gts.nbh.res$P.value[ind])), scan.bin.mr$lod[ind], col='blue', pch=16, ylim=c(0,20),xlim=c(0,10))
 points( -log(as.numeric(gts.nbh.sen$P.value[ind])),scan.bin.mr$lod[ind], col='red', pch=16)
 ##points(-log(as.numeric(gts.nbh$P.value[ind])),scan.bin.mr$lod[ind],  col='green', pch=16)

 text(2550,15,i)
}
dev.off()



scan.nbh <- scanone(nbh, method = "mr", model = "binary", pheno.col = 4)
scan.elr <- scanone(elr, method = "mr", model = "binary", pheno.col = 4)
scan.new <- scanone(new, method = "mr", model = "binary", pheno.col = 4)


png(paste0('~/public_html/allpop_lodxSeg.png'),height=900,width=900)
par(mfrow=c(5,5))
for (i in 1:24){
 ind.elr <- which(elr.chr==i)
 ind.new <- which(new.chr==i)
 ind.nbh <- which(nbh.chr==i)

 plot(-log10(as.numeric(gts.nbh$P.value[ind.nbh])),scan.nbh$lod[ind.nbh], col='blue', pch=16, ylim=c(0,20),xlim=c(0,5))
 ##points( -log10(as.numeric(gts.nbh.sen$P.value[ind])),scan.bin.mr$lod[ind], col='red', pch=16)
 ##points(-log10(as.numeric(gts.nbh$P.value[ind])),scan.bin.mr$lod[ind],  col='green', pch=16)
 points(-log10(as.numeric(gts.new$P.value[ind.new])),scan.new$lod[ind.new], col='green', pch=16)
  points(-log10(as.numeric(gts.elr$P.value[ind.elr])),scan.elr$lod[ind.elr], col='red', pch=16)

 text(2550,15,i)
}
dev.off()







png(paste0('~/public_html/NBH_lodxSeg.png'),height=900,width=900)
par(mfrow=c(5,5))
for (i in 1:24){

 ind.nbh <- which(nbh.chr==i)

 plot( -log10(as.numeric(gts.nbh.sen$P.value[ind.nbh])),scan.nbh$lod[ind.nbh], col='red', pch=16, ylim=c(0,20),xlim=c(0,5))
 points(-log10(as.numeric(gts.nbh.res$P.value[ind.nbh])),scan.nbh$lod[ind.nbh],  col='blue', pch=16)
 points(-log10(as.numeric(gts.nbh$P.value[ind.nbh])),scan.nbh$lod[ind.nbh], col='green')
 text(2550,15,i)
}
dev.off()


png(paste0('~/public_html/NEW_lodxSeg.png'),height=900,width=900)
par(mfrow=c(5,5))
for (i in 1:24){

 ind.new <- which(new.chr==i)

 plot( -log10(as.numeric(gts.new.sen$P.value[ind.new])),scan.new$lod[ind.new], col='red', pch=16, ylim=c(0,20),xlim=c(0,5))
  #plot( -log10(as.numeric(gts.new$P.value[ind.new])),scan.new$lod[ind.new], col='green', pch=16, ylim=c(0,20),xlim=c(0,5))

 points(-log10(as.numeric(gts.new.res$P.value[ind.new])),scan.new$lod[ind.new],  col='blue', pch=16)
 points( -log10(as.numeric(gts.new$P.value[ind.new])),scan.new$lod[ind.new], col='green')
 text(3,10,i)
}
dev.off()








png(paste0('~/public_html/ELR_lodxSeg.png'),height=900,width=900)
par(mfrow=c(5,5))
for (i in 1:24){

 ind.elr <- which(elr.chr==i)

 plot( -log10(as.numeric(gts.elr.sen$P.value[ind.elr])),scan.elr$lod[ind.elr], col='red', pch=16, ylim=c(0,20),xlim=c(0,5))
 points(-log10(as.numeric(gts.elr.res$P.value[ind.elr])),scan.elr$lod[ind.elr],  col='blue', pch=16)
 points(-log10(as.numeric(gts.elr$P.value[ind.elr])),scan.elr$lod[ind.elr], col='green')
 text(2550,15,i)
}
dev.off()











cross.par <- pull.markers(cross.par, markernames(cross))

par1 <- cbind(
          sapply(1:24,function(X){ table(pull.geno(cross.par,chr=X)[1,])['1'] } ),
          sapply(1:24,function(X){ table(pull.geno(cross.par,chr=X)[1,])['2'] } ),
          sapply(1:24,function(X){ table(pull.geno(cross.par,chr=X)[1,])['3'] } )
)

par2 <- cbind(
          sapply(1:24,function(X){ table(pull.geno(cross.par,chr=X)[2,])['1'] } ),
          sapply(1:24,function(X){ table(pull.geno(cross.par,chr=X)[2,])['2'] } ),
          sapply(1:24,function(X){ table(pull.geno(cross.par,chr=X)[2,])['3'] } )
)


cbind(1:24,par1,par2)


par2 <- sapply(1:24,function(X){ table(pull.geno(cross.par,chr=X)[2,]) } )


lapply(par1,"[[",1)


summary(scan.new)



a <- find.marker(new, 18, 84.4)
png("/home/jmiller1/public_html/NEW18_pxg.png")
plotPXG(new, a, pheno.col = 4, jitter = 1.5, infer = F, pch = 19, main = a)
dev.off()


a <- find.marker(new, 2, 25.6)
png("/home/jmiller1/public_html/NEW2_pxg.png")
plotPXG(new, a, pheno.col = 4, jitter = 1.5, infer = F, pch = 19, main = a)
dev.off()



a <- find.marker(nbh, 2, 79.34)
png("/home/jmiller1/public_html/NBH2_pxg.png")
plotPXG(nbh, a, pheno.col = 4, jitter = 1.5, infer = F, pch = 19, main = a)
dev.off()

a <- find.marker(nbh, 18, 26.75)
png("/home/jmiller1/public_html/NBH18_pxg.png")
plotPXG(nbh, a, pheno.col = 4, jitter = 1.5, infer = F, pch = 19, main = a)
dev.off()








nbh.sg <- sim.geno(nbh)


png("/home/jmiller1/public_html/nbh_effect_scan.png",width=10000)
effectscan(nbh.sg, pheno.col=5, get.se=T, draw=TRUE,
           gap=10, mtick="triangle",
           add.legend=TRUE, alternate.chrid=T)

dev.off()



nbh.rf <- markerlrt(subset(nbh,chr=2))
png("/home/jmiller1/public_html/markerlrt2.png",width=1000)
plotRF(nbh.rf)
dev.off()
