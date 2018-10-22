## Determine which individuals to keep

source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')

cross <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
  format='csvr', geno=c(1:3),estimate.map=FALSE)

## Pull names from plinkfile
path <- file.path(indpops,paste(pop,'.ped',sep=''))
popname <- system(paste('cut -f1 -d\' \'',path),intern = TRUE)
indname <- system(paste('cut -f2 -d\' \'',path),intern = TRUE)
cross$pheno$ID <- paste(popname,indname,sep='_')

### Parent Cross
cross.pars <- subset(cross,ind=is.na(cross$pheno$Pheno))
### Drop parents before filter
cross <- subset(cross,ind=!is.na(cross$pheno$Pheno))
gt.cross.par <- geno.table(cross)
miss <- 8
cutoff <- 1.0e-06
cross <- drop.markers(cross,rownames(gt.cross.par[gt.cross.par$missing>miss,]))
cross <- drop.markers(cross,rownames(gt.cross.par[gt.cross.par$P.value<cutoff,]))

### Drop
drop <- markernames(cross.pars)[!markernames(cross.pars) %in% markernames(cross)]
cross.pars  <- drop.markers(cross.pars,drop)
cross <- c(cross.pars,cross)

## Calculate matrix
rela <- comparegeno(cross)
colnames(rela) <- cross$pheno$ID
rownames(rela) <- cross$pheno$ID
rela[rela==NaN] <- NA
diag(rela) <- NA
rela <- rela[rowSums(is.na(rela)) < nind(cross),colSums(is.na(rela)) < nind(cross)]

main = paste(pop,'kinship before filter (proportion of shared genotypes, 0-1)')
png(file.path(popdir,paste(pop,'_kinship_heatmap_before_filter.png',sep='')))
heatmap(rela,symm=T,main=main,sub="red 0 to 1.0 yellow/white" )
dev.off()

#Drop dup inds
wh <- which(rela > 0.75, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh2 <- wh[,2]
names(wh2) <- rownames(rela)[wh2]
wh <- wh[,1]
numiss <- nmissing(cross)
in.drop <- sapply(1:length(wh),function(X){
  gdn <- which.max(c(numiss[wh[X]],numiss[wh2[X]]))
  return(names(c(wh[X],wh2[X]))[gdn])
  }
)
### Drop Simiar ind
cross <- subset(cross, ind=!cross$pheno$ID %in% in.drop)
### Drop ind with greater than 50% missing data
cross <- subset(cross, ind=nmissing(cross)/sum(nmar(cross)) < 0.50)

## Calculate matrix
rela <- comparegeno(cross)
colnames(rela) <- cross$pheno$ID
rownames(rela) <- cross$pheno$ID
rela[rela==NaN] <- NA
diag(rela) <- NA
rela <- rela[rowSums(is.na(rela)) < nind(cross),colSums(is.na(rela)) < nind(cross)]

main = paste(pop,'kinship before mapping (proportion of shared genotypes, 0-1)')

png(file.path(popdir,paste(pop,'_kinship_histogram_to.map.png',sep='')))
hist(rela,main=main,sub=sum(nmar(cross)))
dev.off()

png(file.path(popdir,paste(pop,'_kinship_heatmap_to.map.png',sep='')))
heatmap(rela,symm=T,main=main,sub="red 0 to 1.0 yellow/white" )
dev.off()

fileConn <- file(file.path(popdir,'kinship.keep.ind.txt'))
writeLines(cross$pheno$ID, fileConn)
close(fileConn)
