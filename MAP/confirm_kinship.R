### DATA
ids <- read.table('~/Dropbox/QTL_Paper/ELR.rel.id')
n.len <- file.info('~/Dropbox/QTL_Paper/ELR.rel.bin')$size
dis <- readBin('~/Dropbox/QTL_Paper/ELR.rel.bin',what='double',n=n.len)
con <- file(file.path('~/Dropbox/QTL_Paper/Rough Figures/Kinship Analysis/ELR.kinship.keep.ind.txt'), open='r')
### DATA
ids <- read.table('~/Dropbox/QTL_Paper/NBH.rel.id')
n.len <- file.info('~/Dropbox/QTL_Paper/NBH.rel.bin')$size
dis <- readBin('~/Dropbox/QTL_Paper/NBH.rel.bin',what='double',n=n.len)
con <- file(file.path('~/Dropbox/QTL_Paper/Rough Figures/Kinship Analysis/NBH.kinship.keep.ind.txt'), open='r')
#######
### DATA
ids <- read.table('~/Dropbox/QTL_Paper/NEW.rel.id')
n.len <- file.info('~/Dropbox/QTL_Paper/NEW.rel.bin')$size
dis <- readBin('~/Dropbox/QTL_Paper/NEW.rel.bin',what='double',n=n.len)
con <- file(file.path('~/Dropbox/QTL_Paper/Rough Figures/Kinship Analysis/NEW.kinship.keep.ind.txt'), open='r')
#######

keepers <- readLines(con)
close(con)

### TRANSFORM
dis <- as.numeric(dis)
dis <- matrix(dis, nrow=length(ids$V1),ncol=length(ids$V1))
colnames(dis) <- ids$V2
rownames(dis) <- ids$V2

nam <- unlist(sapply(strsplit(keepers,'_'),'[[',2))
nam <- rownames(dis) %in% nam
dis <- dis[nam,nam]
dis[dis==NaN] <- NA
diag(dis) <- 0
d <- dist(dis)
#### MDS COLORS
cols <- as.numeric(rownames(dis))
cols[!is.na(cols)] <- 'green'
cols[is.na(cols)] <- 'red'
names(cols) <- rownames(dis)
### MDS PLOT
fit <- cmdscale(as.dist(1-dis),eig=TRUE, k=2)
fit <- cmdscale(as.dist(sqrt(1-dis)),eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]
png('~/Dropbox/QTL_Paper/Rough Figures/Kinship Analysis/ELR.MDS.png')
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
  main="Metric	MDS",	type="n")
text(x, y, labels = row.names(dis),col=cols, cex=2)
dev.off()
## qGraph
heatmap(dis,symm=T,main=NA,labRow=NULL,na.rm=T,cexRow=1,labCol=NULL,key.title=NA,key.xlab=NA,key.ylab=NA,
  srtCol=45,cexCol=1,dendrogram="row", margins = c(15,10),trace="none",keysize=.5)


##




qgraph(dist_mi, layout='spring', vsize=3)
qgraph(dis, layout='spring', vsize=3)
qgraph(dis,layout = L, directed = TRUE)

qgraph(d, layout='spring', vsize=3)

pca <- prcomp(dis, scale=T)
plot(pca$x, pch=20, col="blue", type="n")
text(pca$x, rownames(pca$x), cex=0.8)

m <- dis
mds <- cmdscale(as.dist(1-m))
plot(mds,pch=20,col=cols)
