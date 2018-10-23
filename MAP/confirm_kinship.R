## Libraries
packs <- c('qtl','foreach','doParallel')
lapply(packs, require, character.only = TRUE)
require(qtl2,lib.loc='/share/apps/rmodules')

## For plotting
marker_dens <- list()

### parents ###
cross.pars <- read.cross.jm(file=file.path(indpops,paste(pop,'.parents.csvr',sep='')),
                format='csvr', geno=c(1:3),estimate.map=FALSE)

gt.par <- geno.table(cross.pars)
par.confirm.marks <- rownames(gt.par[which(gt.par$AA==1 & gt.par$BB==1),])

## read in the QTL cross
cross.18 <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
                format='csvr', geno=c(1:3),estimate.map=FALSE)

cross.pars <- subset(cross.pars, chr=c(1:24))

cross.18 <- c(cross.pars,cross.18)

## Map each QTL chro independently
allbut <- c(1:24)[-c(1,2,3,18)]
subset.qtl <- chrnames(cross.18)[!chrnames(cross.18) %in% allbut]
cross.18 <- subset(cross.18, chr=subset.qtl)

## Specific to 'outname'
if(mapped.only==T){
cross.18 <- drop.markers(cross.18,markernames(cross.18)[grep('NW',markernames(cross.18))])
}

#### Filter Conservative
print('Dropping markers with more than 5 genotypes missing')
cross.18 <- drop.missing(cross.18,7)
gt.missing <- geno.table(cross.18)

### invariants
gt.cross.par <- geno.table(cross.18)
cross.18 <- drop.markers(cross.18,rownames(gt.cross.par[gt.cross.par$P.value<cutoff,]))
gt.pval <- geno.table(cross.18)

cross.18 <- formLinkageGroups(cross.18, max.rf=0.1, reorgMarkers=TRUE)
cross.18 <- subset(cross.18, chr=c(1:10))


cross2 <- convert2cross2(cross.18)
map <- insert_pseudomarkers(cross2$gmap, step=1)
pr <- calc_genoprob(cross2, map, err=ers, cores=slurmcore)
pr <- clean_genoprob(pr)
apr <- genoprob_to_alleleprob(pr)

rela <- calc_kinship(pr, type = "overall", omit_x = FALSE,
  use_allele_probs = F, quiet = TRUE, cores = 12)

save.image('kinship.rsave')

############

packs <- c('qtl','foreach','doParallel','gplots')
lapply(packs, require, character.only = TRUE)
require(qtl2,lib.loc='/share/apps/rmodules')
#load('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/kinship.rsave')
load('~/Dropbox/QTL_Paper/Kinship.analysis.Rsave')
load('~/Dropbox/QTL_Paper/plink.NBH.Rsave')
NBH.kinship.Rsave
plot(1,1)
#heatmap(rela,Rowv = NA,Colv = NA)

M <- rela$"24"
#M[lower.tri(M)] <- NA
#heatmap.2(M,scale="row", dendrogram="row",symm = T)
heatmap(M,symm = T,Rowv=NULL,scale="row")

hilo <- lapply(rela,function(M){
a <- names(table(which(M == min(M), arr.ind = TRUE)))
#b <- table(which(M == max(M), arr.ind = TRUE))
#return(list(a,b))
return(a)
}
)


hilo <- lapply(rela,function(M){
return(min(M))
}
)

hilo <- lapply(rela,function(M){
return(table(which(M == min(M), arr.ind = TRUE)))
}
)





#heatmap(M,Rowv = NA,Colv = NA)

## indv 28,47,84 look odd
## which are the least related among fak linkage groups?  Should be parents
table(unlist(lapply(rela,function(mm){
  names(table(rownames(which(mm == min(mm), arr.ind = TRUE))))
  })))
ela <- calc_kinship(pr, type = "chr", omit_x = FALSE,
  use_allele_probs = T, quiet = TRUE, cores = 12)


####use plink (captured NBH parents ok)
ids <- read.table('~/Dropbox/QTL_Paper/NBH.rel.id')
n.len <- file.info('~/Dropbox/QTL_Paper/NBH.rel.bin')$size
dis <- readBin('~/Dropbox/QTL_Paper/NBH.rel.bin',what='double',n=n.len)

### ELR dis[c('ER1124F','BI1124M'),] dis['ER1124F','BI1124M']
packs <- c('qtl','foreach','doParallel','gplots','RColorBrewer','qtl2')
lapply(packs, require, character.only = TRUE)
ids <- read.table('~/Dropbox/QTL_Paper/ELR.rel.id')
n.len <- file.info('~/Dropbox/QTL_Paper/ELR.rel.bin')$size
dis <- readBin('~/Dropbox/QTL_Paper/ELR.rel.bin',what='double',n=n.len)
dis <- as.numeric(dis)
dis <- matrix(dis, nrow=length(ids$V1),ncol=length(ids$V1))
colnames(dis) <- ids$V2
rownames(dis) <- ids$V2
dis[dis==NaN] <- NA
diag(dis) <- NA

#dis[lower.tri(dis)] <- NA
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
my_palette <- colorRampPalette(c("red","green","red"))(n = 299)

col_breaks = c(seq(-0.5,0.0,length=100),seq(0.01,0.5,length=100),seq(0.51,2,length=100))

my_palette <- colorRampPalette(c("green","red"))(n = 299)
col_breaks = c(seq(-0.5,0.5,length=150),seq(0.51,1,length=150))


heatmap(dis,Colv=NA,Rowv=NA,col=my_palette,scale='none',symm =T,breaks=col_breaks)



ord <- names(sort(dis['BI1124M',],decreasing = T))
heatmap(dis,Colv=ord,col=hmcol,Rowv=NA)



heatmap(dis,symm = T,Rowv=NULL,scale="row",col=hmcol)
heatmap(lower.tri(dis))

ord <- names(sort(dis['NBH1M',],decreasing = T))


ord <- order(dis['BI1124M',],decreasing = F)
heatmap(dis,symm = F,Colv=ord,Rowv=NA,scale="row",col=hmcol)

plot(na.omit(dis))

dis.dist <- dist(dis)
try <- isoMDS(dis.dist,eig=TRUE, k=2)
x <- try$points[,1]
y <- try$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
  main="Metric	MDS",	type="n")
text(x, y, labels = rownames(try$points), cex=.7)

which(dis==max(dis,na.rm =T),arr.ind = TRUE)
which(dis==min(dis,na.rm =T),arr.ind = TRUE)

hist(dis)




sort(table(which(rela>.8,arr.ind = TRUE)))

hist.data = hist(dis[which(dis>0.4)], plot=F)
hist.data = hist(dis)
#hist.data$counts = log(hist.data$counts, 10)
#hist.data$counts[hist.data$counts==-Inf] <- 0
plot(hist.data)

plot(mydata_hist$count, log="y", type='h', lwd=10, lend=2)
table(rownames(which(dis>0.5,arr.ind = TRUE)))
sort(colSums(abs(dis)))

ord <- order(dis['ER1124F',])
plot(dis['ER1124F',ord])
points(dis['BI1124M',ord],col='red')


dis[c('ER1124F','BI1124M'),]

##In ELR, IND10869 is more distant from parents than they are from one another





allbut <- rownames(dis)[-(which(rownames(dis)=='5528'))]
dis <- dis[allbut,allbut]

diag(dis) <- rowMeans(dis)



#####3 In RQTL #####
rela <- comparegeno(cross.18)

save.image('kinship.rsave')

wh <- which(rela < 0.15, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]

g <- pull.geno(cross)

table(g[wh[,1],], g[wh[,2],])


table(g[214,], g[216,])
table(g[238,], g[288,])


> for(i in 1:nrow(wh)) {
+ tozero <- !is.na(g[wh[i,1],]) & !is.na(g[wh[i,2],]) & g[wh[i,1],] != g[wh[i,2],]
+ mapthis$geno[[1]]$data[wh[i,1],tozero] <- NA +}
mapthis <- subset(mapthis, ind=-wh[,2])




rf <- pull.rf(cross.18)
lod <- pull.rf(cross.18, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")
