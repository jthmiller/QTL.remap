## Determine which individuals to keep
### SAMPLE 10869 in ELR is really odd


source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')
#source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/debug.R')

cross <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
  format='csvr', geno=c(1:3),estimate.map=FALSE)

# Colors for the heat map
my_palette <- colorRampPalette(c("white", "blue", "white"))(n = 74)
col_breaks = c(seq(0,0.3,length=30),seq(0.31,0.41,length=5),seq(0.42,1,length=40))

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
miss <- 10
cutoff <- 1.0e-06
cross <- drop.markers(cross,rownames(gt.cross.par[gt.cross.par$missing>miss,]))
cross <- drop.markers(cross,rownames(gt.cross.par[gt.cross.par$P.value<cutoff,]))
cross.max <- subset(cross, ind=nmissing(cross)<median(nmissing(cross)))

### Drop
drop <- markernames(cross.pars)[!markernames(cross.pars) %in% markernames(cross)]
cross.pars  <- drop.markers(cross.pars,drop)
cross.max <- c(cross.pars,cross.max)
cross <- c(cross.pars,cross)

rela <- rels(cross)
## Calculate matrix
rela <- rels(cross.max)
name <- paste(pop,'_kinship_heatmap_before_filter.pdf',sep='')
main <- paste(pop,'kinship before filter (proportion of shared genotypes, 0-1)')
feet(rela,name,main)

#Drop dup inds
wh <- which(rela > 0.75, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],] ## just take larger index of all
wh2 <- wh[,2]
names(wh2) <- rownames(rela)[wh2]
wh <- wh[,1]
numiss <- nmissing(cross)
in.drop <- sapply(1:length(wh),function(X){
  gdn <- which.max(c(numiss[wh[X]],numiss[wh2[X]]))
  return(names(c(wh[X],wh2[X]))[gdn])
  }
)

### Lowdata
cross.min <- subset(cross, ind=nmissing(cross)>median(nmissing(cross)) | is.na(cross$pheno$Pheno))
### Drop Simiar ind
cross <- subset(cross, ind=(!gsub(paste(pop,'_',sep=''),'',cross$pheno$ID) %in% in.drop) | is.na(cross$pheno$Pheno))
### Drop ind with greater than 50% missing data
cross <- subset(cross, ind=nmissing(cross)/sum(nmar(cross)) < 0.50 | is.na(cross$pheno$Pheno))
cross.max <- subset(cross, ind=nmissing(cross)<median(nmissing(cross)) | is.na(cross$pheno$Pheno))

## Calculate matrix
rela <- rels(cross.max)
name <- paste(pop,'_kinship_histogram_to.map.pdf',sep='')
main <- paste(pop,'kinship before mapping (proportion of shared genotypes, 0-1)')
feet(rela,name,main)

## Lowdata
rela <- rels(cross.min)
name <- paste(pop,'_kinship_histogram_low_data.map.pdf',sep='')
main <- paste(pop,'kinship low data (proportion of shared genotypes, 0-1)')
feet(rela,name,main)

fileConn <- file(file.path(popdir,'kinship.keep.ind.txt'))
writeLines(cross$pheno$ID, fileConn)
close(fileConn)
