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


chr5; 25 max, 25x49 min
chr6; 25X88 max ,min 25,88,91
chr7; 3 15 71 , min 3 15 71 72
chr



rela <- calc_kinship(pr, type = "chr", omit_x = FALSE,
  use_allele_probs = T, quiet = TRUE, cores = 12)


use plink (captured NBH parents ok)
ids <- read.table('ELR.rel.id')
dis <- readBin('ELR.rel.bin',what='numeric',n=9216)
dis <- matrix(dis, nrow=89,ncol=89)
colnames(dis) <- ids$V2
rownames(dis) <- ids$V2

dis[dis==NaN] <- 0
