## Libraries
packs <- c('qtl','foreach','doParallel')
lapply(packs, require, character.only = TRUE)
## Load a couple fixed rQTL functions
source(file.path(basedir,'rQTL/scripts/QTL_remap/custom_rQTL_functions.R'))

marker.density <- function(cross,gt.1){
  index <- gsub('\\:.*','',rownames(gt.1))==X
  x <- order(as.numeric(gsub('.*\\:','',rownames(gt.1)))[index])
  y <- sort(as.numeric(gsub('.*\\:','',rownames(gt.1)))[index])
  return(list(cor=cor(x,y)^2,pos=y))
}
fix.pheno <- function(cross){

  cross$pheno$ID <- paste("ind", 1:nind(cross), sep="")
  cross$pheno[which(cross$pheno[,1]<2),1] <- 0
  cross$pheno[which(cross$pheno[,1]>1),1] <- 1
  cross <- subset(cross, ind=(!is.na(cross$pheno$Pheno))) ## drop g.parents from main set
  return(cross)
}
drop.missing <- function(cross,M){
  gt <- geno.table(cross)
  todrop <- rownames(gt[which(gt$missing>M),])
  paste(length(todrop),'markers dropped')
  marker_dens[['drop.missing']][['before']] <<- marker.density(cross,gt)
  cross <- drop.markers(cross,unlist(todrop))
  gt[!rownames(gt) %in% unlist(todrop),]
  marker_dens[['drop.missing']]['aft'] <<- marker.density(cross,gt[!rownames(gt) %in% unlist(todrop),])
  return(cross)
}
distort <- function(cross,p){
  gt <- geno.table(cross)
  todrop <- rownames(gt[gt$P.value < p,])
  paste(length(todrop),'markers dropped')
  marker_dens[['distort']][['before']] <<- marker.density(cross,gt)
  cross <- drop.markers(cross,unlist(todrop))
  gt[!rownames(gt) %in% unlist(todrop),]
  marker_dens[['distort']]['aft'] <<- marker.density(cross,gt[!rownames(gt) %in% unlist(todrop),])

  return(cross)
}
distort.18 <- function(cross,p){
  gt <- geno.table(cross)
  todrop <- rownames(gt[gt$P.value < p,])
  todrop <- todrop[!todrop %in% tokeep]
  paste(length(todrop),'markers dropped')
  marker_dens[['distort.18']][['before']] <<- marker.density(cross,gt)
  cross <- drop.markers(cross,unlist(todrop))
  gt[!rownames(gt) %in% unlist(todrop),]
  marker_dens[['distort.18']]['aft'] <<- marker.density(cross,gt[!rownames(gt) %in% unlist(todrop),])

  return(cross)
}
drop.missing.18 <- function(cross,missing){
  M <- nind(cross.18)-round(nind(cross.18)*missing)
  gt <- geno.table(cross)
  todrop <- rownames(gt[which(gt$missing>M),])
  todrop <- todrop[!todrop %in% tokeep]
  paste(length(todrop),'markers dropped')
  marker_dens[['drop.missing.18']][['before']] <<- marker.density(cross,gt)
  cross <- drop.markers(cross,unlist(todrop))
  gt[!rownames(gt) %in% unlist(todrop),]
  marker_dens[['drop.missing.18']]['aft'] <<- marker.density(cross,gt[!rownames(gt) %in% unlist(todrop),])
  return(cross)
}
drop.mark <- function(crossZ,y){
  for (t in 1:y){
     dropone <- droponemarker(crossZ, chr=X, error.prob=0.002, maxit=5)
     badmar <- rownames(summary(dropone, lod.column=1))
     crossZ <- drop.markers(crossZ, badmar)
   }
    return(crossZ)
}
keepQTL <- function(Z,i){
  pos.m <- test.QTLs$str[qtl.index]
  gt.1 <- geno.table(i)
  chr <- gsub('\\:.*','',rownames(gt.1))==X
  gt.1 <- gt.1[chr,]
  pos <- as.numeric(gsub('.*\\:','',rownames(gt.1)))
  markerVec <- row.names(gt.1[order(abs(pos.m-pos)) < 10,])
  return(markerVec)
}
dropone.par <- function(cross,p=0.05,chr=X,map.function = 'kosambi',
  maxit=2,sex.sp = F,verbose=F,parallel=T){
  y <- round(sum(nmar(cross))*p) #get rid of 5% of markers
  cross.drops <- parallel.droponemarker(cross.18,
        chr=X,map.function = 'kosambi',maxit=2,
        sex.sp = F,verbose=F,parallel=T)
  todrop <- rownames(cross.drops)[1:y]
  cross <- drop.markers(cross.18,unlist(todrop))
  return(cross)
}
