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
  before <- marker.density(cross,gt)
  cross <- drop.markers(cross,unlist(todrop))
  after <- marker.density(cross,gt[!rownames(gt) %in% unlist(todrop),])
  marker_dens[['drop.missing']] <<- list(before,after)
  return(cross)
}
distort <- function(cross,p){
  gt <- geno.table(cross)
  todrop <- rownames(gt[gt$P.value < p,])
  paste(length(todrop),'markers dropped')
  before <- marker.density(cross,gt)
  cross <- drop.markers(cross,unlist(todrop))
  after <- marker.density(cross,gt[!rownames(gt) %in% unlist(todrop),])
  marker_dens[['distort']] <<- list(before,after)

  return(cross)
}
distort.18 <- function(cross,p){
  gt <- geno.table(cross)
  todrop <- rownames(gt[gt$P.value < p,])
  todrop <- todrop[!todrop %in% tokeep]
  paste(length(todrop),'markers dropped')
  before <- marker.density(cross,gt)
  cross <- drop.markers(cross,unlist(todrop))
  after <- marker.density(cross,gt[!rownames(gt) %in% unlist(todrop),])
  marker_dens[['distort.18']] <<- list(before,after)

  return(cross)
}
drop.missing.18 <- function(cross,missing){
  M <- nind(cross.18)-round(nind(cross.18)*missing)
  gt <- geno.table(cross)
  todrop <- rownames(gt[which(gt$missing>M),])
  todrop <- todrop[!todrop %in% tokeep]
  paste(length(todrop),'markers dropped')
  before <- marker.density(cross,gt)
  cross <- drop.markers(cross,unlist(todrop))
  after <- marker.density(cross,gt[!rownames(gt) %in% unlist(todrop),])
  marker_dens[['drop.missing.18']] <<- list(before,after)
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
dropone.par <- function(cross,p,chr,map.function = 'kosambi',
  maxit,sex.sp = F,verbose=F,parallel=T,error.prob = 0.03){
  y <- round(sum(nmar(cross))*p)
  ### p = percent of longest markers to drop
  cross.drops <- parallel.droponemarker(cross,
        chr, map.function = 'kosambi',maxit=12,
        sex.sp = F,verbose=F,parallel=T)
  index.lod <- 

  cross.drops <-
  ### Positive value in Ldif = decrease in length
  ### Positive value in LOD = increase in ocerall lod
  todrop <- rownames(cross.drops)[1:y]
  cross <- drop.markers(cross.18,unlist(todrop))
  return(cross)
}
marker.warning <- function(cross=cross.18){
  print(paste('Starting markers mapped =',
    sum(markernames(cross) %in% markernames(cross.18,chr=X))))

  print(paste('Starting markers un-mapped =',
    sum(!markernames(cross) %in% markernames(cross.18,chr=X))))

}
er.rate <- function(cross){
  loglik <- err <- c(0.005, 0.01, 0.015,
     0.02,0.025, 0.03, 0.04, 0.05)
      registerDoParallel(slurmcore)
      hoods <- foreach(i=seq(along=err),
        .inorder=T,.packages = "qtl") %dopar% {
        tempmap <- est.map(cross, error.prob=err[i],maxit=1000)
        return(attr(tempmap[[1]], "loglik"))
        #loglik[i] <- attr(tempmap[[1]], "loglik")
      }
      return(err[which.max(abs(unlist(hoods)))])
}

all.crossed <- function(X=X,i=pop){
  read.cross(format='csv',file=X,geno=c('AA','AB','BB'),
  alleles=c("A","B"))
}

drop.errlod <- function(cross,lod=lod,ers=ers){
  print(paste('cut below erlod >',lod))
  mapthis <- calc.errorlod(cross, error.prob=ers)
  toperr <- top.errorlod(mapthis, cutoff=lod)
  if (toperr)
  step <- sort(as.numeric(names(table(round(toperr$errorlod)))), decreasing=T)
  sapply(step,function(Z){
    mapthis <- calc.errorlod(cross, error.prob=ers)
    toperr <- top.errorlod(mapthis, cutoff=Z)
    apply(toperr,1, function(D){
      D <- as.character(D)
      cross$geno[[D[1]]]$data[cross$pheno$ID==D[2], D[3]] <<- NA
      }
    )
    }
  )
  return(cross)
}
reconst <- function(X,pop,out){
  temp <- file.path(basedir,'rQTL',pop,paste('REMAPS/temp.',X,sep=''))
  myfiles <- lapply(temp, function(X,pop=pop){
    all.crossed(X)
    }
  )
  ID <- myfiles[[1]]$pheno$ID
  pheno <- myfiles[[1]]$pheno$Pheno
  sex <- myfiles[[1]]$pheno$sex

  map <- unlist(sapply(seq(along=myfiles),
    function(i){myfiles[[i]]$geno[[1]]$map}))

  chr <-  unlist(sapply(seq(along=myfiles),
    function(i){sapply(1:nmar(myfiles[[i]]),
      function(Z)chrnames(myfiles[[i]]))}))

  registerDoParallel(slurmcore)

  cross <- foreach(i=seq(along=myfiles),
    .combine=cbind,.packages = "qtl") %dopar% {
      marks <- colnames(myfiles[[i]]$geno[[1]]$data)
      data <- myfiles[[i]]$geno[[1]]$data
      colnames(data) <- marks
      data
  }
  chr <- c('','','',chr)
  map <- c('','','',map[colnames(cross)])
  cross <- cbind(pheno,sex,ID=as.character(ID),cross)
  cross <- rbind(colnames(cross),chr,map,cross)

  write.table(cross,file=out,
      col.names=F,row.names=F,quote=F,sep=',')

  return(read.cross.jm(file=out,format='csv',
    geno=c(1:3),estimate.map=FALSE))
}
