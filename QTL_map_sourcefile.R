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
dropone.par <- function(cross,p,chr,map.function = 'kosambi',length.imp = 1, LOD.imp = 0,
  maxit,sex.sp = F,verbose=F,parallel=T,error.prob = 0.03){
  y <- round(sum(nmar(cross))*p)
  ### p = percent of longest markers to drop
  cross.drops <- parallel.droponemarker(cross,
        chr, map.function = 'kosambi',maxit=12,
        sex.sp = F,verbose=F,parallel=T)

  Len <- quantile(as.numeric(cross.drops$Ldiff),p)
  Lod <- quantile(as.numeric(cross.drops$LOD),p)

  index.lod <- rownames(cross.drops[which(as.numeric(cross.drops$LOD) > Lod & as.numeric(cross.drops$LOD) > LOD.imp),])
  index.ldif <- rownames(cross.drops[which(as.numeric(cross.drops$Ldiff) > Len & as.numeric(cross.drops$Ldiff) > length.imp),])
  drops <- unique(c(index.lod,index.ldif))
  if (length(drops)>0){
    cross <- drop.markers(cross.18,unlist(drops))
    print(paste('dropping',cross.drops[drops,1]))
    print(paste('dropping',cross.drops[drops,3]))
    print(paste('dropping',cross.drops[drops,4]))
  } else {
    print('no drops made')
  }
  ### Positive value in Ldif = decrease in length
  ### Positive value in LOD = increase in ocerall lod
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
drop.errlod <- function(cross,lod=lod,ers=ers){
  print(paste('remove genotypes with erlod >',lod))
  mapthis <- calc.errorlod(cross, error.prob=ers)
  toperr <- top.errorlod(mapthis, cutoff=lod)
  dropped <- 0
  if (length(toperr[,1]) > 0){
    step <- sort(as.numeric(names(table(floor(toperr$errorlod)))), decreasing=T)
    while (!sum(step)==0) {
      sapply(step,function(Z){
        apply(toperr[which(toperr$errorlod>Z),],1, function(marks){
          marks <- as.character(marks)
          cross$geno[[marks[1]]]$data[cross$pheno$ID==marks[2], marks[3]] <<- NA
          dropped <<- dropped + 1
          }
        )
        mapthis <<- calc.errorlod(cross, error.prob=ers)
        toperr2 <<- top.errorlod(mapthis, cutoff=lod)
          if (!is.null(toperr2)){ step <<- sort(as.numeric(names(table(floor(toperr2$errorlod)))), decreasing=T)
          } else { step <<- 0 }
        }
      )
    }
  }
  print(paste('done...dropped',dropped,'genotypes'))
  return(cross)
}
all.crossed <- function(X){
  read.cross(format='csv',file=X,geno=c('AA','AB','BB'),
  alleles=c("A","B"))
}

reconst <- function(X,pop,dir){
  temp <- file.path(basedir,'rQTL',pop,paste('REMAPS/temp.',X,sep=''))
  myfiles <- lapply(temp, function(tocross){
    all.crossed(tocross)
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

  write.table(cross,file=file.path(dir,'tempout'),
      col.names=F,row.names=F,quote=F,sep=',')

  return(read.cross.jm(file=file.path(dir,'tempout'),format='csv',
    geno=c(1:3),estimate.map=FALSE))
}
