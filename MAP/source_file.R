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
dropone.par <- function(cross,chr,prop=0.025,map.function = c("haldane",
    "kosambi", "c-f", "morgan"),length.imp = 1, LOD.imp = 0,
  maxit=1,sex.sp = F,verbose=F,parallel=T,error.prob = 0.03,cores=slurmcore)

  {
  y <- round(sum(nmar(cross))*prop)
  ### p = percent of longest markers to drop
  print('starting parallel.droponemarker')
  cross.drops <- parallel.droponemarker(cross,chr,maxit,cores,map.function='kosambi')

  Len <- quantile(as.numeric(cross.drops$Ldiff),prop)
  Lod <- quantile(as.numeric(cross.drops$LOD),prop)

  index.lod <- rownames(cross.drops[which(as.numeric(cross.drops$LOD) > Lod & as.numeric(cross.drops$LOD) > LOD.imp),])
  index.ldif <- rownames(cross.drops[which(as.numeric(cross.drops$Ldiff) > Len & as.numeric(cross.drops$Ldiff) > length.imp),])
  drops <- unique(c(index.lod,index.ldif))
  if (length(drops)>0){
    cross <- drop.markers(cross.18,unlist(drops))
    print(paste('dropping',cross.drops[drops,1],cross.drops[drops,3],cross.drops[drops,4]))
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
er.rate <- function(cross,slurmcore){
  loglik <- err <- c(0.005, 0.01, 0.015,
     0.02,0.025, 0.03, 0.04, 0.05)
      registerDoParallel(slurmcore)
      hoods <- foreach(i=seq(along=err),
        .inorder=T,.packages = "qtl") %dopar% {
        tempmap <- est.map(cross, error.prob=err[i],maxit=10000)
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
markersInInterval <- function(cross, chr, min, max) {

 names(which(pull.map(cross=cross, chr=chr)[[chr]] < max &

             pull.map(cross=cross, chr=chr)[[chr]] > min))

}
singleMarkerInInterval <- function(cross, chr, min, max) {

  tmp <- markersInInterval(cross,chr,min,max)

  val <- ifelse(sum(!is.na(tmp) == 1 | is.na(tmp) == 2), tmp[!is.na(tmp)], FALSE)

}
removeDoubleXO <- function(cross, chr, verbose=TRUE) {

  if (!missing(chr))

      chr <- matchchr(chr, names(cross$geno))

  else chr <- names(cross$geno)

  for (ch in chr) {  # loop through all linkage groups

    if (verbose)

      cat("Starting linkage group",ch,"\n")



    # find all recombination events in cross. xo is a list spanning

    # individuals.  each element of xo is either NULL or a numeric

    # vector with estimated crossover locations



    xo <- locateXO(cross, ch)



    # initialize some variables to keep track of the number

    # of changes



    total_removed <- 0

    tot_genotypes <- length(cross$geno[[ch]]$data[,]) -

                          sum(is.na(cross$geno[[ch]]$data[,]))



    for (ind in 1:length(xo)) {            # loop through individuals

      if (length(xo[[ind]]) <= 1)

        next  # skip individuals with one or fewer recombination events

      # walk along the linkage groups recombination events

      for (location in 3:length(xo[[ind]])-1) {



        # Determine if there is a single marker between each recombination

        # event and the next (location -1 thru location)



        sMar <- singleMarkerInInterval(cross,ch,xo[[ind]][location-2],xo[[ind]][location+1])

        if (sMar!=FALSE) { # if there are double recombination events

          oldValue <- cross$geno[[ch]]$data[ind,sMar] # original genotype call

          cross$geno[[ch]]$data[ind,sMar] <- NA # assign genotype to NA

          total_removed <- total_removed + 1 # count removed genotypes

          if (verbose>1)

            cat("individual", ind, "marker", sMar, "oldvalue=", oldValue,

                "newvalue = ", cross$geno[[ch]]$data[ind,sMar], "\n")

        }

      }

    }

    if (verbose) {

      cat("Removed",total_removed,"of",tot_genotypes,

          "genotyped markers on linkage group",ch,"\n")

    }

  }

  cross

}
read.cross.jm <- function (format = c("csv", "csvr", "csvs", "csvsr", "mm", "qtx",
    "qtlcart", "gary", "karl", "mapqtl", "tidy"), dir = "", file,
    genfile, mapfile, phefile, chridfile, mnamesfile, pnamesfile,
    na.strings = c("-", "NA"), genotypes = c("A", "H", "B", "D",
        "C"), alleles = c("A", "B"), estimate.map = FALSE, convertXdata = TRUE,
    error.prob = 1e-04, map.function = c("haldane", "kosambi",
        "c-f", "morgan"), BC.gen = 0, F.gen = 0, crosstype, ...)
{
    if (format == "csvrs") {
        format <- "csvsr"
        warning("Assuming you mean 'csvsr' rather than 'csvrs'.\n")
    }
    format <- match.arg(format)
    if (format == "csv" || format == "csvr") {
        cross <- read.cross.csv(dir, file, na.strings, genotypes,
            estimate.map, rotate = (format == "csvr"), ...)
    }
    else if (format == "csvs" || format == "csvsr") {
        if (missing(phefile) && !missing(file) && !missing(genfile)) {
            phefile <- genfile
            genfile <- file
        }
        else if (missing(genfile) && !missing(file) && !missing(phefile)) {
            genfile <- file
        }
        cross <- read.cross.csvs(dir, genfile, phefile, na.strings,
            genotypes, estimate.map=FALSE, rotate = (format == "csvsr"),
            ...)
    }
    else if (format == "qtx") {
        cross <- read.cross.qtx(dir, file, estimate.map)
    }
    else if (format == "qtlcart") {
        if (missing(mapfile) && !missing(genfile))
            mapfile <- genfile
        cross <- read.cross.qtlcart(dir, file, mapfile)
    }
    else if (format == "karl") {
        if (missing(genfile))
            genfile <- "gen.txt"
        if (missing(mapfile))
            mapfile <- "map.txt"
        if (missing(phefile))
            phefile <- "phe.txt"
        cross <- read.cross.karl(dir, genfile, mapfile, phefile)
    }
    else if (format == "mm") {
        if (missing(mapfile) && !missing(genfile))
            mapfile <- genfile
        cross <- read.cross.mm(dir, file, mapfile, estimate.map)
    }
    else if (format == "gary") {
        if (missing(genfile))
            genfile <- "geno.dat"
        if (missing(mnamesfile))
            mnamesfile <- "mnames.txt"
        if (missing(chridfile))
            chridfile <- "chrid.dat"
        if (missing(phefile))
            phefile <- "pheno.dat"
        if (missing(pnamesfile))
            pnamesfile <- "pnames.txt"
        if (missing(mapfile))
            mapfile <- "markerpos.txt"
        cross <- read.cross.gary(dir, genfile, mnamesfile, chridfile,
            phefile, pnamesfile, mapfile, estimate.map, na.strings)
    }
    else if (format == "mapqtl") {
        cross <- read.cross.mq(dir = dir, locfile = genfile,
            mapfile = mapfile, quafile = phefile)
    }
    else if (format == "tidy") {
        if (!missing(file) && !missing(genfile) && !missing(mapfile) &&
            missing(phefile)) {
            phefile <- mapfile
            mapfile <- genfile
            genfile <- file
        }
        if (missing(genfile))
            genfile <- "gen.csv"
        if (missing(phefile))
            phefile <- "phe.csv"
        if (missing(mapfile))
            mapfile <- "map.csv"
        cross <- read.cross.tidy(dir = dir, genfile = genfile,
            phefile = phefile, mapfile = mapfile, na.strings = na.strings,
            genotypes = genotypes)
    }
    estimate.map <- cross[[2]]
    cross <- cross[[1]]
    chrnam <- names(cross$geno)
    if (all(regexpr("^[Cc][Hh][Rr]", chrnam) > 0)) {
        chrnam <- substr(chrnam, 4, nchar(chrnam))
        if (all(regexpr("^[Oo][Mm][Oo][Ss][Oo][Mm][Ee]", chrnam) >
            0))
            chrnam <- substr(chrnam, 8, nchar(chrnam))
    }
    if (sum(chrnam == "x") > 0)
        chrnam[chrnam == "x"] <- "X"
    names(cross$geno) <- chrnam
    for (i in 1:length(cross$geno)) if (names(cross$geno)[i] ==
        "X")
        class(cross$geno[[i]]) <- "X"
    chrtype <- sapply(cross$geno, class)
    if (any(chrtype == "X") && convertXdata) {
        if (class(cross)[1] == "bc")
            cross <- fixXgeno.bc(cross)
        if (class(cross)[1] == "f2") {
            if (missing(alleles))
                alleles <- c("A", "B")
            cross <- fixXgeno.f2(cross, alleles)
        }
    }
    cross <- read.cross.bcsft(cross = cross, BC.gen = BC.gen,
        F.gen = F.gen, ...)
    if (estimate.map) {
        cat(" --Estimating genetic map\n")
        map.function <- match.arg(map.function)
        newmap <- est.map(cross, error.prob = error.prob, map.function = map.function)
        cross <- replace.map(cross, newmap)
    }
    for (i in 1:nchr(cross)) storage.mode(cross$geno[[i]]$data) <- "integer"
    if (class(cross)[1] != "4way") {
        if (length(alleles) > 2) {
            warning("length of arg alleles should be 2")
            alleles <- alleles[1:2]
        }
        if (length(alleles) < 2)
            stop("length of arg alleles should be 2")
    }
    else {
        if (missing(alleles))
            alleles <- c("A", "B", "C", "D")
        if (length(alleles) > 4) {
            warning("length of arg alleles should be 4 for a 4-way cross")
            alleles <- alleles[1:4]
        }
        if (length(alleles) < 4)
            stop("length of arg alleles should be 4 for a 4-way cross")
    }
    if (any(nchar(alleles)) != 1) {
        warning("Each item in arg alleles should be a single character")
        alleles <- substr(alleles, 1, 1)
    }
    attr(cross, "alleles") <- alleles
    type <- class(cross)[1]
    if (!missing(crosstype)) {
        if (crosstype == "risib")
            cross <- convert2risib(cross)
        else if (crosstype == "riself")
            cross <- convert2riself(cross)
        else class(cross)[1] <- crosstype
    }
    #summary(cross)
    #cat(" --Cross type:", class(cross)[1], "\n")
    cross
}
parallel.droponemarker <- function (cross, chr, error.prob=0.03, map.function = c("haldane",
    "kosambi", "c-f", "morgan"), m = 0, p = 0, maxit = 2, cores=slurmcore,
    tol = 1e-06, sex.sp = FALSE, verbose = TRUE , parallel=T)
{
    if (!("cross" %in% class(cross)))
        stop("Input must have class \"cross\".")
    if (!missing(chr))
        cross <- subset(cross, chr = chr)
    if (any(nmar(cross) < 3)) {
        if (all(nmar(cross) < 3))
            stop("No chromosomes with at least three markers\n")
        todrop <- names(cross$geno)[nmar(cross) < 3]
        tokeep <- names(cross$geno)[nmar(cross) > 2]
        warning("Dropping chr with <3 markers: ", paste(todrop,
            collapse = ", "))
        cross <- subset(cross, chr = tokeep)
    }


    map.function <- match.arg(map.function)
    if (verbose)
        cat(" -Re-estimating map\n")
    origmap <- qtl:::est.map(cross, error.prob = 0.03, map.function = map.function,
        maxit = maxit, tol = tol, sex.sp = sex.sp,m = 0, p = 0)
    cat(" Done Re-estimating map\n")
    cross <- replace.map(cross, origmap)
    origmaptab <- pull.map(cross, as.table = TRUE)
    origmaptab <- cbind(origmaptab, LOD = rep(NA, nrow(origmaptab)))
    if (is.matrix(origmap[[1]])) {
        origmaptab <- cbind(origmaptab, Ldiff.female = rep(NA,
            nrow(origmaptab)), Ldiff.male = rep(NA, nrow(origmaptab)))
        sex.sp <- TRUE
    } else {
        origmaptab <- cbind(origmaptab, Ldiff = rep(NA, nrow(origmaptab)))
        sex.sp <- FALSE
    }

    for (i in names(cross$geno)) {
        if (sex.sp) {
            Lf <- diff(range(origmap[[i]][1, ]))
            Lm <- diff(range(origmap[[i]][2, ]))
        } else {
         L <- diff(range(origmap[[i]]))
        }

        if (verbose){cat(" -Chromosome", i, "\n")}

        mnames <- markernames(cross, chr = i)
        temp <- subset(cross, chr = i)

        if (parallel) {
              print('starting parallel droponemarker')
              registerDoParallel(cores)
              lod.dif <- foreach(j=seq(along=mnames),
                .inorder=T,.combine='rbind',.packages = "qtl") %dopar% {


                if (verbose > 1) cat(" ---Marker", j, "of", length(mnames), "\n")

                if (sex.sp) {
                  origmaptab[mnames[j], 4] <- -(attr(origmap[[i]],
                    "loglik") - markerll - attr(newmap[[1]], "loglik"))/log(10)
                  origmaptab[mnames[j], 5] <- Lf - diff(range(newmap[[1]][1,
                    ]))
                  origmaptab[mnames[j], 6] <- Lm - diff(range(newmap[[1]][2,
                    ]))
                }

                markerll <- qtl:::markerloglik(cross, mnames[j], error.prob)

                newmap <- qtl:::est.map(drop.markers(temp, mnames[j]),
                  error.prob = 0.03, map.function='kosambi', m=0, p=0,
                  maxit=maxit, tol=tol,sex.sp=FALSE)


                markit <- mnames[j]
                k <- -(attr(origmap[[i]], "loglik") - markerll - attr(newmap[[1]], "loglik"))/log(10)
                Z <- L - diff(range(newmap[[1]]))

                N <- cbind(markit,k,Z)

                return(N)
            }
              rownames(lod.dif) <- lod.dif[,1]
              origmaptab[mnames,'LOD'] <- lod.dif[mnames,2]
              origmaptab[mnames,'Ldiff'] <- lod.dif[mnames,3]

          } else { print('use rqtl if multi cpus not avail')}

      }
      print('done with parallel on all chrs')
      class(origmaptab) <- c("scanone", "data.frame")
      origmaptab$chr <- factor(origmaptab$chr, levels = unique(origmaptab$chr))
      origmaptab
}
environment(read.cross.jm) <- asNamespace('qtl')
environment(parallel.droponemarker) <- asNamespace('qtl')
