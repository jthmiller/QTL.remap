marker.density <- function(cross,gt.1){
  index <- gsub('\\:.*','',rownames(gt.1))==X
  x <- order(as.numeric(gsub('.*\\:','',rownames(gt.1)))[index])
  y <- sort(as.numeric(gsub('.*\\:','',rownames(gt.1)))[index])
  return(list(cor=cor(x,y)^2,pos=y))
}
fix.pheno <- function(cross){
  cross$pheno$ID <- paste("ind", 1:nind(cross), sep="")
  cross$pheno$Pheno_05 <- cross$pheno$Pheno
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
dropone.par <- function(cross,chr, drop.its=1, ...)
  {
  ### Wrapper for parallel.droponemarker... see arguments there ##
  print('starting parallel.droponemarker')
  for (i in 1:drop.its) {
      cross.drops <- parallel.droponemarker(cross,chr, ...)
      drops <- unique(rownames(cross.drops[c(which.max(cross.drops$Ldiff),which.max(cross.drops$LOD)),]))
      cross <- drop.markers(cross,drops)
    }
  return(cross)
}
marker.warning <- function(cross=cross.18){
  print(paste('Starting markers mapped =',
    sum(markernames(cross) %in% markernames(cross.18,chr=X))))

  print(paste('Starting markers un-mapped =',
    sum(!markernames(cross) %in% markernames(cross.18,chr=X))))

}
er.rate <- function(cross,cpus,maxit){
  loglik <- err <- c(0.0025,0.005,0.0075,0.01,0.015,0.02)
      registerDoParallel(cpus)
      hoods <- foreach(i=seq(along=err),.combine=c,
        .inorder=T,.packages = "qtl") %dopar% {
        tempmap <- est.map(cross, error.prob=err[i],maxit)
        return(sum(sapply(tempmap,attr,"loglik")))
      }
      lod <- (hoods - max(hoods))/log(10)
      png(file.path(popdir,paste(X,'_error.png',sep='')))
      plot(err, lod, xlab="Genotyping error rate", xlim=c(0,0.05),
        ylab=expression(paste(log[10], " likelihood")))
      dev.off()
      return(err[which.max(lod)])
}
drop.errlod <- function(cross,cutoff,error.prob){

  print(paste('remove genotypes with erlod >',cutoff))
  mapthis <- calc.errorlod(cross)
  toperr <- top.errorlod(mapthis)
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
        mapthis <<- calc.errorlod(cross)
        toperr2 <<- top.errorlod(mapthis)
          if (!is.null(toperr2)){ step <<- sort(as.numeric(names(table(floor(toperr2$errorlod)))), decreasing=T)
          } else { step <<- 0 }
        }
      )
    }
  }
  print(paste('done...dropped',dropped,'genotypes'))
  return(cross)
}
all.crossed <- function(X,a){
    read.cross(format='csv',file=X,geno=c('AA','AB','BB'),
    alleles=c("A","B"))
}
reconst <- function(X,pop,temp.dir,a, ... ){
  temp <- file.path(basedir,'rQTL',pop,paste('REMAPS/temp.',X,sep=''))
  if (a==1){
    myfiles <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
                format='csvr', geno=c(1:3),estimate.map=FALSE)
    gt.cross.par <- geno.table(myfiles)
    myfiles <- drop.markers(myfiles,rownames(gt.cross.par[gt.cross.par$missing<=miss,]))
    myfiles <- drop.markers(myfiles,rownames(gt.cross.par[gt.cross.par$P.value<cutoff,]))
    pheno <- myfiles[[1]]$pheno$Pheno
    sex <- myfiles[[1]]$pheno$sex
    ID <- myfiles[[1]]$pheno$ID
    #pheno_05 <- myfiles[[1]]$pheno$Pheno

  } else if (!a==1) {
    myfiles <- lapply(temp, function(tocross){
      all.crossed(tocross,a=2)
      }
    )
    ID <- myfiles[[1]]$pheno$ID
    pheno <- myfiles[[1]]$pheno$Pheno
    sex <- myfiles[[1]]$pheno$sex
    gt <- myfiles[[1]]$pheno$gt
  }

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
  chr <- c('','','','',chr)
  map <- c('','','','',map[colnames(cross)])
  cross <- cbind(pheno,sex,ID=as.character(ID),gt,cross)
  cross <- rbind(colnames(cross),chr,map,cross)

  write.table(cross,file=file.path(temp.dir,'tempout'),
      col.names=F,row.names=F,quote=F,sep=',')

  return(read.cross.jm(file=file.path(temp.dir,'tempout'),format='csv',
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
    tol = 1e-06, sex.sp = FALSE, verbose = F , parallel=T)

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
      class(origmaptab) <- c("scanone", "data.frame")
      origmaptab$chr <- factor(origmaptab$chr, levels = unique(origmaptab$chr))
      origmaptab
}
plot.draws <- function (x, chr, reorder = FALSE, main = "Genotype data", alternate.chrid = FALSE,...){
    cross <- x
    if (!any(class(cross) == "cross"))
        stop("Input should have class \"cross\".")
    if (!missing(chr))
        cross <- subset(cross, chr = chr)
    type <- class(cross)[1]
    if (type == "bc" || type == "f2") {
        chrtype <- sapply(cross$geno, class)
        if (any(chrtype == "X")) {
            for (i in which(chrtype == "X")) cross$geno[[i]]$data <- reviseXdata(type,
                "simple", getsex(cross), geno = cross$geno[[i]]$data,
                cross.attr = attributes(cross))
        }
    }
    Geno <- pull.draws(cross)[,,1]
    maxgeno <- max(Geno, na.rm = TRUE)
    if (type != "4way") {
        thecolors <- c("white", "#E41A1C", "#377EB8", "#4DAF4A",
            "#984EA3", "#FF7F00")
        thebreaks <- seq(-0.5, 5.5, by = 1)
    }
    else {
        if (maxgeno <= 5) {
            thecolors <- c("white", "#E41A1C", "#377EB8", "#4DAF4A",
                "#984EA3", "#FF7F00")
            thebreaks <- seq(-0.5, 5.5, by = 1)
        }
        else {
            thecolors <- c("white", "#8DD3C7", "#FFFFB3", "#BEBADA",
                "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
                "#D9D9D9", "#BC80BD")
            thebreaks <- seq(-0.5, 10.5, by = 1)
        }
    }
    thecolors <- thecolors[1:(maxgeno + 1)]
    thebreaks <- thebreaks[1:(maxgeno + 2)]
    o <- 1:nrow(Geno)
    if (reorder) {
        if (is.numeric(reorder)) {
            if (reorder < 1 || reorder > nphe(cross))
                stop("reorder should be TRUE, FALSE, or an integer between 1 and ",
                  nphe(cross))
            o <- order(cross$pheno[, reorder])
        }
        else {
            wh <- sapply(cross$pheno, is.numeric)
            o <- order(apply(cross$pheno[, wh, drop = FALSE],
                1, sum))
        }
    }
    g <- t(Geno[o, ])
    g[is.na(g)] <- 0
    old.xpd <- par("xpd")
    old.las <- par("las")
    par(xpd = TRUE, las = 1)
    on.exit(par(xpd = old.xpd, las = old.las))
    plot_image_sub <- function(g, ylab = "Individuals", xlab = "Markers",
        col = thecolors, ...) {
        if (length(thebreaks) != length(col) + 1)
            stop("Must have one more break than color\n", "length(breaks) = ",
                length(thebreaks), "\nlength(col) = ", length(col))
        image(1:nrow(g), 1:ncol(g), g, col = col, xlab = xlab,
            ylab = ylab, breaks = thebreaks, ...)
    }
    plot_image_sub(g, ...)
    n.mar <- nmar(cross)
    n.chr <- nchr(cross)
    a <- c(0.5, cumsum(n.mar) + 0.5)
    b <- par("usr")
    segments(a, b[3], a, b[4] + diff(b[3:4]) * 0.02)
    abline(h = 0.5 + c(0, ncol(g)), xpd = FALSE)
    a <- par("usr")
    wh <- cumsum(c(0.5, n.mar))
    x <- 1:n.chr
    for (i in 1:n.chr) x[i] <- mean(wh[i + c(0, 1)])
    thechr <- names(cross$geno)
    if (!alternate.chrid || length(thechr) < 2) {
        for (i in seq(along = x)) axis(side = 3, at = x[i], thechr[i],
            tick = FALSE, line = -0.5)
    }
    else {
        odd <- seq(1, length(x), by = 2)
        even <- seq(2, length(x), by = 2)
        for (i in odd) axis(side = 3, at = x[i], labels = thechr[i],
            line = -0.75, tick = FALSE)
        for (i in even) axis(side = 3, at = x[i], labels = thechr[i],
            line = +0, tick = FALSE)
    }
    title(main = main)
    invisible()
}
return.dropped.markers <- function(){
  if(length(mar <- tokeep[!tokeep %in% markernames(cross.18)])>0){
    for (i in tokeep[!tokeep %in% markernames(cross.18)]){
      cross.18 <<- addmarker(cross.18,chr=which.max(nmar(cross.18)),pos=which(tokeep==i),
      markername=i,genotypes=gi[,i])
    }
  }
}
plot.geno <- function(L,gen.main){
  plot(NULL,xlim=c(min(pos),max(pos)),ylim=c(min(pval),max(pval)),main=gen.main)
  points(as.numeric(gsub(paste(X,':',sep=''),'',rownames(L))), log10(L$P.value),pch=20)
}
plot.geno.2 <- function(L,gen.main,wtf=FALSE){
  plot(NULL,xlim=c(min(pos),max(pos)),ylim=c(min(pval),max(pval)),main=gen.main)
  points(as.numeric(gsub(paste(X,':',sep=''),'',rownames(L))), log10(L$P.value),pch=20)
  if(wtf==TRUE){
    ##points(as.numeric(gsub(paste(X,':',sep=''),'',par.gt)),rep(0,length=length(par.gt)), col='green',pch=20)
    points(as.numeric(gsub(paste(X,':',sep=''),'',par.confirm.marks)),log10(L[par.confirm.marks,8]), col='green',pch=20)
  }
}
hist.geno <- function(gt){
  hist(log10(gt),breaks=100,ylim=c(0,300),xlim=c(-20,0))
}
feet <- function(X,Y,Z){
  pdf(file.path(popdir,Y),width = 26, height = 21)
  heatmap.2(X,symm=T,main=NA,notecol="black",labRow=NULL,na.rm=T,cexRow=2,labCol=NULL,key.title=NA,key.xlab=NA,key.ylab=NA,
  cellnote = round(rela,digits=2),notecex=1.2,srtCol=45,cexCol=2,dendrogram="row", margins = c(15,10),trace="none",keysize=.5,col=my_palette,breaks=col_breaks)
  title(Z,cex.main=3, line = -1)
  dev.off()
}
feet2 <- function(X,Y,Z,dir){
  pdf(file.path(dir,Y),width = 30, height = 25)
  heatmap.2(X,symm=T,main=NA,notecol="black",labRow=NULL,na.rm=T,cexRow=1,labCol=NULL,key.title=NA,key.xlab=NA,key.ylab=NA,
  cellnote = round(rela,digits=2),notecex=0.75,srtCol=45,cexCol=1,dendrogram="row", margins = c(15,10),trace="none",keysize=.5,col=my_palette,breaks=col_breaks)
  title(Z,cex.main=3, line = -1)
  dev.off()
}
mdees <- function(X,Y,Z,dirp){
  diag(X) <- 1
  cols <- unlist(sapply(strsplit(rownames(X),'_'),'[[',1))
  labs <- unlist(sapply(strsplit(rownames(X),'_'),'[[',2))

  fit <- cmdscale(as.dist(1-X),eig=TRUE, k=2)
  x <- fit$points[,1]
  y <- fit$points[,2]
  pdf(file.path(dirp,Y),width = 20, height = 20)
    plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
    main="Metric	MDS",	type="n")
  text(x, y, labels = labs,col=brewer.pal(11,"Spectral")[as.factor(cols)], cex=2)
  dev.off()
}
mdees.single <- function(X,Y,Z,dirp){
  diag(X) <- 1
  cols <- unlist(sapply(strsplit(rownames(X),'_'),'[[',2))
  cols <- as.numeric(rownames(X))
  cols[!is.na(cols)] <- 'green'
  cols[is.na(cols)] <- 'red'
  names(cols) <- rownames(X)
  fit <- cmdscale(as.dist(1-X),eig=TRUE, k=2)
  x <- fit$points[,1]
  y <- fit$points[,2]
  pdf(file.path(dirp,Y),width = 20, height = 20)
    plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
    main="Metric	MDS",	type="n")
  text(x, y, labels = row.names(X),col=cols, cex=2)
  dev.off()
}
mdees.single.IBS <- function(X,Y,Z,dirp,dist=T){
  diag(X) <- 1
  pop <- unlist(sapply(strsplit(rownames(X),'_'),'[[',1))
  popcol <- as.numeric(as.factor(pop))
  samp <- unlist(sapply(strsplit(rownames(X),'_'),'[[',2))
  cols <- as.numeric(samp)
  index <- is.na(cols)

  cols <- brewer.pal(8,"Dark2")[popcol]
  names(cols) <- rownames(X)
  ofs <- rownames(X)[!index]
  parso <- rownames(X)[index]
  if (dist==F){
  fit <- cmdscale(as.dist(1-X),eig=TRUE, k=2)
  } else { fit <- cmdscale(as.dist(X),eig=TRUE, k=2)}
  x <- fit$points[ofs,1]
  y <- fit$points[ofs,2]
  A <- fit$points[parso,1]
  B <- fit$points[parso,2]
  raw <- range(fit$points)
  pdf(file.path(dirp,Y),width = 20, height = 20)

    plot(raw, raw, xlab="Coordinate 1", ylab="Coordinate 2",
    main="Metric	MDS",type="n" )

  points(x, y,pch=19,col=cols[names(x)],cex=3)
  text(A, B, labels = names(A),col=cols[names(A)], cex=2)
  legend(0.06, 0.04, legend=unique(pop),pch=19,
       col=brewer.pal(8,"Dark2")[unique(popcol)],cex=2)
  dev.off()
}
newt <- function(X,Y,Z,dir){
  cols <- unlist(sapply(strsplit(rownames(X),'_'),'[[',1))
  labs <- unlist(sapply(strsplit(rownames(X),'_'),'[[',2))
  names(cols) <- rownames(X)
  pdf(file.path(dir,Y),width = 20, height = 20)
  qgraph(X, layout='spring', vsize=3,label.cex=1,groups=cols,labels=labs)
  dev.off()
}
rels <- function(X){
  bela <- comparegeno(X)
  #colnames(bela) <- gsub(paste(pop,'_',sep=''),'',X$pheno$ID)
  #rownames(bela) <- gsub(paste(pop,'_',sep=''),'',X$pheno$ID)
  colnames(bela) <- X$pheno$ID
  rownames(bela) <- X$pheno$ID

  bela[bela==NaN] <- NA
  diag(bela) <- NA
  bela <- bela[rowSums(is.na(bela)) < nind(X),colSums(is.na(bela)) < nind(X)]
  return(bela)
}
mega.cross <- function(pope){
  pop <- pope
  basedir <- '/home/jmiller1/QTL_Map_Raw/popgen'
  indpops <- file.path(basedir,'plinkfiles/ind.pops')

  cross <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
    format='csvr', geno=c(1:3),estimate.map=FALSE)
  path <- file.path(indpops,paste(pop,'.ped',sep=''))
  popname <- system(paste('cut -f1 -d\' \'',path),intern = TRUE)
  indname <- system(paste('cut -f2 -d\' \'',path),intern = TRUE)
  cross$pheno$ID <- paste(popname,indname,sep='_')
  cross <- subset(cross, ind=nmissing(cross)<nmissing(cross)[1:(round(sum(nind(cross))/2,digits=0))] | is.na(cross$pheno$Pheno))
  return(cross)
}
get.par.snp.number <- function(X){
  sapply(1:24,function(X){
    length(which(gt.pars$chr==X & gt.pars$AA==1 & gt.pars$BB==1))/length(which(gt.pars$chr==X & gt.pars$missing==0 & gt.pars$AA<2 & gt.pars$BB<2))
    }
  )
}
repRipple.jm<-function(cross, chr = NULL, window = 5,method = "countxo", verbose = T,re.est.map=F,
                    map.function = "kosambi", sex.sp=F, clean1st = FALSE, ripVerb = TRUE, ...){
  loadNamespace("qtl")
  if(clean1st) cross<-clean(cross)
  if(is.null(chr)){
    chr<-chrnames(cross)
  }
  for(j in chr){
    if(length(markernames(cross, chr = j))>=3){
      if(verbose) cat(j,"...")
      new.xo<-0
      orig.xo<-1
      while(new.xo<orig.xo){
        mars<-lapply(chrnames(cross), function(x) markernames(cross,chr = x))
        names(mars)<-chrnames(cross)
        verb<-ifelse(new.xo == 0 & ripVerb , TRUE, FALSE)

        s<-summary(ripple(cross, chr = j, window = window,
                          method = method, verbose = verb))
        best<-as.numeric(which.min(s[,ncol(s)])[1])
        ord<-s[best,-ncol(s)]
        orig.xo<-s[1,ncol(s)]
        new.xo<- s[best,ncol(s)]
        if(verbose){
          if(orig.xo == new.xo){
            cat("no reduction in XOs found\n")
          }else{
            cat("orig n XO = ", orig.xo, "new n XO = ",new.xo,"\n")
          }
        }
        mars[[as.character(j)]]<-mars[[as.character(j)]][ord]
        cross<-qtlTools::newLG(cross = cross, markerList = mars)
        names(cross$geno) <- j
      }
    }
  }

  if(verbose & re.est.map==T){ cat("final map estimation")
  map<-est.map(cross, map.function = map.function, sex.sp=sex.sp,...)
  cross<-replace.map(cross, map)
  }
  return(cross)
}


qb.BayesFactor.jm <- function(qbObject,
                           items = c("nqtl","pattern","chrom","pairs"),
                           cutoff.pattern = ifelse(epistasis, 0.25, 0.5),
                           cutoff.pairs = 1, nmax = 15,
                           epistasis = TRUE, ...)
{
  qtlbim:::qb.exists(qbObject)

  assess <- list()

  if(any(pmatch(tolower(items), "nqtl", nomatch = 0)))
    assess$nqtl <- qtlbim:::qb.numqtl(qbObject, ...)

  if(any(pmatch(tolower(items), "pattern", nomatch = 0)))
    assess$pattern <- qb.pattern.jm(qbObject, cutoff.pattern, nmax, epistasis, ...)

  if(any(pmatch(tolower(items), "chrom", nomatch = 0)))
    assess$chrom <- qtlbim:::qb.chrom(qbObject, ...)

  if(any(pmatch(tolower(items), "pairs", nomatch = 0)))
    assess$pairs <- qtlbim:::qb.pairs(qbObject, cutoff.pairs, nmax, ...)

  class(assess) <- c("qb.BayesFactor", "list")
  assess
}



###
qb.pattern.jm <- function(qbObject, cutoff = 1, nmax = 15, epistasis = TRUE, ...)
{
  iterdiag <- qtlbim:::qb.get(qbObject, "iterdiag", ...)
  n.iter <- qtlbim:::qb.niter(qbObject)

  mainloci <- qtlbim:::qb.get(qbObject, "mainloci", ...)
  pairloci <- qtlbim:::qb.get(qbObject, "pairloci", ...)
  pattern <- qtlbim:::qb.makepattern(qbObject, epistasis, iterdiag = iterdiag,
                            mainloci = mainloci, pairloci = pairloci)

  posterior <- rev(sort(table(pattern)))
  posterior <- posterior / sum(posterior)
  tmp <- posterior >= cutoff / 100
  if(sum(tmp))
    posterior <- posterior[tmp]
  else {
    cat("warning: pattern posterior cutoff", cutoff,
        "is too large and is ignored\n",
        "posterior range is", range(posterior), "\n")
  }
  if(length(posterior) > nmax)
    posterior <- posterior[1:nmax]
  ucount <- match(names(posterior), pattern)

  ## prior for pattern
  rng <- max(iterdiag$nqtl)
  pr <- qtlbim:::qb.prior(qbObject, 0:rng)
  bf <- posterior
  map <- pull.map(qtlbim:::qb.cross(qbObject, genoprob = FALSE))
  chrlen <- unlist(lapply(map, max))
  nchrom <- length(chrlen)
  chrlen <- chrlen / sum(chrlen)

  fact <- rep(1, rng)
  for(i in 2:(rng+1))
    fact[i] <- fact[i-1] * i

  ## New plan. Use only subset of mainloci matching ucount's.
  ## bundle table and for loop into one.

  ## Set up prior using null.
  prior <- rep(pr[1], length(posterior))
  names(prior) <- names(posterior)

  ## Find prior proportional to qb.prior and lengths of chromosomes.
  ## Use factorial adjustments for multiple linked QTL.
  tmpfn <- function(x, chrlen, fact) {
    ct <- c(table(x))
    prod(chrlen[names(ct)] ^ ct) * fact[sum(ct)] / prod(fact[ct])
  }

  ## Subset on non-null iterations corresponding to posterior.
  sub.post <- iterdiag$niter[ucount]
  is.depen <- qtlbim:::qb.get(qbObject, "depen")
  if(epistasis & is.depen)
    tmp <- !is.na(match(pairloci$niter, sub.post))
  sub.post <- !is.na(match(mainloci$niter, sub.post))

  ## Adjust for epistatic effects.
  ## For now only considering hierarchical model case.
  ## That is epistasis only if main effects.
  ## In order to handle epistasis with 1 or no main effect
  ## we would need to compute probabilities for each iteration.
  ## This will be a lot more work!
  if(epistasis) {
    if(is.depen) {
      ## Epistatic priors depend on number of main loci.
      ## Not correct yet when c1 or c0 > 0.
      ## Find number of possible epistatic pairs of main loci.
      tbl.main <- c(table(mainloci[sub.post, "niter"]))
      tmp2 <- tbl.main * (tbl.main - 1) / 2

      tmp <- c(table(pairloci[tmp,"niter"]))
      prop <- qtlbim:::qb.get(qbObject, "prop")
      if(sum(prop[-1]) > 0)
        warning("Bayes factor computations assume all QTL are main (may not be correct)")

      ## Epistasis with two main QTL.
      tmp2[names(tmp)] <- (prop[1] ^ tmp) *
        ((1 - prop[1]) ^ (tmp2[names(tmp)] - tmp))
      tmp2[is.na(match(names(tmp2), names(tmp)))] <-
        (1 - prop[1]) ^ tmp2[is.na(match(names(tmp2), names(tmp)))]
      tbl.main <- pr[1 + tbl.main] * tmp2
    }
    else {
      ## Independent prior for epistasis.
      ## Product of prior for main QTL * prior for epis-only QTL.

      ## Find mainloci entries matching patterns with high posterior.
      ## Kludge to catch iterations with 0 QTL.
      tmp <- rep(0, nrow(iterdiag))
      names(tmp) <- iterdiag$niter
      tmp2 <- c(table(mainloci[, "niter"]))
      tmp[names(tmp2)] <- tmp2
      tmp2 <- rep(!is.na(match(pattern, names(posterior))), tmp)

      ## Sum up main loci variance components and tally those not zero.
      tmp <- mainloci[tmp2, "varadd"]
      is.bc <- (qtlbim:::qb.cross.class(qbObject) == "bc")
      if(!is.bc)
        tmp <- tmp + mainloci[tmp2, "vardom"]
      tmp <- tapply(tmp, mainloci[tmp2, "niter"], function(x) sum(x > 0))
      ## Table number of total QTL per iteration.
      tbl.main <- c(table(mainloci[tmp2, "niter"]))

      ## Product of priors for main and epis-only QTL.
      main.nqtl <- qtlbim:::qb.get(qbObject, "main.nqtl")
      mean.nqtl <- qtlbim:::qb.get(qbObject, "mean.nqtl")
      ## Average product of poisson probabilities by pattern.
      pr.main <- qtlbim:::qb.prior(qbObject, seq(0, max(tmp)), main.nqtl)
      pr <- qtlbim:::qb.prior(qbObject, seq(0, max(tbl.main - tmp)), mean.nqtl - main.nqtl)
      tbl.main <- tapply(pr.main[1 + tmp] * pr[1 + tbl.main - tmp],
                         pattern[names(tbl.main)], mean)
      ## Rearrange in order to match calculations below.
      tbl.main <- tbl.main[pattern[as.character(sort(iterdiag$niter[ucount]))]]
      names(tbl.main) <- sort(iterdiag$niter[ucount])
      tbl.main <- tbl.main[!is.na(tbl.main)]
    }
  }
  else
    tbl.main <- pr[1 + c(table(mainloci[sub.post, "niter"]))]

  geno.names <- names(map)
  tmp <- c(tapply(ordered(geno.names[mainloci[sub.post, "chrom"]], geno.names),
                  mainloci[sub.post, "niter"],
                  tmpfn, chrlen, fact))
  tmp <- tmp * tbl.main

  prior[pattern[names(tmp)]] <- tmp

  ## Bayes factor ratio = rescaled version of posterior / prior.
  bf <- posterior / prior
  ## rescale bf so smallest value is 1 (avoid NA, 0/0)
  minbf <- bf[!is.na(bf)]
  if(length(minbf))
    bf <- bf / min(minbf)

  ## bfse = approximate Monte Carlo standard error for bf
  ## (actually binomial error)
  ## note that this is rescaled since bf[1] forced to be 1
  tmp <- names(posterior)
  nqtl <- sapply(strsplit(tmp, ","),length)
  names(nqtl) <- tmp
  bf <- as.data.frame(bf)$Freq
  names(bf) <- tmp
  posterior  <- as.data.frame(posterior )$Freq
  names(posterior ) <- tmp

  data.frame(nqtl = nqtl,
             posterior = posterior, prior = prior, bf = bf,
             bfse = sqrt((1 - posterior) / (posterior * n.iter)) * bf)
}
cnv.popgen <- function(cross2, popgen, top) {
  nbhmap <- convert2cross2(cross2)
  nbhmap$pmap <- nbhmap$gmap
  nbhmap$pmap <- lapply(nbhmap$pmap, function(X) {
    return(as.numeric(gsub("[0-9]+:", "", names(X))))
  })
  for (i in 1:24) {
    names(nbhmap$pmap[[i]]) <- names(nbhmap$gmap[[i]])
  }

  popgen$chrom <- gsub("chr", "", popgen$chrom)
  colnames(popgen)[1] <- "chr"
  pop.list <- split(popgen, popgen$chr)
  pop.list <- pop.list[as.character(c(1:24))]
  popgen.pos1 <- lapply(pop.list, "[[", 2)
  popgen.pos2 <- lapply(pop.list, "[[", 3)
  popgen.pos1 <- lapply(popgen.pos1, as.numeric)
  popgen.pos2 <- lapply(popgen.pos2, as.numeric)
  pos1 <- interp_map(popgen.pos1, nbhmap$pmap, nbhmap$gmap)
  pos2 <- interp_map(popgen.pos2, nbhmap$pmap, nbhmap$gmap)
  pop.list <- mapply(cbind, pop.list, pos1 = pos1, SIMPLIFY = FALSE)
  pop.list <- mapply(cbind, pop.list, pos = pos2, SIMPLIFY = FALSE)
  fraps <- ldply(pop.list, data.frame, .id = NULL)
  fraps$chr <- factor(fraps$chr, levels = as.character(1:24))
  fraps <- fraps[which(fraps$rank < top), ]
  fraps$lod <- rescale(fraps$rank, to = c(20, 0))
  fraps$sz <- rescale(fraps$rank, to = c(4, 1))
  return(fraps)
}
cnv.premap <- function(cross2, tomap) {
  nbhmap <- convert2cross2(cross2)
  nbhmap$pmap <- nbhmap$gmap
  nbhmap$pmap <- lapply(nbhmap$pmap, function(X) {
    return(as.numeric(gsub("[0-9]+:", "", names(X))))
  })
  for (i in 1:24) {
    names(nbhmap$pmap[[i]]) <- names(nbhmap$gmap[[i]])
  }
  map <- map2table(pull.map(tomap))
  map$pos <- gsub("^[[:digit:]]+:", "", rownames(map))
  map.list <- split(map, map$chr)
  map.list <- map.list[as.character(c(1:24))]
  map.pos <- lapply(map.list, "[[", 2)
  map.pos <- lapply(map.pos, as.numeric)
  pos <- interp_map(map.pos, nbhmap$pmap, nbhmap$gmap)
  map.list <- mapply(cbind, map.list, pos = pos, SIMPLIFY = FALSE)
  fraps <- ldply(map.list, data.frame, .id = NULL)
  fraps$chr <- factor(fraps$chr, levels = as.character(1:24))
  rownames(fraps) <- rownames(map)
  for (i in 1:24) {
    ord <- order(fraps$pos.1[which(fraps$chr == i)])
    tomap <- switch.order(tomap, chr = i, ord, error.prob = 0.1, map.function = "kosambi",
      maxit = 1, tol = 0.1, sex.sp = F)
  }
  fraps$pos <- fraps$pos.1

  fraps <- fraps[, -3]
  rownames(fraps) <- rownames(map2table(pull.map(tomap)))
  map <- table2map(fraps)
  brp.newmap <- replacemap(tomap, map)
  brp.newmap

}
cnv.ahrs <- function(cross2, AHRdf, EXP) {
  mapo <- convert2cross2(cross2)
  mapo$pmap <- mapo$gmap

  mapo$pmap <- lapply(mapo$pmap, function(X) {
    return(as.numeric(gsub("[0-9]+:", "", names(X))))
  })
  for (i in 1:24) {
    names(mapo$pmap[[i]]) <- names(mapo$gmap[[i]])
  }
  AHR.list <- split(AHRdf, AHRdf$chrom)
  AHR.pos1 <- lapply(AHR.list, "[[", 2)
  AHR.pos2 <- lapply(AHR.list, "[[", 3)
  AHR.pos1 <- lapply(AHR.pos1, as.numeric)
  AHR.pos2 <- lapply(AHR.pos2, as.numeric)
  pos1 <- interp_map(AHR.pos1, mapo$pmap, mapo$gmap)
  pos2 <- interp_map(AHR.pos2, mapo$pmap, mapo$gmap)
  AHR.list <- mapply(cbind, AHR.list, pos1 = pos1, SIMPLIFY = FALSE)
  AHR.list <- mapply(cbind, AHR.list, pos = pos2, SIMPLIFY = FALSE)
  ah.gens <- ldply(AHR.list, data.frame, .id = NULL)
  ah.gens$chrom <- factor(ah.gens$chrom, levels = as.character(1:24))
  colnames(ah.gens)[1] <- "chr"
  # ah.gens$lod <- rep_len(c(-1:-10), length(ah.gens[,1]))
  ah.gens$lod <- 0
  if (EXP == F) {
    ah.gens <- ah.gens[which(!ah.gens$gene == "EXPRESSED"), ]
  }
  ah.gens <- ah.gens[!grepl("*many*", ah.gens$gene), ]
  ah.gens <- ah.gens[!grepl("*many*", ah.gens$gene), ]
  return(ah.gens)
}
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}
getone <- function(mar1,mar2,phenocol1,phenocol2){

     return(
       data.frame(po1 = am[,mar1],
          po2 = am[,mar2],
          pt1=pt[,phenocol1],
          pt2=pt[,phenocol2])
        )
}
getem <- function(mar1,mar2,phenocol1,phenocol2){

     return(
       data.frame(po1 = am[,mar1],
          po2 = am[,mar2],
          tb=paste(tra[am[,mar1]],tra[am[,mar2]],sep=''),
          pt1=pt[,phenocol1],
          pt2=pt[,phenocol2])
        )
}
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}
cleanGeno_jm <-
function (cross, chr, maxdist = 2.5, maxmark = 2, verbose = TRUE)
{
    if (!missing(chr))
        cleaned <- subset(cross, chr = chr)
    else cleaned <- cross
    thechr <- names(cleaned$geno)
    totdrop <- 0
    maxmaxdist <- max(maxdist)
    for (i in thechr) {
        xoloc <- locateXO(cleaned, chr = i, full.info = TRUE)
        nxo <- sapply(xoloc, function(a) if (is.matrix(a))
            return(nrow(a))
        else return(0))
        g <- pull.geno(cleaned, chr = i)
        ndrop <- 0
        for (j in which(nxo > 1)) {

            maxd <- xoloc[[j]][-1, "right"] - xoloc[[j]][-nrow(xoloc[[j]]),
                "left"]
            wh <- maxd <= maxmaxdist
            if (any(wh)) {
                for (k in which(wh)) {
                  nt <- sum(!is.na(g[j, (xoloc[[j]][k, "ileft"] +
                    1):(xoloc[[j]][k + 1, "iright"] - 1)]))


                  if (nt > 0 && any(nt <= maxmark & maxd[k] <
                    maxdist)) {
                      mks <- cleaned$geno[[i]]$data[ j , (xoloc[[j]][k,"ileft"] + 1) : (xoloc[[j]][k + 1, "iright"] - 1) ]

                      if (any(mks==1|3)){
                       cleaned$geno[[i]]$data[j, (xoloc[[j]][k,
                       "ileft"] + 1):(xoloc[[j]][k + 1, "iright"] -
                       1)] <- ifelse(mks==1|3,2,NA)
                      }else{
                       cleaned$geno[[i]]$data[j, (xoloc[[j]][k,
                       "ileft"] + 1):(xoloc[[j]][k + 1, "iright"] -
                       1)] <- NA
                      }
                    ndrop <- ndrop + nt
                    totdrop <- totdrop + nt
                  }
                }
            }
        }
        if (verbose && ndrop > 0) {
            totgen <- sum(ntyped(subset(cross, chr = i)))
            cat(" ---Dropping ", ndrop, " genotypes (out of ",
                totgen, ") on chr ", i, "\n", sep = "")
        }
    }
    if (verbose && nchr(cleaned) > 1 && totdrop > 0) {
        totgen <- sum(ntyped(subset(cross, chr = thechr)))
        cat(" ---Dropped ", totdrop, " genotypes (out of ", totgen,
            ") in total\n", sep = "")
    }
    for (i in names(cleaned$geno)) cross$geno[[i]] <- cleaned$geno[[i]]
    cross
}

cleanGeno_jm_2 <-
function (cross, chr, maxdist = 2.5, maxmark = 2, verbose = TRUE)
{
    if (!missing(chr))
        cleaned <- subset(cross, chr = chr)
    else cleaned <- cross
    thechr <- names(cleaned$geno)
    totdrop <- 0
    maxmaxdist <- max(maxdist)
    for (i in thechr) {
        xoloc <- locateXO(cleaned, chr = i, full.info = TRUE)
        nxo <- sapply(xoloc, function(a) if (is.matrix(a))
            return(nrow(a))
        else return(0))
        g <- pull.geno(cleaned, chr = i)
        ndrop <- 0
        for (j in which(nxo > 1)) {

            maxd <- xoloc[[j]][-1, "right"] - xoloc[[j]][-nrow(xoloc[[j]]),
                "left"]
            wh <- maxd <= maxmaxdist
            if (any(wh)) {
                for (k in which(wh)) {
                  nt <- sum(!is.na(g[j, (xoloc[[j]][k, "ileft"] +
                    1):(xoloc[[j]][k + 1, "iright"] - 1)]))


                  if (nt > 0 && any(nt <= maxmark & maxd[k] <
                    maxdist)) {

                       cleaned$geno[[i]]$data[j, (xoloc[[j]][k,
                       "ileft"] + 1):(xoloc[[j]][k + 1, "iright"] -
                       1)] <- NA
                    ndrop <- ndrop + nt
                    totdrop <- totdrop + nt
                  }
                }
            }
        }
        if (verbose && ndrop > 0) {
            totgen <- sum(ntyped(subset(cross, chr = i)))
            cat(" ---Dropping ", ndrop, " genotypes (out of ",
                totgen, ") on chr ", i, "\n", sep = "")
        }
    }
    if (verbose && nchr(cleaned) > 1 && totdrop > 0) {
        totgen <- sum(ntyped(subset(cross, chr = thechr)))
        cat(" ---Dropped ", totdrop, " genotypes (out of ", totgen,
            ") in total\n", sep = "")
    }
    for (i in names(cleaned$geno)) cross$geno[[i]] <- cleaned$geno[[i]]
    cross
}
environment(plot.draws) <- asNamespace('qtl')
environment(read.cross.jm) <- asNamespace('qtl')
environment(parallel.droponemarker) <- asNamespace('qtl')
assignInNamespace
