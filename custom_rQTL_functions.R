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


parallel.droponemarker <- function (cross, chr, error.prob = 1e-04, map.function = c("haldane",
    "kosambi", "c-f", "morgan"), m = 0, p = 0, maxit = 5, cores=slurmcore
    tol = 1e-06, sex.sp = TRUE, verbose = TRUE,parallel=T )
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
    origmap <- est.map(cross, error.prob = error.prob, map.function = map.function,
        m = m, p = p, maxit = maxit, tol = tol, sex.sp = sex.sp)
    cross <- replace.map(cross, origmap)
    origmaptab <- pull.map(cross, as.table = TRUE)
    origmaptab <- cbind(origmaptab, LOD = rep(NA, nrow(origmaptab)))
    if (is.matrix(origmap[[1]])) {
        origmaptab <- cbind(origmaptab, Ldiff.female = rep(NA,
            nrow(origmaptab)), Ldiff.male = rep(NA, nrow(origmaptab)))
        sexsp <- TRUE
    }
    else {
        origmaptab <- cbind(origmaptab, Ldiff = rep(NA, nrow(origmaptab)))
        sexsp <- FALSE
    }
### Start of par.drop.markers

    for (i in names(cross$geno)) {
        if (sexsp) {
            Lf <- diff(range(origmap[[i]][1, ]))
            Lm <- diff(range(origmap[[i]][2, ]))
        }
        else L <- diff(range(origmap[[i]]))
        if (verbose)
            cat(" -Chromosome", i, "\n")
        mnames <- markernames(cross, chr = i)
        temp <- subset(cross, chr = i)
        for (j in seq(along = mnames)) {
            if (verbose > 1)
                cat(" ---Marker", j, "of", length(mnames), "\n")
            markerll <- qtl:::markerloglik(cross, mnames[j], error.prob)
            newmap <- est.map(drop.markers(temp, mnames[j]),
                error.prob = error.prob, map.function = map.function,
                m = m, p = p, maxit = maxit, tol = tol, sex.sp = sex.sp)
            if (sexsp) {
                origmaptab[mnames[j], 4] <- -(attr(origmap[[i]],
                  "loglik") - markerll - attr(newmap[[1]], "loglik"))/log(10)
                origmaptab[mnames[j], 5] <- Lf - diff(range(newmap[[1]][1,
                  ]))
                origmaptab[mnames[j], 6] <- Lm - diff(range(newmap[[1]][2,
                  ]))
            }

            if (parallel) {
                #cores=detectCores()
                #cl <- makeCluster(cores[1])
                #registerDoParallel(cl)
                registerDoParallel(slurmcore)
                lod.dif <- foreach(j=seq(along=mnames),
                  .inorder=T,.combine='rbind',.packages = "qtl") %dopar% {

                    markerll <- qtl:::markerloglik(cross, mnames[j], error.prob)

                    newmap <- qtl:::est.map(qtl:::drop.markers(temp, mnames[j]), error.prob=error.prob,
                      map.function=map.function, m=m, p=p, maxit=maxit, tol=tol,
                      sex.sp=sex.sp)

                    N <- cbind(mnames[j],
                      (attr(origmap[[i]], "loglik") - markerll - attr(newmap[[1]], "loglik"))/log(10),
                      L - diff(range(newmap[[1]])))

                    return(N)
                  }
                names(lod.dif) <- lod.dif[,1]
                origmaptab[mnames,'LOD'] <- lod.dif[,2]
                origmaptab[mnames,'Ldiff'] <- lod.dif[,3]
            } else {
                origmaptab[mnames[j], 3] <- -(attr(origmap[[i]],
                  "loglik") - markerll - attr(newmap[[1]], "loglik"))/log(10)
                origmaptab[mnames[j], 4] <- L - diff(range(newmap[[1]]))
            }
        }
    }
    class(origmaptab) <- c("scanone", "data.frame")
    origmaptab$chr <- factor(origmaptab$chr, levels = unique(origmaptab$chr))
    origmaptab
}
environment(read.cross.jm) <- asNamespace('qtl')
environment(parallel.droponemarker) <- asNamespace('qtl')
