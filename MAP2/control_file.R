## Directories
basedir <- "/home/jmiller1/QTL_Map_Raw/popgen"
plotdir <- file.path(basedir, "rQTL/plots")
indpops <- file.path(basedir, "plinkfiles/ind.pops")
popdir <- file.path(basedir, "rQTL", pop, "REMAPS")
qtldir <- file.path(basedir, "rQTL/remap_out")
errfile <- file.path(qtldir, "genotyping_error_rate.txt")
dirso <- "/home/jmiller1/QTL_Map_Raw/popgen/rQTL/data/"

## Funtions for processing rQTL map data
source(file.path(basedir, "rQTL/scripts/QTL_remap/MAP/source_file.R"))
# source(file.path(basedir, 'rQTL/scripts/QTL_remap/QTL/model_source_file.R'))

## Libraries
flib <- "/share/apps/rmodules"
fpacks <- c("devtools", "httr", "RColorBrewer", "qtl")
lapply(fpacks, require, character.only = TRUE, lib.loc = flib)

mylib <- "/home/jmiller1/R/x86_64-pc-linux-gnu-library/3.5"
mpacks <- c("qtl", "foreach", "qtl2", "qtlTools", "doParallel", "plyr")
lapply(mpacks, require, character.only = TRUE, lib.loc = mylib)

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if (trace)
      cat(nm, ":")
    source(file.path(path, nm), ...)
    if (trace)
      cat("\n")
  }
}

sourceDir("doParallel/R")

dropByDropone<-function(cross, droponeRes,
                        endMarkerThresh = 12, midMarkerThresh = 4,
                        which.map = NULL,
                        map.function = "kosambi", re.est.map = F, sex.sp=F,
                        ...){

  if(class(cross)[1] == "4way"){
    if(is.null(which.map)){
      which.map<-"Ldiff.female"
    }
    ldCol = which(colnames(droponeRes)==which.map)
    map<-lapply(pull.map(cross), colnames)
  }else{
    ldCol = which(colnames(droponeRes)=="Ldiff")
    map<-lapply(pull.map(cross), names)
  }

  endMars<-unlist(lapply(map, function(x) c(x[1],x[length(x)])))

  todrop<-rownames(
    droponeRes[droponeRes[,ldCol]>midMarkerThresh,])
  todrop<-todrop[!todrop %in% endMars]
  tail2drop<-rownames(
    droponeRes[rownames(droponeRes) %in% endMars &
                 droponeRes[,ldCol]>endMarkerThresh,])

  cross<-drop.markers(cross, c(todrop,tail2drop))

  if(re.est.map){
    map<-est.map(cross, map.function = "kosambi", sex.sp=F, ...)
    cross<-replace.map(cross, map)
  }

  return(cross)
}

refine_maps <- function(X){

 tmp <- orderMarkers(X, window=4,verbose=FALSE,
                 use.ripple=FALSE, error.prob=0.01, sex.sp=FALSE,
                 map.function="kosambi",maxit=1000, tol=1e-3)

 tmp <- calc.errorlod(tmp, err=0.01)

 tmp_map <-  est.map(tmp, error.prob=0.01,
              map.function="kosambi",
              maxit=1000, tol=1e-4, sex.sp=FALSE,
              verbose=FALSE, n.cluster=6)

 tmp <- qtl:::replace.map(tmp,tmp_map)

 drp1 <- droponemarker(tmp, error.prob=0.01,
                    map.function="kosambi",
                    maxit=1000, tol=1e-3, sex.sp=FALSE,
                    verbose=FALSE)

 dropByDropone(tmp, drp1, endMarkerThresh = 20,
  midMarkerThresh = 20, map.function = "kosambi",
  re.est.map = T, error.prob=0.01,maxit=1000, tol=1e-4, sex.sp=FALSE,
  verbose=FALSE, n.cluster=6)

}
