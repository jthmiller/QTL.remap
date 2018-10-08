## Data informations

# This is a data using selective genotyping in an F2 population. There are four 150-cM
# chromosomes, each has 11 equally spaced markers. Only consider additive effects in the QTL
# model.

## Input Data

# Data format:
# 1. Save the data of trait values and marker genotypes as a example.csv file with column names.
#    The first column is the trait values, others are corresponding marker genotypes.
#    Each row shows the data of an individual.
# 2. Marker genotype: 2 for AA, 1 for Aa, 0 for aa ,and NA for ungenotyped individuals.

source('~/Dropbox/QTL_Paper/CODE/QTL/local_debug.R')

#data<-read.csv("~/Dropbox/QTL_Paper/Lee_et_all/example.csv")
cross.18 <- read.cross(file=file.path('~/Dropbox/QTL_Paper/DATA/tempout.csv'),format='csv',
  geno=c(1:3),estimate.map=FALSE)

gts <- c(2,1,0)
names(gts) <- c('1','2','3')
tr.geno <- function(gtvec){
  gts[as.character(gtvec)]
  }

try <- sapply(1:24,function(Z){
  cross.18$geno[[as.character(Z)]]$"data" <<- apply(cross.18$geno[[as.character(Z)]]$"data",2,tr.geno)
  }
)

### Set a phenotype column for strata to indicate genotyped individuals
cross.18$pheno$stata <- gsub('[0-9]','',cross.18$pheno$ID)

## Give the chromosome number.
# n.chrom: number of chromosomes
##n.chrom<-4
n.chrom <- nchr(cross.18)

## Give the marker number in chromosomes.
# info.chrom: a vector of number of markers on each chromosome
###info.chrom<-c(11,11,11,11)
info.chrom <- nmar(cross.18)
data <- unlist(do.call(cbind, lapply(cross.18$geno,'[[',1)))

## Trait values
y <- cross.18$pheno$pheno_05
data <- cbind(y,data) ## Not nec. but keeping with Lee structure
ngy.ind<-is.na(data[,ncol(data)])
ys<-y[ngy.ind==0]
yu<-y[ngy.ind==1]

## Give the length of marker intervals on each chromosome (in Morgan).
## marker.map[[i]]: the length of marker intervals in the ith chromosome

#marker.map<-vector("list",n.chrom)
#marker.map[[1]]<-c(rep(0.15,10))
#marker.map[[2]]<-c(rep(0.15,10))
#marker.map[[3]]<-c(rep(0.15,10))
#marker.map[[4]]<-c(rep(0.15,10))

marker.map <- vector("list",24)
marker.map <- lapply(1:n.chrom,function(Z){
    a <- cross.18$geno[[as.character(Z)]]$"map"
    b <- length(cross.18$geno[[as.character(Z)]]$"map")
    d <- c(0,cross.18$geno[[as.character(Z)]]$"map")[-b]
    return(a-d)
  }
)

## Marker genotypes on each chromosome.
# marker.gty[[i]]: marker genotypes in the ith chromosome


mar.gty<-data[ngy.ind==0,-1]
marker.gty <- vector("list",n.chrom)
marker.gty[[1]] <- mar.gty[,1:info.chrom[1]]
for(infc in 2: n.chrom){
  marker.gty[[infc]]<-mar.gty[,(sum(info.chrom[1:(infc-1)])+1):sum(info.chrom[1:infc])]
  }

### Translate marker.gty to correct gt code

## Specify the chromosome region you would like to search for QTLs and the QTLs you would like to
## condition in in the model for the serach.
# For example: If the detected QTLs at (1, 3, 0.01), (1, 5, 0.12) and (4, 9, 0) are used for
#              conditioning in the model in the search of chromosomes 2 and 3 (one QTL at a time)
#              for QTLs. Here (1, 3, 0.01) indicates the position is at 0.01 Morgan away from the
#              left marker of the third marker interval in the first chromosome.
## aka
## chr, marker_number, morgan_away
## (1, 3, 0.01)

## Give the QTL number.

nQTL<-2

## Genetic design matrix.
D.matrix<-Design(nQTL,domi=F,epis=F)

## Give the positions of detected QTLs.
# Note that the order is sort according to chromosome and marker.

#posi<-vector("list",nQTL)
#posi[[1]]<-c(1,3,0.01)
#posi[[2]]<-c(1,5,0.12)

posi<-vector("list",nQTL)
posi[[1]]<-c(2,280,1.395633)
posi[[2]]<-c(18,25,5.834301)

## Give the walking speed along chromsome in Morgan.

sp<-0.01

## Give the search informations.
# sear.chrom: the information for chromosome
# sear.mar: the information for marker interval



#sear.chrom <- c(1,3:17,19:24)
sear.mar <- markernames(cross.18,chr=8)
sear.chrom <- c(8)
## Give model type
# model: model=1 indicates using the proposed method in Lee et al. (2014);
#        model=2 indicates using the truncated model

model<-1
# If model=2, users need to give the values of truncation points (tL and tR).

wi<-0
n.search<-2*10*16

result<-vector("list",n.search)







  for(sear.chrom in 8){
    for(sear.mar in 1:(info.chrom[sear.chrom]-1)){
      mar.len<-marker.map[[c(sear.chrom,sear.mar)]]
      mov.seq<-seq(0,mar.len,sp)
      for(sear.int in 1:length(mov.seq)){
        wi<-wi+1
        posi[[3]]<-c(sear.chrom,sear.mar,mov.seq[sear.int])
        if(model==1)
          {result[[wi]]<-QTL.estim.p(nQTL, posi, ys, yu, D.matrix, marker.map, marker.gty)}
        if(model==2)
          {result[[wi]]<-QTL.estim.t(nQTL, posi, ys, tL, tR, D.matrix, marker.map, marker.gty)}
        }
      }
    }

## Summary the result

sear.posi<-matrix(0,n.search,3*nQTL)
sear.LRT<-rep(0,n.search)
sear.parameter<-matrix(0,n.search,ncol(D.matrix)+2)
for(L.ind in 1:n.search){
  s.posi<-c()
  for(n.po in 1:nQTL){
    s.posi<-c(s.posi,result[[c(L.ind,1,n.po)]])
    }
  sear.posi[L.ind,]<-s.posi
  sear.LRT[L.ind]<-result[[c(L.ind,2)]]
  sear.parameter[L.ind,]<-result[[c(L.ind,3)]]
  }
summ.result<-cbind(sear.posi,sear.LRT,sear.parameter)
QTLs<-result[[which.max(sear.LRT)]]
