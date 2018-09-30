

#data<-read.csv("~/Dropbox/QTL_Paper/Lee_et_all/example.csv")


cross.18 <- read.cross.jm(file=file.path(popdir,'tempout'),format='csv',
  geno=c(1:3),estimate.map=FALSE)

## Trait values

y <- cross.18$pheno$Pheno_05
ngy.ind<-is.na(data[,ncol(data)])
ys<-y[ngy.ind==0]
yu<-y[ngy.ind==1]


## Give the chromosome number.
# n.chrom: number of chromosomes

##n.chrom<-4
n.chrom<-24

## Give the marker number in chromosomes.
# info.chrom: a vector of number of markers on each chromosome

###info.chrom<-c(11,11,11,11)
info.chrom <- nmar(cross.18)

## Give the length of marker intervals on each chromosome (in Morgan).
# marker.map[[i]]: the length of marker intervals in the ith chromosome

#marker.map<-vector("list",n.chrom)
#marker.map[[1]]<-c(rep(0.15,10))
#marker.map[[2]]<-c(rep(0.15,10))
marker.map[[3]]<-c(rep(0.15,10))
marker.map[[4]]<-c(rep(0.15,10))

marker.map<-vector("list",1)
marker.map <- lapply(1:n.chrom,function(Z){
  return(cross.18$geno[[as.character(Z)]]$'map' - c(0,cross.18$geno[[as.character(Z)]]$'map')[-length(cross.18$geno[[as.character(Z)]]$'map')])
  }
)

## Marker genotypes on each chromosome.
# marker.gty[[i]]: marker genotypes in the ith chromosome


##data <- cross.18$geno[['1']]$'data'

marker.map <- lapply(1:n.chrom,function(Z){
  return(cross.18$geno[[as.character(Z)]]$
)





mar.gty <- data[ngy.ind==0,]  ###remove phenotype col
marker.gty <- vector("list",n.chrom)
marker.gty[[1]] <- mar.gty[,1:info.chrom[1]]

for(infc in 2: n.chrom){
  marker.gty[[infc]]<-mar.gty[,(sum(info.chrom[1:(infc-1)])+1):sum(info.chrom[1:infc])]
  }



mar.gty<-data[ngy.ind==0,-1]
marker.gty<-vector("list",n.chrom)
marker.gty[[1]]<-mar.gty[,1:info.chrom[1]]
for(infc in 2: n.chrom){
  marker.gty[[infc]]<-mar.gty[,(sum(info.chrom[1:(infc-1)])+1):sum(info.chrom[1:infc])]
  }





## Specify the chromosome region you would like to search for QTLs and the QTLs you would like to
## condition in in the model for the serach.
# For example: If the detected QTLs at (1, 3, 0.01), (1, 5, 0.12) and (4, 9, 0) are used for
#              conditioning in the model in the search of chromosomes 2 and 3 (one QTL at a time)
#              for QTLs. Here (1, 3, 0.01) indicates the position is at 0.01 Morgan away from the
#              left marker of the third marker interval in the first chromosome.

## Give the QTL number.

nQTL<-4

## Genetic design matrix.
D.matrix<-Design(nQTL,domi=F,epis=F)

## Give the positions of detected QTLs.
# Note that the order is sort according to chromosome and marker.

posi<-vector("list",nQTL)
posi[[1]]<-c(1,3,0.01)
posi[[2]]<-c(1,5,0.12)
posi[[4]]<-c(4,9,0)

## Give the walking speed along chromsome in Morgan.

sp<-0.01

## Give the search informations.
# sear.chrom: the information for chromosome
# sear.mar: the information for marker interval

## Give model type
# model: model=1 indicates using the proposed method in Lee et al. (2014);
#        model=2 indicates using the truncated model

model<-1
# If model=2, users need to give the values of truncation points (tL and tR).

wi<-0
n.search<-2*10*16
result<-vector("list",n.search)
  for(sear.chrom in 2:3){
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
