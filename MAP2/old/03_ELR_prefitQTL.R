#!/bin/R
### Map QTLs 1 of 3
#debug.cross <- T
#source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library('qtl')
pop <- 'ELR'
source("/home/jmiller1/QTL_Map_Raw/ELR_final_map/CODE/control_file.R")
mpath <- '/home/jmiller1/QTL_Map_Raw/ELR_final_map'

################################################################################
## put chromosomes together
################################################################################
filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_all_mark_',i,'_tsp')
write.cross(cross,chr=i,filestem=filename,format="csv")

file_list <- list.files(mpath, 'ELR_all_mark_.*tsp.csv')

chr <- gsub("ELR_gts_CHR",'',file_list)
chr <- as.numeric(gsub("_downsmpl_reordered.csv",'',chr))

elr <- lapply(file_list,function(X){ read.cross(file=X,format = "csv", dir=mpath, genotypes=c("AA","AB","BB"), alleles=c("A","B"),estimate.map = FALSE)})

gnos <- lapply(elr,function(X){
 data.frame(X[[1]][[1]][['data']],stringsAsFactors=F)
})
gnos <- do.call(cbind,gnos)
gnos <- cbind(elr[[1]]$pheno,gnos)
gnos$ID <- as.character(gnos$ID)

m_names <- unlist(sapply(elr,function(X){
 markernames(X)
}))

colnames(gnos) <- c('Pheno','sex','ID','bin','pheno_norm',m_names)
rownames(gnos) <- elr[[1]]$pheno$ID

map <- c(colnames(elr[[1]]$pheno),unname(unlist(sapply(elr,pull.map))))
chr <- c(colnames(elr[[1]]$pheno),gsub(":.*","",m_names))
info <- c(colnames(elr[[1]]$pheno),m_names)
headers <- rbind(info,chr,map)
colnames(headers) <- headers[1,]
headers[2:3,1:5] <- ''

headers.u <- unname(data.frame(headers,row.names=NULL,stringsAsFactors=FALSE))
gnos.u <- unname(data.frame(lapply(gnos, as.character),row.names=NULL,stringsAsFactors=FALSE))
colnames(headers.u) <- colnames(gnos.u) <- headers.u[1,]
to_write <- rbind(headers.u,gnos.u)

################################################################################
write.table(to_write,file.path(mpath,'elr_all_mapped.tsp.csv'),sep=',',row.names=F,quote=F,col.names = F)
################################################################################

fl <- file.path(mpath,'elr_all_mapped.tsp.csv')

cross <- read.cross(
 file = fl,
 format = "csv", genotypes=c("1","2","3"),
 estimate.map = FALSE
)

cross <- sim.geno(cross)
################################################################################
#### DOWNSAMPLE CHROMOSOMES

ints <- seq(0, 100, by=0.1)
f <- function(a,x,b) abs(b - length(pickMarkerSubset(pull.map(cross)[[as.character(x)]],a)))
mrks <- sapply(ints,f,x=1,b=500)

f(a,x,b)

x <- 1
a <- 0
a <- 0.0000000000001
b <- 500

f  <- function(x) ifelse(x > -1, ifelse(x < 4, exp(-1/abs(x - 1)), 10), 10)
optimize(f, c(-4, 20), tol = 0.1)   # doesn't see the minimum
optimize(f, c(-7, 20))   # ok

fp <- function(x) { print(x); f(x) }
optimize(fp, c(-4, 20))   # doesn't see the minimum
optimize(fp, c(-7, 20))   # ok

optimize(f, c(1e-25,2), x = 1, b = 500)

fp <- function(a,x,b) { print(a); f(a,x,b) }
optimize(fp, c(1e-25,5), x = 1, b = 500)


f(0.000000000000000001,1,500)





f <- function (x, a) (x - a)^2
xmin <- optimize(f, c(0, 1), tol = 0.0001, a = 1/3)



x <- 1
a <- 1


interval(

f <- function (x, a) (x - a)^2

f <- function(x,a) abs(x-a)

optimize(f, c(1, 1000), tol = 0.1, a = 500)


xmin <- optimize(f, c(0, 1), tol = 0.0001, a = 1/3)

seq(0,1,by=0.0001)


x <- 1
a <- 0.3444185
a <- 0.0000000000001
b <- 500

ints <- seq(0.00000000001,1,by=0.0001)
f <- function(a,x,b) abs(b - length(pickMarkerSubset(pull.map(cross)[[as.character(x)]],ints[a])))


f <- function(a,x,b) length(pickMarkerSubset(pull.map(cross)[[as.character(x)]],a))
optimize(f, c(1,500), x = 1, b = 500)


f <- function(x) print(x)
optimize(f, c(1,500), x = 1, b = 500)


optimize(f, c(1,500), x = 1, b = 500)

shouls be 300

f(1,1,300)

length(pickMarkerSubset(pull.map(cross)[[as.character(x)]],0.04)

optimize(f, c(0.0000000000001, 10), tol = 0.0000000001, x = 1, b = 843)

f  <- function(x) ifelse(x > -1, ifelse(x < 4, exp(-1/abs(x - 1)), 10), 10)
fp <- function(x) { print(x); f(x) }

png(paste0('~/public_html/optimize.png'))

plot(f, -2,5, ylim = 0:1, col = 2)
plot(f, -2,5, ylim = 0:1, col = 2)
plot(fp, 1e-25,5, x = 1, b = 500)
dev.off()

optimize(fp, c(-4, 20), tol = 0.1)   # doesn't see the minimum
optimize(fp, c(-7, 20))   # ok

optimize(f, c(-7, 20))



abs(b - length(pickMarkerSubset(pull.map(cross)[[as.character(x)]],a)))


pickMarkerSubset(pull.map(cross)[[as.character(x)]],0.3444185)



subs <- sapply(1:24, function(X){

 if (nmar(cross) > 500){
  dis <- max(pull.map(cross)[[as.character(X)]])/1000
  pickMarkerSubset(pull.map(cross)[[as.character(X)]],dis)
 } else {
  markernames(cross, chr=X)
 }

})

filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_dwns_CHR',i,'_tsp')
write.cross(cross.d,chr=i,filestem=filename,format="csv")

  png(paste0('~/public_html/ELR_RF_concord',i,'_tsp.png'))
    plotRF(cross.d)
  dev.off()

################################################################################
cross_map <-  est.map(cross, error.prob=0.025,map.function="kosambi",maxit=100000,
  tol=1e-7, sex.sp=FALSE, verbose=FALSE, n.cluster=10)

cross_map <- shiftmap(cross_map, offset=0)

cross <- qtl:::replace.map(cross,cross_map)

dwnsmpl <- pickMarkerSubset(pull.map(cross)[[1]],100)

filename <- paste0('/home/jmiller1/QTL_Map_Raw/ELR_final_map/ELR_gts_ordered_TSP')
write.cross(cross,filestem=filename,format="csv")
################################################################################

png(paste0('~/public_html/ELR_RF_remap.png'))
plotRF(cross)
dev.off()
