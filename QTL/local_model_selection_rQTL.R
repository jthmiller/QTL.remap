cross.18 <- read.cross(file=file.path('~/Dropbox/QTL_Paper/DATA/tempout.csv'),format='csv',
  geno=c(1:3),estimate.map=FALSE)

dups <- findDupMarkers(cross.18, exact.only=FALSE, adjacent.only=FALSE)
### remove markers that are exactly the same.
cross.18 <- drop.markers(cross.18,unlist(dups))

### Ids for strata permutations
cross.18$pheno$stata <- gsub('[0-9]','',cross.18$pheno$ID)

### fix Phenos (not needed in final)
pheno.all <- phen <- read.table('/Users/jeffreymiller/Dropbox/QTL_Paper/Metadata/ALL_phenotype_Dist.txt',header=T)
phen$Pheno_05 <- phen$pheno_all
index <- which(phen$pop_all==pop)
count.pheno <- sapply(0:5, function(pt){
      pt <- as.character(pt)
      total <- sum(phen$Pheno_05[index]==pt)
      incross <- sum(cross.18$pheno$Pheno_05==pt)
      return(as.numeric(total-incross))
  }
)
names(count.pheno) <- as.character(0:5)
count.pheno <- count.pheno[!is.na(count.pheno)]
count.pheno <- rep.int(names(count.pheno), times=as.numeric(count.pheno))
trsl.bin <- c(0,0,0,1,1,1)
names(trsl.bin) <- as.character(0:5)
phenos <- data.frame(count.pheno,Pheno=as.numeric(trsl.bin[as.character(count.pheno)]))
rownames(phenos) <- paste('NG',1:length(count.pheno),sep='')
ind.inx <- grep('NG',cross.18$pheno$ID)
repl <- phenos[as.character(cross.18$pheno$ID[grep('NG',cross.18$pheno$ID)]),1]
cross.18$pheno$pheno_05[ind.inx] <- as.character(repl)
cross.18$pheno$pheno_05 <- as.numeric(cross.18$pheno$pheno_05)
cross.18$pheno$pheno <- as.numeric(cross.18$pheno$pheno_05 >= 3) ### Split pheotypes

### If M.Map position is used
#x <- 1
#order(as.numeric(gsub(paste(x,':',sep=''),'',markernames(cross.18,chr=1))),decreasing=F)
#set new order of the chr
#restimate the map


### Error rate and genoprobs
ers <- 0.002
cross.18 <- calc.genoprob(cross.18, step=1,error.prob=ers,map.function='kosambi')

#### Nonparametric scan
scan.np.em <- scanone(cross.18, method='em', model="np", pheno.col=2, maxit=5000)
#perms.np.em <- scanone(cross.18, model="np", pheno.col=2, n.perm=500, method='em',perm.strata=cross.18$pheno$stata)

### 2 part model
scan.2p.em <- scanone(cross.18, method='em', model="2part", pheno.col=1, maxit=5000)
#perms.2p.em <- scanone(cross.18, model="2part",maxit=500,n.perm=500,perm.strata=cross.18$pheno$stata,pheno.col=2)

### Binary model with EM
scan.bin.em <- scanone(cross.18, method="em",model='binary', pheno.col=1, maxit=5000)
#perms.bin.em <- scanone(cross.18, method="em",model='binary',maxit=500,n.perm=500,perm.strata=cross.18$pheno$stata,pheno.col=1)

### Binary model with marker regression
scan.bin.mr <- scanone(cross.18, method="mr",model='binary',pheno.col=1)
#perms.bin.mr <- scanone(cross.18, method="mr",model='binary', n.perm=500, perm.strata=cross.18$pheno$stata)

### Binary model with haley knott (heavy inflation of LOD)
scan.bin.hk <- scanone(cross.18, method="hk",model='binary',pheno.col=1)
#perms.bin.em <- scanone(cross.18, method="hk",model='binary',maxit=500,n.perm=500,perm.strata=cross.18$pheno$stata,pheno.col=1)

#### Imputations (better for selective genotyping than HK)
Z <- transformPheno(cross.18, pheno.col=2, transf=nqrank)
Z <- sim.geno(Z,error.prob=0.002)
qtl.uns <- makeqtl(Z, chr=norm.qtl$chr, pos=norm.qtl$pos)
full <- stepwiseqtl(Z,additive.only=T,method="imp",pheno.col=2, scan.pairs=T)

### Normal scan on transformed phenotype
scan.norm.em <- scanone(Z, method="em",model='normal',maxit=500, pheno.col=2)
#perms.norm.em <- scanone(Z, method="em",model='normal',maxit=500,n.perm=500,perm.strata=cross.18$pheno$stata,pheno.col=2)

### Normal scan on transformed phenotype w/Extended haley knott (better for selective/missing genos at non-gt'ed ind)
scan.norm.ehk <- scanone(Z, method="ehk",model='normal',maxit=500, pheno.col=2)
#perms.norm.ehk <- scanone(Z, method="ehk",model='normal',maxit=500,n.perm=500,perm.strata=cross.18$pheno$stata,pheno.col=2)

### Normal scan on transformed phenotype fast haley knott (not robust to missing data. LOD inflation)
scan.norm.hk <- scanone(Z, method="hk",model='normal',maxit=500, pheno.col=2)
#perms.norm.hk <- scanone(Z, method="hk",model='normal',maxit=500,n.perm=500,perm.strata=cross.18$pheno$stata,pheno.col=2)

### Imputation on transformed phenotype
rorororororororororororororororo <- scanone(Z, method="imp",model='normal',maxit=500, pheno.col=2)
perms.norm.imp <- scanone(Z, method="imp",model='normal',maxit=500,n.perm=500,perm.strata=cross.18$pheno$stata,pheno.col=2)


### Imputation on transformed phenotype
indy <- Z$pheno$ID[grep('ind',Z$pheno$ID)]
cross.gt <- subset(Z,ind=as.character(indy))
scan.norm.imp.gt <- scanone(cross.gt, method="imp",model='normal',maxit=500, pheno.col=2)

### Imputation on un-transformed scores
scan.norm.imp.05 <- scanone(cross.18, method="imp",model='normal',maxit=500, pheno.col=2)

## EM scan without ungenotyped
scan.norm.em.gt <- scanone(cross.gt, method="em",model='normal',maxit=500, pheno.col=2)




### Multi QTL Models ###
### New permutation to run analysis without chr 2
perms.norm.imp.2 <- scanone(Z, method="imp",model='normal',chr=-2, maxit=500,n.perm=500,perm.strata=cross.18$pheno$stata,pheno.col=2)

## make a cofactor vector of QTL genotypes to consider in the scan
mar <- find.marker(Z, chr=2, pos=92.6)
g <- pull.geno(fill.geno(Z))[,mar]
g <- cbind(as.numeric(g==1), as.numeric(g==2))

scan.norm.imp.2ad <- scanone(Z.sg, method="imp",model='normal',maxit=500, pheno.col=2,addcovar=g)
scan.norm.imp.2in <- scanone(Z.sg, method="imp",model='normal',maxit=500, pheno.col=2,addcovar=g,intcovar=g)
plot(scan.norm.imp.2in - scan.norm.imp.2ad, ylab="interaction LOD score")

### Fit full model QTL stepwise and manually with imputation and transformed data
th <- summary(perms.norm.imp)[1,]
norm.qtl <- summary(scan.norm.imp, perms=perms.norm.imp, alpha=0.05)
full.sw <- fitqtl(Z, qtl=qtl.sw, get.ests=FALSE, method="imp",pheno.col=2)
full.man <- fitqtl(Z, qtl=qtl.uns, get.ests=FALSE, method="imp",pheno.col=2)








mar <- find.marker(Z.sg, chr=norm.qtl$chr, pos=norm.qtl$pos)
g <- pull.geno(Z.sg)[,mar]
We will use three marker covariates, and window sizes of 20 and 40 cM,
as well as infinity (meaning the entire length of the chromosome).
out.cim.40 <- cim(hyper, n.marcovar=3, window=40)
out.cim.inf <- cim(Z.sg, n.marcovar=3, window=Inf,pheno.col=2,method="imp", error.prob=0.002)

# plot(scan.2p.em,chr=qs, lty=1,col="blue", add=TRUE)

"scan.bin.em"
"scan.bin.hk"
"scan.bin.mr"

"scan.2p.em"
"scan.np.em"

"scan.norm.ehk"
"scan.norm.em"
"scan.norm.hk"
"scan.norm.imp"

#### Plots
#qs <- c(1,2,8,9,18,22,24)
qs <- c(1:24)

## Binary vs nonparametric vs transformed
plot(scan.np.em,chr=qs, lty=1,col="black")
plot(scan.bin.mr,chr=qs, lty=1,col="black", add=TRUE)
#plot(scan.2p.em,chr=qs, lty=1,col="blue", add=TRUE)
plot(scan.bin.em,chr=qs, lty=1,col="green", add=TRUE)


## Violate normal assumptions, transform
plot(scan.norm.em,chr=qs, lty=1,col="red", main='Violate normal assumptions, transform')
plot(scan.norm.ehk, chr=qs, lty=1,col="blue", add=TRUE)
plot(scan.norm.em.tf,chr=qs, lty=1,col="black", add=TRUE)
plot(scan.norm.imp,chr=qs, lty=1,col="purple", add=TRUE)

## Normal vs non-normal plot (EM-normal vs EM-binary vs EM-NP)
plot(scan.norm.em,chr=qs, lty=1,col="red", main='Violate normal assumptions, all EM')
plot(scan.bin.em,chr=qs, lty=1,col="green", add=TRUE)
plot(scan.np.em,chr=qs, lty=1,col="black",add=TRUE)
## Binary appears to behave differently- norm and nonparamtric are similar.


## Plot all binary
plot(scan.bin.hk,chr=qs, lty=1,col="red", main='Haley-Knott Inflation of LOD')
plot(scan.bin.mr,chr=qs, lty=1,col="black",add=TRUE)
plot(scan.bin.em,chr=qs, lty=1,col="green", add=TRUE)
plot(scan.norm.imp,chr=qs, lty=1,col="purple", add=TRUE)
## EM,HK,MR behave almost the same
## IMP diverges from all of the above, likely because it is dealing with missing GTs (EM, MR, and HK probably drop individuals)

## All EM
plot(scan.norm.em,chr=qs, lty=1,col="red", main='EM')
plot(scan.bin.em,chr=qs, lty=1,col="green", add=TRUE)
plot(scan.np.em,chr=qs, lty=1,col="black", add=TRUE)
plot(scan.norm.em.gt,chr=qs, lty=1,col="black", add=TRUE)
## Norm and bin em are similar except in a spot in chr3
## NP tracks just under norm, NP greater emphasized at large effect QTLs
## NP vs bin, NP tracks but is much inflated
## EM does seem to consider ungenotyped ind.. but far inflated)

### Non parametric and transformed IMP
plot(scan.norm.imp,chr=qs, lty=1,col="purple")
plot(scan.np.em,chr=qs, lty=1,col="black", add=TRUE)
### It seems that NP test inflates and doesn't deal with missing GTs

## Transformed Normal
plot(scan.bin.mr,chr=qs, lty=1,col="black")
plot(scan.norm.ehk, chr=qs, lty=1,col="red", add=TRUE)
plot(scan.norm.em,chr=qs, lty=1,col="blue", add=TRUE)
plot(scan.norm.imp,chr=qs, lty=1,col="purple", add=TRUE)
### EHK/Mr are similar but MR spikes where EHK doesnt (Due to missing GT info?)
### EHK/IMP are similar, except LOD is inflated at QTL peaks.
## EM/EHK results same with transformed data
### IMP tracks others, but with lower LOD

## IMP (compare with and without ungenotyped ind)
plot(scan.norm.imp.gt, chr=qs, lty=1,col="black")
plot(scan.norm.imp,chr=qs, lty=1,col="purple", add=TRUE)

## Imp with un-transformed and transformed
plot(scan.norm.imp.05,chr=qs, lty=1,col="black")
plot(scan.norm.imp,chr=qs, lty=1,col="purple", add=TRUE)
## Seems that I have more power at minor QTL- bu otherwise similar











plot(scan.norm.ehk, chr=qs, lty=1,col="red")
plot(scan.bin.mr,chr=qs, lty=1,col="purple", add=TRUE)
plot(scan.bin.em,chr=qs, lty=1,col="green", add=TRUE)
plot(scan.np.em,chr=qs, lty=1,col="black", add=TRUE)
plot(scan.norm.em.tf,chr=qs, lty=1,col="red",add=TRUE)
plot(scan.norm.imp,chr=qs, lty=1,col="red",add=TRUE)
abline(h=summary(perms.norm.em)[1,],col="red")
abline(h=summary(perms.bin.em)[1,],col="green")
abline(h=summary(perms.np.em)[1,],col="black")
abline(h=summary(perms.2p.em)[1,],col="blue")
abline(h=summary(perms.bin.mr)[1,],col="purple")







## test
maug <- fill.geno(cross.18,error.prob=0.001)
dups <- findDupMarkers(maug, exact.only=FALSE, adjacent.only=FALSE)
### remove markers that are exactly the same.
maug <- drop.markers(maug,unlist(dups))
scan.mqm <- scanone(maug, method='em', model="np", pheno.col=2, maxit=500)
#maug <- mqmaugment(cross.18, minprob=1.0)
mqm <- mqmscan(maug)
max(mqm)
## add that marker to the model
multitoset <- find.marker(maug, chr=2, pos=106.95)
setcofactors <- mqmsetcofactors(maug, cofactors=multitoset)
### scan with chr2 QTL cofactor
mqm_co1 <- mqmscan(maug, setcofactors)
mqmscanall(cross, multicore=TRUE, n.clusters = 1,batchsize=10,cofactors=NULL, ...)
## plot
par(mfrow = c(2,1))
plot(mqmgetmodel(mqm_co1))
plot(mqm_co1)
