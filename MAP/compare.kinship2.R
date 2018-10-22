source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')

load(paste('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/',pop,'/REMAPS/QTLmap.Rsave',sep=''))

cross <- subset(cross.18,ind=cross.18$pheno$stata=='ind')
cross <- dropSimilarMarkers(cross, rf.threshold = 0.002)
cross <- est.rf(cross)
cross <- est.map(cross,error.prob=ers,map.function="kosambi",maxit=1000)
zero.map <- shiftmap(pull.map(cross))
cross <- replace.map(cross, zero.map)
cross <- removeDoubleXO(cross, verbose=T)
cross <- est.rf(cross)
cross <- repRipple(cross, error.prob=ers, map.function="kosambi",window = 6)

cross.pars <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
        format='csvr', geno=c(1:3),estimate.map=FALSE)

cross.pars <- subset(cross,ind=is.na(cross$pheno$Pheno))
cross.pars$pheno$ID <-  c('par1','par2')
drop <- markernames(cross.pars)[!markernames(cross.pars) %in% markernames(cross.18)]
cross.pars  <- drop.markers(cross.pars,drop)
cross.pars <- c(cross.pars,cross.18)

save.image(paste('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/',pop,'/REMAPS/kinship_QTLmap.Rsave',sep=''))

cross.18 <- orderMarkers(cross.Z,window=5,use.ripple=T,
  error.prob=0.05, map.function='kosambi',sex.sp=F,maxit=2000,tol=1e-2)

  rela <- comparegeno(cross.18)
  colnames(rela) <- cross.b$pheno$ID
  rownames(rela) <- cross.b$pheno$ID
  rela[rela==NaN] <- NA
  diag(rela) <- NA

  hist(rela)
  ##resize
  heatmap(rela,symm=T)

save.image(paste('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/',pop,'/REMAPS/kinship_QTLmap.Rsave',sep=''))



#drop <- markernames(cross.18)[!markernames(cross.18) %in% markernames(cross)]
#cross.18 <- drop.markers(cross.18,drop)
#cross.Z <- removeDoubleXO(cross.Z, verbose=T)
#cross.18 <- orderMarkers(cross.Z,window=5,use.ripple=T,
#  error.prob=0.05, map.function='kosambi',sex.sp=F,maxit=2000,tol=1e-2)
#cross.18 <- repRipple(cross.18, error.prob=ers, map.function="kosambi",window = 6)
#POS.map.18 <- est.map(cross.18,error.prob=ers,map.function="kosambi",maxit=1000)
#cross.18 <- replace.map(cross.18, POS.map.18)
#cross.18 <- jittermap(cross.18, amount=1e-6)
#cross.18 <- sim.geno(cross.18,error.prob=ers)
#cross.pars <- read.cross.jm(file=file.path(indpops,paste(pop,'.unphased.f2.csvr',sep='')),
#        format='csvr', geno=c(1:3),estimate.map=FALSE)
#cross.pars <- subset(cross,ind=is.na(cross$pheno$Pheno))
#cross.pars$pheno$ID <-  c('par1','par2')
#drop <- markernames(cross.pars)[!markernames(cross.pars) %in% markernames(cross.18)]
#cross.pars  <- drop.markers(cross.pars,drop)
#cross.pars <- c(cross.pars,cross.18)
