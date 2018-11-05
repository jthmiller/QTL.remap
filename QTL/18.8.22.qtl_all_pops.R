#!/bin/bash
#debug.cross<- T
source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R')
library(ggplot2)
library("ggridges")
library('devtools')

NBH <- new.env()
load('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NBH/REMAPS/QTLmap.Rsave', envir=NBH)
ELR <- new.env()
load('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/ELR/REMAPS/QTLmap.Rsave', envir=ELR)
NEW <- new.env()
load('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/NEW/REMAPS/QTLmap.Rsave', envir=NEW)

ELR$cross.18$pheno$ID_pop <- 'ELR'
NEW$cross.18$pheno$ID_pop <- 'NEW'
NBH$cross.18$pheno$ID_pop <- 'NBH'

#### Sim and reduce to grid ####
NBH.x <- qtl::clean(NBH$cross.18)
NEW.x <- qtl::clean(NEW$cross.18)
ELR.x <- qtl::clean(ELR$cross.18)

nbh.grid <- sim.geno(NBH.x, n.draws=250, step=5, off.end=0,
  error.prob=0.01,map.function="kosambi",stepwidth="fixed")
new.grid <- sim.geno(NEW.x, n.draws=250, step=5, off.end=0,
  error.prob=0.01,map.function="kosambi",stepwidth="fixed")
elr.grid <- sim.geno(ELR.x, n.draws=250, step=5, off.end=0,
  error.prob=0.01,map.function="kosambi",stepwidth="fixed")

nbh.grid <- reduce2grid(nbh.grid)
new.grid <- reduce2grid(new.grid)
elr.grid <- reduce2grid(elr.grid)

scan.NBH <- scanone(nbh.grid, method="imp",model='normal',pheno.col=6)
scan.NEW <- scanone(new.grid, method="imp",model='normal',pheno.col=6)
scan.ELR <- scanone(elr.grid, method="imp",model='normal',pheno.col=6)

melted.nbh <- data.frame(pop='NBH',chr=scan.NBH$chr,pos=scan.NBH$pos,lod=scan.NBH$lod)
melted.new <- data.frame(pop='NEW',chr=scan.NEW$chr,pos=scan.NEW$pos,lod=scan.NEW$lod)
melted.elr <- data.frame(pop='ELR',chr=scan.ELR$chr,pos=scan.ELR$pos,lod=scan.ELR$lod)
melted <- rbind(melted.nbh,melted.new,melted.elr)
pals <- brewer.pal(11,"Set3")
#melted$pop <- factor(melted$pop, levels=c('NBH','NEW','ELR'))

### Saved to home
save.image(file.path('~/NEW.NBH.ELR.QTLmap.Rsave'))
#load(file.path('~/NEW.NBH.ELR.QTLmap.Rsave'))


phenosN <- c(NBH$cross.18$pheno$pheno_norm,ELR$cross.18$pheno$pheno_norm,NEW$cross.18$pheno$pheno_norm)
phenos5 <- c(NBH$cross.18$pheno$pheno_05,ELR$cross.18$pheno$pheno_05,NEW$cross.18$pheno$pheno_05)
names.pop <- c(NBH$cross.18$pheno$ID_pop,ELR$cross.18$pheno$ID_pop,NEW$cross.18$pheno$ID_pop)
phenos5 <- data.frame(pop=as.factor(names.pop),phen=as.numeric(phenos5))
phenosN <- data.frame(pop=as.factor(names.pop),phen=as.numeric(phenosN))

#Density ridges
png('/home/jmiller1/public_html/phenotypes.png',width = 600)
ggplot(phenosN, aes(phen,pop,height=..density..,group=pop,fill=as.factor(pop),color=as.factor(pop))) +
  geom_density_ridges(alpha=0.5,scale=1.3, bandwidth=.5)+
     xlim(-1, 6) +
     xlab("Phenotype Score") +
     ylab("Density") +
     theme(text = element_text(size=20),axis.text.y = element_blank())
dev.off()

png('/home/jmiller1/public_html/phenotypes.png',width = 600)
ggplot(phenos5, aes(phen,pop,height=..density..,group=pop,fill=as.factor(pop),color=as.factor(pop)),show.legend = F) +
  geom_density_ridges(alpha=0.5,scale=1.3, bandwidth=.5)+
     xlim(0, 5) +
     xlab("Phenotype Score") +
     ylab("Density") +
     scale_alpha(guide = 'none')
     scale_fill_discrete(name = "Populations") +
     theme(text = element_text(size=20),axis.text.y = element_blank())
dev.off()





pheno.all <- phen <- read.table('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/metadata/ALL_phenotype_Dist.txt',header=T)
phen$Pheno_05 <- phen$pheno_all
index <- which(phen$pop_all=='BRP')
















https://cran.rstudio.com/web/packages/mapfuser/vignettes/mapfuser.html


png('/home/jmiller1/public_html/ggplot2.qtl.png',width = 3000)
ggplot(melted,
  aes(y=as.factor(pop),x=pos,height=lod,fill=pop)) +
  geom_ridgeline(stat='identity',alpha=0.5,scale=0.25) +
  facet_wrap(~as.factor(chr),scales='free_x',nrow = 1, ncol = 24,strip.position = NULL) +
  #scale_x_continuous('pos') +
  #scale_y_continuous('lod') +
  scale_y_discrete('lod')+
  theme(axis.text=element_text(size=10))
dev.off()

#Phenotypes
phenos5 <- c(NBH$cross.18$pheno$pheno_05,ELR$cross.18$pheno$pheno_05,NEW$cross.18$pheno$pheno_05)
names.pop <- c(NBH$cross.18$pheno$ID_pop,ELR$cross.18$pheno$ID_pop,NEW$cross.18$pheno$ID_pop)
phenos <- data.frame(pop=as.factor(names.pop),phen=factor(phenos5, ordered = TRUE))



png('/home/jmiller1/public_html/phenotypes.png',width = 600)
ggplot(phenos, aes(phen,fill=pop)) +
  geom_histogram(aes(y = ..density..), breaks=seq(0,5),position="dodge")
dev.off()


png('/home/jmiller1/public_html/phenotypes.png',width = 600)
ggplot(phenos, aes(phen,fill=as.factor(pop),color=as.factor(pop))) +
  geom_density(aes(y = ..count..),binwidth=1,position="dodge")+
  scale_x_discrete(breaks=seq(1,5))
dev.off()


png('/home/jmiller1/public_html/phenotypes.png',width = 600)
ggplot(phenos, aes(phen,fill=as.factor(pop),color=as.factor(pop))) +
  geom_density(aes(y = ..density..), alpha=0.5)+
  scale_x_discrete(breaks=seq(1,5))
dev.off()


png('/home/jmiller1/public_html/phenotypes.png',width = 600)
ggplot(phenos, aes(phen,pop,height=..density..,group=pop,fill=as.factor(pop),color=as.factor(pop))) +
  geom_density_ridges2(alpha=0.5,scale=1, bandwidth=1)+
     xlim(0, 5)
dev.off()


png('/home/jmiller1/public_html/phenotypes.png',width = 600)
ggplot(phenos, aes(phen,pop,height=..density..,group=pop,fill=as.factor(pop),color=as.factor(pop))) +
  geom_density_ridges(alpha=0.5,scale=1.3, bandwidth=.5)+
     xlim(0, 5) +
     xlab("Phenotype Score") +
     ylab("Density") +
     theme(text = element_text(size=20))
dev.off()


phenos5 <- c(NBH$cross.18$pheno$pheno_norm,ELR$cross.18$pheno$pheno_norm,NEW$cross.18$pheno$pheno_norm)
names.pop <- c(NBH$cross.18$pheno$ID_pop,ELR$cross.18$pheno$ID_pop,NEW$cross.18$pheno$ID_pop)
phenos <- data.frame(pop=as.factor(names.pop),phen=as.numeric(phenos5))









png('/home/jmiller1/public_html/phenotypes.png',width = 300)
ggplot(phenos, aes(x=phen, colour=pop)) + geom_density()
dev.off()

png('/home/jmiller1/public_html/phenotypes.png',width = 3000)
ggplot(phenos, aes(as.numeric(phenos$phen),fill=pop)) +
  geom_bar(stat = 'count',position="dodge",width = 1) +
  scale_x_discrete(limits = as.numeric(phenos$phen))
dev.off()

ggplot(data=phenos, aes(x=factor(dose), y=length, fill=supp)) +
    geom_bar(stat="identity", position=position_dodge())

ggplot(phenos) + geom_density(aes(x = yield, color = site))

png('/home/jmiller1/public_html/phenotypes.png',width = 3000)
ggplot(phenos, aes(phen)) +
  stat_density(aes(fill=factor(pop)), alpha=0.8)
dev.off()

p <-ggplot(mydata, aes(months, values))
p +geom_bar(stat = "identity", aes(fill = type), position = "dodge") +
  xlab("Months") + ylab("Count") +
  ggtitle("Chickens & Eggs") +
  theme_bw()


  png('/home/jmiller1/public_html/phenotypes.png',width = 3000)

  qplot(phenos, data = phenos$phen, geom = "density",
      color = pop, linetype = pop)
dev.off()



png('/home/jmiller1/public_html/phenotypes.png',width = 3000)
ggplot(phenos, aes(phen)) +
  geom_bar(aes(fill=pop),stat="count", position = "dodge",
  bins=6, bin = 1,
   col="black",
   size=.1)
dev.off()



png('/home/jmiller1/public_html/phenotypes.png',width = 3000)
ggplot(phenos, aes(phen, colour = pop)) +
  geom_freqpoly(binwidth = 1,stat="count")
  dev.off()




  #Phenotypes
  phenos5 <- c(NBH$cross.18$pheno$pheno_05,ELR$cross.18$pheno$pheno_05,NEW$cross.18$pheno$pheno_05)
  names.pop <- c(NBH$cross.18$pheno$ID_pop,ELR$cross.18$pheno$ID_pop,NEW$cross.18$pheno$ID_pop)
  phenos <- data.frame(pop=as.factor(names.pop),phen=as.numeric(phenos5))



tabso <- table(phenos)

png('/home/jmiller1/public_html/phenotypes.png',width = 3000)
 ggplot(tabso, aes(x = phen,y=pop, color=pop))+
 geom_bar(stat="count", position = "dodge")
  dev.off()










check.ELR <- scanone(ELR$cross.18, method="imp",model='normal',pheno.col=6)

cm.len <- cbind(
  NBH=summaryMap(NBH$cross.18)$length,
  NEW=summaryMap(NEW$cross.18)$length,
  ELR=summaryMap(ELR$cross.18)$length)

cm.len.max <- apply(cm.len,1,max)
cm.len <- cm.len.max-cm.len[,1:2]
cm.len <- cm.len[1:24,]


a <- unlist(lapply(strsplit(markernames(NEW$cross.18),':',fixed=T),'[[',1))
b <- unlist(lapply(strsplit(markernames(NEW$cross.18),':',fixed=T),'[[',2))

new.pos <- data.frame(chr=a,pos=b)
new.map <- pull.map(NEW$cross.18)
nbh.map <- pull.map(NBH$cross.18)
try <- interpPositions(new.pos, new.map, nbh.map)


NBH.x <- qtl::clean(NBH$cross.18)
NEW.x <- qtl::clean(NEW$cross.18)
ELR.x <- qtl::clean(ELR$cross.18)

cross <- c(NBH.x,NEW.x,ELR.x)

cross <- sim.geno(cross, n.draws=20, step=1, off.end=10,
  error.prob=0.001,map.function="kosambi",stepwidth="fixed")

cross.nbh <- subset(cross,ind=cross$pheno$cross1==1)
cross.new <- subset(cross,ind=cross$pheno$cross2==1)
cross.elr <- subset(cross,ind=cross$pheno$cross3==1)

scan.norm.imp.NBH <- scanone(cross.nbh, method="imp",model='normal',pheno.col=6)
scan.norm.imp.NEW <- scanone(cross.new, method="imp",model='normal',pheno.col=6)
scan.norm.imp.ELR <- scanone(cross.elr, method="imp",model='normal',pheno.col=6)
### Saved to home
#save.image(file.path('~/NEW.NBH.ELR.QTLmap.Rsave'))
#load(file.path('~/NEW.NBH.ELR.QTLmap.Rsave'))






# Read a table with positions to interpolate and/or extrapolate
predict_file <- list.files(fpath,
                       pattern = "BaySha_physical",
                       full.names = T)
to_predict <- read.table(predict_file, sep = ",", header = T)
# Predict, accessible under MF.obj$predictions
MF.obj <- predict(MF.obj, to_predict)



MF.obj <-  genphys_fit(MF.obj, type = "consensus", z = 5)
plot(MF.obj, which = "mareymap", maps = "consensus")











melted.nbh <- data.frame(pop='NBH',chr=scan.norm.imp.NBH$chr,pos=scan.norm.imp.NBH$pos,lod=scan.norm.imp.NBH$lod)
melted.new <- data.frame(pop='NEW',chr=scan.norm.imp.NEW$chr,pos=scan.norm.imp.NEW$pos,lod=scan.norm.imp.NEW$lod)
melted <- rbind(melted.nbh,melted.new)
pals <- brewer.pal(11,"Set3")
melted$pop <- factor(melted$pop, levels=c('NBH','NEW'))
########### works #######
png('/home/jmiller1/public_html/ggplot2.qtl.png',width = 3000)
ggplot(melted,
  aes(y=as.factor(pop),x=pos,height=lod,fill=pop)) +
  geom_ridgeline(stat='identity',alpha=0.5,scale=0.25) +
  facet_wrap(~as.factor(chr),scales='free_x',nrow = 1, ncol = 24,strip.position = NULL) +
  #scale_x_continuous('pos') +
  #scale_y_continuous('lod') +
  scale_y_discrete('lod')+
  theme(axis.text=element_text(size=10))
dev.off()





#scan.norm.imp.NBH <- scanone(comp, method="imp",model='normal',pheno.col=6,ind.noqtl=comp$pheno$cross1==1)
#scan.norm.imp.NEW <- scanone(comp, method="imp",model='normal',pheno.col=6,ind.noqtl=comp$pheno$cross2==1)
#scan.norm.imp.cov <- scanone(comp, method="imp",model='normal',pheno.col=6,addcovar=comp$pheno$cross1)
#NBH.perm <- NBH$perms.norm.imp
#NEW.perm <- NEW$perms.norm.imp
#th.NBH <- summary(NBH.perm)[1,]
#th.NEW <- summary(NEW.perm)[1,]
#norm.qtl.NBH <- summary(scan.norm.imp.NBH, perms=NBH.perm, alpha=0.05)
#norm.qtl.NEW <- summary(scan.norm.imp.NEW, perms=NEW.perm, alpha=0.05)
#qtl.uns.NBH <- makeqtl(cross, chr=norm.qtl.NBH$chr, pos=norm.qtl.NBH$pos,what="draws")
#qtl.uns.NEW <- makeqtl(cross, chr=norm.qtl.NEW$chr, pos=norm.qtl.NEW$pos,what="draws")
#chrs.model <- unique(qtl.uns.NEW$chr, qtl.uns.NEW$chr)


scan.all <- c(scan.norm.imp.NBH,scan.norm.imp.NEW)

png('~/NEW_NBH.png',width = 2000)
plot(scan.norm.imp.NBH,scan.norm.imp.NEW, chr=1:24, lodcolumn=1,
  incl.markers=FALSE,lty=1, col=c("black","blue"), lwd=3,
  add=FALSE, gap=25, mtick = "line",bandcol="gray70",
  show.marker.names=FALSE, alternate.chrid=T,
  type="l", cex=1, pch=1, bg="grey",
  bgrect=NULL)
dev.off()

png('~/NBH_NEW_ALL.png',width = 1000)
plot(scan.norm.imp.NBH, scan.norm.imp.NEW, chr=1:24, lodcolumn=1,
  incl.markers=FALSE,lty=1, col=c("violetred","black"), lwd=2,
  add=FALSE, gap=25, mtick = "line",bandcol="gray70",
  show.marker.names=FALSE, alternate.chrid=T,
  type="l", cex=1, pch=1, bg="grey",
  bgrect=NULL)
dev.off()


png('~/NBH_RF_2.png')
plotRF(cross.nbh, chr = 2)
dev.off()


## Positionmaps instead
for (X in 1:24){
  ord <- order(as.numeric(gsub(paste(X,":",sep=''),'',markernames(NBH$cross.18,chr=X))))
  NBH$cross.18 <- switch.order(NBH$cross.18, chr=X, ord, error.prob=0.01,
    map.function="kosambi",maxit=100, tol=1e-6, sex.sp=F)
  POS.map.18 <- est.map(NBH$cross.18,error.prob=0.01,map.function="kosambi", chr=X,maxit=10)
  NBH$cross.18 <- replace.map(NBH$cross.18, POS.map.18)
}

f

save.image(file.path('~/NEW.NBH.QTLmap.Rsave'))
#perms.norm.imp <- scanone(comp, method="imp",model='normal',n.perm=5,perm.strata=cross.18$pheno$stata,pheno.col=6, n.cluster=slurmcore)















png('~/ggplot2.qtl.png')
ggplot(melted,
  aes(y=lod,x=lod,height=lod)) +
  geom_ridgeline(alpha=0.5) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0, 0))+
  theme(axis.text=element_text(size=5))
dev.off()



png('/home/jmiller1/public_html/ggplot2.qtl.png',width = 1000)
ggplot(melted,
  aes(y=as.factor(pop),x=as.factor(chr),height=lod,width=pos)) +
  geom_ridgeline(stat='identity',alpha=0.5,scale=1) +
  scale_x_discrete() +
  scale_y_discrete()+
  theme(axis.text=element_text(size=10))
dev.off()

pals <- brewer.pal(11,"Set3")
melted$pop <- factor(melted$pop, levels=c('NBH','NEW'))
########### works #######
png('/home/jmiller1/public_html/ggplot2.qtl.png',width = 3000)
ggplot(melted,
  aes(y=as.factor(pop),x=pos,height=lod,fill=pop)) +
  geom_ridgeline(stat='identity',alpha=0.5,scale=0.25) +
  facet_wrap(~as.factor(chr),scales='free_x',nrow = 1, ncol = 24,strip.position = NULL) +
  #scale_x_continuous('pos') +
  #scale_y_continuous('lod') +
  scale_y_discrete('lod')+
  theme(axis.text=element_text(size=10))
dev.off()

melted2 <- melted[which(melted$chr==2),]

png('/home/jmiller1/public_html/ggplot2.qtl.png',width = 300)
ggplot(melted2,
  aes(y=as.factor(pop),x=pos,height=lod)) +
  geom_ridgeline(stat='identity',alpha=0.5,scale=0.25) +
  theme(axis.text=element_text(size=10))
dev.off()

png('/home/jmiller1/public_html/ggplot2.qtl.png',width = 300)
ggplot(melted2,
  aes(y=as.factor(pop),x=pos,height=lod)) +
  stat_smooth(method = lm, formula = y ~ poly(x,2), aes(colour = 'polynomial'), se= FALSE) +
  geom_ridgeline(stat='identity',alpha=0.5,scale=0.25) +
  theme(axis.text=element_text(size=10))
dev.off()



png('/home/jmiller1/public_html/ggplot2.qtl.png',width = 300)
ggplot(melted2,
  aes(y=lod,x=pos,colour=pop,height=lod)) +
  stat_smooth(method = lm, formula = y ~ poly(x,2), aes(colour = 'polynomial'), se= FALSE) +
  geom_ridgeline(stat='identity',alpha=0.5,scale=0.25) +
  theme(axis.text=element_text(size=10))
dev.off()





ggplot(melted, aes(x = lod, y = Species)) + geom_density_ridges2()





trym <- combineMap(NBH$cross.18,NEW$cross.18, keep.all = FALSE, merge.by = "marker" )
tryg <- combineMap(NBH$cross.18,NEW$cross.18, keep.all = TRUE, merge.by = "genotypes")
test3 <- mstmap.cross(cross.18, chr= 1:24, id = "ID", bychr = TRUE, suffix = "numeric",
     anchor = FALSE, dist.fun = "kosambi", objective.fun = "COUNT",
     p.value = 1e-06, noMap.dist = 15, noMap.size = 0, miss.thresh = .75,
     mvest.bc = FALSE, detectBadData = FALSE, return.imputed = FALSE,
     trace = FALSE)

test4 <- mstmap(test3, bychr = TRUE, id = "ID", dist.fun = "kosambi", p.value = 2, trace = FALSE,miss.thresh = .75)



BC.gen=2, F.gen=3


try <- reduce2grid(NBH$cross.18)


hyper <- calc.genoprob(hyper, step=2)

try <- argmax.geno(NBH$cross.18, step=0, off.end=0, error.prob=0.002,map.function="kosambi",stepwidth="fixed")
#try <- pull.argmaxgeno(try)
try <- reduce2grid(try)



hypersub <- reduce2grid(hyper)
## Not run: out <- scanone(hypersub)
plot(out, incl.markers=FALSE)
## End(Not run)
