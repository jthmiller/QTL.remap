#!/bin/bash
#source('/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/debug.R')
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

check.ELR <- scanone(ELR$cross.18, method="imp",model='normal',pheno.col=6)

cm.len <- cbind(
  NBH=summaryMap(NBH$cross.18)$length,
  NEW=summaryMap(NEW$cross.18)$length,
  ELR=summaryMap(ELR$cross.18)$length)

cm.len.max <- apply(cm.len,1,max)
cm.len <- cm.len.max-cm.len[,1:2]
cm.len <- cm.len[1:24,]

ELR$cross.18$ID <- paste('ELR',ELR$cross.18$ID,sep='_')
NBH$cross.18$ID <- paste('NBH',NBH$cross.18$ID,sep='_')
NEW$cross.18$ID <- paste('NEW',NEW$cross.18$ID,sep='_')

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

#### Sim and reduce to grid ####
NBH.x <- qtl::clean(NBH$cross.18)
NEW.x <- qtl::clean(NEW$cross.18)
ELR.x <- qtl::clean(ELR$cross.18)


nbh.grid <- sim.geno(NBH.x, n.draws=20, step=5, off.end=0,
  error.prob=0.01,map.function="kosambi",stepwidth="fixed")
new.grid <- sim.geno(NEW.x, n.draws=20, step=5, off.end=0,
  error.prob=0.01,map.function="kosambi",stepwidth="fixed")
elr.grid <- sim.geno(ELR.x, n.draws=20, step=5, off.end=0,
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
