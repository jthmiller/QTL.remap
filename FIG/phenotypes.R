#!/bin/R
##install.packages("RColorBrewer")

library(RColorBrewer)
library(ggplot2)
pheno.dist <- read.csv("~/Dropbox/QTL_Paper/METADATA/ALL_pheno-ALL_phenotype.csv", header=T)
str(pheno.dist)

###
BRP <- read.table('~/Desktop/bj/BRP.id.txt', sep=' ', stringsAsFactors = F)
BRP <- paste(BRP$V1,BRP$V2,sep='_')
NEW <- read.table('~/Desktop/bj/NEW.id.txt', sep=' ', stringsAsFactors = F)
NEW <- paste(NEW$V1,NEW$V2,sep='_')
ELR <- read.table('~/Desktop/bj/ELR.id.txt', sep=' ', stringsAsFactors = F)
ELR <- paste(ELR$V1,ELR$V2,sep='_')
NBH <- read.table('~/Desktop/bj/NBH.id.txt', sep=' ', stringsAsFactors = F)
NBH <- paste(NBH$V1,NBH$V2,sep='_')

ids <- paste(pheno.dist$pop_all,pheno.dist$Sample,sep = '_')
pheno.dist$GT_NG_ALT <- 'NG'
pheno.dist$GT_NG_ALT[which(ids %in% c(NBH,ELR,BRP,NEW))] <- 'GT'

write.csv(pheno.dist,'~/Dropbox/QTL_Paper/METADATA/PhenoDist.csv')

###

##########################################################################################
##########################################################################################
##Important stuff
##getting the proportions to put in a new csv by hand
brp.gt <- data.frame(pop='BRP',table(pheno.dist[which(pheno.dist$pop_all=="BRP" & pheno.dist$GT_NG_ALT=='GT'),]$pheno_all))
brp <- data.frame(pop='BRP',table(pheno.dist[which(pheno.dist$pop_all=="BRP"),]$pheno_all))
colnames(brp.gt) <- c('pop','phenotype','gt')
colnames(brp) <- c('pop','phenotype','ngt')

NBH.gt <- data.frame(pop='NBH',table(pheno.dist[which(pheno.dist$pop_all=="NBH" & pheno.dist$GT_NG_ALT=='GT'),]$pheno_all))
NBH <- data.frame(pop='NBH',table(pheno.dist[which(pheno.dist$pop_all=="NBH"),]$pheno_all))
colnames(NBH.gt) <- c('pop','phenotype','gt')
colnames(NBH) <- c('pop','phenotype','ngt')

ELR.gt <- data.frame(pop='ELR',table(pheno.dist[which(pheno.dist$pop_all=="ELR" & pheno.dist$GT_NG_ALT=='GT'),]$pheno_all))
ELR <- data.frame(pop='ELR',table(pheno.dist[which(pheno.dist$pop_all=="ELR"),]$pheno_all))
colnames(ELR.gt) <- c('pop','phenotype','gt')
colnames(ELR) <- c('pop','phenotype','ngt')

NEW.gt <- data.frame(pop='NEW',table(pheno.dist[which(pheno.dist$pop_all=="NEW" & pheno.dist$GT_NG_ALT=='GT'),]$pheno_all))
NEW <- data.frame(pop='NEW',table(pheno.dist[which(pheno.dist$pop_all=="NEW"),]$pheno_all))
colnames(NEW.gt) <- c('pop','phenotype','gt')
colnames(NEW) <- c('pop','phenotype','ngt')

df1 <- merge(brp,brp.gt,all=T)
df2 <- merge(NEW,NEW.gt,all=T)
df3 <- merge(ELR,ELR.gt,all=T)
df4 <- merge(NBH,NBH.gt,all=T)
df <- rbind(df1,df2,df3,df4)
            
##keeping colors consistent
all.pops <- c("NBH", "BRP", "ELR", "NEW")
popcol <- brewer.pal(8, "Paired")[c(2, 4, 6, 8)]
names(popcol) <- all.pops

df$pop <-factor(df$pop,levels=c("NBH","BRP","NEW","ELR"))
str(df)
df$phenotype <- as.numeric(as.character(df$phenotype))
##pretty graph
p <- ggplot(df, aes(x=phenotype, y=ngt,width = 1 ))+
  geom_bar(aes(color=pop), stat="identity")+
  geom_bar(aes(x=phenotype, y=gt, fill=pop), stat="identity")+
  facet_grid(pop~.)+
  xlab("Deformity Index")+
  ylab("Number of Embryos")+
  scale_x_continuous(breaks=seq(0,5,by=1),expand=c(0,0.1))+
  theme_bw()+
  scale_fill_manual(values=popcol)+
  scale_color_manual(values=popcol)+
  scale_y_continuous(breaks=c(25,50,75,100))+
  theme(legend.position="none")+
  theme(axis.text.y = element_text(face='bold', size = 12),
        axis.title = element_text(face='bold', size = 12),
        axis.text.x = element_text(face='bold', size = 12),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text.y = element_text(size=12, face="bold"))


g <- ggplot_gtable(ggplot_build(p))
strip_both <- which(grepl('strip-', g$layout$name))
fills <-as.vector(popcol[c(1,2,4,3)])
k <- 1
for (i in strip_both) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)



