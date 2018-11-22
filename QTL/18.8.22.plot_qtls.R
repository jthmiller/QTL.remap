#!/bin/R
### first run combine pops for multi-pop cross objects

debug.cross <- T
source("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/scripts/QTL_remap/MAP/control_file.R")
library("ggridges")

load(file.path('~/NEW.NBH.ELR.QTLmap.Rsave'))
# load(file.path('~/NEW.NBH.ELR.QTLmap.Rsave'))

## Lod on 2
png("/home/jmiller1/public_html/chr2.lod.qtl.png", width = 600)
plotRF(cross.NEW,chr=2)
dev.off()

# All chr not scaled
png("/home/jmiller1/public_html/ggplot2.qtl.png", width = 3000)
ggplot(melted, aes(y = as.factor(pop), x = pos, height = lod, fill = pop)) + ylab(NULL) +
  geom_ridgeline(stat = "identity", scale = 0.2) + facet_wrap(~as.factor(chr),
  scales = "free_x", nrow = 1, ncol = 24, strip.position = NULL) + scale_y_discrete("lod") +
  scale_fill_manual(name = "Populations", guide = F, values = popcol) + theme(axis.text = element_text(size = 10),
  axis.text.y = element_blank())
dev.off()

# All chr Not scaled either
png("/home/jmiller1/public_html/ggplot2.qtl.png", width = 3000)
ggplot(melted, aes(y = as.factor(pop), x = pos, height = lod, fill = pop)) + geom_ridgeline(stat = "identity",
  alpha = 0.5, scale = 0.25) + facet_wrap(~as.factor(chr), scales = "free_x", nrow = 1,
  ncol = 24, strip.position = NULL) + # scale_x_continuous('pos') + scale_y_continuous('lod') +
scale_y_discrete("lod") + theme(axis.text = element_text(size = 10))
dev.off()

# gene expression
png("/home/jmiller1/public_html/ggplot2.expr.qtl.png", width = 3000)
ggplot(melting, aes(x = pos1,y=-logp.DNT, height = logp.DNT)) +
geom_point( aes(color=as.factor(chem)), alpha=0.8, size=1.3) +
facet_wrap(~as.factor(chr), scales = "free_x", nrow = 1,
  ncol = 24, strip.position = NULL) + # scale_x_continuous('pos') + scale_y_continuous('lod') +

scale_x_continuous(label = chr) +
scale_y_continuous(expand = c(0, 0) ))    # remove space between plot area and x axis
dev.off()


# Phenotype figure
##### For Phenotype plots Phenotype plots Phenotype density
pheno.all <- phen <- read.csv("/home/jmiller1/QTL_Map_Raw/popgen/rQTL/metadata/ALL_phenotype_Dist.csv",
  header = T)
phen <- phen[!is.na(phen$pheno_all), ]
phen$Pheno_05 <- phen$pheno_all
# index <- which(phen$pop_all == 'BRP')
phen$pop_all <- factor(phen$pop_all, levels = c("NBH", "BP", "NEW", "ELR"))

png("/home/jmiller1/public_html/phenotypes.png", width = 600)
ggplot(phen, aes(Pheno_05, y = pop_all, fill = pop_all, group = pop_all, height = ..density..,
  show.legend = F, order = as.numeric(pop_all))) +
  geom_density_ridges(scale = 1.33, bandwidth = 0.5) + xlab("Phenotype Score") + ylab("Density") + scale_alpha(guide = "none") +
  scale_y_discrete(limits = rev(levels(phen$pop_all)), expand = c(0.01, 0.01)) +
  scale_x_continuous(expand = c(0.02, 0.02), limits = c(0, 5)) + scale_fill_manual(name = "Populations",
  guide = "legend", values = popcol) + theme(text = element_text(size = 20), axis.text.y = element_blank(),
  panel.background = element_rect(fill = "white", colour = "black", size = 0.1,
    linetype = 1), )
dev.off()


#### Phenotype by genotype plots
outdir <-  '~/QTL_Map_Raw/popgen/rQTL/plots'
plot_pxg





### line for gene expression (each pcb/pah), pcb popgen outliers, AHR gene locations












#
# guess_phase - turn imputed genotypes into phased genotypes along chromosomes
# sim_geno - multiple imputations of underlying genotypes given marker data
# predict_snpgeno - predict SNP genotypes in a multiparent population from
# inferred genotypes plus founder strains’ SNP alleles.

perms.unstrat <- scan1perm(pr, cross2$pheno, model = "binary", cores = 12, n_perm = 200)

par(mar = c(4.1, 4.1, 0.6, 0.6))
plot_pxg(nbh2, iron$pheno[, "liver"], ylab = "Liver phenotype")


pxg1 <- find_marker(nbh2, chr = "2", 198)

png("/home/jmiller1/public_html/Q1.png", width = 300)
plotPXG(nbh2, pheno.col = 2)
