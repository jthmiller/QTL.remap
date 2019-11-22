## NBH

### CHECK UNMAPPED AHR PATHWAY GENES
unmapped <- read.cross.jm(file = file.path(indpops, paste0(pop, ".um.unmapped.f2.csvr")),
format = "csvr", geno = c(1:3), estimate.map = FALSE)


################################################################################
### Pull names from plinkfile
path <- file.path(indpops, paste(pop, ".ped", sep = ""))
popname <- system(paste("cut -f1 -d' '", path), intern = TRUE)
indname <- system(paste("cut -f2 -d' '", path), intern = TRUE)
unmapped$pheno$ID <- paste(popname, indname, sep = "_")
################################################################################

#### PHENO #####################################################################
unmapped$pheno$bin <- ifelse(unmapped$pheno$Pheno > 2, 1 , 0)
unmapped$pheno$pheno_norm <- round(nqrank(unmapped$pheno$Pheno))
################################################################################

arnt2 <- pull.markers(unmapped,'NW_012225110.1:86344')
unmapped <- subset(unmapped,chr=c('NW_012234311.1'))

#### Pvalue and Missing ##############################################
toss.missing <- c("NBH_5525","NBH_6177")
gt <- geno.table(subset(unmapped, ind=!unmapped$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F')))
bfixA <- rownames(gt[which(gt$P.value > 0.0001 & gt$missing < 5),])

###### FILTER #######################################################
unmapped <- pull.markers(unmapped,bfixA)
unmapped <- subset(unmapped,ind=!unmapped$pheno$ID %in% c(toss.missing,'NBH_NBH1M','NBH_NBH1F'))
################################################################################
### no-qtls
scan.norm.mr <- scanone(unmapped, method = "mr", model = "normal", pheno.col = 5)
scan.arnt2 <- scanone(arnt2, method = "mr", model = "normal", pheno.col = 5)
################################################################################

summary(scan.arnt2)


save.image(file.path(mpath,'temp_single_scans.nbh.rsave'))

################################################################################

summary(sc2_normal_imp, thresholds,what="best",
             perms=sc2_normal_imp_perms, alphas=c(0.05,0.1), lodcolumn=1,
             pvalues=FALSE, allpairs=TRUE)

summary(sc2_normal_imp, thresholds,what="full",
             perms=sc2_normal_imp_perms, alphas=c(0.05,0.1), lodcolumn=1,
             pvalues=FALSE, allpairs=TRUE)

summary(sc2_normal_imp, thresholds,what="add",
             perms=sc2_normal_imp_perms, alphas=c(0.05,0.1), lodcolumn=1,
             pvalues=FALSE, allpairs=TRUE)

summary(sc2_normal_imp, thresholds,what="int",
             perms=sc2_normal_imp_perms, alphas=c(0.05,0.1), lodcolumn=1,
             pvalues=FALSE, allpairs=TRUE)

sm <- summary(sc2_normal_imp, thresholds,what=c("best", "full", "add", "int"),
             perms=sc2_normal_imp_perms, alphas=c(0.05,0.1), lodcolumn=1,
             pvalues=FALSE, allpairs=TRUE)

##########################################################################

NBH, Is it arnt?    chr8 16483144 16483822  ARNT
##c8.loc88,8:6520005
pos8 <- as.numeric(gsub("8:","",markernames(cross,chr=8)))


arnt <- 16483144
mark_close_to_arnt <- which.min(abs(pos8 - arnt))
mark_close_to_arnt <- markernames(cross,chr=8)[mark_close_to_arnt]



which(pos8 == 6520005)

png(paste0('~/public_html/NBH_gts8.png'),height=2500,width=4500)
 plotGeno(cross, chr=i, cex=2)
dev.off()


################################################################################
## ex formula=y~Q1*Q4+Q2+Q3

qtl <- makeqtl(cross, chr=c(2,18), pos=c(80,27))
fitted <- fitqtl(cross,qtl=qtl, formula=y~Q1*Q2, pheno.col=5, model='normal',get.ests=T, maxit=10000)

qtl.2 <- makeqtl(cross, chr=c(2), pos=c(80),  what="draws")
out.a.2 <- addqtl(cross, qtl=qtl.2, formula=y~Q1+Q2, method="imp",model="binary",pheno.col=4)
out.i.2 <- addqtl(cross, qtl=qtl.2, formula=y~Q1*Q2, method="imp",model="binary",pheno.col=4)
## indicates possible interaction on 13 and 18

qtl.8 <- makeqtl(cross, chr=c(8), pos=c(88),  what="draws")
out.i.8 <- addqtl(cross, qtl=qtl.8, formula=y~Q1*Q2, method="imp",model="normal",pheno.col=4)
out.a.8 <- addqtl(cross, qtl=qtl.8, formula=y~Q1+Q2, method="imp",model="normal",pheno.col=4)

qtl.18 <- makeqtl(cross, chr=c(18), pos=c(26.7),  what="draws")
out.a.18 <- addqtl(cross, qtl=qtl.18, formula=y~Q1+Q2, method="imp",model="binary",pheno.col=4)
out.i.18 <- addqtl(cross, qtl=qtl.18, formula=y~Q1*Q2, method="imp",model="binary",pheno.col=4)
################################################################################

################################################################################
## manual at stepwise positions
qtl.18 <- makeqtl(cross, chr=c(18), pos=c(27),  what="draws")
out.i.18.b <- addqtl(cross, qtl=qtl.18, formula=y~Q1*Q2, method="imp",model="normal",pheno.col=5)
out.a.18.b <- addqtl(cross, qtl=qtl.18, formula=y~Q1+Q2, method="imp",model="normal",pheno.col=5)

qtl.2 <- makeqtl(cross, chr=c(2), pos=c(80),  what="draws")
out.i.2.b <- addqtl(cross, qtl=qtl.2, formula=y~Q1*Q2, method="imp",model="normal",pheno.col=5)
out.a.2.b <- addqtl(cross, qtl=qtl.2, formula=y~Q1+Q2, method="imp",model="normal",pheno.col=5)

qtl.8 <- makeqtl(cross, chr=c(8), pos=c(67),  what="draws")
out.i.8.b <- addqtl(cross, qtl=qtl.8, formula=y~Q1*Q2, method="imp",model="normal",pheno.col=5)
out.a.8.b <- addqtl(cross, qtl=qtl.8, formula=y~Q1+Q2, method="imp",model="normal",pheno.col=5)
################################################################################

####################################################################################
## Multi/interacting QTL

################################################################################
qtl <- makeqtl(cross, chr=c(2,8,18), pos=c(80,67,27),  what="draws")
fitted <- fitqtl(cross,qtl=qtl, formula=y~Q1+Q2+Q3, pheno.col=4)

qtl <- makeqtl(cross, chr=c(2,8,18), pos=c(80,67,27),  what="draws")
fitted <- fitqtl(cross, qtl=qtl,  formula=y~Q1+Q2+Q3, pheno.col=5, method="imp",model="normal")
fitted.noint.add <- addqtl(cross, qtl=qtl,  formula=y~Q1+Q2+Q3, pheno.col=5, method="imp",model="normal")

qtl <- makeqtl(cross, chr=c(2,18), pos=c(80,27))
fitted.int.only <- fitqtl(cross, qtl=qtl,  formula=y~Q1+Q2+Q1:Q2, pheno.col=5, method="imp",model="normal")
out.i <- addqtl(cross, qtl=qtl, formula=y~Q1+Q2+Q3, pheno.col=5, method="imp",model="normal")

qtl <- makeqtl(cross, chr=c(2,18), pos=c(80,27))
out.ap13 <- addqtl(cross, qtl=qtl, formula = y~Q1 + Q1*Q3 + Q2,  model='normal', method = "imp", pheno.col = 5)
out.ap23 <- addqtl(cross, qtl=qtl, formula = y~Q1 + Q2 + Q2*Q3,  model='normal', method = "imp", pheno.col = 5)

qtl <- makeqtl(cross, chr=c(2,3,13,18,19), pos=c(80,11.18,38.09,27,41.3))
out.ap13 <- fitqtl(cross, qtl=qtl, formula = y~Q1 + Q1*Q2 + Q1*Q3 + Q4 + Q4*Q5,  model='normal', method = "imp", pheno.col = 5)

### 19 interacts w 2 (from picture)
### 18 and 23
####################################################################################
