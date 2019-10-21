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

qtl8 is near
241698-
8080728
6754602 - 6036644
1628405 - 1389278



                  8755941  8804728   154311   154443   154615  8414376
9016170  9016269  9016366  9016106  8960809  8941638  8804828  8804713
8708306  8553028  8414541  8414523  8414600  8414352   241698   241704
 241845   602168   602184   602367   602357   602223   602582   683297
 683369   683337   683360  1046214  1046553  1046230  1046506  1046612
1046324  1046647  1126640  8028077  8028158  8027976  8027911  8027806
8027830  8015744  1127062  1195402  1134997  1195417  1195640  1201903
1195454  1195370  1127109  1126769  1126708  1126656  1126738  8015594
7598274  7598560  7598637  7598213  7598479  7226687  7370622  7370792
6839989  6840336  6845477  6901975  6845454  6845816  6856802  6902045
6965615  6965661  6983706  7370806  7370853  7370928  8205513  8205468
8203472  8203629  8205462  8116921  8205241  8116914  8116703  8116891
8116721  8116840  8116622  8099036  8098709  8098919  8088477  8088497
8098939  8099127  8099086  8116683  8099042  8098881  8088462  8080764
8088448  8088356  8088429  8080728  6754602  6754585  6754636  6754573
6749298  6744037  6749644  6744378  6743962  6740741  6714294  6714061
6714285  6714243  6714048  6659176  6630644  6520269  6520005  6576731
6576898  6576618  6576451  6069866  6069553  6069645  6069531  6036774
6069612  6069670  6069650  6069598  6036730  6036659  6036644  1628405

 CODE/
 MAP2/
 MAP2/


################################################################################
## ex formula=y~Q1*Q4+Q2+Q3

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
