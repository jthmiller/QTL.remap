### NBH
a <- find.marker(cross.18, 2, 180)
b <- find.marker(cross.18, 18, 110)


tb <- getem(a, b, 5, 1)
tgc <- summarySE(tb, measurevar = "pt1", groupvars = c("po1", "po2"))
tb$pt3 <- as.numeric(as.factor(tb$pt2))
set.seed(666)
p <- ggplot(tb, aes(x = po1, y = pt1))
p <- p + geom_point(aes(x = po1 - 0.2, y = pt1, group = po2, color = po2), size = 3.5, 
  position = position_jitterdodge(dodge.width = 0.25, jitter.height = 0.15, jitter.width = 0)) + 
  
scale_y_continuous(breaks = c(0, 1), limits = c(-0.25, 1.25), labels = NULL) + scale_x_continuous(breaks = c(1, 
  2, 3), limits = c(0.66, 3.33), labels = c(`1` = "AA", `2` = "AB", `3` = "BB")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor.y = element_line(colour = "white", 
    size = 0.75), panel.grid.minor.x = element_line(colour = "white", size = 0.75), 
    axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), 
    legend.position = "none", axis.text.x = element_text(face = "bold", color = "black", 
      size = 16)) + geom_vline(xintercept = c(1.5, 2.5), colour = "white", 
  size = 0.75)
# p <- p + geom_errorbar(data=tgc, aes(ymin=pt1-se, ymax=pt1+se,color=po2),
# width=.1,size=1.3)
p <- p + stat_summary(aes(y = pt1, group = po2, color = po2), fun.y = mean, geom = "line", 
  size = 1.25)
p <- p + stat_summary(aes(y = pt1, group = po2, color = po2), size = 10, pch = 18, 
  fun.y = mean, geom = "point", )
pdf("/home/jmiller1/public_html/NBH_2_18.pdf")
# pdf('/home/jmiller1/public_html/NEW_2_18.pdf')
# pdf('/home/jmiller1/public_html/ELR_8_13.pdf')
p
dev.off()

tb1 <- getone(a, b, 5, 1)
tgc1 <- summarySE(tb1, measurevar = "pt1", groupvars = c("po1"))
p <- ggplot(tb1, aes(x = po1, y = pt1))
p <- p + geom_jitter(aes(x = po1 - 0.2, y = pt1, color = po2, size = 1.5), height = 0.1, 
  width = 0.075)
p <- p + scale_y_continuous(breaks = c(0, 1), limits = c(-0.25, 1.25), labels = NULL) + 
  # scale_color_manual(values = po1, limits = c('1' ,'2' ,'3'))+
scale_x_continuous(breaks = c(1, 2, 3), limits = c(0.66, 3.33), labels = c(`1` = "AA", 
  `2` = "AB", `3` = "BB")) + theme(panel.grid.major = element_blank(), panel.grid.minor.y = element_line(colour = "white", 
  size = 0.75), panel.grid.minor.x = element_line(colour = "white", size = 0.75), 
  axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), 
  legend.position = "none", axis.text.x = element_text(face = "bold", color = "black", 
    size = 16)) + geom_vline(xintercept = c(1.5, 2.5), colour = "white", size = 0.75)
p <- p + stat_summary(aes(y = pt1), fun.y = mean, geom = "line", size = 1.25)
# p <- p + geom_errorbar(data=tgc1, aes(ymin=pt1-se, ymax=pt1+se,
# width=.1,size=1.3))

# p <- p + stat_summary(aes(y = pt1),cex=1.2,pch=2, fun.y=mean, geom='point')
pdf("/home/jmiller1/public_html/NBH_2.pdf")
# pdf('/home/jmiller1/public_html/NEW_2.pdf')
p
dev.off()

tb1 <- getone(b, a, 5, 1)
tgc1 <- summarySE(tb1, measurevar = "pt1", groupvars = c("po1"))
p <- ggplot(tb1, aes(x = po1, y = pt1))
p <- p + geom_jitter(aes(x = po1 - 0.2, y = pt1, color = po1, size = 1.5), height = 0.1, 
  width = 0.075)
p <- p + scale_y_continuous(breaks = c(0, 1), limits = c(-0.25, 1.25), labels = NULL) + 
  # scale_color_manual(values = po1, limits = c('1' ,'2' ,'3'))+
scale_x_continuous(breaks = c(1, 2, 3), limits = c(0.66, 3.33), labels = c(`1` = "AA", 
  `2` = "AB", `3` = "BB")) + theme(panel.grid.major = element_blank(), panel.grid.minor.y = element_line(colour = "white", 
  size = 0.75), panel.grid.minor.x = element_line(colour = "white", size = 0.75), 
  axis.title.x = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), 
  legend.position = "none", axis.text.x = element_text(face = "bold", color = "black", 
    size = 16)) + geom_vline(xintercept = c(1.5, 2.5), colour = "white", size = 0.75)
p <- p + stat_summary(aes(y = pt1), fun.y = mean, geom = "line", size = 1.25)
# p <- p + geom_errorbar(data=tgc1, aes(ymin=pt1-se, ymax=pt1+se,
# width=.1,size=1.3)) p <- p + stat_summary(aes(y = pt1),cex=1.2,pch=2,
# fun.y=mean, geom='point')
pdf("/home/jmiller1/public_html/NBH_18.pdf")
# pdf('/home/jmiller1/public_html/NEW_18.pdf')
p
dev.off()
