library(ggplot2)
infolder  <- "./CVresults/"
outfolder <- "./material"
mets <- c("exhaustive_old", "random", "spatial", "intensity","heteroscedastic", "modelbased", "knndm","nndm")
colnms <- c("method", "variate", "design", "number", "RMSE")
outtab <- data.frame(matrix(NA, 0, 5))
names(outtab) <- colnms

for(m in mets){
  p <- file.path(infolder, m)
  f_ins <- list.files(p, glob2rx("???_*.Rdata"))
  for(f_in in f_ins){
    lchar <- nchar(f_in)
    variate <- substr(f_in, 1, 3)
    design <- substr(f_in, 5, lchar-9)
    number <- as.numeric(substr(f_in, lchar-8, lchar-6))
    load(file.path(p, f_in))
    RMSE   <- mean(RMSE)
    if(m == "nndm") MEC = 0
    else MEC <- mean(MEC)
    
    newrow <- data.frame(method = m, variate = variate, design = design,
                         number = number, RMSE = RMSE, MEC=MEC)
    outtab <- rbind(outtab, newrow)
  }
}

write.csv(outtab, file.path(outfolder, "outtab100_ranger.csv"))


outtab <- read.csv(file.path(outfolder, "outtab100_ranger.csv"))
outtab$methodID <- with(outtab, ifelse(method=="exhaustive_old",0,1))

outtab[outtab$methodID==0,]
# relative RMSE & MEC
outtab$rRMSE=outtab$rMEC <- NA

# 86, 88, 93
numbers=c(1:85, 87, 89:92, 94:100)

for(variate in c("AGB", "OCS")){
  for(design in unique(outtab$design)){
    for(number in numbers){
      idx1 <- which(outtab$design == design & outtab$variate == variate & 
                      outtab$number == number)
      idx2 <- which(outtab$design == design & outtab$variate == variate & 
                      outtab$methodID == 0 & outtab$number == number)
      
      outtab$rRMSE[idx1] <- 100 * (outtab$RMSE[idx1] - outtab$RMSE[idx2])/
        outtab$RMSE[idx2]
      outtab$rMEC[idx1] <- 100 * (outtab$MEC[idx1] - outtab$MEC[idx2])/
        outtab$MEC[idx2]
    }
  }
}


outtab$method <- factor(outtab$method, levels=c("exhaustive_old", "random", "spatial", "intensity",
                                                "heteroscedastic", "modelbased", "knndm","nndm"))
outtab$design <- factor(outtab$design, levels=c("simpleRandom", "regular", "clusterMedium", 
                                                "clusterStrong", "clusterGapped"))
outtab$variate <- as.factor(outtab$variate)
xlabs <- c("conventional", "spatial", "intensity", 
           "kNNDM", "NNDM")
collabs <- c("SRS","syst","clustMed","clustStr","clustGap")
cols <- c("white", "brown", "pink", "orange", "lightgreen")

(rmse <- ggplot(outtab[!outtab$method%in%c("exhaustive_old", "heteroscedastic", "modelbased"),], 
                aes(x=method, y=rRMSE, fill=design)) +
    geom_boxplot(linetype = "dashed", outlier.shape = NA) +
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) +
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) +
    stat_summary(fun.y=mean, geom="point", shape=1, size=1.5, stroke = 1, colour="black",
                 position = position_dodge2(width = 0.75,   
                                            preserve = "single")) +
    geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5),  color="grey60", linetype="dashed", alpha=0.6) +
    geom_hline(yintercept = 0, linetype="solid", alpha=0.6) +
    scale_fill_manual(values=cols, labels=collabs) +
    scale_x_discrete(labels=xlabs) +
    scale_y_continuous(limits=c(-60,60), breaks=seq(-60,60,20)) +
    ylab("relative RMSE [%]") + 
    xlab("") +
    theme_classic() +
    theme(strip.text=element_text(hjust = 0.55,  # hjust = 0.5 centers the title
                                  size = 14,
                                  face = "bold"),
          panel.border = element_rect(linetype = "solid",
                                      colour = "black", fill = "NA", size = 0.5),
          legend.position = NaN) +
    facet_wrap(~variate)  )

ggsave("material/rmse_ranger.pdf")


(mec <- ggplot(outtab[!outtab$method%in%c("exhaustive_old", "heteroscedastic", "modelbased"),], aes(x=method, y=rMEC, fill=design)) +
    geom_boxplot(linetype = "dashed", outlier.shape = NA) +
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..)) +
    stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..)) +
    stat_summary(fun.y=mean, geom="point", shape=1, size=1.5, colour="black",
                 position = position_dodge2(width = 0.75,   
                                            preserve = "single")) +
    geom_vline(xintercept=c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5),  color="grey60", linetype="dashed", alpha=0.6) +
    geom_hline(yintercept = 0, linetype="solid", alpha=0.6) +
    scale_fill_manual(values=cols, labels=collabs) +
    scale_x_discrete(labels=xlabs) +
    scale_y_continuous(limits=c(-40,40), breaks=seq(-40,40,20)) +
    scale_x_discrete(labels=xlabs) + 
    xlab("") +
    ylab("relative MEC [%]") + 
    theme_classic() +
    theme(strip.text=element_text(hjust = 0.5,  # hjust = 0.5 centers the title
                                  size = 14,
                                  face = "bold"),
          panel.border = element_rect(linetype = "solid",
                                      colour = "black", fill = "NA", size = 0.5),legend.box.just = "left",
          legend.justification = c(0,0),
          legend.position = c(0,0),
          legend.text.align = 0,
          legend.title = element_blank(),
          legend.key.size = unit(0.8, "lines"),
          legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')) +
    facet_wrap(~variate))

bpl <- gridExtra::grid.arrange(rmse, mec)

ggsave(file.path(outfolder, "comp_both_wo.pdf"), bpl, height = unit(8, "cm"), width=unit(12, "cm"))
