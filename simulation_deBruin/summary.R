library(ggplot2)
infolder  <- "./CVresults/"
outfolder <- "./material"
mets <- c("exhaustive_test", "random_caret_new", "nndm_caret","knndm_caret")
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
    
    if(m %in% c("knndm_caret", "nndm_caret")) {
      time_mod <- mean(time_mod)
      time_alg <- mean(time)
      time_overall <- time_mod + time_alg
    } else {
      time_mod=time_alg=time_overall <- NaN
    }
      
    newrow <- data.frame(method = m, variate = variate, design = design,
                         number = number, RMSE = RMSE,
                         time_mod=time_mod, time_alg=time_alg, time_overall=time_overall)
    outtab <- rbind(outtab, newrow)
  }
}

write.csv(outtab, file.path(outfolder, "outtab100.csv"))


outtab <- read.csv(file.path(outfolder, "outtab100.csv"))
outtab$methodID <- with(outtab, ifelse(method=="exhaustive_test",0,1))
                                       
                                       
# relative RMSE & MEC
outtab$rRMSE <- NA
#outtab$rMEC  <- NA
library(dplyr)
outtab <- outtab |> 
  group_by(method, design, variate) |> 
  slice(1:50)


numbers=1:50
#18
 

for(variate in c("AGB", "OCS")){
  for(design in unique(outtab$design)){
    for(number in numbers){
      idx1 <- which(outtab$design == design & outtab$variate == variate & 
                      outtab$number == number)
      idx2 <- which(outtab$design == design & outtab$variate == variate & 
                      outtab$methodID == 0 & outtab$number == number)
      
      outtab$rRMSE[idx1] <- 100 * (outtab$RMSE[idx1] - outtab$RMSE[idx2])/
        outtab$RMSE[idx2]
      #outtab$rMEC[idx1] <- 100 * (outtab$MEC[idx1] - outtab$MEC[idx2])/
      #  outtab$MEC[idx2]
    }
  }
}

outtab$method <- factor(outtab$method, levels=c("exhaustive_test", "random_caret_new",
                                                "nndm_caret","knndm_caret"))
outtab$design <- factor(outtab$design, levels=c("simpleRandom", "regular", "clusterMedium", 
                                                "clusterStrong", "clusterGapped"))
outtab$variate <- as.factor(outtab$variate)
xlabs <- c("conventional", "kNNDM", "modelbased")
collabs <- c("SRS","syst","clustMed","clustStr","clustGap")
cols <- c("white", "brown", "pink", "orange", "lightgreen")
outtab[!outtab$method%in%c("exhaustive_test")&outtab$variate=="OCS",]

(rmse <- ggplot(outtab[outtab$method%in%c("random_caret_new", "knndm_caret",
                                          "nndm_caret", "spatial"),], 
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
  #scale_x_discrete(labels=xlabs) +
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

ggsave("material/rmse.pdf")


(mec <- ggplot(outtab[outtab$methodID!=0,], aes(x=method, y=rMEC, fill=design)) +
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

ggsave(file.path(outfolder, "comp_global.pdf"), bpl, height = unit(8, "cm"), width=unit(12, "cm"))


## density
# global validation results in under/overestimation of map error in random - regualar cases
(pl <- ggplot(data=outtab[outtab$method %in% c("nndm_caret", "knndm_caret"),], aes(fill=method)) +
    geom_density(aes(y=rRMSE), alpha=0.5) +
    geom_hline(yintercept = 0) +
    scale_y_continuous(limits=c(-60,60)) +
    xlab("") +
    ylab("relative RMSE [%]") + 
    theme_bw() +
    facet_wrap(~design))
ggsave(paste0(outfolder, "rrmse_density10.pdf"),pl)


# time
time_stats_all <- outtab |> 
  group_by(method, design) |> 
  summarise(mean_alg=mean(time_alg),
            mean_mod=mean(time_mod),
            sd_overall=sd(time_overall)) |> 
  na.omit()


tp <- ggplot(time_stats_all, aes(x=method)) +
  geom_col(aes(y=mean_alg+mean_mod, fill="model"), width=0.65) +
  geom_col(aes(y=mean_alg, fill="algorithm"), width=0.65) +
  geom_errorbar(aes(ymin = mean_alg+mean_mod - sd_overall, 
                    ymax = mean_alg+mean_mod + sd_overall), width = 0.2) +
  xlab("") + ylab("time [s]") + labs(fill="mean time for ") +
  theme_gray(base_size=16) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size=9)) +
  scale_x_discrete(labels=c("NNDM LOOCV", "kNNDM CV")) +
  scale_y_log10() +
  facet_wrap(~design, nrow=1)

tp
ggsave("material/time.pdf", tp, height=5, width=10)





# ffs ---------------------------------------------------------------------


read.csv("CVresults/ffs_knndm/ffs_knndm/AGB_clusterGapped004.csv")
warnings()

