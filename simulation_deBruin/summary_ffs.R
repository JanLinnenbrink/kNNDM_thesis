library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

t <- do.call(rbind,list.files("CVresults/ffs_knndm/", full.names = TRUE) |> 
          lapply(read.csv) |>
          lapply(function(x)x[,1:16]))

tr <- do.call(rbind,list.files("CVresults/ffs_rand/", full.names = TRUE) |> 
                lapply(read.csv))


t$CV <- "knndm"
tr$CV <- "random"
aoa_df <- t[,c("variate", "smpl", "inAOA", "RMSE_diff_aoa", "RMSE_diff", "CV")] |> 
  rbind(tr[,c("variate", "smpl", "inAOA", "RMSE_diff_aoa", "RMSE_diff", "CV")]) |>
  mutate(design = case_when(smpl == "regular" ~ "syst",
                            smpl == "simpleRandom" ~ "SRS",
                            smpl == "clusterMedium" ~ "clustMed",
                            smpl == "clusterStrong" ~ "clustStr",
                            smpl == "clusterGapped" ~ "clustGap"),
         design = as.factor(design),
         design = fct_relevel(design, c("SRS", "syst", "clustMed","clustStr","clustGap"))) |> 
  pivot_longer(c(RMSE_diff, RMSE_diff_aoa))


(diff_RMSE <- ggplot(aoa_df, aes(x=design, fill=name)) +
  geom_boxplot(aes(y=value), width=0.45, position = position_dodge(.65)) +
  geom_hline(yintercept = 0, linetype="solid", alpha=0.6) +
  scale_fill_discrete(name="", labels=c("whole area", "inside AOA")) +
  ylab("relative RMSE [%]") + 
  xlab("") +
    scale_y_continuous(limits=c(-60,60), breaks=seq(-60,60,20)) +
  facet_grid(CV~variate) +
  theme_bw() +
    theme(legend.position = "right",
          axis.text = element_text(size=8)))

(diff_AOA <- ggplot(aoa_df, aes(x=design, fill=CV)) +
  geom_boxplot(aes(y=inAOA), width=0.45, position = position_dodge(.65)) +
  geom_hline(yintercept = 100) +
    scale_y_continuous(limits = c(0,100)) +
    scale_fill_discrete(name="", labels=c("kNNDM", "NNDM")) +
  facet_grid(""~variate) +
    ylab("cells inside AOA [%]") + 
    xlab("") +
  theme_bw()+ 
    theme(legend.position = "right",
          axis.text = element_text(size=8)))

#ggsave("material/rmse_aoa_diff.pdf", diff_RMSE, width=6.5, height = 4.5)
#ggsave("material/aoa.pdf", diff_AOA, width=6.5, height = 3)


pl <- cowplot::plot_grid(diff_RMSE, diff_AOA, nrow = 2, rel_heights = c(1,0.6),
                   align = "v", labels = c("a)", "b)"))


ggsave("material/rmse_aoa_db.pdf", pl, width=7, height = 6.5)


ggplot(tr) +
  geom_histogram(aes(x=xy_imp)) +
  facet_wrap(~smpl)
   

ggplot(tr) +
  geom_histogram(aes(x=RMSE_gains_ffs)) +
  geom_vline(aes(xintercept=0)) +
  facet_wrap(~smpl)            


ggplot(tr) +
  geom_histogram(aes(x=AOA_gains_ffs)) +
  geom_vline(aes(xintercept=0)) +
  scale_x_continuous(limits = c(-20,20)) +
  facet_wrap(~smpl) 

names(tr)
(diff_RMSE <- ggplot(tr, aes(x=smpl)) +
    geom_boxplot(aes(y=RMSE_diff_caoa_aoa), width=0.45, position = position_dodge(.65)) +
    geom_hline(yintercept = 0, linetype="solid", alpha=0.6) +
    ylab("relative RMSE [%]") + 
    xlab("") +
    scale_y_continuous(limits=c(-60,60), breaks=seq(-60,60,20)) +
    facet_grid(CV~variate) +
    theme_bw() +
    theme(legend.position = "right",
          axis.text = element_text(size=8)))
