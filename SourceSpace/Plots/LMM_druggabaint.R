setwd("/Users/alistairperry/Documents/Cambridge/Project/Analysis/LFPs/newmaxfilter_icafixes/LFPs_COH_wrad_allMEGch_wfids/C_Plots/DrugAnalysis/MRS")

library("tidyr")

library("lmerTest")

library("ggplot2")

library("dplyr")

library("afex")

library("BayesFactor")

library("bayestestR")

library("interactions")


FigOutDir  <- "/Users/alistairperry/Documents/Cambridge/Project/Analysis/LFPs/newmaxfilter_icafixes/LFPs_COH_wrad_allMEGch_wfids/C_Plots/DrugAnalysis/MRS/Scatt/"
dir.create(FigOutDir)


FigOutPrefix <- "allsubjs"


# Right hemi only

ROInames <- c("RIFG", "RSTG", "RAUD")

Labnames <- c("R IFG", "R STG", "R AUD")


# width and height variables for saved plots - leave it as default for now

w = 6
h = 4



for (i in seq(1,3)) {

  
datfname <- paste("adv_ssst_newmaxf_fixICA_wfids_nodipcor_nop10p19_FIX_LFPs_MRSIFGgabacorassoc_meanMMN3_DrugDiff_Pats_fullsample_stats_LMMtableforR_", ROInames[i], ".txt", sep="")   

LFPdata_wide <- read.table(datfname, header = TRUE)

LFPdata_long <- gather(LFPdata_wide, drug, MMN, PLA:MEM, factor_key = TRUE)


LFPdata_long$Subject <- factor(LFPdata_long$Subject)


LFPdata_long$MMN_z <- scale(LFPdata_long$MMN)

LFPdata_long$GABA_z <- scale(LFPdata_long$GABA)


contrasts(LFPdata_long$drug) <- contr.sum

mixed.lmer <- lmer(MMN_z ~ GABA_z*drug + (1|Subject), data = LFPdata_long)

lme_res <- anova(mixed.lmer)


#Write table - should write to current dir

outfname <- paste("LMEres_druggabaint_", ROInames[i],".txt", sep="")
write.table(lme_res, file = outfname)


#Change Plots

LFPdata_wide$DrugDiff <- LFPdata_wide$PLA-LFPdata_wide$MEM

data <- LFPdata_long
diff_data <- LFPdata_wide

# panel A
MMN_drug_lc_plot <- ggplot(data = data, aes(x = GABA, y = MMN)) +
  
  # gray vertical lines to connect observations from same subject
  geom_segment(data = diff_data, aes(x = GABA, y = PLA,
                                     xend = GABA, yend = MEM),
               size = 0.75, colour = "grey50") +
  # dark gray arrows pointing from placebo to atomoxetine
  geom_segment(data = diff_data %>%
                 # only show arrows for participants with meaningful drug effect
                 dplyr::filter(abs(DrugDiff) > .005),
               aes(x = GABA, y = PLA,
                   xend = GABA, yend = MEM),
               arrow = arrow(length = unit(0.175, "inches"), angle = 30),
               lineend = "round", linejoin = "mitre",
               size = .75, colour = "grey50") +
  
  # dots, filled by drug condition
  geom_point(aes(fill = drug), shape = 21, colour = "black", size = 2.5, stroke = 1.2) +
  # trend lines separately for each drug condition
  geom_smooth(aes(colour = drug), method = "lm", se = FALSE, size = 1.25) +
  # specify drug condition colours
  scale_fill_manual(
    limits = c("MEM", "PLA"),
    values = c(`MEM` = "#E41A1C", `PLA` = "#377EB8"),
    labels = c(`MEM` = "MEM", `PLA` = "PLA")
  ) +
  scale_colour_manual(
    limits = c("MEM", "PLA"),
    values = c(`MEM` = "#E41A1C", `PLA` = "#377EB8"),
    labels = c(`MEM` = "MEM", `PLA` = "PLA")
  ) +
  # final tweaks
  # ylims here will best fit R AUD
  coord_cartesian(xlim = c(1.25, 2.75), ylim = c(-0.31, 0.11), expand = FALSE) +
  labs(x = "GABA rIFG (cor.)",
       y = paste("MMN", Labnames[i], sep = " ")) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 15, colour = "black"),
    #panel.grid.major = ggplot2::element_line(colour = "#e0e0e0", size = 0.25),
    #panel.grid.minor = ggplot2::element_blank(),
    #axis.ticks = ggplot2::element_line(color = "black", size = 17, face = "bold"),
    #axis.line = ggplot2::element_line(size = 1, colour = "#cccccc"),
    axis.text = ggplot2::element_text(color = "black", size = 17, face = "bold"),
    axis.title = ggplot2::element_text(color = "black", size = 19, face = "bold"),
    #plot.margin = ggplot2::margin(t = 10, b = 10, l = 10, r = 15)   
  ) 


MMN_drug_lc_plot

figfname<-paste(FigOutDir, "LFP_DrugGABAInt_", ROInames[i], ".tiff", sep="")

ggsave(figfname, width = w, height = h)


#Bayes Analysis

# `generalTestBF` computes the Bayes Factor for all possible restrictions of the below model formula.
# The random intercept by subject is always included in all model variants.
# Note the change in formula syntax for random effects.
mixed.lmer_BF_alt <- BayesFactor::generalTestBF(
  MMN_z ~ GABA_z*drug + Subject, data = LFPdata_long,
  whichRandom = "Subject", neverExclude = "Subject")
# For comparison, compute the Bayes Factor of the null model, with only a random intercept by subject
mixed.lmer_BF_null <- BayesFactor::lmBF(
  MMN_z ~ Subject, data = LFPdata_long, whichRandom = "Subject")
# The final Bayes Factors for the alternative models are relative to the null model
mixed.lmer_BF_final <- mixed.lmer_BF_alt / mixed.lmer_BF_null
# Compute Bayes Factor for each *effect*, i.e. the Inclusion Bayes Factors

bayestestR::bayesfactor_inclusion(mixed.lmer_BF_final, match_models = TRUE)


#Save output to file

outfname <- paste("LMEres_druggabaint_", ROInames[i], "_BF.txt")
write.table(bayestestR::bayesfactor_inclusion(mixed.lmer_BF_final, match_models = TRUE), file = outfname) 


#Correlations
#Make label sizes bigger for conference

diff_data$Diag <- as.factor(diff_data$Diag)

fig <- ggplot(diff_data, aes(GABA, DrugDiff, col = Diag)) + 
  
  geom_point(size = 4, color = "red", fill = "red", aes(shape=Diag)) +
  
  scale_shape_manual(values=c(22, 24)) +
  
  #scale_color_manual(values = c('#E69F00', "#00BA38")) +
  
  geom_smooth(method = lm, size = 2, color = "black") + 
  
  
  labs(y = expression(bold(paste(Delta, ' MMN (PLA-MEM)')))) + 
  
  theme_classic() + 
  
  theme(axis.text.x = element_text(color = "black", size = 17, face = "bold"), axis.text.y = element_text(color = "black", size = 17, face = "bold"), axis.title.x = element_text(color = "black", size = 19, face = "bold"), axis.title.y = element_text(color = "black", size = 19, face = "bold")) +
  
  xlab("GABA rIFG (cor.)") +
  
  
  #Turn caption off
  
  theme(legend.position = "none")


fig

figfname<-paste(FigOutDir, "LFP_DrugGABAInt_", ROInames[i], "_scat.tiff", sep="")

ggsave(figfname, width = w, height = h)



#Now with subgroup fits

diff_data$Diag <- as.factor(diff_data$Diag)

fig <- ggplot(diff_data, aes(x = GABA, y = DrugDiff, col = Diag)) + 
  
  geom_point(size = 4) + 
  
  scale_color_manual(values = c('#E69F00', "#00BA38")) +
  
  geom_smooth(method = lm, size = 2, color = "black") + 
  
  geom_smooth(method = lm, size = 2, fill = NA) + 
  
  labs(y = expression(bold(paste(Delta, ' MMN (PLA-MEM)')))) + 
  
  theme_classic() + 
  
  theme(axis.text.x = element_text(color = "black", size = 17, face = "bold"), axis.text.y = element_text(color = "black", size = 17, face = "bold"), axis.title.x = element_text(color = "black", size = 19, face = "bold"), axis.title.y = element_text(color = "black", size = 19, face = "bold")) +
  
  xlab("GABA rIFG (cor.)") +
  
  
  #Turn caption off
  
  theme(legend.position = "none")


fig

figfname<-paste(FigOutDir, "LFP_DrugGABAInt_", ROInames[i], "_scat_wsubgrpfit.tiff", sep="")

ggsave(figfname, width = w, height = h)


}


# #Supplementary: Mediating effect of Age ---------------------------------

#For info on plotting continous moderator, see: https://cran.r-project.org/web/packages/interactions/vignettes/interactions.html#plotting-observed-data


for (i in seq(1,3)) {
  
  
  #Setup and load ROI text file
  
  datfname <- paste("adv_ssst_newmaxf_fixICA_wfids_nodipcor_nop10p19_FIX_LFPs_MRSIFGgabacorassoc_meanMMN3_DrugDiff_Pats_fullsample_stats_LMMtableforR_", ROInames[i], ".txt", sep="")   
  
  LFPdata_wide <- read.table(datfname, header = TRUE)
  
  diff_data <- LFPdata_wide
  

diff_data$Diag <- as.factor(diff_data$Diag)


#Linear model fitting interaction

fiti <- lm(Diff ~ GABA * Age, data = diff_data)


#Main interaction plot

fig <- interact_plot(fiti, pred = GABA, modx = Age, plot.points = TRUE, x.label = "GABA rIFG (cor.)", y.label = expression(bold(paste(Delta, ' MMN (PLA-MEM) R AUD'))), colors = "red", point.size = 4) + 
  
  
  #Manually edit text labels
  
  theme_classic() +
  
  theme(axis.text.x = element_text(color = "black", size = 17, face = "bold"), axis.text.y = element_text(color = "black", size = 17, face = "bold"), axis.title.x = element_text(color = "black", size = 19, face = "bold"), axis.title.y = element_text(color = "black", size = 19, face = "bold")) +
  
  
  #Legend features
  
  theme(legend.key.size = unit(1, 'cm')) +
  
  theme(legend.text = element_text(size=12, color = "black", face = "bold"), legend.title = element_text(size=14, color = "black", face = "bold"))


fig

figfname<-paste(FigOutDir, "LFP_DrugGABAInt_", ROInames[i], "_scat_AgeMod.tiff", sep="")

ggsave(figfname, width = w, height = h)


}


#Now E/I balance, just for RAUD

LFPdata_wide <- read.table("adv_ssst_newmaxf_fixICA_wfids_nodipcor_remo_nop10p19_LFPs_MRSIFGglugabaratassoc_meanMMN3_DrugDiff_Pats_fullsample_stats_LMMtableforR_RAud.txt", header = TRUE) 

LFPdata_wide$Diag <- as.factor(LFPdata_wide$Diag)


fig <- ggplot(LFPdata_wide, aes(GluGABArat, Diff, color = GABA)) + 
  
  geom_smooth(method = lm, size = 2, color = "black") +
  
  geom_point(size = 4, aes(color = GABA, shape=Diag)) +
               
  scale_shape_manual(values=c(15, 17)) +
  
  scale_colour_gradient(low = "#fc9272", high = "#cb181d") +
  
  #scale_colour_brewer(palette = "Reds") +
  
  
  labs(y = expression(bold(paste(Delta, ' MMN (PLA-MEM) R AUD')))) + 
  
  theme_classic() + 
  
  theme(axis.text.x = element_text(color = "black", size = 17, face = "bold"), axis.text.y = element_text(color = "black", size = 17, face = "bold"), axis.title.x = element_text(color = "black", size = 19, face = "bold"), axis.title.y = element_text(color = "black", size = 19, face = "bold")) +
  
  xlab("Glu/GABA rIFG") +
  

  #Turn caption off
  
  theme(legend.position = "none")


  fig

  figfname<-paste(FigOutDir, "LFP_DrugEIInt_", "RAUD", "_scat.tiff", sep="")

  ggsave(figfname, width = w, height = h)


  
  #Now including controls
  
  #ANCOVA
  
  
  LFPdata_wide <- read.table("adv_ssst_newmaxf_fixICA_wfids_nodipcor_remo_nop10p19_LFPs_MRSIFGglugabaratassoc_meanMMN3_DrugDiff_ConsandPats_fullsample_stats_LMMtableforR_RAud.txt", header = TRUE)
  
  LFPdata_wide$Group <- as.factor(LFPdata_wide$Group)
  
  
  EIDrugDiffGrpInt <- aov(Diff ~ GluGABArat*Group, data = LFPdata_wide)
  
  
  sink("ANCOVA_drugglugabaratGroupInt_RAUD.txt")
  
  summary(EIDrugDiffGrpInt)
  
  sink()
  
  
  
  #Lastly group means
  
  LFPdata_wide <- read.table("adv_ssst_newmaxf_fixICA_wfids_nodipcor_remo_nop10p19_LFPs_MRSIFGglugabaratassoc_ConsandPats.txt", header = TRUE)
  
  
  
  CN <- subset(LFPdata_wide, Group == 1)
  
  Pat <- subset(LFPdata_wide, Group == 2)
  

  
  CNdat <- CN$GluGABArat
  Patdat <- Pat$GluGABArat
  
  ncnt <- length(CNdat)
  npat <- length(Patdat)
  
  totalscans <- ncnt+npat
  
  
  x1 <- as.integer(totalscans)
  x1[1:ncnt] <- rep(1, each=ncnt)
  x1[(ncnt+1):totalscans] <- rep(2, each=npat)
  
  
  d <- data.frame(y = c(CNdat, Patdat), x = x1, PatSubgrp = LFPdata_wide$PatSubgrp)
  
  
  #Start plotting
  
  set.seed(321)
  
  #xj <- jitter(x1, amount = .09)
  xj <- jitter(x1)
  
  d$xj <- xj
  
  d$x <- as.factor(d$x)
  
  d$PatSubgrp <- as.factor(d$PatSubgrp)
  
  my_sum <- d %>%
    group_by(x) %>% 
    dplyr::summarise( 
      n=n(),
      mean=mean(y),
      sd=sd(y)
    ) %>%
    mutate( se=sd/sqrt(n))
  
  
  my_sum_box <- d %>% 
    group_by(x) %>% 
    dplyr::summarise(
      y0 = quantile(y, 0.05), 
      y25 = quantile(y, 0.25), 
      y50 = mean(y), 
      y75 = quantile(y, 0.75), 
      y100 = quantile(y, 0.95))
  
  #Set axis limits
  
  #y_lim_min <- -2
  #y_lim_max <- 0.5
  
  #coord_cartesian(ylim=c(y_lim_min, y_lim_max)) 
  
  
  #Group means
  
  #Start
  
  
  f1 <- ggplot(my_sum, aes(fill=x, y=mean, x=x)) + 
    
    
    
    geom_bar(position="dodge", stat="identity", alpha = 0.7, width = 0.75, color = "black", size = 0.75) +
    
    
    
    scale_fill_manual(values = c("Blue", "Red")) +
    
    geom_errorbar(aes(x=x, ymin=mean, ymax=mean+se), width=0.4, colour="black", alpha=0.9, size=0.75) +
    
    
    geom_point(data = d %>% filter(x =="1"), aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .6) +
    
    geom_point(data = d %>% filter(x =="2"), aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .6) + 
    
    scale_x_discrete(breaks=c(1,2), labels=c("CON", "FTLD")) +
    
    ylab("Glu/GABA rIFG") +
    
    xlab("") +
    
    theme_classic() +
    
    
    
    theme(plot.title = element_text(face="bold", size = 18, hjust = 0.5, color = "black"), axis.title.y = element_text(face="bold", size = 19)) + 
    
    theme(axis.text.x = element_text(face="bold", size = 19, color = "black"), axis.text.y = element_text(face="bold", size = 18, color = "black")) +
    
    
    #Turn caption off
    
    theme(legend.position = "none")
  
  
  
  f1
  
  
  figfname<-paste(FigOutDir, "GroupGluGABAratmean.tiff", sep="")
  
  ggsave(figfname, width = w, height = h)
  
  
  #Boxplot
  
  
  f1_box <- ggplot(my_sum_box, aes(fill=x, x=x, y = y50)) +
    
    #stat_boxplot(geom ='errorbar', width = 0.6) +
    
    geom_boxplot(aes(ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100),
                 position="dodge", stat="identity", alpha = 0.7, width = 0.75, color = "black", size = 0.75) + 
    
    scale_fill_manual(values = c("Blue", "Red")) + 
    
    geom_point(data = d %>% filter(x =="1"), aes(x = xj, y = y), color = 'black', size = 3, alpha = .8) +
    
    geom_point(data = d %>% filter(x =="2"), aes(x = xj, y = y, shape = PatSubgrp), color = 'black', fill = 'black', size = 3, alpha = .8) + 
    
    
    scale_shape_manual(values=c(22, 24)) +
    
    #scale_color_manual(values = c('#E69F00', "#00BA38")) +
    
    scale_x_discrete(breaks=c(1,2), labels=c("CON", "bvFTD/PSP")) +
    
    ylab("Glu/GABA rIFG") +
    
    xlab("") +
    
    theme_classic() +
    
    
    
    theme(plot.title = element_text(face="bold", size = 18, hjust = 0.5, color = "black"), axis.title.y = element_text(face="bold", size = 19)) + 
    
    theme(axis.text.x = element_text(face="bold", size = 19, color = "black"), axis.text.y = element_text(face="bold", size = 18, color = "black")) +
    
    
    #Turn caption off
    
    theme(legend.position = "none")
  
  
  
  f1_box
  
  
  figfname<-paste(FigOutDir, "GroupGluGABAratmean_box.tiff", sep="")
  
  ggsave(figfname, width = w, height = h)
  
  
  
  #Rep6
  
  #Correlations

  
  for (i in seq(1,3)) {
    
    
    datfname <- paste("adv_ssst_newmaxf_fixICA_wfids_nodipcor_nop10p19_FIX_LFPs_MRSIFGgabacorassoc_meanMMN_DrugDiff_Pats_fullsample_stats_LMMtableforR_", ROInames[i], ".txt", sep="")   
    
    
    #Not sure why of the multiple steps but oh well..
    
    
    LFPdata_wide <- read.table(datfname, header = TRUE)
    
    LFPdata_wide$DrugDiff <- LFPdata_wide$PLA-LFPdata_wide$MEM
    
    diff_data <- LFPdata_wide
  
    diff_data$Diag <- as.factor(diff_data$Diag)
  
 
   fig <- ggplot(diff_data, aes(GABA, DrugDiff, col = Diag)) + 
    
    geom_point(size = 4) + 
    
    scale_color_manual(values = c('#E69F00', "#00BA38")) +
    
    geom_smooth(method = lm, size = 2, color = "black") + 
    
    
    
    theme_classic() + 
    
    theme(axis.text.x = element_text(color = "black", size = 17, face = "bold"), axis.text.y = element_text(color = "black", size = 17, face = "bold"), axis.title.x = element_text(color = "black", size = 19, face = "bold"), axis.title.y = element_text(color = "black", size = 19, face = "bold")) +
    
    xlab("GABA rIFG (cor.)") +
    
    labs(y = expression(bold(paste(Delta, ' MMN (PLA-MEM)')))) + 
    
    #Turn caption off
    
    theme(legend.position = "none")
  
  
  fig
  
  figfname<-paste(FigOutDir, "LFP_DrugGABAInt_", ROInames[i], "rep6_scat.tiff", sep="")
  
  ggsave(figfname, width = w, height = h)
  
  
  }
  
  
  
  i <- 1
  
  datfname <- paste("adv_ssst_newmaxf_fixICA_wfids_nodipcor_nop10p19_FIX_LFPs_MRSIFGgabacorassoc_meanMMN1_DrugDiff_Pats_fullsample_stats_LMMtableforR_", ROInames[i], ".txt", sep="")   
  
  
  #Not sure why of the multiple steps but oh well..
  
  
  LFPdata_wide <- read.table(datfname, header = TRUE)
  
  LFPdata_wide$DrugDiff <- LFPdata_wide$PLA-LFPdata_wide$MEM
  
  diff_data <- LFPdata_wide
  
  diff_data$Diag <- as.factor(diff_data$Diag)
  
  
  fig <- ggplot(diff_data, aes(GABA, DrugDiff, col = Diag)) + 
    
    geom_point(size = 4) + 
    
    scale_color_manual(values = c('blue', "red")) +
    
    geom_smooth(method = lm, size = 2, color = "black") + 
    
    
    labs(y = expression(bold(paste(Delta, ' MMN (PLA-MEM)')))) + 
    
    theme_classic() + 
    
    theme(axis.text.x = element_text(color = "black", size = 17, face = "bold"), axis.text.y = element_text(color = "black", size = 17, face = "bold"), axis.title.x = element_text(color = "black", size = 19, face = "bold"), axis.title.y = element_text(color = "black", size = 19, face = "bold")) +
    
    xlab("GABA rIFG (cor.)") 
    
    
    #Turn caption off
    
    #theme(legend.position = "none")
  
  
  fig
  
  figfname<-paste(FigOutDir, "LFP_DrugGABAInt_", ROInames[i], "rep1_scat_fkforcols.tiff", sep="")
  
  ggsave(figfname, width = w, height = h)