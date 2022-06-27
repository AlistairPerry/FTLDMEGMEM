setwd("/Users/alistairperry/Documents/Cambridge/Project/Analysis/LFPs/newmaxfilter_icafixes/LFPs_COH_wrad_allMEGch_wfids/C_Plots/DrugAnalysis/MRS")

library("tidyr")

library("lmerTest")

library("ggplot2")

library("dplyr")

library("afex")

library("BayesFactor")

library("bayestestR")


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

  
datfname <- paste("adv_ssst_newmaxf_fixICA_wfids_nodipcor_nop10p19_FIX_LFPs_MRSIFGgluassoc_meanMMN3_DrugDiff_Pats_fullsample_stats_LMMtableforR_", ROInames[i], ".txt", sep="")   

LFPdata_wide <- read.table(datfname, header = TRUE)

LFPdata_long <- gather(LFPdata_wide, drug, MMN, PLA:MEM, factor_key = TRUE)


LFPdata_long$Subject <- factor(LFPdata_long$Subject)


LFPdata_long$MMN_z <- scale(LFPdata_long$MMN)

LFPdata_long$Glu_z <- scale(LFPdata_long$Glu)


contrasts(LFPdata_long$drug) <- contr.sum

mixed.lmer <- lmer(MMN_z ~ Glu_z*drug + (1|Subject), data = LFPdata_long)

lme_res <- anova(mixed.lmer)


#Write table - should write to current dir

outfname <- paste("LMEres_druggluint_", ROInames[i],".txt", sep="")
write.table(lme_res, file = outfname)


#Change Plots

LFPdata_wide$DrugDiff <- LFPdata_wide$PLA-LFPdata_wide$MEM

data <- LFPdata_long
diff_data <- LFPdata_wide

# panel A
MMN_drug_lc_plot <- ggplot(data = data, aes(x = Glu, y = MMN)) +
  
  # gray vertical lines to connect observations from same subject
  geom_segment(data = diff_data, aes(x = Glu, y = PLA,
                                     xend = Glu, yend = MEM),
               size = 0.75, colour = "grey50") +
  # dark gray arrows pointing from placebo to atomoxetine
  geom_segment(data = diff_data %>%
                 # only show arrows for participants with meaningful drug effect
                 dplyr::filter(abs(DrugDiff) > .005),
               aes(x = Glu, y = PLA,
                   xend = Glu, yend = MEM),
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
  coord_cartesian(xlim = c(1.25, 2.75), ylim = c(-0.4, 0.20), expand = FALSE) +
  labs(x = "Glu rIFG (cor.)",
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

figfname<-paste(FigOutDir, "LFP_DrugGluInt_", ROInames[i], ".tiff", sep="")

ggsave(figfname, width = w, height = h)


#Bayes Analysis

# `generalTestBF` computes the Bayes Factor for all possible restrictions of the below model formula.
# The random intercept by subject is always included in all model variants.
# Note the change in formula syntax for random effects.
mixed.lmer_BF_alt <- BayesFactor::generalTestBF(
  MMN_z ~ Glu_z*drug + Subject, data = LFPdata_long,
  whichRandom = "Subject", neverExclude = "Subject")
# For comparison, compute the Bayes Factor of the null model, with only a random intercept by subject
mixed.lmer_BF_null <- BayesFactor::lmBF(
  MMN_z ~ Subject, data = LFPdata_long, whichRandom = "Subject")
# The final Bayes Factors for the alternative models are relative to the null model
mixed.lmer_BF_final <- mixed.lmer_BF_alt / mixed.lmer_BF_null
# Compute Bayes Factor for each *effect*, i.e. the Inclusion Bayes Factors

bayestestR::bayesfactor_inclusion(mixed.lmer_BF_final, match_models = TRUE)


#Save output to file

outfname <- paste("LMEres_druggluint_", ROInames[i], "_BF.txt")
write.table(bayestestR::bayesfactor_inclusion(mixed.lmer_BF_final, match_models = TRUE), file = outfname) 


#Correlations
#Make label sizes bigger for conference

diff_data$Diag <- as.factor(diff_data$Diag)

fig <- ggplot(diff_data, aes(Glu, DrugDiff, col = Diag)) + 
  
  geom_point(size = 4, color = "red", fill = "red", aes(shape=Diag)) +
  
  scale_shape_manual(values=c(22, 24)) +
  
  #scale_color_manual(values = c('#E69F00', "#00BA38")) +
  
  geom_smooth(method = lm, size = 2, color = "black") + 
  
  
  labs(y = expression(bold(paste(Delta, ' MMN (PLA-MEM)')))) + 
  
  theme_classic() + 
  
  theme(axis.text.x = element_text(color = "black", size = 17, face = "bold"), axis.text.y = element_text(color = "black", size = 17, face = "bold"), axis.title.x = element_text(color = "black", size = 19, face = "bold"), axis.title.y = element_text(color = "black", size = 19, face = "bold")) +
  
  xlab("Glu rIFG (cor.)") +
  
  
  #Turn caption off
  
  theme(legend.position = "none")


fig

figfname<-paste(FigOutDir, "LFP_DrugGluInt_", ROInames[i], "_scat.tiff", sep="")

ggsave(figfname, width = w, height = h)


}
  
  
  
  #Rep6
  
  #Correlations

  
  for (i in seq(1,3)) {
    
    
    datfname <- paste("adv_ssst_newmaxf_fixICA_wfids_nodipcor_nop10p19_FIX_LFPs_MRSIFGglucorassoc_meanMMN_DrugDiff_Pats_fullsample_stats_LMMtableforR_", ROInames[i], ".txt", sep="")   
    
    
    #Not sure why of the multiple steps but oh well..
    
    
    LFPdata_wide <- read.table(datfname, header = TRUE)
    
    LFPdata_wide$DrugDiff <- LFPdata_wide$PLA-LFPdata_wide$MEM
    
    diff_data <- LFPdata_wide
  
    diff_data$Diag <- as.factor(diff_data$Diag)
  
 
   fig <- ggplot(diff_data, aes(Glu, DrugDiff, col = Diag)) + 
    
    geom_point(size = 4) + 
    
    scale_color_manual(values = c('#E69F00', "#00BA38")) +
    
    geom_smooth(method = lm, size = 2, color = "black") + 
    
    
    
    theme_classic() + 
    
    theme(axis.text.x = element_text(color = "black", size = 17, face = "bold"), axis.text.y = element_text(color = "black", size = 17, face = "bold"), axis.title.x = element_text(color = "black", size = 19, face = "bold"), axis.title.y = element_text(color = "black", size = 19, face = "bold")) +
    
    xlab("Glu rIFG (cor.)") +
    
    labs(y = expression(bold(paste(Delta, ' MMN (PLA-MEM)')))) + 
    
    #Turn caption off
    
    theme(legend.position = "none")
  
  
  fig
  
  figfname<-paste(FigOutDir, "LFP_DrugGluInt_", ROInames[i], "rep6_scat.tiff", sep="")
  
  ggsave(figfname, width = w, height = h)
  
  
  }
  
  