
# Setup -------------------------------------------------------------------

# Setup working dir and paths

setwd("/Users/alistairperry/Documents/Cambridge/Project/Analysis/LFPs/newmaxfilter_icafixes/LFPs_COH_wrad_allMEGch_wfids/C_Plots/DrugAnalysis/MRS")

library("tidyr")

library("lmerTest")

library("ggplot2")

library("dplyr")

library("afex")

library("BayesFactor")

library("bayestestR")

library("interactions")


# Output directory and prefix

FigOutDir  <- "/Users/alistairperry/Documents/Cambridge/Project/Analysis/LFPs/newmaxfilter_icafixes/LFPs_COH_wrad_allMEGch_wfids/C_Plots/DrugAnalysis/MRS/Scatt/"
dir.create(FigOutDir)


FigOutPrefix <- "allsubjs"


# Width and height variables for saved plots 

w = 6
h = 4


# Regions to load and analyse

ROInames <- c("RIFG", "RSTG", "RAUD")

Labnames <- c("R IFG", "R STG", "R AUD")



# Correlation Plots - Drug Difference by GABA -----------------------------

# Figures 4A

for (i in seq(1,length(ROInames))) {

  
# Load data
    
datfname <- paste("adv_ssst_newmaxf_fixICA_wfids_nodipcor_nop10p19_FIX_LFPs_MRSIFGgabacorassoc_meanMMN3_DrugDiff_Pats_fullsample_stats_LMMtableforR_", ROInames[i], ".txt", sep="")   

LFPdata_wide <- read.table(datfname, header = TRUE)


# Calculate drug difference in MMN

LFPdata_wide$DrugDiff <- LFPdata_wide$PLA-LFPdata_wide$MEM


# Correlations

LFPdata_wide$Diag <- as.factor(LFPdata_wide$Diag)


# Start scatter

fig <- ggplot(LFPdata_wide, aes(GABA, DrugDiff, col = Diag)) + 
  
  
  # Individual data points and aes
  
  geom_point(size = 4, color = "red", fill = "red", aes(shape=Diag)) +
  
  # Scale shape by disease group (Diag)
  
  scale_shape_manual(values=c(22, 24)) +
  
  # Insert linear fit 
  
  geom_smooth(method = lm, size = 2, color = "black") + 
  
  # Axis labels
  
  labs(y = expression(bold(paste(Delta, ' MMN (PLA-MEM)')))) + 
  
  theme_classic() + 
  
  theme(axis.text.x = element_text(color = "black", size = 17, face = "bold"), axis.text.y = element_text(color = "black", size = 17, face = "bold"), axis.title.x = element_text(color = "black", size = 19, face = "bold"), axis.title.y = element_text(color = "black", size = 19, face = "bold")) +
  
  xlab("GABA rIFG (cor.)") +

  #Lastly, turn caption off
  
  theme(legend.position = "none")


  #Produce figure to window

  fig


  #And save to disk

  figfname<-paste(FigOutDir, "LFP_DrugGABAInt_", ROInames[i], "_scat.tiff", sep="")

  ggsave(figfname, width = w, height = h)

}

# Loop ended


# Drug Session Change Plots -----------------------------------------------

# Figure 4B

# Start loop again


for (i in seq(1,length(ROInames))) {
  
  
  # Load data, and setup long version
  
  datfname <- paste("adv_ssst_newmaxf_fixICA_wfids_nodipcor_nop10p19_FIX_LFPs_MRSIFGgabacorassoc_meanMMN3_DrugDiff_Pats_fullsample_stats_LMMtableforR_", ROInames[i], ".txt", sep="")   
  
  LFPdata_wide <- read.table(datfname, header = TRUE)
  
  LFPdata_long <- gather(LFPdata_wide, drug, MMN, PLA:MEM, factor_key = TRUE)
  
  LFPdata_long$Subject <- factor(LFPdata_long$Subject)
  
  
  # Calculate drug difference in MMN
  
  LFPdata_wide$DrugDiff <- LFPdata_wide$PLA-LFPdata_wide$MEM
  
  data <- LFPdata_long
  diff_data <- LFPdata_wide
  
  # Start panel
  MMN_drug_lc_plot <- ggplot(data = data, aes(x = GABA, y = MMN)) +
    
  # Add all panel features  
    
    # gray vertical lines to connect observations from same subject
    geom_segment(data = diff_data, aes(x = GABA, y = PLA,
                                       xend = GABA, yend = MEM),
                 size = 0.75, colour = "grey50") +
    
    # dark gray arrows pointing from placebo to drug
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
      axis.text = ggplot2::element_text(color = "black", size = 17, face = "bold"),
      axis.title = ggplot2::element_text(color = "black", size = 19, face = "bold"),
    ) 
  
  
  MMN_drug_lc_plot
  
  figfname<-paste(FigOutDir, "LFP_DrugGABAInt_", ROInames[i], ".tiff", sep="")
  
  ggsave(figfname, width = w, height = h)
 
   
}
  
# Loop finished


# Glu/GABA ratio group differences  ------------------------------------

# (Figure 4C)

LFPdata_wide <- read.table("adv_ssst_newmaxf_fixICA_wfids_nodipcor_remo_nop10p19_LFPs_MRSIFGglugabaratassoc_ConsandPats.txt", header = TRUE)


# Split data into groups

CN <- subset(LFPdata_wide, Group == 1)

Pat <- subset(LFPdata_wide, Group == 2)

CNdat <- CN$GluGABArat
Patdat <- Pat$GluGABArat


# Setup for jittering points

ncnt <- length(CNdat)
npat <- length(Patdat)

totalscans <- ncnt+npat

x1 <- as.integer(totalscans)
x1[1:ncnt] <- rep(1, each=ncnt)
x1[(ncnt+1):totalscans] <- rep(2, each=npat)

d <- data.frame(y = c(CNdat, Patdat), x = x1, PatSubgrp = LFPdata_wide$PatSubgrp)


# Now jitter

set.seed(321)

xj <- jitter(x1)


# Calc means for box plots

d$xj <- xj

d$x <- as.factor(d$x)

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


# Start box plot

d$PatSubgrp <- as.factor(d$PatSubgrp)

f1_box <- ggplot(my_sum_box, aes(fill=x, x=x, y = y50)) +
  
  
  #Start box
  
  geom_boxplot(aes(ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100),
               position="dodge", stat="identity", alpha = 0.7, width = 0.75, color = "black", size = 0.75) + 
  
  
  # Set box colour
  
  scale_fill_manual(values = c("Blue", "Red")) + 
  
  
  # Individual data point aesthetics
  
  geom_point(data = d %>% filter(x =="1"), aes(x = xj, y = y), color = 'black', size = 3, alpha = .8) +
  
  geom_point(data = d %>% filter(x =="2"), aes(x = xj, y = y, shape = PatSubgrp), color = 'black', fill = 'black', size = 3, alpha = .8) + 
  
  scale_shape_manual(values=c(22, 24)) +
  
  
  # Set labels and font style
  
  scale_x_discrete(breaks=c(1,2), labels=c("CON", "bvFTD/PSP")) +
  
  ylab("Glu/GABA rIFG") +
  
  xlab("") +
  
  theme_classic() +
  
  theme(plot.title = element_text(face="bold", size = 18, hjust = 0.5, color = "black"), axis.title.y = element_text(face="bold", size = 19)) + 
  
  theme(axis.text.x = element_text(face="bold", size = 19, color = "black"), axis.text.y = element_text(face="bold", size = 18, color = "black")) +
  
  
  #Turn caption off
  
  theme(legend.position = "none")


  #Produce plot in window

  f1_box

  #Now save
  
  figfname<-paste(FigOutDir, "GroupGluGABAratmean_box.tiff", sep="")

  ggsave(figfname, width = w, height = h)


  

# Glu/GABA and Drug Diff Interaction --------------------------------------

# Figure 4D
  
# Only for R AUD    
  

# Read data  
    
LFPdata_wide <- read.table("adv_ssst_newmaxf_fixICA_wfids_nodipcor_remo_nop10p19_LFPs_MRSIFGglugabaratassoc_meanMMN3_DrugDiff_Pats_fullsample_stats_LMMtableforR_RAud.txt", header = TRUE) 

LFPdata_wide$Diag <- as.factor(LFPdata_wide$Diag)


# Start scatter plot

fig <- ggplot(LFPdata_wide, aes(GluGABArat, Diff, color = GABA)) + 
  
  
  #Regression line
  
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
  
  

  
  
  #Supplementary: Moderating effect of Age ---------------------------------
  
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
  
#Finished