
# Setup -------------------------------------------------------------------

# Setup working directory and paths

setwd("/Users/alistairperry/Documents/Cambridge/Project/Analysis/LFPs/newmaxfilter_icafixes/LFPs_COH_wrad_allMEGch_wfids/C_Plots/DrugAnalysis")

# The following packages are required

library("plyr")
library("lattice")
library("ggplot2")
library("dplyr")
library("readr")
library("rmarkdown")
library("Rmisc")
library("devtools")
library("gghalves")
library("readxl")

# width and height variables for saved plots - leave it as default for now

w = 4
h = 3

# Set figure output dir

FigOutDir <- "/Users/alistairperry/Documents/Cambridge/Project/Analysis/LFPs/newmaxfilter_icafixes/LFPs_COH_wrad_allMEGch_wfids/C_Plots/DrugAnalysis/MeanPlots/"

dir.create(FigOutDir)


# Text and ROI names for reading, writing files

basefname = "commconpatplac_remo_wbothsess_noc22_LFP_"

ROInames = c("RIFG", "RSTG", "RAUD")

Labnames = c("R IFG", "R STG", "R AUD")



# Drug Group Interaction Setup --------------------------------------------------------

# Loop through each region in for loop

# Pickup text file


for (i in seq(1,3)) {
  
datfname <- paste(basefname, ROInames[i], "_MMNmean_rep3_ConsandPats_DrugInt_fullsample_forJASP.txt", sep="")


# Read file

MMNtab <- read.delim(datfname)


# Split into cons and pats

CN <- subset(MMNtab, Group == 1)

Pat <- subset(MMNtab, Group == 2)


# And also drug session

CN1dat <- CN$meanMMNcol_PLA
CN2dat <- CN$meanMMNcol_MEM
Pat1dat <- Pat$meanMMNcol_PLA
Pat2dat <- Pat$meanMMNcol_MEM


# Setup metadata for plotting

ncnt <- length(CN1dat)
npat <- length(Pat1dat)

ncntscans <- ncnt*2
npatscans <- npat*2
totalscans <- ncntscans+npatscans

x1 <- as.integer(totalscans)
x1[1:ncntscans] <- rep(c(1,2), each=ncnt)
x1[(ncntscans+1):totalscans] <- rep(c(3,4), each=npat)

cntids <- rep(seq(1,ncnt),2)
patids <- rep(seq(ncnt+1,ncnt+npat),2)


# Now compile into appropriate dataframe

d <- data.frame(y = c(CN1dat, CN2dat, Pat1dat, Pat2dat), x = x1, ID_Order = as.factor(c(cntids, patids)))

y_lim_min <- min(d$y)
y_lim_max <- max(d$y)


# Descriptive stats for calculating and connecting group means

# First do the descriptives

group <- c(1,2,1,2)

score_mean_1 <- mean(CN1dat)
score_mean_2 <- mean(CN2dat)
score_mean_3 <- mean(Pat1dat)
score_mean_4 <- mean(Pat2dat)

score_mean <- c(score_mean_1, score_mean_2, score_mean_3, score_mean_4)


#### get confidence intervals
score_ci_1 <- CI(CN1dat, ci = 0.95) #i.e CN1
score_ci_2 <- CI(CN2dat, ci = 0.95) #CN2
score_ci_3 <- CI(Pat1dat, ci = 0.95) #PAT1
score_ci_4 <- CI(Pat2dat, ci = 0.95) #PAT2


# CI's around mean
ci <- c((score_ci_1[1] - score_ci_1[3]), (score_ci_2[1] - score_ci_2[3]), (score_ci_3[1] - score_ci_3[3]), (score_ci_4[1] - score_ci_4[3]))

summary_df <- data.frame(group, score_mean, ci)


# First we must again define the x-coordinates of the means.
x_tick_means_x <- c(1.15, 1.85) #same as above


# Jitter to space datapoints in plot

set.seed(321)

x1fk <- numeric(ncnt+npat)

x1fk[1:ncntscans] <- rep(x_tick_means_x, each=ncnt)
x1fk[(ncntscans+1):totalscans] <- rep(x_tick_means_x, each=npat)

xj <- jitter(x1fk, amount = .09)

d$xj <- xj

d$x_tick_means_x <- as.factor(x_tick_means_x)


# ggplot: Drug Group Interaction ------------------------------------------


f1 <- ggplot(data = d, aes(y = y)) +
  
  #Add geom_() objects
  
  #All data points
  
  geom_point(data = d %>% filter(x =="3"), aes(x = xj), color = 'red', size = 1.5,alpha = .6) +
  geom_point(data = d %>% filter(x =="4"), aes(x = xj), color = 'red', size = 1.5,alpha = .6) +
  
  geom_point(data = d %>% filter(x =="1"), aes(x = xj), color = 'blue', size = 1.5, alpha = .6) +
  geom_point(data = d %>% filter(x =="2"), aes(x = xj), color = 'blue', size = 1.5, alpha = .6) +
  
  
  #Lines connecting each individual across drug condition
  
  geom_line(data = d %>% filter(x =="3" | x =="4"), aes(x = xj, group = ID_Order), color = 'red', alpha = 0.2) +
  
  geom_line(data = d %>% filter(x =="1" | x =="2"), aes(x = xj, group = ID_Order), color = 'blue', alpha = 0.2) +
  
 
  #Axis labels
  
  scale_x_continuous(breaks=x_tick_means_x, labels=c("PLA", "MEM"))+
  
  theme_classic()+
  
  ggtitle(Labnames[i]) + 
  
  xlab("Drug Session") +
  
  ylab("Mean MMN (125-175ms)") +
  
  
  
  #Font and size of axis labels
  
  theme(plot.title = element_text(face="bold", size = 14, hjust = 0.5, color = "black"), axis.title.y = element_text(face="bold", size = 11), axis.title.x = element_text(face="bold", size = 12, color = "black")) + 
  
  theme(axis.text.x = element_text(face="bold", size = 12, color = "black"), axis.text.y = element_text(face="bold", size = 11, color = "black")) +
  
  
  #Turn caption off
  
  theme(legend.position = "right") +
  
  coord_cartesian(ylim=c(y_lim_min, y_lim_max)) +
  
  
  #Add thick lines connecting the means of groups as function of drug condition 
  
  geom_line(data = summary_df[3:4,], aes(x = x_tick_means_x, y = score_mean),
            color = 'red', size = 1.5, alpha = 1) + 
  
  geom_line(data = summary_df[1:2,], aes(x = x_tick_means_x, y = score_mean),
            color = 'blue', size = 1.5, alpha = 1) +
  
  
  #And error bars (again both groups as function of drug condition)
  geom_errorbar(data = d %>% filter(x=="1"), aes(x = 1, y = score_mean[1], ymin = score_mean[1]-ci[1], ymax = score_mean[1]+ci[1]),
                position = position_nudge(0.15),
                color = "blue", width = 0.10, size = 0.8, alpha = 0.8) +
  
  geom_point(data = d %>% filter(x=="1"), aes(x = 1, y = score_mean[1]),
             position = position_nudge(x = 0.15), shape = 21, color = "black", fill = "blue",  size = 3, alpha = 0.8) +
  
  
  geom_errorbar(data = d %>% filter(x=="3"), aes(x = 1, y = score_mean[3], ymin = score_mean[3]-ci[3], ymax = score_mean[3]+ci[3]),
                position = position_nudge(0.15),
                color = "red", width = 0.10, size = 0.8, alpha = 0.8) +
  
  geom_point(data = d %>% filter(x=="3"), aes(x = 1, y = score_mean[3]),
             position = position_nudge(x = 0.15), shape = 21, color = "black", fill = "red",  size = 3, alpha = 0.8) +
  
  
  geom_errorbar(data = d %>% filter(x=="2"), aes(x = 2, y = score_mean[2], ymin = score_mean[2]-ci[2], ymax = score_mean[2]+ci[2]),
                position = position_nudge(-0.15),
                color = "blue", width = 0.10, size = 0.8, alpha = 0.8) +
  
  geom_point(data = d %>% filter(x=="2"), aes(x = 2, y = score_mean[2]),
             position = position_nudge(x = -0.15), shape = 21, color = "black", fill = "blue",  size = 3, alpha = 0.8) +
  
  
  geom_errorbar(data = d %>% filter(x=="4"), aes(x = 2, y = score_mean[4], ymin = score_mean[4]-ci[4], ymax = score_mean[4]+ci[4]),
                position = position_nudge(-0.15),
                color = "red", width = 0.10, size = 0.8, alpha = 0.8) +
  
  geom_point(data = d %>% filter(x=="4"), aes(x = 2, y = score_mean[4]),
             position = position_nudge(x = -0.15), shape = 21, color = "black", fill = "red",  size = 3, alpha = 0.8)


# Write final plot out

f1


# And save

figfname<-paste(FigOutDir, '/LFP_MMN3_ConPatDrugInt_wfig_', ROInames[i], '.tiff', sep="")
ggsave(figfname, width = w, height = h)


} # Loop ends


# 3 Groups: Drug Pat Subgroup Interaction Setup ---------------------------------

#rep3

for (i in seq(1,3)) {
  
  datfname <- paste("adv_ssst_newmaxf_fixICA_wfids_nodipcor_nop10p19_LFP_", ROInames[i], "_MMNmean_rep3_ConsandPats_DrugInt_fullsample_forJASP.txt", sep="")
  
  MMNtab <- read.delim(datfname)
  
  
  MMNtab <- read.delim(datfname)
  
  CN <- subset(MMNtab, PatSubGrp == 3)
  
  BV <- subset(MMNtab, PatSubGrp == 1)
  
  PSP <- subset(MMNtab, PatSubGrp == 2)
  
  
  CN1dat <- CN$meanMMNcol_PLA
  CN2dat <- CN$meanMMNcol_MEM
  
  BV1dat <- BV$meanMMNcol_PLA
  BV2dat <- BV$meanMMNcol_MEM
  
  PSP1dat <- PSP$meanMMNcol_PLA
  PSP2dat <- PSP$meanMMNcol_MEM
  
  #Re-create diag ids
  
  ncnt <- length(CN$meanMMNcol_PLA)
  nbv <- length(BV$meanMMNcol_PLA)
  npsp <- length(PSP$meanMMNcol_PLA)
  
  totalscans <- ncnt+nbv+npsp
  
  
  cntids <- rep(seq(1,ncnt),2)
  bvids <- rep(seq(ncnt+1,ncnt+nbv),2)
  pspids <- rep(seq((ncnt+nbv+1),totalscans),2)
  
  
  ncntscans <- ncnt*2
  nbvscans <- nbv*2
  npspscans <- npsp*2
  totalscans <- ncntscans+nbvscans+npspscans
  
  
  x1 <- as.integer(totalscans)
  x1[1:ncntscans] <- rep(c(1,2), each=ncnt)
  x1[(ncntscans+1):(ncntscans+nbvscans)] <- rep(c(3,4), each=nbv)
  x1[(ncntscans+nbvscans+1):totalscans] <- rep(c(5,6), each=npsp)
  
  
  set.seed(321)
  
  x1fk <- numeric(totalscans)
  
  x1fk[1:ncntscans] <- rep(c(1,2), each=ncnt)
  x1fk[(ncntscans+1):(ncntscans+nbvscans)] <- rep(c(1,2), each=nbv)
  x1fk[(ncntscans+nbvscans+1):totalscans] <- rep(c(1,2), each=npsp)
  
  xj <- jitter(x1fk, amount = .09)
  
  
  d <- data.frame(y = c(CN1dat, CN2dat, BV1dat, BV2dat, PSP1dat, PSP2dat), x = x1, xj = xj, ID_Order = as.factor(c(cntids, bvids, pspids)))
  
  
  #y_lim_min <- min(d$y)
  #y_lim_max <- max(d$y)
  
  
  #Descriptive stats for calculating and connecting group means
  
  #First do the descriptives
  
  group <- c(1,2,1,2,1,2)
  
  score_mean_1 <- mean(CN1dat)
  score_mean_2 <- mean(CN2dat)
  score_mean_3 <- mean(BV1dat)
  score_mean_4 <- mean(BV2dat)
  score_mean_5 <- mean(PSP1dat)
  score_mean_6 <- mean(PSP2dat)
  
  score_mean <- c(score_mean_1, score_mean_2, score_mean_3, score_mean_4, score_mean_5, score_mean_6)
  
  
  #cis
  
  #### get confidence intervals
  score_ci_1 <- CI(CN1dat, ci = 0.95) #CN1
  score_ci_2 <- CI(CN2dat, ci = 0.95) #CN1
  score_ci_3 <- CI(BV1dat, ci = 0.95) #CN1
  score_ci_4 <- CI(BV2dat, ci = 0.95) #CN1
  score_ci_5 <- CI(PSP1dat, ci = 0.95) #CN1
  score_ci_6 <- CI(PSP2dat, ci = 0.95) #CN1
  
  # CI's around mean
  ci <- c((score_ci_1[1] - score_ci_1[3]), (score_ci_2[1] - score_ci_2[3]), (score_ci_3[1] - score_ci_3[3]), (score_ci_4[1] - score_ci_4[3]), (score_ci_5[1] - score_ci_5[3]), (score_ci_6[1] - score_ci_6[3]))
  
  
  summary_df <- data.frame(group, score_mean, ci)
  
  
  
  f1 <- ggplot(data = d, aes(y = y)) +
    
    
    #Add geom_() objects
    
    #PSP
    
    geom_point(data = d %>% filter(x =="5"), aes(x = xj), color = '#00BA38', size = 1.5,alpha = .6) +
    geom_point(data = d %>% filter(x =="6"), aes(x = xj), color = '#00BA38', size = 1.5,alpha = .6) +
    
    
    #BV
    
    geom_point(data = d %>% filter(x =="3"), aes(x = xj), color = '#E69F00', size = 1.5,alpha = .6) +
    geom_point(data = d %>% filter(x =="4"), aes(x = xj), color = '#E69F00', size = 1.5,alpha = .6) +
    
    
    #CN
    geom_point(data = d %>% filter(x =="1"), aes(x = xj), color = 'blue', size = 1.5, alpha = .6) +
    geom_point(data = d %>% filter(x =="2"), aes(x = xj), color = 'blue', size = 1.5, alpha = .6) +
    
    
    #Lines connecting each individual
    
    geom_line(data = d %>% filter(x =="5" | x =="6"), aes(x = xj, group = ID_Order), color = '#00BA38', alpha = 0.2) +
    
    geom_line(data = d %>% filter(x =="3" | x =="4"), aes(x = xj, group = ID_Order), color = '#E69F00', alpha = 0.2) +
    
    geom_line(data = d %>% filter(x =="1" | x =="2"), aes(x = xj, group = ID_Order), color = 'blue', alpha = 0.2) +
    
    
    #Aesthetic stuff
    
    scale_x_continuous(breaks=c(1,2), labels=c("PLA", "MEM"))+
    
    theme_classic()+
    
    ggtitle(Labnames[i]) + 
    
    xlab("Drug Session") +
    
    ylab("Mean MMN (125-175ms)") +
    
    
    
    theme(plot.title = element_text(face="bold", size = 14, hjust = 0.5, color = "black"), axis.title.y = element_text(face="bold", size = 11), axis.title.x = element_text(face="bold", size = 12, color = "black")) + 
    
    theme(axis.text.x = element_text(face="bold", size = 12, color = "black"), axis.text.y = element_text(face="bold", size = 11, color = "black")) +
    
    
    #Turn caption off
    
    
    #Add lines connecting the two means
    
    #PSP
    
    geom_line(data = summary_df[5:6,], aes(x = c(1,2), y = score_mean),
              color = '#00BA38', size = 1.5, alpha = 1) +
    
    #BV
    
    geom_line(data = summary_df[3:4,], aes(x = c(1,2), y = score_mean),
              color = '#E69F00', size = 1.5, alpha = 1) + 
    
    #CN
    
    geom_line(data = summary_df[1:2,], aes(x = c(1,2), y = score_mean),
              color = 'blue', size = 1.5, alpha = 1) +
    
    
    # means with shape = 16, filled circle
    
    #Placebo
    geom_errorbar(data = d %>% filter(x=="1"), aes(x = 1, y = score_mean[1], ymin = score_mean[1]-ci[1], ymax = score_mean[1]+ci[1]),
                  position = position_nudge(0),
                  color = "blue", width = 0.10, size = 0.8, alpha = 0.8) +
    
    geom_point(data = d %>% filter(x=="1"), aes(x = 1, y = score_mean[1]),
               position = position_nudge(0), shape = 21, color = "black", fill = "blue",  size = 3, alpha = 0.8) +
    
    
    geom_errorbar(data = d %>% filter(x=="3"), aes(x = 1, y = score_mean[3], ymin = score_mean[3]-ci[3], ymax = score_mean[3]+ci[3]),
                  position = position_nudge(0),
                  color = "#E69F00", width = 0.10, size = 0.8, alpha = 0.8) +
    
    geom_point(data = d %>% filter(x=="3"), aes(x = 1, y = score_mean[3]),
               position = position_nudge(0), shape = 21, color = "black", fill = "#E69F00",  size = 3, alpha = 0.8) +
    
    
    geom_errorbar(data = d %>% filter(x=="5"), aes(x = 1, y = score_mean[5], ymin = score_mean[5]-ci[5], ymax = score_mean[5]+ci[5]),
                  position = position_nudge(0),
                  color = "#00BA38", width = 0.10, size = 0.8, alpha = 0.8) +
    
    geom_point(data = d %>% filter(x=="5"), aes(x = 1, y = score_mean[5]),
               position = position_nudge(0), shape = 21, color = "black", fill = "#00BA38",  size = 3, alpha = 0.8) +
    
    
    #Drug
    geom_errorbar(data = d %>% filter(x=="2"), aes(x = 2, y = score_mean[2], ymin = score_mean[2]-ci[2], ymax = score_mean[2]+ci[2]),
                  position = position_nudge(0),
                  color = "blue", width = 0.10, size = 0.8, alpha = 0.8) +
    
    geom_point(data = d %>% filter(x=="2"), aes(x = 2, y = score_mean[2]),
               position = position_nudge(0), shape = 21, color = "black", fill = "blue",  size = 3, alpha = 0.8) +
    
    
    geom_errorbar(data = d %>% filter(x=="4"), aes(x = 2, y = score_mean[4], ymin = score_mean[4]-ci[4], ymax = score_mean[4]+ci[4]),
                  position = position_nudge(0),
                  color = "#E69F00", width = 0.10, size = 0.8, alpha = 0.8) +
    
    geom_point(data = d %>% filter(x=="4"), aes(x = 2, y = score_mean[4]),
               position = position_nudge(0), shape = 21, color = "black", fill = "#E69F00",  size = 3, alpha = 0.8) +
    
    
    geom_errorbar(data = d %>% filter(x=="6"), aes(x = 2, y = score_mean[6], ymin = score_mean[6]-ci[6], ymax = score_mean[6]+ci[6]),
                  position = position_nudge(0),
                  color = "#00BA38", width = 0.10, size = 0.8, alpha = 0.8) +
    
    geom_point(data = d %>% filter(x=="6"), aes(x = 2, y = score_mean[6]),
               position = position_nudge(0), shape = 21, color = "black", fill = "#00BA38",  size = 3, alpha = 0.8)
  
  f1
  
  figfname<-paste(FigOutDir, '/LFP_MMN3_ConPatSubGrpDrugInt_', ROInames[i], '.tiff', sep="")
  ggsave(figfname, width = w, height = h)
  
  
}