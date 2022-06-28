setwd("/Users/alistairperry/Documents/Cambridge/Project/Analysis/LFPs/newmaxfilter_icafixes/LFPs_COH_wrad_allMEGch_wfids/C_Plots/Plac")


#Load packages
library("plyr")
library("lattice")
library("ggplot2")
library("dplyr")
library("readr")
library("rmarkdown")
library("Rmisc")
library("devtools")
library("gghalves")
library("ggsignif")
library("readxl")


# width and height variables for saved plots - leave it as default for now

w = 6
h = 4


# Set fig output dir

FigOutDir <- "/Users/alistairperry/Documents/Cambridge/Project/Analysis/LFPs/newmaxfilter_icafixes/LFPs_COH_wrad_allMEGch_wfids/C_Plots/Plac/MeanPlots/"
dir.create(FigOutDir)



# MMN

Locfname <- c("RIFG","RSTG","RAUD")

LabNames <- c("R IFG","R STG","R AUD")


datfname <- paste("commconpatplac_remo_wbothsess_noc22_LFPs_MMN3diffmean_ConsandPats_forJASP_fullsample.txt")    


MMNtab <- read.delim(datfname)

CN <- subset(MMNtab, Group == 1)

Pat <- subset(MMNtab, Group == 2)


#RIFG

  CNdat <- CN[,4]
  Patdat <- Pat[,4]
  
  ncnt <- length(CNdat)
  npat <- length(Patdat)
  
  totalscans <- ncnt+npat
  
  
  x1 <- as.integer(totalscans)
  x1[1:ncnt] <- rep(1, each=ncnt)
  x1[(ncnt+1):totalscans] <- rep(2, each=npat)
  
  
  d <- data.frame(y = c(CNdat, Patdat), x = x1)
  
  
  #Start plotting
  
  set.seed(321)
  
  #xj <- jitter(x1, amount = .09)
  xj <- jitter(x1)
  
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
  
  
  #Set axis limits
  
  #y_lim_min <- -2
  #y_lim_max <- 0.5
  
  #coord_cartesian(ylim=c(y_lim_min, y_lim_max)) 
  
  
  #Start
  
  
  f1 <- ggplot(my_sum, aes(fill=x, y=mean, x=x)) + 
    
    
    
    geom_bar(position="dodge", stat="identity", alpha = 0.7, width = 0.75, color = "black", size = 0.75) +
    
    
    
    scale_fill_manual(values = c("Blue", "Red")) +
    
      geom_errorbar(aes(x=x, ymin=mean, ymax=mean+se), width=0.4, colour="black", alpha=0.9, size=0.75) +

  
  geom_point(data = d %>% filter(x =="1"), aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .6) +
    
    geom_point(data = d %>% filter(x =="2"), aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .6) + 
    
    scale_x_discrete(breaks=c(1,2), labels=c("CON", "bvFTD/PSP")) +
    
    ylab("Mean MMN (125-175ms)") +
    
    xlab("") +
    
    theme_classic() +
    
    ggtitle(LabNames[1]) + 
    
    theme(plot.title = element_text(face="bold", size = 18, hjust = 0.5, color = "black"), axis.title.y = element_text(face="bold", size = 19)) + 
    
    theme(axis.text.x = element_text(face="bold", size = 19, color = "black"), axis.text.y = element_text(face="bold", size = 18, color = "black")) +
    
    
    #Turn caption off
    
    theme(legend.position = "none")
  
  #coord_cartesian(ylim=c(y_lim_min, y_lim_max)) 
  
  
  f1
  
  
  figfname <- paste(FigOutDir, 'ConsandPats_meanMMN_rep3_',Locfname[1],'.tiff', sep="")
  
  ggsave(figfname, width = 6, height = 4)
  
  
  #Now boxplot
  

  
  f1_box <- ggplot(my_sum_box, aes(fill=x, x=x, y = y50)) +
    
    #stat_boxplot(geom ='errorbar', width = 0.6) +
    
    geom_boxplot(aes(ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100),
                 position="dodge", stat="identity", alpha = 0.7, width = 0.75, color = "black", size = 0.75) + 
    
    scale_fill_manual(values = c("Blue", "Red")) + 
    
    geom_point(data = d %>% filter(x =="1"), aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .8) +
    
    geom_point(data = d %>% filter(x =="2"), aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .8) + 
    
    scale_x_discrete(breaks=c(1,2), labels=c("CON", "bvFTD/PSP")) +
    
    ylab("Mean MMN (125-175ms)") +
    
    xlab("") +
    
    theme_classic() +
    
    #ggtitle(LabNames[1]) + 
    
    theme(plot.title = element_text(face="bold", size = 18, hjust = 0.5, color = "black"), axis.title.y = element_text(face="bold", size = 19)) + 
    
    theme(axis.text.x = element_text(face="bold", size = 19, color = "black"), axis.text.y = element_text(face="bold", size = 18, color = "black")) +
    
    
    #Turn caption off
    
    theme(legend.position = "none")
    
  
  f1_box
  
  figfname <- paste(FigOutDir, 'ConsandPats_meanMMN_rep3_',Locfname[1],'_box.tiff', sep="")
  
  ggsave(figfname, width = 6, height = 4)
  
  
#RSTG
  
  CNdat <- CN[,5]
  Patdat <- Pat[,5]
  
  ncnt <- length(CNdat)
  npat <- length(Patdat)
  
  totalscans <- ncnt+npat
  
  
  x1 <- as.integer(totalscans)
  x1[1:ncnt] <- rep(1, each=ncnt)
  x1[(ncnt+1):totalscans] <- rep(2, each=npat)
  
  
  d <- data.frame(y = c(CNdat, Patdat), x = x1)
  
  
  #Start plotting
  
  set.seed(321)
  
  #xj <- jitter(x1, amount = .09)
  xj <- jitter(x1)
  
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
  
  
  #Set axis limits
  
  #y_lim_min <- -2
  #y_lim_max <- 0.5
  
  #coord_cartesian(ylim=c(y_lim_min, y_lim_max)) 
  
  
  #Start
  
  
  f1 <- ggplot(my_sum, aes(fill=x, y=mean, x=x)) + 
    
    
    
    geom_bar(position="dodge", stat="identity", alpha = 0.7, width = 0.75, color = "black", size = 0.75) +
    
    
    
    scale_fill_manual(values = c("Blue", "Red")) +
    
    geom_errorbar(aes(x=x, ymin=mean-se, ymax=mean), width=0.4, colour="black", alpha=0.9, size=0.75) +
    
    
    geom_point(data = d %>% filter(x =="1"), aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .6) +
    
    geom_point(data = d %>% filter(x =="2"), aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .6) + 
    
    scale_x_discrete(breaks=c(1,2), labels=c("CON", "bvFTD/PSP")) +
    
    ylab("Mean MMN (125-175ms)") +
    
    xlab("") +
    
    theme_classic() +
    
    ggtitle(LabNames[2]) + 
    
    theme(plot.title = element_text(face="bold", size = 18, hjust = 0.5, color = "black"), axis.title.y = element_text(face="bold", size = 19)) + 
    
    theme(axis.text.x = element_text(face="bold", size = 19, color = "black"), axis.text.y = element_text(face="bold", size = 18, color = "black")) +
    
    
    #Turn caption off
    
    theme(legend.position = "none")
  
  #coord_cartesian(ylim=c(y_lim_min, y_lim_max)) 
  
  
  f1
  
  
  figfname <- paste(FigOutDir, 'ConsandPats_meanMMN_rep3_',Locfname[2],'.tiff', sep="")
  
  ggsave(figfname, width = 6, height = 4)
  
  
  #Now boxplot
  
  
  
  f1_box <- ggplot(my_sum_box, aes(fill=x, x=x, y = y50)) +
    
    #stat_boxplot(geom ='errorbar', width = 0.6) +
    
    geom_boxplot(aes(ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100),
                 position="dodge", stat="identity", alpha = 0.7, width = 0.75, color = "black", size = 0.75) + 
    
    scale_fill_manual(values = c("Blue", "Red")) + 
    
    geom_point(data = d %>% filter(x =="1"), aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .8) +
    
    geom_point(data = d %>% filter(x =="2"), aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .8) + 
    
    scale_x_discrete(breaks=c(1,2), labels=c("CON", "bvFTD/PSP")) +
    
    ylab("Mean MMN (125-175ms)") +
    
    xlab("") +
    
    theme_classic() +
    
    #ggtitle(LabNames[1]) + 
    
    theme(plot.title = element_text(face="bold", size = 18, hjust = 0.5, color = "black"), axis.title.y = element_text(face="bold", size = 19)) + 
    
    theme(axis.text.x = element_text(face="bold", size = 19, color = "black"), axis.text.y = element_text(face="bold", size = 18, color = "black")) +
    
    
    #Turn caption off
    
    theme(legend.position = "none")
  
  
  f1_box
  
  figfname <- paste(FigOutDir, 'ConsandPats_meanMMN_rep3_',Locfname[2],'_box.tiff', sep="")
  
  ggsave(figfname, width = 6, height = 4)
  
  
  #RAUD
  
  
  CNdat <- CN[,6]
  Patdat <- Pat[,6]
  
  ncnt <- length(CNdat)
  npat <- length(Patdat)
  
  totalscans <- ncnt+npat
  
  
  x1 <- as.integer(totalscans)
  x1[1:ncnt] <- rep(1, each=ncnt)
  x1[(ncnt+1):totalscans] <- rep(2, each=npat)
  
  
  d <- data.frame(y = c(CNdat, Patdat), x = x1)
  
  
  #Start plotting
  
  set.seed(321)
  
  #xj <- jitter(x1, amount = .09)
  xj <- jitter(x1)
  
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
  
  #Set axis limits
  
  #y_lim_min <- -2
  #y_lim_max <- 0.5
  
  #coord_cartesian(ylim=c(y_lim_min, y_lim_max)) 
  
  
  #Start
  
  
  f1 <- ggplot(my_sum, aes(fill=x, y=mean, x=x)) + 
    
    
    
    geom_bar(position="dodge", stat="identity", alpha = 0.7, width = 0.75, color = "black", size = 0.75) +
    
    
    
    scale_fill_manual(values = c("Blue", "Red")) +
    
    geom_errorbar(aes(x=x, ymin=mean-se, ymax=mean), width=0.4, colour="black", alpha=0.9, size=0.75) +
    
    
    geom_point(data = d %>% filter(x =="1"), aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .6) +
    
    geom_point(data = d %>% filter(x =="2"), aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .6) + 
    
    scale_x_discrete(breaks=c(1,2), labels=c("CON", "bvFTD/PSP")) +
    
    ylab("Mean MMN (125-175ms)") +
    
    xlab("") +
    
    theme_classic() +
    
    ggtitle(LabNames[3]) + 
    
    theme(plot.title = element_text(face="bold", size = 18, hjust = 0.5, color = "black"), axis.title.y = element_text(face="bold", size = 19)) + 
    
    theme(axis.text.x = element_text(face="bold", size = 19, color = "black"), axis.text.y = element_text(face="bold", size = 18, color = "black")) +
    
    
    #Turn caption off
    
    theme(legend.position = "none")
  
  #coord_cartesian(ylim=c(y_lim_min, y_lim_max)) 
  
  
  f1
  
  
  figfname <- paste(FigOutDir, 'ConsandPats_meanMMN_rep3_',Locfname[3],'.tiff', sep="")
  
  ggsave(figfname, width = 6, height = 4)
  
  
  
  #Boxplot
  
  f1_box <- ggplot(my_sum_box, aes(fill=x, x=x, y = y50)) +
    
    #stat_boxplot(geom ='errorbar', width = 0.6) +
    
    geom_boxplot(aes(ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100),
                 position="dodge", stat="identity", alpha = 0.7, width = 0.75, color = "black", size = 0.75) + 
    
    scale_fill_manual(values = c("Blue", "Red")) + 
    
    geom_point(data = d %>% filter(x =="1"), aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .8) +
    
    geom_point(data = d %>% filter(x =="2"), aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .8) + 
    
    scale_x_discrete(breaks=c(1,2), labels=c("CON", "bvFTD/PSP")) +
    
    ylab("Mean MMN (125-175ms)") +
    
    xlab("") +
    
    theme_classic() +
    
    #ggtitle(LabNames[1]) + 
    
    theme(plot.title = element_text(face="bold", size = 18, hjust = 0.5, color = "black"), axis.title.y = element_text(face="bold", size = 19)) + 
    
    theme(axis.text.x = element_text(face="bold", size = 19, color = "black"), axis.text.y = element_text(face="bold", size = 18, color = "black")) +
    
    
    #Turn caption off
    
    theme(legend.position = "none")
  
  
  f1_box
  
  figfname <- paste(FigOutDir, 'ConsandPats_meanMMN_rep3_',Locfname[3],'_box.tiff', sep="")
  
  ggsave(figfname, width = 6, height = 4)
  
  
  
  #Boxplot
  
  f1_box <- ggplot(my_sum_box, aes(fill=x, x=x, y = y50)) +
    
    #stat_boxplot(geom ='errorbar', width = 0.6) +
    
    geom_boxplot(aes(ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100),
                 position="dodge", stat="identity", alpha = 0.7, width = 0.75, color = "black", size = 0.75) + 
    
    scale_fill_manual(values = c("Blue", "Red")) + 
    
    geom_point(data = d %>% filter(x =="1"), aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .8) +
    
    geom_point(data = d %>% filter(x =="2"), aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .8) + 
    
    scale_x_discrete(breaks=c(1,2), labels=c("CON", "bvFTD/PSP")) +
    
    ylab("Mean MMN (125-175ms)") +
    
    xlab("") +
    
    theme_classic() +
    
    #ggtitle(LabNames[1]) + 
    
    theme(plot.title = element_text(face="bold", size = 18, hjust = 0.5, color = "black"), axis.title.y = element_text(face="bold", size = 19)) + 
    
    theme(axis.text.x = element_text(face="bold", size = 19, color = "black"), axis.text.y = element_text(face="bold", size = 18, color = "black")) +
    
    
    #Turn caption off
    
    theme(legend.position = "none")
  
  
  f1_box
  
  figfname <- paste(FigOutDir, 'ConsandPats_meanMMN_',Locfname[3],'_box.tiff', sep="")
  
  ggsave(figfname, width = 6, height = 4)
  
  
  f1_box <- ggplot(my_sum_box, aes(fill=x, x=x, y = y50)) +
    
    #stat_boxplot(geom ='errorbar', width = 0.6) +
    
    geom_boxplot(aes(ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100),
                 position="dodge", stat="identity", alpha = 0.7, width = 0.75, color = "black", size = 0.75) + 
    
    scale_fill_manual(values = c("Blue", "Red")) + 
    
    geom_point(data = d %>% filter(x =="1"), aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .8) +
    
    geom_point(data = d %>% filter(x =="2"), aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .8) + 
    
    scale_x_discrete(breaks=c(1,2), labels=c("CON", "bvFTD/PSP")) +
    
    ylab("Mean MMN (125-175ms)") +
    
    xlab("") +
    
    theme_classic() +
    
    #ggtitle(LabNames[1]) + 
    
    theme(plot.title = element_text(face="bold", size = 18, hjust = 0.5, color = "black"), axis.title.y = element_text(face="bold", size = 19)) + 
    
    theme(axis.text.x = element_text(face="bold", size = 19, color = "black"), axis.text.y = element_text(face="bold", size = 18, color = "black"))
    
    
    #Turn caption off
    
    #theme(legend.position = "none")
  
  
  f1_box
  
  figfname <- paste(FigOutDir, 'ConsandPats_meanMMN_',Locfname[3],'_box_wleg.tiff', sep="")
  
  ggsave(figfname, width = 6, height = 4)
  
  

# 3 Groups ----------------------------------------------------------------

  
  # width and height variables for saved plots - leave it as default for now
  
  w = 7
  h = 4
  
  
  CN <- subset(MMNtab, PatSubgrp == 3)
  
  BV <- subset(MMNtab, PatSubgrp == 1)
  
  PSP <- subset(MMNtab, PatSubgrp == 2)
  
  
  
  #Re-create diag ids
  
  ncnt <- length(CN$meanMMNcol_1)
  nbv <- length(BV$meanMMNcol_1)
  npsp <- length(PSP$meanMMNcol_1)
  
  totalscans <- ncnt+nbv+npsp
  
  
  x1 <- as.integer(totalscans)
  x1[1:ncnt] <- rep(1, each=ncnt)
  x1[(ncnt+1):(ncnt+nbv)] <- rep(2, each=nbv)
  x1[(ncnt+nbv+1):totalscans] <- rep(3, each=npsp)
  
  
  set.seed(321)
  
  #xj <- jitter(x1, amount = .09)
  xj <- jitter(x1)
  
  
  #RAUD
  
  
  for (i in seq(4,6)) {
    
    
    d <- data.frame(y = c(CN[,i], BV[,i], PSP[,i]), x = x1, xj = xj)
    
    
    #Ensure factors
    
    
    d$x <- as.factor(d$x)
    
    
    my_sum_box <- d %>% 
      group_by(x) %>% 
      summarise(
        y0 = quantile(y, 0.05), 
        y25 = quantile(y, 0.25), 
        y50 = mean(y), 
        y75 = quantile(y, 0.75), 
        y100 = quantile(y, 0.95))
    
    
    
    #Boxplot
    
    f1_box <- ggplot(my_sum_box, aes(fill=x, x=x, y = y50)) +
      
      #stat_boxplot(geom ='errorbar', width = 0.6) +
      
      geom_boxplot(aes(ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100),
                   position="dodge", stat="identity", alpha = 0.7, width = 0.75, color = "black", size = 0.75) + 
      
      scale_fill_manual(values = c("Blue", "Red", "Red")) + 
      
      
      geom_point(data = d, aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .8) +
      
      
      scale_x_discrete(breaks=c(1,2,3), labels=c("CON", "bvFTD", "PSP")) +
      
      
      ylab("Mean MMN (125-175ms)") +
      
      xlab("") +
      
      theme_classic() +
      
      ggtitle(LabNames[i-3]) + 
      
      
      theme(plot.title = element_text(face="bold", size = 18, hjust = 0.5, color = "black"), axis.title.y = element_text(face="bold", size = 19)) + 
      
      theme(axis.text.x = element_text(face="bold", size = 19, color = "black"), axis.text.y = element_text(face="bold", size = 18, color = "black")) +
      
      
      #Turn caption off
      
      theme(legend.position = "none")
    
    
    f1_box
    
    figfname <- paste(FigOutDir, 'ConsandPatSubGrp_meanMMN_rep3_',Locfname[i-3],'_box.tiff', sep="")
    
    ggsave(figfname, width = w, height = h)
    
    
  }
  
  
  #For RIFG - add significance 
  
  d <- data.frame(y = c(CN[,4], BV[,4], PSP[,4]), x = x1, xj = xj)
  
  
  #Ensure factors
  
  
  d$x <- as.factor(d$x)
  
  
  my_sum_box <- d %>% 
    group_by(x) %>% 
    summarise(
      y0 = quantile(y, 0.05), 
      y25 = quantile(y, 0.25), 
      y50 = mean(y), 
      y75 = quantile(y, 0.75), 
      y100 = quantile(y, 0.95))
  
  
  f1_box <- ggplot(my_sum_box, aes(fill=x, x=x, y = y50)) +
    
    #stat_boxplot(geom ='errorbar', width = 0.6) +
    
    geom_boxplot(aes(ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100),
                 position="dodge", stat="identity", alpha = 0.7, width = 0.75, color = "black", size = 0.75) + 
    
    scale_fill_manual(values = c("Blue", "Red", "Red")) + 
    
    
    
    
    geom_point(data = d, aes(x = xj, y = y), color = 'black', size = 2.5, alpha = .8) +
    
    
    scale_x_discrete(breaks=c(1,2,3), labels=c("CON", "bvFTD", "PSP")) +
    
    geom_signif(data = d, aes(x = x, y = y), comparisons = list(c(2, 3)), 
                annotation = c("*")) +
    
    ylab("Mean MMN (125-175ms)") +
    
    xlab("") +
    
    theme_classic() +
    
    ggtitle(LabNames[1]) + 
    
    
    theme(plot.title = element_text(face="bold", size = 18, hjust = 0.5, color = "black"), axis.title.y = element_text(face="bold", size = 19)) + 
    
    theme(axis.text.x = element_text(face="bold", size = 19, color = "black"), axis.text.y = element_text(face="bold", size = 18, color = "black")) +
    
    
    #Turn caption off
    
    theme(legend.position = "none")
  
  
  f1_box
  
  figfname <- paste(FigOutDir, 'ConsandPatSubGrp_meanMMN_rep3_',Locfname[1],'_box_wsig.tiff', sep="")
  
  ggsave(figfname, width = w, height = h)
  