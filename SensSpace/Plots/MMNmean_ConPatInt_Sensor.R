setwd("/Users/alistairperry/Documents/Cambridge/Project/Analysis/Sensor/newmaxfilter_icafixes/ERF/C_Plots/Plac")


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

library("readxl")


# width and height variables for saved plots - leave it as default for now

w = 6
h = 4


# Set fig output dir

FigOutDir <- "/Users/alistairperry/Documents/Cambridge/Project/Analysis/Sensor/newmaxfilter_icafixes/ERF/C_Plots/Plac/MeanPlots/"

dir.create(FigOutDir)



# MMN

Locfname <- c("GRAD","Frontal GRAD","Temporal GRAD")

LabNames <- c("All Gradiometers","Frontal","Temporal")


for (i in seq(1,3)) {
  
  
datfname <- paste("adv_ssst_newmaxf_fixICA_normALLtrls_automeg_noc22_nonorm_remo_bothsess_newsensors_ERF_MMNdiffmean_rep3_ConsandPats_", Locfname[i], "_fullsample.txt", sep="")    


MMNtab <- read.delim(datfname, header = TRUE)


CN <- subset(MMNtab, Group == 1)

Pat <- subset(MMNtab, Group == 2)


CNdat <- CN$meanMMNcol
Patdat <-Pat$meanMMNcol

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
  
  ggtitle(LabNames[i]) + 
  
  theme(plot.title = element_text(face="bold", size = 18, hjust = 0.5, color = "black"), axis.title.y = element_text(face="bold", size = 19)) + 
  
  theme(axis.text.x = element_text(face="bold", size = 19, color = "black"), axis.text.y = element_text(face="bold", size = 18, color = "black")) +
  
  
  #Turn caption off
  
  theme(legend.position = "none")

#coord_cartesian(ylim=c(y_lim_min, y_lim_max)) 


f1


figfname <- paste(FigOutDir, 'ConsandPats_meanMMN_rep3_',Locfname[i],'.tiff', sep="")

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
  
  ggtitle(LabNames[i]) + 
  
  theme(plot.title = element_text(face="bold", size = 20, hjust = 0.5, color = "black"), axis.title.y = element_text(face="bold", size = 21), axis.title.x = element_text(face="bold", size = 21)) + 
  
  theme(axis.text.x = element_text(face="bold", size = 21, color = "black"), axis.text.y = element_text(face="bold", size = 20, color = "black")) +
  
  
  #Turn caption off
  
  theme(legend.position = "none")


f1_box

figfname <- paste(FigOutDir, 'ConsandPats_meanMMN_rep3_',Locfname[i],'_box.tiff', sep="")

ggsave(figfname, width = 6, height = 4)


}

