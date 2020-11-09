##########################################################
##                                                      ##   
## Author: Sherry Lin                                   ##
## Student ID: u6890762                                 ##
## Project: evolution of bee viruses                    ##
## Purpose: plotting overall viral load and prevalence  ##
## Version: 18                                          ##
## Date: 10/11/20                                       ##
##                                                      ##
##########################################################

library(ggplot2) #ggplot()+geom_boxplot(), geom_bar(), geom_line, geom_point()
library(dplyr) #group_by(): groups a dataframe based on certain fields / summarise()
library(ggpubr) #ggarrange(): make multiple ggplots into one / annotate_figure()
library(reshape2)

rm(list = ls())

# import data
# let 'Geocluster' become orderd factor(from south to north)
data <- read.csv("new_data.csv")

data$Geocluster <- factor(data$Geocluster, level=c("WEL","BOR","TAM","VER"))
data$year.1 <- as.factor(data$year.1)
data$year.2 <- factor(data$year.2, level=c('0','<1',"<2","<3","<5","<7"))
data$year.3 <- factor(data$year.3, level=c('0','<1','<2','>2'))
data$V.arrival <- as.factor(data$V.arrival)

# plotting viral population over time
## by TPM
name_full <- c('Acute bee paralysis virus','Apis mellifera filamentous virus',
               'Apis rhabdovirus 2','Black queen cell virus', 'Chronic bee paralysis virus',
               'Deformed wing virus-A','Deformed wing virus-C','Israeli acute paralysis virus',
               'Kashmir bee virus','Lake Sinai virus','Slow bee paralysis virus','Sacbrood virus',
               'Varroa destructor virus 1','Varroa destructor virus 2','Varroa destructor virus 3',
               'Varroa orthomyxovirus-1','Varroa Tymo-like virus')

# BOXPLOT: y=tpm, x=year #================================================================================
## Time since first Varroa detection: -8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,5,7
for(i in 25:41){
  p1.1<-ggplot(data=data[-which(is.na(data$year.1)),], aes_string(x='year.1',y=names(data)[i],fill='V.arrival'))+
    geom_boxplot(outlier.size = 0.8)+
    scale_y_sqrt()+
    facet_grid(rows=vars(Geocluster))+
    scale_fill_manual(values = c('white','grey'))+
    theme_classic()+
    theme(legend.position="none",plot.margin=unit(c(-1.8,1,-0.5,-0.5),"lines"))+
    labs(title = "",y="",x="")
    scale_x_discrete(labels=c('<8',"<7",'<6','<5','<4','<3','<2','<1','0','<1','<2','<3','<5','<7'))
    
  
  p1.2<-ggplot(data=data[-which(is.na(data$year.1)),], aes_string(x='year.1',y=names(data)[i],fill='V.arrival'))+
    geom_boxplot(outlier.size = 0.8)+
    scale_y_sqrt()+
    scale_fill_manual(values = c('white','grey'))+
    theme_classic()+
    theme(plot.margin=unit(c(1,2.2,0,-0.6),"lines"),legend.position = 'none')+
    labs(y="",x="")
    scale_x_discrete(labels=c('<8',"<7",'<6','<5','<4','<3','<2','<1','0','<1','<2','<3','<5','<7'))
    
  p1 <- ggarrange(p1.2,p1.1,nrow=2,ncol=1,heights=c(0.3, 0.7))
  p1 <- annotate_figure(p1,top=text_grob(paste(name_full[i-24])),
                        left = text_grob("Transcript per million",rot = 90),
                        bottom=text_grob("Time since first Varroa detection (years)"))
  print(p1)
  
}

#-----------------------------------------------------------------------------------------------------
## exposure time groups: 0,<1,<2,<3,<5,<7
for(i in 25:41){
  p1.1 <- ggplot(data, aes_string(x='year.2',y=names(data)[i],fill='V.arrival'))+
    geom_boxplot(outlier.shape = NA)+
    scale_y_sqrt()+
    facet_grid(rows=vars(Geocluster))+
    scale_fill_manual(values = c('white','grey'))+
    theme_classic()+
    theme(legend.position="none",plot.margin=unit(c(-1.8,1,-0.5,-0.5),"lines"))+
    labs(title = "",y="",x="")
  
  p1.2<-ggplot(data, aes_string(x='year.2',y=names(data)[i],fill='V.arrival'))+
    geom_boxplot(outlier.shape = NA)+
    scale_y_sqrt()+
    scale_fill_manual(values = c('white','grey'))+
    theme_classic()+
    theme(plot.margin=unit(c(1,2.2,0,-0.6),"lines"),legend.position = 'none')+
    labs(y="",x="")
  p1 <- ggarrange(p1.2,p1.1,nrow=2,ncol=1,heights=c(0.3, 0.7))
  p1 <- annotate_figure(p1,top=text_grob(paste(name_full[i-24])),
                        left = text_grob("Transcript per million",rot = 90),
                        bottom=text_grob("Varroa exposure time (years)"))
  print(p1)
  
}

#------------------------------------------------------------------------------------------------------
## exposure time groups: 0,<1,<2,>2
for(i in 25:41){
  p1.1 <- ggplot(data, aes_string(x='year.3',y=names(data)[i],fill='V.arrival'))+
    geom_boxplot(outlier.shape = NA)+
    scale_y_sqrt()+
    facet_grid(rows=vars(Geocluster))+
    scale_fill_manual(values = c('white','grey'))+
    theme_classic()+
    theme(legend.position="none",plot.margin=unit(c(-1.8,1,-0.5,-0.5),"lines"))+
    labs(title = "",y="",x="")
  
  p1.2<-ggplot(data, aes_string(x='year.3',y=names(data)[i],fill='V.arrival'))+
    geom_boxplot(outlier.shape = NA)+
    scale_y_sqrt()+
    scale_fill_manual(values = c('white','grey'))+
    theme_classic()+
    theme(plot.margin=unit(c(1,2.2,0,-0.6),"lines"),legend.position = "none")+
    labs(y="",x="")
  
  p1 <- ggarrange(p1.2,p1.1,nrow=2,ncol=1,heights=c(0.3, 0.7))
  p1 <- annotate_figure(p1,top=text_grob(paste(name_full[i-24])),
                        left = text_grob("Transcript per million",rot = 90),
                        bottom=text_grob("Varroa exposure time (years)"))
  print(p1)
  
}

#--------------------------------------------------------------------------------------------------------
## OVERALL
data_melt <- melt(data[24:41], id=c('V.arrival'))

ggplot(data_melt,aes(x=variable,y=value,fill=as.factor(V.arrival)))+
  geom_boxplot(outlier.shape = NA)+
  labs(x='',y='Transcript per million',title='Overall viral load')+
  scale_y_sqrt()+
  scale_fill_manual(labels=c('Before Varroa','After Varroa'),values=c('white','grey'))+
  theme_classic()+
  theme(legend.position = c(0.9,0.9),
        legend.background = element_rect(linetype="solid",colour ="darkgrey"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size=11))

for (i in 25:41) {
  box <- ggplot(data,aes_string(x='year.3',y=names(data)[i],fill='V.arrival'))+
    geom_boxplot(outlier.shape = NA)+
    labs(title='Overall viral load',subtitle=name_full[i-24],x='Varroa exposure time (years)',y='Transcript per million',fill='')+
    scale_y_sqrt()+
    scale_fill_manual(values=c('white','grey'))+
    theme_classic()+
    theme(legend.position = 'none')
  
  print(box)
  
}

#======================================================================================================
## prevalence #========================================================================================

## Time since first Varroa detection: -8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,5,7
data.1 <- data[,c(11,19,24,25:41)] %>% group_by(year.1,V.arrival,Geocluster) %>% summarise(ABPV=length(which(ABPV>0))/length(ABPV),
                                                                             Amfv=length(which(Amfv>0))/length(Amfv),
                                                                             Arv2=length(which(Arv2>0))/length(Arv2),
                                                                             BQCV=length(which(BQCV>0))/length(BQCV),
                                                                             CBPV=length(which(CBPV>0))/length(CBPV),
                                                                             DWV_A=length(which(DWV_A>0))/length(DWV_A),
                                                                             DWV_C=length(which(DWV_C>0))/length(DWV_C),
                                                                             IAPV=length(which(IAPV>0))/length(IAPV),
                                                                             KBV=length(which(KBV>0))/length(KBV),
                                                                             LSV=length(which(LSV>0))/length(LSV),
                                                                             SBPV=length(which(SBPV>0))/length(SBPV),
                                                                             SBV=length(which(SBV>0))/length(SBV),
                                                                             DWV_B=length(which(DWV_B>0))/length(DWV_B),
                                                                             VDV2=length(which(VDV2>0))/length(VDV2),
                                                                             VDV3=length(which(VDV3>0))/length(VDV3),
                                                                             Vov1=length(which(Vov1>0))/length(Vov1),
                                                                             VTV=length(which(VTV>0))/length(VTV))
data.1[,4:20] <- data.1[,4:20]*100

data.2 <- data[,c(11,24,25:41)] %>% group_by(year.1,V.arrival) %>% summarise(ABPV=length(which(ABPV>0))/length(ABPV),
                                                                                           Amfv=length(which(Amfv>0))/length(Amfv),
                                                                                           Arv2=length(which(Arv2>0))/length(Arv2),
                                                                                           BQCV=length(which(BQCV>0))/length(BQCV),
                                                                                           CBPV=length(which(CBPV>0))/length(CBPV),
                                                                                           DWV_A=length(which(DWV_A>0))/length(DWV_A),
                                                                                           DWV_C=length(which(DWV_C>0))/length(DWV_C),
                                                                                           IAPV=length(which(IAPV>0))/length(IAPV),
                                                                                           KBV=length(which(KBV>0))/length(KBV),
                                                                                           LSV=length(which(LSV>0))/length(LSV),
                                                                                           SBPV=length(which(SBPV>0))/length(SBPV),
                                                                                           SBV=length(which(SBV>0))/length(SBV),
                                                                                           DWV_B=length(which(DWV_B>0))/length(DWV_B),
                                                                                           VDV2=length(which(VDV2>0))/length(VDV2),
                                                                                           VDV3=length(which(VDV3>0))/length(VDV3),
                                                                                           Vov1=length(which(Vov1>0))/length(Vov1),
                                                                                           VTV=length(which(VTV>0))/length(VTV))
data.2[,3:19] <- data.2[,3:19]*100

for (i in 4:20){
  a<-ggplot(data.1[-which(is.na(data.1$year.1)),], aes_string(x='year.1', y=names(data.1)[i], fill='V.arrival'))+
    geom_bar(stat="identity",color='black')+
    ylim(c(0,100))+
    facet_grid(rows=vars(Geocluster))+
    labs(x="",y="")+
    scale_fill_manual(values=c('white','grey'))+
    theme_classic()+
    theme(legend.position = "none",plot.margin=unit(c(0,0.5,-0.3,0),"lines"),
          axis.text.y=element_text(size=8))+
    scale_x_discrete(labels=c('<8',"<7",'<6','<5','<4','<3','<2','<1','0','<1','<2','<3','<5','<7'))
  
  b<-ggplot(data.2[-which(is.na(data.2$year.1)),], aes_string(x='year.1', y=names(data.2)[i-1],fill='V.arrival'))+
    geom_bar(stat="identity",color='black')+
    ylim(c(0,100))+
    labs(x="",y="")+
    scale_fill_manual(values=c('white','grey'))+
    theme_classic()+
    theme(plot.margin=unit(c(1,1.6,-0.3,-0.1),"lines"),
          axis.text.y=element_text(size=8),legend.position = 'none')+
    scale_x_discrete(labels=c('<8',"<7",'<6','<5','<4','<3','<2','<1','0','<1','<2','<3','<5','<7'))
  
  c <- ggarrange(b,a,nrow=2,ncol=1,heights=c(0.3, 0.7))
  c <- annotate_figure(c,top=text_grob(paste("Prevalence of ",name_full[i-3]),size=12),
                       bottom=text_grob("Time since first Varroa detection (years)"),
                       left = text_grob("Pathogen prevalence (% positive colonies)",rot = 90),)
  print(c)
  
}

#------------------------------------------------------------------------------------------
## exposure time groups: 0,<1,<2,<3,<5,<7
data.3 <- data[,c(12,19,24,25:41)] %>% group_by(year.2,V.arrival,Geocluster) %>% summarise(ABPV=length(which(ABPV>0))/length(ABPV),
                                                                                           Amfv=length(which(Amfv>0))/length(Amfv),
                                                                                           Arv2=length(which(Arv2>0))/length(Arv2),
                                                                                           BQCV=length(which(BQCV>0))/length(BQCV),
                                                                                           CBPV=length(which(CBPV>0))/length(CBPV),
                                                                                           DWV_A=length(which(DWV_A>0))/length(DWV_A),
                                                                                           DWV_C=length(which(DWV_C>0))/length(DWV_C),
                                                                                           IAPV=length(which(IAPV>0))/length(IAPV),
                                                                                           KBV=length(which(KBV>0))/length(KBV),
                                                                                           LSV=length(which(LSV>0))/length(LSV),
                                                                                           SBPV=length(which(SBPV>0))/length(SBPV),
                                                                                           SBV=length(which(SBV>0))/length(SBV),
                                                                                           DWV_B=length(which(DWV_B>0))/length(DWV_B),
                                                                                           VDV2=length(which(VDV2>0))/length(VDV2),
                                                                                           VDV3=length(which(VDV3>0))/length(VDV3),
                                                                                           Vov1=length(which(Vov1>0))/length(Vov1),
                                                                                           VTV=length(which(VTV>0))/length(VTV))
data.3[,4:20] <- data.3[,4:20]*100

data.4 <- data[,c(12,24,25:41)] %>% group_by(year.2,V.arrival) %>% summarise(ABPV=length(which(ABPV>0))/length(ABPV),
                                                                             Amfv=length(which(Amfv>0))/length(Amfv),
                                                                             Arv2=length(which(Arv2>0))/length(Arv2),
                                                                             BQCV=length(which(BQCV>0))/length(BQCV),
                                                                             CBPV=length(which(CBPV>0))/length(CBPV),
                                                                             DWV_A=length(which(DWV_A>0))/length(DWV_A),
                                                                             DWV_C=length(which(DWV_C>0))/length(DWV_C),
                                                                             IAPV=length(which(IAPV>0))/length(IAPV),
                                                                             KBV=length(which(KBV>0))/length(KBV),
                                                                             LSV=length(which(LSV>0))/length(LSV),
                                                                             SBPV=length(which(SBPV>0))/length(SBPV),
                                                                             SBV=length(which(SBV>0))/length(SBV),
                                                                             DWV_B=length(which(DWV_B>0))/length(DWV_B),
                                                                             VDV2=length(which(VDV2>0))/length(VDV2),
                                                                             VDV3=length(which(VDV3>0))/length(VDV3),
                                                                             Vov1=length(which(Vov1>0))/length(Vov1),
                                                                             VTV=length(which(VTV>0))/length(VTV))
data.4[,3:19] <- data.4[,3:19]*100

for (i in 4:20){
  a<-ggplot(data.3, aes_string(x='year.2', y=names(data.3)[i], fill='V.arrival'))+
    geom_bar(stat="identity",color='black')+
    ylim(c(0,100))+
    facet_grid(rows=vars(Geocluster))+
    labs(x="",y="")+
    scale_fill_manual(values=c('white','grey'))+
    theme_classic()+
    theme(legend.position = "none",plot.margin=unit(c(0,0.5,-0.3,0),"lines"),
          axis.text.y=element_text(size=8))
  
  b<-ggplot(data.4, aes_string(x='year.2', y=names(data.4)[i-1], fill='V.arrival'))+
    geom_bar(stat="identity",color='black')+
    ylim(c(0,100))+
    labs(x="",y="")+
    scale_fill_manual(values=c('white','grey'))+
    theme_classic()+
    theme(legend.position='none',
          plot.margin=unit(c(1,1.6,-0.3,-0.1),"lines"),
          axis.text.y=element_text(size=8))
  
  c <- ggarrange(b,a,nrow=2,ncol=1,heights=c(0.3, 0.7))
  c <- annotate_figure(c,top=text_grob(paste("Prevalence of ",name_full[i-3]),size=12),
                       bottom=text_grob("Varroa exposure time (years)"),
                       left = text_grob("Pathogen prevalence (% positive colonies)",rot = 90),)
  print(c)
  
}

#-----------------------------------------------------------------------------------------
## exposure time groups: 0,<1,<2,>2
data.5 <- data[,c(13,19,24,25:41)] %>% group_by(year.3,V.arrival,Geocluster) %>% summarise(ABPV=length(which(ABPV>0))/length(ABPV),
                                                                                           Amfv=length(which(Amfv>0))/length(Amfv),
                                                                                           Arv2=length(which(Arv2>0))/length(Arv2),
                                                                                           BQCV=length(which(BQCV>0))/length(BQCV),
                                                                                           CBPV=length(which(CBPV>0))/length(CBPV),
                                                                                           DWV_A=length(which(DWV_A>0))/length(DWV_A),
                                                                                           DWV_C=length(which(DWV_C>0))/length(DWV_C),
                                                                                           IAPV=length(which(IAPV>0))/length(IAPV),
                                                                                           KBV=length(which(KBV>0))/length(KBV),
                                                                                           LSV=length(which(LSV>0))/length(LSV),
                                                                                           SBPV=length(which(SBPV>0))/length(SBPV),
                                                                                           SBV=length(which(SBV>0))/length(SBV),
                                                                                           DWV_B=length(which(DWV_B>0))/length(DWV_B),
                                                                                           VDV2=length(which(VDV2>0))/length(VDV2),
                                                                                           VDV3=length(which(VDV3>0))/length(VDV3),
                                                                                           Vov1=length(which(Vov1>0))/length(Vov1),
                                                                                           VTV=length(which(VTV>0))/length(VTV))
data.5[,4:20] <- data.5[,4:20]*100

data.6 <- data[,c(13,24,25:41)] %>% group_by(year.3,V.arrival) %>% summarise(ABPV=length(which(ABPV>0))/length(ABPV),
                                                                             Amfv=length(which(Amfv>0))/length(Amfv),
                                                                             Arv2=length(which(Arv2>0))/length(Arv2),
                                                                             BQCV=length(which(BQCV>0))/length(BQCV),
                                                                             CBPV=length(which(CBPV>0))/length(CBPV),
                                                                             DWV_A=length(which(DWV_A>0))/length(DWV_A),
                                                                             DWV_C=length(which(DWV_C>0))/length(DWV_C),
                                                                             IAPV=length(which(IAPV>0))/length(IAPV),
                                                                             KBV=length(which(KBV>0))/length(KBV),
                                                                             LSV=length(which(LSV>0))/length(LSV),
                                                                             SBPV=length(which(SBPV>0))/length(SBPV),
                                                                             SBV=length(which(SBV>0))/length(SBV),
                                                                             DWV_B=length(which(DWV_B>0))/length(DWV_B),
                                                                             VDV2=length(which(VDV2>0))/length(VDV2),
                                                                             VDV3=length(which(VDV3>0))/length(VDV3),
                                                                             Vov1=length(which(Vov1>0))/length(Vov1),
                                                                             VTV=length(which(VTV>0))/length(VTV))
data.6[,3:19] <- data.6[,3:19]*100

for (i in 4:20){
  a<-ggplot(data.5, aes_string(x='year.3', y=names(data.5)[i], fill='V.arrival'))+
    geom_bar(stat="identity",color='black')+
    ylim(c(0,100))+
    facet_grid(rows=vars(Geocluster))+
    labs(x="",y="")+
    scale_fill_manual(values=c('white','grey'))+
    theme_classic()+
    theme(legend.position = "none",plot.margin=unit(c(0,0.5,-0.3,0),"lines"),
          axis.text.y=element_text(size=8))
  
  b<-ggplot(data.6, aes_string(x='year.3', y=names(data.6)[i-1],fill='V.arrival'))+
    geom_bar(stat="identity",color='black')+
    ylim(c(0,100))+
    labs(x="",y="")+
    scale_fill_manual(values=c('white','grey'))+
    theme_classic()+
    theme(plot.margin=unit(c(1,1.6,-0.3,-0.1),"lines"),
          axis.text.y=element_text(size=8),legend.position = 'none')
  
  c <- ggarrange(b,a,nrow=2,ncol=1,heights=c(0.3, 0.7))
  c <- annotate_figure(c,top=text_grob(paste("Prevalence of ",name_full[i-3]),size=12),
                       bottom=text_grob("Varroa exposure time (years)"),
                       left = text_grob("Pathogen prevalence (% positive colonies)",rot = 90),)
  print(c)
  
}

#-----------------------------------------------------------------------------------------------
## OVERALL
data.7 <- data[,c(24,25:41)] %>% group_by(V.arrival) %>% summarise(ABPV=length(which(ABPV>0))/length(ABPV),
                                                                             Amfv=length(which(Amfv>0))/length(Amfv),
                                                                             Arv2=length(which(Arv2>0))/length(Arv2),
                                                                             BQCV=length(which(BQCV>0))/length(BQCV),
                                                                             CBPV=length(which(CBPV>0))/length(CBPV),
                                                                             DWV_A=length(which(DWV_A>0))/length(DWV_A),
                                                                             DWV_C=length(which(DWV_C>0))/length(DWV_C),
                                                                             IAPV=length(which(IAPV>0))/length(IAPV),
                                                                             KBV=length(which(KBV>0))/length(KBV),
                                                                             LSV=length(which(LSV>0))/length(LSV),
                                                                             SBPV=length(which(SBPV>0))/length(SBPV),
                                                                             SBV=length(which(SBV>0))/length(SBV),
                                                                             DWV_B=length(which(DWV_B>0))/length(DWV_B),
                                                                             VDV2=length(which(VDV2>0))/length(VDV2),
                                                                             VDV3=length(which(VDV3>0))/length(VDV3),
                                                                             Vov1=length(which(Vov1>0))/length(Vov1),
                                                                             VTV=length(which(VTV>0))/length(VTV))
data.7[,2:18] <- data.7[,2:18]*100


data.7_melt <- melt(data.7, id=c('V.arrival'))

ggplot(data.7_melt,aes(x=variable,y=value,fill=V.arrival))+
  geom_bar(stat='identity',position=position_dodge(),color='black')+
  labs(x='',y='Pathogen prevalence (% positive colonies)')+
  scale_fill_manual(labels=c('Before Varroa','After Varroa'),values=c('white','grey'))+
  ylim(c(0,100))+
  theme_classic()+
  theme(legend.position = c(0.9,0.9),
        legend.background = element_rect(linetype="solid",colour ="darkgrey"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.text.x = element_text(size=11))

for (i in 3:19) {
  a <-ggplot(data.6,aes_string(x='year.3',y=names(data.6)[i],fill='V.arrival'))+
    geom_bar(stat='identity',color='black')+
    scale_fill_manual(values = c('white','grey'))+
    labs(title=name_full[i-2],x='Varroa exposure time (years)',y='Pathogen prevalence (% positive colonies)')+
    ylim(c(0,100))+
    theme_classic()+
    theme(legend.position = 'none',axis.text.x = element_text(size=11))
  print(a)
}

