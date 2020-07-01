##############################################################################
##                                                                          ##   
## Author: Sherry Lin                                                       ##
## Student ID: u6890762                                                     ##
## Project: evolution of bee viruses                                        ##
## Purpose: Varroa presence, using different data to see when Varroa had    ##
##          arrived at various sites                                        ##
## Version: 15                                                              ##
## Date: 01/07/20                                                           ##
##                                                                          ##
##############################################################################


library(ggplot2)
library(ggsci)  #scale_color_jco() / scale_fill_jco(): color in Journal of Clinical Oncology
library(ggpubr) #ggarrange()

rm(list = ls())

# Big Dataset Honeybees US
data.big_c <- read.csv("new_BigDataset_Varroa_YearCountry.csv")
data.big_c <- data.big_c[which(data.big_c$Varroa=="YES"), ]
data.big_c$Country <- ordered(data.big_c$Country, level = c("BELIZE","MEXICO","USA","Unknown"))

data.big_s <- read.csv("new_BigDataset_Varroa_YearState.csv")
data.big_s <- data.big_s[which(data.big_s$Varroa=="YES"), ]
data.big_s$State <- ordered(data.big_s$State, level = c("BELIZE","CAMPECHE","VERACRUZ","QUINTANA ROO","TAMAULIPAS","TEXAS","Unknown"))

data.big_g <- read.csv("new_BigDataset_Varroa_YearGeo.csv")
data.big_g <- data.big_g[which(data.big_g$Varroa=="YES"), ]
data.big_g$Geocluster <- ordered(data.big_g$Geocluster, level = c("CAN","BEL","CAM","VER","TAM","BOR","WEL","HOU","SAT","DAL","Unknown"))

p1 <-ggplot(data.big_c, aes(x=Year, y=percent, fill=Country))+ 
  geom_bar(stat = "identity")+
  facet_grid(.~Country)+
  labs(title="Country",x="",y="")+
  geom_text(aes(label=percent),position = position_dodge(0.8),
            vjust = -0.3, size = 3.5)+
  scale_fill_jco()+
  theme_light()

p2 <-ggplot(data.big_s, aes(x=Year, y=percent, fill=State))+ 
  geom_bar(stat = "identity")+
  facet_grid(.~State)+
  labs(title="State",x="",y="")+
  geom_text(aes(label=percent),position = position_dodge(0.8),
            vjust = -0.3, size = 3.5)+
  scale_fill_jco()+
  theme_light()

p3 <-ggplot(data.big_g, aes(x=Year, y=percent, fill=Geocluster))+ 
  geom_bar(stat = "identity")+
  facet_grid(.~Geocluster)+
  labs(title="Geocluster",x="",y="")+
  geom_text(aes(label=percent),position = position_dodge(0.8),
            vjust = -0.3, size = 3.5)+
  scale_fill_jco()+
  theme_light()

pVarroa.big <-ggarrange(ggarrange(p1,p2,nrow=1,ncol=2,labels=c("A","B"),font.label=list(size=11)),
                    p3,nrow=2,ncol=1,labels=c("","C"),font.label=list(size=11))

pVarroa.big <- annotate_figure(pVarroa.big,top=text_grob("Varroa Presence",size=15),
                           bottom=text_grob("Year"),
                           left = text_grob("Percentage (%)",rot = 90),)
pVarroa.big

####################################################################################################
# meta data
data.meta_c <- read.csv("new_meta_Varroa_YearCountry.csv")
data.meta_c <- data.meta_c[which(data.meta_c$Varroa=="YES"), ]
data.meta_c$Country <- ordered(data.meta_c$Country, level = c("BELIZE","MEXICO","USA"))

data.meta_s <- read.csv("new_meta_Varroa_YearState.csv")
data.meta_s <- data.meta_s[which(data.meta_s$Varroa=="YES"), ]
data.meta_s$State <- ordered(data.meta_s$State, level = c("BELIZE", "CAMPECHE","VERACRUZ","QUINTANA ROO","TAMAULIPAS","TEXAS"))

data.meta_g <- read.csv("new_meta_Varroa_YearGeo.csv")
data.meta_g <- data.meta_g[which(data.meta_g$Varroa=="YES"), ]
data.meta_g$Geocluster <- ordered(data.meta_g$Geocluster, level = c("CAN","BEL","CAM","VER","TAM","BOR","WEL","HOU","SAT"))


p4 <-ggplot(data.meta_c, aes(x=Year, y=percent, fill=Country))+ 
  geom_bar(stat = "identity")+
  facet_grid(.~Country)+
  labs(title="Country",x="",y="")+
  geom_text(aes(label=percent),position = position_dodge(0.8),
            vjust = -0.3, size = 3.5)+
  scale_fill_jco()+
  theme_light()

p5 <-ggplot(data.meta_s, aes(x=Year, y=percent, fill=State))+ 
  geom_bar(stat = "identity")+
  facet_grid(.~State)+
  labs(title="State",x="",y="")+
  geom_text(aes(label=percent),position = position_dodge(0.8),
            vjust = -0.3, size = 3.5)+
  scale_fill_jco()+
  theme_light()

p6 <-ggplot(data.meta_g, aes(x=Year, y=percent, fill=Geocluster))+ 
  geom_bar(stat = "identity")+
  facet_grid(.~Geocluster)+
  labs(title="Geocluster",x="",y="")+
  geom_text(aes(label=percent),position = position_dodge(0.8),
            vjust = -0.3, size = 3.5)+
  scale_fill_jco()+
  theme_light()

pVarroa.meta <-ggarrange(ggarrange(p4,p5,nrow=1,ncol=2,labels=c("A","B"),font.label=list(size=11)),
                        p6,nrow=2,ncol=1,labels=c("","C"),font.label=list(size=11))

pVarroa.meta <- annotate_figure(pVarroa.meta,top=text_grob("Varroa Presence",size=15),
                               bottom=text_grob("Year"),
                               left = text_grob("Percentage (%)",rot = 90),)
pVarroa.meta

###############################################################################################
# sequenced information
data_c <- read.csv("new_data_Varroa_YearCountry.csv")
data_c <- data_c[which(data_c$Varroa=="YES"), ]
data_c$Country <- ordered(data_c$Country, level = c("Unknown","MEXICO","USA"))

data_s <- read.csv("new_data_Varroa_YearState.csv")
data_s <- data_s[which(data_s$Varroa=="YES"), ]
data_s$State <- ordered(data_s$State, level = c("VERACRUZ","TAMAULIPAS","TEXAS"))

data_g <- read.csv("new_data_Varroa_YearGeo.csv")
data_g <- data_g[which(data_g$Varroa=="YES"), ]
data_g$Geocluster <- ordered(data_g$Geocluster, level = c("VER","TAM","BOR","WEL"))

p7 <-ggplot(data_c, aes(x=Year, y=percent, fill=Country))+ 
  geom_bar(stat = "identity")+
  facet_grid(.~Country)+
  labs(title="Country",x="",y="")+
  geom_text(aes(label=percent),position = position_dodge(0.8),
            vjust = -0.3, size = 3.5)+
  scale_fill_jco()+
  theme_light()

p8 <-ggplot(data_s, aes(x=Year, y=percent, fill=State))+ 
  geom_bar(stat = "identity")+
  facet_grid(.~State)+
  labs(title="State",x="",y="")+
  geom_text(aes(label=percent),position = position_dodge(0.8),
            vjust = -0.3, size = 3.5)+
  scale_fill_jco()+
  theme_light()

p9 <-ggplot(data_g, aes(x=Year, y=percent, fill=Geocluster))+ 
  geom_bar(stat = "identity")+
  facet_grid(.~Geocluster)+
  labs(title="Geocluster",x="",y="")+
  geom_text(aes(label=percent),position = position_dodge(0.8),
            vjust = -0.3, size = 3.5)+
  scale_fill_jco()+
  theme_light()

pvarroa <-ggarrange(p7,p8,p9,nrow=2,ncol=2,labels=c("A","B","C"),font.label=list(size=11))

pvarroa <- annotate_figure(pvarroa,top=text_grob("Varroa Presence",size=15),
                           bottom=text_grob("Year"),
                           left = text_grob("Percentage (%)",rot = 90),)
pvarroa













