###########################################################################################
##                                                                                       ##   
## Author: Sherry Lin                                                                    ##
## Student ID: u6890762                                                                  ##
## Project: evolution of bee viruses                                                     ##
## Purpose:    ##
## Version: 8                                                                            ##
## Date: 01/07/20                                                                        ##
##                                                                                       ##
###########################################################################################

library(ggplot2)
library(ggpubr)

rm(list = ls())

data.1 <- read.csv("new_data_mean_Year.csv")


data.2 <- read.csv("new_data_mean_YearGeo.csv")
data.2$Geocluster <- ordered(data.2$Geocluster, level = c("WEL","BOR","TAM","VER"))

data.3 <- read.csv("new_data_Varroa_Year.csv")
data.3 <- data.3[which(data.3$Varroa=="YES"), ]

data.4 <- read.csv("new_data_Varroa_YearGeo.csv")
data.4 <- data.4[which(data.4$Varroa=="YES"), ]
data.4$Geocluster <- ordered(data.4$Geocluster, level = c("WEL","BOR","TAM","VER"))

name_list <- colnames(data.1[2:17])

name_full <- c('Acute bee paralysis virus','Apis mellifera filamentous virus',
               'Apis rhabdovirus 2','Black queen cell virus', 'Chronic bee paralysis virus',
               'Deformed wing virus','Israeli acute paralysis virus','Kashmir bee virus',
               'Lake Sinai virus','Slow bee paralysis virus','Sacbrood virus',
               'Varroa destructor virus 1','Varroa destructor virus 2','Varroa destructor virus 3',
               'Varroa orthomyxovirus-1','Varroa Tymo-like virus')

for (i in 1:16) {
  p1 <-ggplot()+
    geom_line(data=data.1,size=1,aes_string(x='Year',y=name_list[i], color='"Viral TPH"'))+
    geom_line(data=data.1,size=1,aes(x=Year,y=Afr.Extent, color="Africanization (%)"))+
    geom_line(data=data.3,size=1,aes(x=Year,y=percent,color="Presece of Varroa (%)"))+
    geom_point(data=data.1,aes_string(x='Year',y=name_list[i],shape='"Viral TPH"'))+
    geom_point(data=data.1,aes(x=Year,y=Afr.Extent,shape="Africanization (%)"))+
    geom_point(data=data.3,aes(x=Year,y=percent,shape="Presece of Varroa (%)"))+
    theme(legend.position ="none",plot.margin=unit(c(1,2.2,-1,0.8),"lines"))+
    labs(subtitle = "Overview",x="",y="")+
    scale_y_sqrt()+
    scale_color_manual(values=c("#0b5db5","#737373","#EAB509"))
  
  p2 <-ggplot()+
    geom_line(data=data.2,size=1,aes_string(x='Year',y=name_list[i],color='"Viral TPH"'))+
    geom_line(data=data.2,size=1,aes(x=Year,y=Afr.Extent, color="Afrricanization (%)"))+
    geom_line(data=data.4,size=1,aes(x=Year,y=percent,color="Presece of Varroa (%)"))+
    geom_point(data=data.2,aes_string(x='Year',y=name_list[i],shape='"Viral TPH"'))+
    geom_point(data=data.2,aes(x=Year,y=Afr.Extent,shape="Africanization (%)"))+
    geom_point(data=data.4,aes(x=Year,y=percent,shape="Presece of Varroa (%)"))+
    facet_grid(rows = vars(Geocluster))+
    theme(legend.position = "none",plot.margin=unit(c(-0.7,1,0,0.5),"lines"))+
    labs(subtitle = "",x="Year",y="")+
    scale_y_sqrt()+
    scale_color_manual(values=c("#0b5db5","#737373","#EAB509"))

  p <- ggarrange(p1,p2,nrow=2,ncol=1,heights=c(0.3, 0.7),common.legend = T,legend = "bottom")
  p <- annotate_figure(p,top=text_grob(paste(name_full[i])))
  print(p)

}



  



