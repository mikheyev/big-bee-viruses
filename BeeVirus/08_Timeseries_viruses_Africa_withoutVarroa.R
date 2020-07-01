library(dplyr)
library(ggplot2)
library(ggpubr)
library(GGally)

rm(list = ls())

data <- read.csv("new_data.csv")

data_NoV <- data[which(data$V.arrival==0), ]
data_NoV$Geocluster <- ordered(data_NoV$Geocluster, level = c("WEL","BOR","TAM","VER"))


data_NoV.a <- group_by(data_NoV, Geocluster, Year) %>% summarise(ABPV=mean(ABPV, na.rm=T ), AM=mean(AM, na.rm=T), 
                                                                 AR=mean(AR, na.rm=T), BQCV=mean(BQCV, na.rm=T),
                                                                 CBPV=mean(CBPV,na.rm=T),DWV=mean(DWV, na.rm=T), 
                                                                 IAPV=mean(IAPV, na.rm=T),KBV=mean(KBV, na.rm=T),
                                                                 LSV=mean(LSV, na.rm=T),SBPV=mean(SBPV, na.rm=T),
                                                                 SV=mean(SV, na.rm=T),VDV1=mean(VDV1, na.rm=T),  
                                                                 VDV2=mean(VDV2, na.rm=T),VDV3=mean(VDV3, na.rm=T),    
                                                                 VO1=mean(VO1, na.rm=T),VTV=mean(VTV, na.rm=T),
                                                                 Afr.Extent = mean(Afr.Extent, na.rm = T))


data_NoV.aa <- group_by(data_NoV, Year) %>% summarise(ABPV=mean(ABPV, na.rm=T ), AM=mean(AM, na.rm=T), 
                                                      AR=mean(AR, na.rm=T), BQCV=mean(BQCV, na.rm=T),
                                                      CBPV=mean(CBPV,na.rm=T),DWV=mean(DWV, na.rm=T), 
                                                      IAPV=mean(IAPV, na.rm=T),KBV=mean(KBV, na.rm=T),
                                                      LSV=mean(LSV, na.rm=T),SBPV=mean(SBPV, na.rm=T),
                                                      SV=mean(SV, na.rm=T),VDV1=mean(VDV1, na.rm=T),  
                                                      VDV2=mean(VDV2, na.rm=T),VDV3=mean(VDV3, na.rm=T),    
                                                      VO1=mean(VO1, na.rm=T),VTV=mean(VTV, na.rm=T),
                                                      Afr.Extent = mean(Afr.Extent, na.rm = T))


name_list <- colnames(data_NoV.a[3:18])
name_full <- c('Acute bee paralysis virus','Apis mellifera filamentous virus',
               'Apis rhabdovirus 2','Black queen cell virus', 'Chronic bee paralysis virus',
               'Deformed wing virus','Israeli acute paralysis virus','Kashmir bee virus',
               'Lake Sinai virus','Slow bee paralysis virus','Sacbrood virus',
               'Varroa destructor virus 1','Varroa destructor virus 2','Varroa destructor virus 3',
               'Varroa orthomyxovirus-1','Varroa Tymo-like virus')


for(i in 1:16){
  p1.1 <- ggplot(data_NoV.a, aes(x=Year))+
    geom_line(size=1,aes_string(y=name_list[i], colour='"Virus"'))+
    geom_line(size=1,aes(y=Afr.Extent,color="Africanization (%)"))+
    geom_point(aes_string(y=name_list[i],shape='"Viral TPH"'))+
    geom_point(aes(y=Afr.Extent, shape="Aricanization (%)"))+
    scale_y_sqrt()+
    facet_grid(rows=vars(Geocluster))+
    theme(legend.position = "none",plot.margin=unit(c(-1.8,1,0,0.5),"lines"))+
    labs(title = "",y="",x="Year")+scale_color_jco()
  
  p1.2<-ggplot(data_NoV.aa, aes(x=Year))+
    geom_line(size=1,aes_string(y=name_list[i], colour='"Viral TPH"'))+
    geom_line(size=1,aes(y=Afr.Extent,color="Africanization (%)"))+
    geom_point(aes_string(y=name_list[i], shape='"Viral TPH"'))+
    geom_point(aes(y=Afr.Extent, shape="Aricanization (%)"))+
    scale_y_sqrt()+
    theme(legend.position="none",plot.margin=unit(c(1,2.1,0,0.5),"lines"))+
    labs(subtitle = "Overview",y="",x="")+scale_color_jco()
  p1 <- ggarrange(p1.2,p1.1,nrow=2,ncol=1,heights=c(0.3, 0.7),common.legend = T, legend = "bottom")
  p1 <- annotate_figure(p1,top=text_grob(paste(name_full[i])))
  print(p1)
  
}

