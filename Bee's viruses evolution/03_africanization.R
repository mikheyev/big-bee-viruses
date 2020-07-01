###########################################################################################
##                                                                                       ##   
## Author: Sherry Lin                                                                    ##
## Student ID: u6890762                                                                  ##
## Project: evolution of bee viruses                                                     ##
## Purpose: plotting how the bee genetic background (africanization) change over time    ##
## Version: 14                                                                           ##
## Date: 01/07/20                                                                        ##
##                                                                                       ##
###########################################################################################

library(ggplot2) #ggplot()+geom_boxplot(), geom_bar(), geom_line, geom_point()
library(ggsci) #scale_color_jco() / scale_fill_jco(): color in Journal of Clinical Oncology
library(ggpubr) #ggarrange(): make multiple ggplots into one / annotate_figure()
library(dplyr) #group_by(): groups a dataframe based on certain fields / summarise()


rm(list = ls())

africanized <- read.table("00_k2.txt")                                                                                      
names(africanized)[1] <- "ID"                                                                                            
africanized$ID <- as.character(africanized$ID)                                                                           
africanized$ID[759] <- "334183"    

data <- read.csv("seq_info.csv")
names(data)[1] <- "ID"
data <- merge(data, africanized, by = "ID")
data$Geocluster <- ordered(data$Geocluster, level=c("WEL","BOR","TAM","VER")) 

# figure out which of the columns most likely corresponds to higher levels of African introgression
E1 <- ggplot(data,aes(x=Year,y=V2,group=Year,color=Geocluster))+
  geom_boxplot()+ geom_point()+ facet_grid(rows=vars(Geocluster))+scale_color_jco()
E2 <- ggplot(data,aes(x=Year,y=V3,group=Year,color=Geocluster))+
  geom_boxplot()+ geom_point()+ facet_grid(rows=vars(Geocluster))+scale_color_jco()

# according to E1 & E2, Extent2 most likely corresponds to higher levels of African introgression                                                                                                                        
data.1 <- read.csv("new_data_mean_YearGeo.csv") 
data.1$Geocluster <- ordered(data.1$Geocluster, level=c("WEL","BOR","TAM","VER"))
data.2 <- read.csv("new_data_mean_Year.csv")

## BAR CHART
p1 <- ggplot(data.1, aes(x=Year, y=Afr.Extent, fill=Geocluster))+                                      
  geom_bar(stat = "identity")+                                                                                           
  facet_grid(rows=vars(Geocluster))+
  theme(axis.title = element_blank(),legend.position = "none",plot.margin=unit(c(-0.5,0.5,0.5,0.5),"lines"))+
  scale_fill_jco()                                    
                                                                                                                         
p2 <- ggplot(data.2, aes(x=Year, y=Afr.Extent))+                   
  geom_bar(stat = "identity")+
  labs(subtitle="Overview",x="",y="")+
  theme(plot.margin = unit(c(0,1.7,0,-0.5),"lines"))
                                                                                                                         
                                                                                                                         
pAfr <- ggarrange(p2,p1,nrow=2,heights=c(0.3, 0.7))

pAfr <- annotate_figure(pAfr,top=text_grob("Africanization over time",size=13),
                        bottom=text_grob("Year"),
                        left = text_grob("Average of extent",rot = 90),)
pAfr
                                                                                                                         

## LINE+BAR CHART
p3 <- ggplot(data.1, aes(x=Year, y=Afr.Extent, fill=Geocluster))+                                      
  geom_bar(stat = "identity")+                                                                                           
  facet_grid(rows=vars(Geocluster))+
  theme(axis.title = element_blank())+
  scale_fill_jco()                                    

p4 <- ggplot(data.1, aes(x=Year, y=Afr.Extent, group=Geocluster, color=Geocluster))+                   
  geom_line(size=1)+                                                                                                     
  geom_point()+
  theme(axis.title = element_blank())+
  scale_color_jco()                                  


pAfr.1 <- ggarrange(p4,p3+theme(plot.margin=margin(t=10,l=5,r=12)),nrow=2,heights=c(0.3, 0.7),
                  labels = c("A","B"),common.legend=T,legend="bottom",font.label=list(size=11))
pAfr.1 <- annotate_figure(pAfr.1,top=text_grob("Africanization over time",size=13),
                        bottom=text_grob("Year"),
                        left = text_grob("Average of extent",rot = 90),)
pAfr.1
                                                                     


