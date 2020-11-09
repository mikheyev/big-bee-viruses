###########################################################################################
##                                                                                       ##   
## Author: Sherry Lin                                                                    ##
## Student ID: u6890762                                                                  ##
## Project: evolution of bee viruses                                                     ##
## Purpose: plotting how the bee genetic background (africanization) change over time    ##
## Version: 15                                                                           ##
## Date: 22/09/20                                                                        ##
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

data <- read.csv("00_seq_info.csv")
names(data)[1] <- "ID"
data <- merge(data, africanized, by = "ID")
data$Geocluster <- factor(data$Geocluster, level=c("WEL","BOR","TAM","VER")) 

# figure out which of the columns most likely corresponds to higher levels of African introgression
E1 <- ggplot(data,aes(x=Year,y=V2,group=Year,color=Geocluster))+
  geom_boxplot()+ facet_grid(rows=vars(Geocluster))+scale_color_jco();E1
E2 <- ggplot(data,aes(x=Year,y=V3,group=Year,color=Geocluster))+
  geom_boxplot()+ facet_grid(rows=vars(Geocluster))+scale_color_jco();E2

# according to E1 & E2, Extent2 most likely corresponds to higher levels of African introgression                                                                                                                        
data <- read.csv('new_data.csv')
data.1 <- data[ ,c(4,19,42)] %>% group_by(Year,Geocluster) %>% summarise_all(mean,na.rm=T)
data.2 <- data[ ,c(4,42)] %>% group_by(Year) %>% summarise_all(mean,na.rm=T)
data.1$Geocluster <- factor(data.1$Geocluster, level=c("WEL","BOR","TAM","VER"))


## BAR CHART
p1 <- ggplot(data.1, aes(x=Year, y=Afr.Extent, fill=Geocluster))+                                      
  geom_bar(stat = "identity")+                                                                                           
  facet_grid(rows=vars(Geocluster))+
  theme(axis.title = element_blank(),legend.position = "none",plot.margin=unit(c(-0.5,0.5,0.5,0.5),"lines"))+
  scale_fill_jco()                                    
                                                                                                                         
p2 <- ggplot(data.2, aes(x=Year, y=Afr.Extent))+                   
  geom_bar(stat = "identity")+
  labs(subtitle='',x="",y="")+
  theme(plot.margin = unit(c(0,1.7,0,-0.5),"lines"))
                                                                                                                         
                                                                                                                         
pAfr <- ggarrange(p2,p1,nrow=2,heights=c(0.3, 0.7))

pAfr <- annotate_figure(pAfr,top=text_grob("Africanization over time",size=13),
                        bottom=text_grob("Year"),
                        left = text_grob("Average of extent (%)",rot = 90),)
pAfr
                                                                                                                         

## LINE+BAR CHART
p3 <- ggplot(data.1, aes(x=Year, y=Afr.Extent, fill=Geocluster))+                                      
  geom_bar(stat = "identity")+
  labs(x='Year',y='')+
  facet_grid(rows=vars(Geocluster))+
  xlim(c(1988,2002))+
  scale_fill_jco()                                    

p4 <- ggplot(data.1, aes(x=Year, y=Afr.Extent, group=Geocluster, color=Geocluster))+                   
  geom_line(size=1)+
  labs(subtitle='')+
  geom_point()+
  xlim(c(1988,2002))+
  theme(axis.title = element_blank(),plot.margin = unit(c(0,2,0,1.3),"lines"))+
  scale_color_jco()                                  


pAfr.1 <- ggarrange(p4,p3+theme(plot.margin=margin(t=10,l=5,r=12)),nrow=2,heights=c(0.25, 0.75),common.legend=T,legend="bottom")
pAfr.1 <- annotate_figure(pAfr.1,top=text_grob("Africanization over time",size=13),
                        left = text_grob("Average of extent (%)",rot = 90),)
pAfr.1

