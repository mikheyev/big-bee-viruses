##########################################################################
##                                                                      ##   
## Author: Sherry Lin                                                   ##
## Student ID: u6890762                                                 ##
## Project: evolution of bee viruses                                    ##
## Purpose: viral community composition over year/Varroa exposure time  ##
## Version: 14                                                          ##
## Date: 10/11/20                                                       ##
##                                                                      ##
##########################################################################

library(dplyr)  #group_by(), summarise()
library(ggplot2)  #theme_set()
library(rnaturalearth)  #ne_countries():returns world country polygons at a specified scale
library(rgeos)
library(scatterpie)  #geom_scatterpie()
library(ggrepel)  #geom_label_repel()


rm(list = ls())

data <- read.csv("new_data.csv")
data.1 <- data[,c(4,19,25:41)] %>% group_by(Year,Geocluster) %>% summarise_all(mean)

data.1$Longitude <- 0
data.1$Latitude <- 0
data.1$Longitude[data.1$Geocluster=="VER"] <- -97.30746  
data.1$Latitude[data.1$Geocluster=="VER"] <- 20.97450
data.1$Longitude[data.1$Geocluster=="TAM"] <- -98.20793
data.1$Latitude[data.1$Geocluster=="TAM"] <- 23.7158
data.1$Longitude[data.1$Geocluster=="BOR"] <- -98.12306
data.1$Latitude[data.1$Geocluster=="BOR"] <- 26.18369
data.1$Longitude[data.1$Geocluster=="WEL"] <- -97.44197
data.1$Latitude[data.1$Geocluster=="WEL"] <- 28.12102 

data.1$radius <- 1

#map setting
theme_set(theme_bw())
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
map <- ggplot(data = world)+
    geom_sf()+
    annotate(geom = "text", x = -90, y = 25, label = "Gulf of Mexico", 
             fontface = "italic", color = "grey22", size = 6)+
    annotate(geom = "text", x = -102.5, y = 22.5, label = "Mexico", color = "darkblue", size = 4)+
    annotate(geom = "text", x = -97.5, y = 32.5, label = "United States", color = "darkblue", size = 4)+
    annotate(geom = "text", x = -87.75, y = 17.5, label = "Bellze", color = "darkblue", size = 4)+
    annotate(geom = "text", x = -90.5, y = 15.5, label = "Guatemala", color = "darkblue", size = 4)+
    coord_sf(xlim = c(-105, -85), ylim = c(7.65, 33.97), expand = FALSE)+
    labs(x="Longitude", y="Latitude")

## map + viral community composition over year
year_list <- c("1988","1989","1990","1991","1992","1993","1994","1995","1996","1997","2001")

for(i in 1:11){
    map_year <- map+
        geom_scatterpie(data=data.1[which(data.1$Year==year_list[i]), ],aes(x=Longitude, y=Latitude, r=radius),
                        cols=colnames(data.1)[3:19],legend_name = "Viruses")+
        scale_fill_manual(values=c("DWV_A"="#1F78B4","KBV"="#2CA02C","CBPV"="#D62728",
                                   "ABPV"="#9467BD","LSV"="#E377C2","BQCV"="#7F7F7F",
                                   "SBV"="#DCBD22","Amfv"="#E77C7C","Arv2"="#5FD35F",
                                   "Vov1"="#BC8B81","VTV"="#C6AEDC","DWV_B"="#17BECF",
                                   "VDV2"="#57A9E2","VDV3"="#FFB574","IAPV"="#FF7F0E",
                                   "SBPV"="#8C564B","DWV_C"="#fee08b"))+
        geom_label_repel(aes(x=Longitude, y=Latitude, label=Geocluster), data=data.1[which(data.1$Year==year_list[i]), ],
                         point.padding=unit(2, "lines"))+
        labs(title=paste("Year: ",year_list[i]))
    print(map_year)
    
}

#-------------------------------------------------------------------------------------------------------------------------
## map + viral community composition over Varroa exposure time
data.2 <- data[,c(13,19,25:41)] %>% group_by(year.3,Geocluster) %>% summarise_all(mean)

data.2$Longitude <- 0
data.2$Latitude <- 0
data.2$Longitude[which(data.2$Geocluster=="VER")] <- -97.30746  
data.2$Latitude[which(data.2$Geocluster=="VER")] <- 20.97450
data.2$Longitude[which(data.2$Geocluster=="TAM")] <- -98.20793
data.2$Latitude[which(data.2$Geocluster=="TAM")] <- 23.7158
data.2$Longitude[which(data.2$Geocluster=="BOR")] <- -98.12306
data.2$Latitude[which(data.2$Geocluster=="BOR")] <- 26.18369
data.2$Longitude[which(data.2$Geocluster=="WEL")] <- -97.44197
data.2$Latitude[which(data.2$Geocluster=="WEL")] <- 28.12102 

data.2$radius <- 1

year.3_list <- c('0','<1','<2','>2')

for(i in 1:4){
    map_year <- map+
        geom_scatterpie(data=data.2[which(data.2$year.3==year.3_list[i]), ],aes(x=Longitude, y=Latitude, r=radius),
                        cols=colnames(data.2)[3:19],legend_name = "Viruses")+
        scale_fill_manual(values=c("DWV_A"="#1F78B4","KBV"="#2CA02C","CBPV"="#D62728",
                                   "ABPV"="#9467BD","LSV"="#E377C2","BQCV"="#7F7F7F",
                                   "SBV"="#DCBD22","Amfv"="#E77C7C","Arv2"="#5FD35F",
                                   "Vov1"="#BC8B81","VTV"="#C6AEDC","DWV_B"="#17BECF",
                                   "VDV2"="#57A9E2","VDV3"="#FFB574","IAPV"="#FF7F0E",
                                   "SBPV"="#8C564B","DWV_C"="#fee08b"))+
        geom_label_repel(aes(x=Longitude, y=Latitude, label=Geocluster), data=data.2[which(data.2$year.3==year.3_list[i]), ],
                         point.padding=unit(2, "lines"))+
        labs(title=paste("Varroa exposure time: ",year.3_list[i], 'years'))
    print(map_year)
    
}

#-----------------------------------------------------------------------------------------------------------------------------------------
## donut chart, time since Varroa first time detected: -8~7 years
rm(list = ls())

data <- read.csv("new_data.csv")
data <- data[ ,c(11,25:41)] %>% group_by(year.1) %>% summarise_all(mean)

data.1 <- as.data.frame(t(data[1:14,]))
data.1 <- data.1[-1, ]
data.1$virus <- rownames(data.1)
rownames(data.1) <- c(1:17)

year_list <- c("-8","-7","-6","-5","-4","-3","-2","-1","0","1","2","3","5","7")

for (i in 1:14){
    a<-ggplot(data.1, aes(x=2, y=data.1[,i], fill=virus))+
        geom_bar(stat='identity',color='white')+
        scale_fill_manual(values=c("DWV_A"="#1F78B4","KBV"="#2CA02C","CBPV"="#D62728",
                                   "ABPV"="#9467BD","LSV"="#E377C2","BQCV"="#7F7F7F",
                                   "SBV"="#DCBD22","Amfv"="#E77C7C","Arv2"="#5FD35F",
                                   "Vov1"="#BC8B81","VTV"="#C6AEDC","DWV_B"="#17BECF",
                                   "VDV2"="#57A9E2","VDV3"="#FFB574","IAPV"="#FF7F0E",
                                   "SBPV"="#8C564B","DWV_C"="#fee08b"))+
        coord_polar(theta='y',start=0)+
        annotate(geom = 'text', x = 0.5, y = 0, label=paste(year_list[i],"years after Varroa arrived",'\n',"Composition of viruses"),size=4)+
        theme_void()+
        xlim(0.5,2.5)
    
    print(a)
}

#==================================================================================================================
## donut chart, Varroa exposure time: 0,1,2,3,5,7
rm(list = ls())

data <- read.csv("new_data.csv", header = TRUE)
data <- data[ ,c(12,25:41)] %>% group_by(year.2) %>% summarise_all(mean)

data.1 <- as.data.frame(t(data))
data.1 <- data.1[-1, ]
data.1$virus <- rownames(data.1)
rownames(data.1) <- c(1:17)
data.1[ ,1:6] <- as.numeric(unlist(data.1[ ,1:6]))

year_list <- c("<1","<2","<3","<5","<7","0")

for (i in 1:6){
    a<-ggplot(data.1, aes(x=2, y=data.1[,i], fill=virus))+
        geom_bar(stat='identity',color='white')+
        scale_fill_manual(values=c("DWV_A"="#1F78B4","KBV"="#2CA02C","CBPV"="#D62728",
                                   "ABPV"="#9467BD","LSV"="#E377C2","BQCV"="#7F7F7F",
                                   "SBV"="#DCBD22","Amfv"="#E77C7C","Arv2"="#5FD35F",
                                   "Vov1"="#BC8B81","VTV"="#C6AEDC","DWV_B"="#17BECF",
                                   "VDV2"="#57A9E2","VDV3"="#FFB574","IAPV"="#FF7F0E",
                                   "SBPV"="#8C564B","DWV_C"="#fee08b"))+
        coord_polar(theta='y',start=0)+
        annotate(geom = 'text', x = 0.5, y = 0, label=paste("Varroa exposure time: ",year_list[i],"years",'\n',"Composition of viruses"),size=4)+
        theme_void()+
        xlim(0.5,2.5)
    
    print(a)
}

#-------------------------------------------------------------------------------------------------------------
## donut chart, Varroa exposure time: 0,<1,<2,>2
rm(list = ls())

data <- read.csv("new_data.csv", header = TRUE)
data <- data[ ,c(13,25:41)] %>% group_by(year.3) %>% summarise_all(mean)

data.1 <- as.data.frame(t(data))
data.1 <- data.1[-1, ]
data.1$virus <- rownames(data.1)
rownames(data.1) <- c(1:17)
data.1[ ,1:4] <- as.numeric(unlist(data.1[ ,1:4]))

year_list <- c("<1","<2",">2","0")

for (i in 1:4){
    a<-ggplot(data.1, aes(x=2, y=data.1[,i], fill=virus))+
        geom_bar(stat='identity',color='white')+
        scale_fill_manual(values=c("DWV_A"="#1F78B4","KBV"="#2CA02C","CBPV"="#D62728",
                                   "ABPV"="#9467BD","LSV"="#E377C2","BQCV"="#7F7F7F",
                                   "SBV"="#DCBD22","Amfv"="#E77C7C","Arv2"="#5FD35F",
                                   "Vov1"="#BC8B81","VTV"="#C6AEDC","DWV_B"="#17BECF",
                                   "VDV2"="#57A9E2","VDV3"="#FFB574","IAPV"="#FF7F0E",
                                   "SBPV"="#8C564B","DWV_C"="#fee08b"))+
        coord_polar(theta='y',start=0)+
        annotate(geom = 'text', x = 0.5, y = 0, label=paste("Composition of viruses",'\n',"Varroa exposure time: ",year_list[i],"years"),size=4.3)+
        theme_void()+
        xlim(0.5,2.5)
    
    print(a)
}
#-------------------------------------------------------------------------------------------------------------
## donut chart, Varroa exposure time: 0,<1,<2,>2 in different Geocluster
rm(list = ls())

data <- read.csv("new_data.csv", header = TRUE)
data <- data[ ,c(13,19,25:41)] %>% group_by(year.3,Geocluster) %>% summarise_all(mean)

data.1 <- as.data.frame(t(data))
data.1 <- data.1[-c(1,2), ]
data.1$virus <- rownames(data.1)
rownames(data.1) <- c(1:17)
data.1[ ,1:11] <- as.numeric(unlist(data.1[ ,1:11]))

year_list <- c('BOR <1','TAM <1','WEL <1','BOR <2','TAM <2','BOR >2','WEL >2','BOR 0','TAM 0','VER 0','WEL 0')

for (i in 1:11){
    a<-ggplot(data.1, aes(x=2, y=data.1[,i], fill=virus))+
        geom_bar(stat='identity',color='white')+
        scale_fill_manual(values=c("DWV_A"="#1F78B4","KBV"="#2CA02C","CBPV"="#D62728",
                                   "ABPV"="#9467BD","LSV"="#E377C2","BQCV"="#7F7F7F",
                                   "SBV"="#DCBD22","Amfv"="#E77C7C","Arv2"="#5FD35F",
                                   "Vov1"="#BC8B81","VTV"="#C6AEDC","DWV_B"="#17BECF",
                                   "VDV2"="#57A9E2","VDV3"="#FFB574","IAPV"="#FF7F0E",
                                   "SBPV"="#8C564B","DWV_C"="#fee08b"))+
        coord_polar(theta='y',start=0)+
        annotate(geom = 'text', x = 0.5, y = 0, label=paste("Composition of viruses",'\n',"Varroa exposure time: ",year_list[i],"years"),size=4)+
        theme_void()+
        xlim(0.5,2.5)
    
    print(a)
}

#---------------------------------------------------------------------------------------------------------------------------------------------
## donut chart, before and after Varroa infested
rm(list = ls())

data <- read.csv("new_data.csv", header = TRUE)
data <- data[ ,24:41] %>% group_by(V.arrival) %>% summarise_all(mean)

data.1 <- data.frame(t(data))
data.1 <- data.1[-1, ]
data.1$virus <- rownames(data.1)
rownames(data.1) <- c(1:17)

year_list <- c("Before","After")

for (i in 1:2){
    a<-ggplot(data.1, aes(x=2, y=data.1[,i], fill=virus))+
        geom_bar(stat='identity',color='white')+
        scale_fill_manual(values=c("DWV_A"="#1F78B4","KBV"="#2CA02C","CBPV"="#D62728",
                                   "ABPV"="#9467BD","LSV"="#E377C2","BQCV"="#7F7F7F",
                                   "SBV"="#DCBD22","Amfv"="#E77C7C","Arv2"="#5FD35F",
                                   "Vov1"="#BC8B81","VTV"="#C6AEDC","DWV_B"="#17BECF",
                                   "VDV2"="#57A9E2","VDV3"="#FFB574","IAPV"="#FF7F0E",
                                   "SBPV"="#8C564B","DWV_C"="#fee08b"))+
        coord_polar(theta='y',start=0)+
        annotate(geom = 'text', x = 0.5, y = 0, label=paste("Composition of viruses",'\n', year_list[i],"Varroa infested"),size=5)+
        theme_void()+
        xlim(0.5,2.5)
    
    print(a)
}


