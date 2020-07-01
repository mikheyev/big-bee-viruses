##########################################################
##                                                      ##   
## Author: Sherry Lin                                   ##
## Student ID: u6890762                                 ##
## Project: evolution of bee viruses                    ##
## Purpose: map+scatterpie->viral population over time  ##
## Version: 12                                          ##
## Date: 01/07/20                                       ##
##                                                      ##
##########################################################

library(dplyr)  #group_by(), summarise()
library(ggplot2)  #theme_set()
library(rnaturalearth)  #ne_countries():returns world country polygons at a specified scale
library(rgeos)
library(scatterpie)  #geom_scatterpie()
library(ggrepel)  #geom_label_repel()


rm(list = ls())

data <- read.csv("new_data_mean_YearGeo.csv")
data <- data[ ,1:18]
data$Geocluster <- ordered(data$Geocluster, level=c("WEL","BOR","TAM","VER"))
data$Longitude <- 0
data$Latitude <- 0
data$Longitude[which(data$Geocluster=="VER")] <- -97.30746  
data$Latitude[which(data$Geocluster=="VER")] <- 20.97450
data$Longitude[which(data$Geocluster=="TAM")] <- -98.20793
data$Latitude[which(data$Geocluster=="TAM")] <- 23.7158
data$Longitude[which(data$Geocluster=="BOR")] <- -98.12306
data$Latitude[which(data$Geocluster=="BOR")] <- 26.18369
data$Longitude[which(data$Geocluster=="WEL")] <- -97.44197
data$Latitude[which(data$Geocluster=="WEL")] <- 28.12102 

data$radius <- 1



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

year_list <- c("1988","1989","1990","1991","1992","1993","1994","1995","1996","1997","2001")

for(i in 1:11){
    map_year <- map+
        geom_scatterpie(data=data[which(data$Year==year_list[i]), ],aes(x=Longitude, y=Latitude, r=radius),
                        cols=colnames(data)[3:18],legend_name = "Viruses")+
        scale_fill_manual(values=c("DWV"="#1F78B4","KBV"="#2CA02C","CBPV"="#D62728","ABPV"="#9467BD","LSV"="#E377C2","BQCV"="#7F7F7F","SV"="#DCBD22","AM"="#E77C7C",
                                   "AR"="#5FD35F","VO1"="#BC8B81","VTV"="#C6AEDC","VDV1"="#17BECF","VDV2"="#57A9E2","VDV3"="#FFB574","IAPV"="#FF7F0E","SBPV"="#8C564B"))+
        geom_label_repel(aes(x=Longitude, y=Latitude, label=Geocluster), data=data[which(data$Year==year_list[i]), ],
                         point.padding=unit(2, "lines"))+
        labs(title=paste("Year: ",year_list[i]))
    print(map_year)
    
}

############################################################################################################################################################

rm(list = ls())

data <- read.csv("new_data_mean_YearGeo.csv")
virus <- colnames(data)[3:18]

data <- t(data)
data <- data[-c(1,2,19), ]
data<-apply(data,2,as.numeric)
data <- as.data.frame(data)
data$virus <- virus

geo_list <- c('1988 BOR','1988 TAM','1989 BOR','1989 TAM','1990 BOR','1990 TAM','1990 VER',
              '1991 BOR','1991 TAM','1992 BOR','1992 TAM','1992 WEL','1993 BOR','1993 TAM',
              '1993 WEL','1994 BOR','1994 WEL','1995 TAM','1995 WEL','1996 BOR','1996 TAM',
              '1996 WEL','1997 WEL','2001 WEL')

for (i in 1:24){
    a<-ggplot(data, aes_string(x=2, y=colnames(data)[i], fill='virus'))+
    geom_bar(stat='identity',color='white')+
    scale_fill_manual(values=c("DWV"="#1F78B4","KBV"="#2CA02C","CBPV"="#D62728","ABPV"="#9467BD","LSV"="#E377C2","BQCV"="#7F7F7F","SV"="#DCBD22","AM"="#E77C7C",
                               "AR"="#5FD35F","VO1"="#BC8B81","VTV"="#C6AEDC","VDV1"="#17BECF","VDV2"="#57A9E2","VDV3"="#FFB574","IAPV"="#FF7F0E","SBPV"="#8C564B"))+
    coord_polar(theta='y',start=0)+
    annotate(geom = 'text', x = 0.5, y = 0, label=paste(geo_list[i],", Composition of viruses"),size=4)+
    theme_void()+
    xlim(0.5,2.5)
    
    print(a)
}

#############################################################################################

rm(list = ls())

data <- read.csv("new_data_mean_Year.csv")

data <- t(data)
data <- as.data.frame(data)
data <- data[-c(1,18), ]
data$virus <- rownames(data)
rownames(data) <- c(1:16)

year_list <- c("1988","1989","1990","1991","1992","1993","1994","1995","1996","1997","2001")

for (i in 1:11){
    b<-ggplot(data, aes_string(x=2, y=colnames(data)[i], fill='virus'))+
    geom_bar(stat='identity',color='white')+
    scale_fill_manual(values=c("DWV"="#1F78B4","KBV"="#2CA02C","CBPV"="#D62728","ABPV"="#9467BD","LSV"="#E377C2","BQCV"="#7F7F7F","SV"="#DCBD22","AM"="#E77C7C",
                               "AR"="#5FD35F","VO1"="#BC8B81","VTV"="#C6AEDC","VDV1"="#17BECF","VDV2"="#57A9E2","VDV3"="#FFB574","IAPV"="#FF7F0E","SBPV"="#8C564B"))+
    coord_polar('y',start=0)+
    annotate(geom = 'text', x = 0.5, y = 0, label=paste(year_list[i],", Composition of viruses"),size=4)+
    theme_void()+
    xlim(0.5, 2.5)
    
    print(b)
}



