#######################################################
##                                                   ##   
## Author: Sherry Lin                                ##
## Student ID: u6890762                              ##
## Project: evolution of bee viruses                 ##
## Purpose: plotting the viral population over time  ##
## Version: 14                                       ##
## Date: 01/07/20                                    ##
##                                                   ##
#######################################################

library(ggcorrplot) #ggcorrplot(): make correalation plot
library(ggplot2) #ggplot()+geom_boxplot(), geom_bar(), geom_line, geom_point()
library(ggsci) #scale_color_jco() / scale_fill_jco(): color in Journal of Clinical Oncology
library(dplyr) #group_by(): groups a dataframe based on certain fields / summarise()
library(ggpubr) #ggarrange(): make multiple ggplots into one / annotate_figure()


rm(list = ls())

# import data
# let 'Geocluster' become orderd factor(from south to north)
data.1 <- read.csv("new_data.csv")
data.2 <- read.csv("new_data_mean_YearGeo.csv")
data.3 <- read.csv("new_data_mean_Year.csv")
data.1$Geocluster <- ordered(data.1$Geocluster, level=c("WEL","BOR","TAM","VER"))
data.2$Geocluster <- ordered(data.2$Geocluster, level=c("WEL","BOR","TAM","VER"))


# plotting viral population over time
## by TPM
name_list <- colnames(data.1[18:33])
name_full <- c('Acute bee paralysis virus','Apis mellifera filamentous virus',
               'Apis rhabdovirus 2','Black queen cell virus', 'Chronic bee paralysis virus',
               'Deformed wing virus','Israeli acute paralysis virus','Kashmir bee virus',
               'Lake Sinai virus','Slow bee paralysis virus','Sacbrood virus',
               'Varroa destructor virus 1','Varroa destructor virus 2','Varroa destructor virus 3',
               'Varroa orthomyxovirus-1','Varroa Tymo-like virus')

# BOXPLOT: y=transcript per hundred, x=year
for(i in 1:16){
  p1.1 <- ggplot(data.1, aes_string(x='Year',y=name_list[i],group='Year',color='Geocluster'))+
    geom_boxplot()+geom_point()+
    scale_y_sqrt()+
    facet_grid(rows=vars(Geocluster))+
    theme(legend.position="none",plot.margin=unit(c(-1.8,1,-0.5,-0.5),"lines"))+
    labs(title = "",y="",x="")+scale_color_jco()
  
  p1.2<-ggplot(data.1, aes_string(x='Year',y=name_list[i],group='Year'))+
    geom_boxplot()+geom_point()+
    scale_y_sqrt()+
    theme(plot.margin=unit(c(1,2.1,0,-0.5),"lines"))+
    labs(subtitle = "Overview",y="",x="")+scale_color_jco()
  p1 <- ggarrange(p1.2,p1.1,nrow=2,ncol=1,heights=c(0.3, 0.7))
  p1 <- annotate_figure(p1,top=text_grob(paste(name_full[i])),
                        left = text_grob("Transcript per hundred",rot = 90),
                        bottom=text_grob("Year"))
  print(p1)
  
}

# BAR CHART: y=average of TPH, x=year
for(i in 1:16){
  p1.1 <- ggplot(data.3, aes_string(x='Year',y=name_list[i],group='Year'))+ 
    geom_bar(stat = "identity")+
    theme(legend.position="none",plot.margin=unit(c(0,2.1,0,-0.5),"lines"))+ 
    scale_y_sqrt()+
    labs(subtitle = "Overview",title="",y="",x="")+scale_color_jco()
  
  p1.2 <- ggplot(data.2, aes_string(x='Year', y=name_list[i], fill='Geocluster'))+ 
    geom_bar(stat = "identity")+
    scale_y_sqrt()+
    facet_grid(rows=vars(Geocluster))+
    theme(legend.position = 'none',plot.margin=unit(c(-1.8,1,-0.5,-0.5), "lines"))+
    labs(title="",y="",x="")+scale_fill_jco()
  
  p <-ggarrange(p1.1,p1.2,nrow=2,ncol=1,heights=c(0.3, 0.7))
  p <-annotate_figure(p,top=text_grob(name_full[i]),
                      left = text_grob("Average of TPH",rot = 90),
                      bottom=text_grob("Year"))
  print(p)
}

#########################################################################################
# plotting viral population over time
## by prevalence

data.4 <- read.csv("new_data_prevalence_YearGeo.csv")
data.5 <- read.csv("new_data_prevalence_Year.csv")

for (i in 1:16){
  a<-ggplot(data.4, aes_string(x='Year', y=name_list[i], fill='Geocluster'))+
    geom_bar(stat="identity")+
    scale_y_sqrt()+
    facet_grid(rows=vars(Geocluster))+
    scale_fill_jco()+
    labs(x="",y="")+
    theme(legend.position = "none",plot.margin=unit(c(0,0.5,-0.3,0),"lines"),
          axis.text.y=element_text(size=8))
  
  b<-ggplot(data.5, aes_string(x='Year', y=name_list[i]))+
    geom_bar(stat="identity")+
    scale_y_sqrt()+
    labs(subtitle="Overview",x="",y="")+
    theme(plot.margin=unit(c(1,1.6,-0.3,-0.1),"lines"),
          axis.text.y=element_text(size=8))
  
  c <- ggarrange(b,a,nrow=2,ncol=1,heights=c(0.3, 0.7))
  c <- annotate_figure(c,top=text_grob(paste("Prevalence of ",name_full[i]),size=12),
                       bottom=text_grob("Year"),
                       left = text_grob("Percentage (%)",rot = 90),)
  print(c)
  
}








