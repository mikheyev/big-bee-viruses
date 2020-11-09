##############################################################################
##                                                                          ##   
## Author: Sherry Lin                                                       ##
## Student ID: u6890762                                                     ##
## Project: evolution of bee viruses                                        ##
## Purpose: Varroa presence, using different data to see when Varroa had    ##
##          arrived at various sites                                        ##
## Version: 16                                                              ##
## Date: 22/09/20                                                           ##
##                                                                          ##
##############################################################################

library(ggplot2)
library(dplyr)
library(ggpubr)

rm(list = ls())

data <- read.csv("new_meta.csv")
data <- data[which(data$Geocluster=='WEL'|data$Geocluster=='BOR'|data$Geocluster=='TAM'|data$Geocluster=="VER"), ]

Sys.setlocale("LC_TIME", "English")
data$Date <- as.Date(paste(data$Year,data$Mo,data$Day,sep='-'), "%Y-%m-%d")

# Date:first Varroa detection
data$Date.V <- '1995-4-11'
data$Date.V[which(data$Geocluster=="BOR")] <- "1992-4-8"
data$Date.V[which(data$Geocluster=="WEL")] <- "1994-4-27"
data$Date.V[which(data$Geocluster=="VER")] <- NA
data$Date.V <- as.Date(data$Date.V)

data$Year.V <- 1995
data$Year.V[which(data$Geocluster=="BOR")] <- 1992
data$Year.V[which(data$Geocluster=="WEL")] <- 1994
data$Year.V[which(data$Geocluster=='VER')] <- NA
data$month <- difftime(data$Date,data$Date.V,units="days")/30
data$quarter <- difftime(data$Date,data$Date.V,units='days')/90
data$year <- difftime(data$Date,data$Date.V,units='days')/365

for (i in 1:2761) {
  if(is.na(data$year[i])&data$Geocluster[i]!='VER'){data$year[i] <- data$Year[i]-data$Year.V[i]}
}


data$year.1 <- NA
for(i in 1:2761){
  if(data$year[i]<0 & !is.na(data$year[i])){data$year.1[i] <- -ceiling(abs(data$year[i]))}
  else if(data$year[i]>=0 & !is.na(data$year[i])){data$year.1[i] <- ceiling(data$year[i])}
}

data$year.2[which(data$year.1==7)] <-'<7'
data$year.2[which(data$year.1==6)] <-'<6'
data$year.2[which(data$year.1==5)] <-'<5'
data$year.2[which(data$year.1==4)] <-'<4'
data$year.2[which(data$year.1==3)] <-'<3'
data$year.2[which(data$year.1==2)] <-'<2'
data$year.2[which(data$year.1==1|data$year.1==0)] <-'<1'
data$year.2[which(data$year.1<0|is.na(data$year.1))] <- '0'

data$year.3[data$year.1==7|data$year.1==6|data$year.1==5|data$year.1==4|data$year.1==3] <- '>2'
data$year.3[data$year.1==2] <- '<2'
data$year.3[data$year.1==1|data$year.1==0] <- '<1'
data$year.3[data$year.1<0|is.na(data$year.1)] <- '0'

#================================================================================
data_all <- data %>%
  group_by(year.1,Varroa) %>%
  summarise(count=length(ID)) %>%
  group_by(year.1) %>%
  mutate(total=sum(count),percent=round(100*count/total,2))

data_g <- data %>% 
  group_by(year.1, Geocluster, Varroa) %>% 
  summarise(count=length(ID)) %>% 
  group_by(year.1, Geocluster) %>% 
  mutate(total=sum(count), percent=round(100*count/total,2))

data_g$Geocluster <- ordered(data_g$Geocluster, level=c("WEL","BOR","TAM","VER"))

#------------------------------------------------------------------------------------
# Varroa_all: NO,YES,NA

a<-ggplot(data_g[-which(is.na(data_g$year.1)),], aes(x=as.factor(year.1), y=percent, fill=Varroa))+
  geom_bar(stat="identity",position=position_dodge(),color='black')+
  facet_grid(rows=vars(Geocluster))+
  scale_fill_manual(values = c('lightgrey','black'),na.value='white')+
  labs(x="",y="",fill='Varroa records')+
  theme_classic()+
  theme(legend.justification=c(0.99,0.02),legend.position=c(0.99,0.02),
        legend.background = element_rect(linetype="solid",colour ="darkgrey"),
        plot.margin=unit(c(0,0.5,-0.3,0),"lines"),
        axis.text.y=element_text(size=8))


b<-ggplot(data_all[-which(is.na(data_all$year.1)),], aes(x=as.factor(year.1), y=percent,fill=Varroa))+
  geom_bar(stat="identity",position=position_dodge(),color='black')+
  labs(x="",y="")+
  scale_fill_manual(values = c('lightgrey','black'),na.value='white')+
  theme_classic()+
  theme(legend.position='none',plot.margin=unit(c(1,1.75,-0.3,-0.1),"lines"),
        axis.text.y=element_text(size=8))

c <- ggarrange(b,a,nrow=2,ncol=1,heights=c(0.3, 0.7))
c <- annotate_figure(c,
                     top=text_grob(paste("Prevalence of Varroa destructor"),size=12),
                     bottom=text_grob("Varroa exposure time (years)"),
                     left = text_grob("Percentage (%)",rot = 90))
print(c)

#===========================================================================================
## exposure time groups: 0,<1,<2,<3,<4,<5,<6,<7
data_all <- data %>%
  group_by(year.2,Varroa) %>%
  summarise(count=length(ID)) %>%
  group_by(year.2) %>%
  mutate(total=sum(count),percent=round(100*count/total,2))

data_all <- ungroup(data_all)
data_all  <- data_all %>% add_row(year.2='0',Varroa='YES',count=0,total=0,percent=0) %>%
  add_row(year.2='<6',Varroa='YES',count=0,total=0,percent=0) %>%
  add_row(year.2='<7',Varroa='YES',count=0,total=0,percent=0)

data_all$year.2 <- ordered(data_all$year.2, levels=c('0','<1','<2','<3','<4','<5','<6','<7'))


data_g <- data %>% 
  group_by(year.2, Geocluster, Varroa) %>% 
  summarise(count=length(ID)) %>% 
  group_by(year.2, Geocluster) %>% 
  mutate(total=sum(count), percent=round(100*count/total,2))

data_g <- ungroup(data_g)
data_g <- data_g %>% add_row(year.2='<2',Geocluster='WEL',Varroa='YES',count=0,total=0,percent=0) %>%
  add_row(year.2='<3',Geocluster='TAM',Varroa='YES',count=0,total=0,percent=0) %>%
  add_row(year.2='<4',Geocluster='TAM',Varroa='YES',count=0,total=0,percent=0) %>%
  add_row(year.2='<4',Geocluster='WEL',Varroa='YES',count=0,total=0,percent=0) %>%
  add_row(year.2='<5',Geocluster='WEL',Varroa='YES',count=0,total=0,percent=0) %>%
  add_row(year.2='<5',Geocluster='TAM',Varroa='YES',count=0,total=0,percent=0) %>%
  add_row(year.2='<6',Geocluster='WEL',Varroa='YES',count=0,total=0,percent=0) %>%
  add_row(year.2='<6',Geocluster='BOR',Varroa='YES',count=0,total=0,percent=0) %>%
  add_row(year.2='<6',Geocluster='TAM',Varroa='YES',count=0,total=0,percent=0) %>%
  add_row(year.2='<7',Geocluster='WEL',Varroa='YES',count=0,total=0,percent=0) %>%
  add_row(year.2='<7',Geocluster='BOR',Varroa='YES',count=0,total=0,percent=0) %>%
  add_row(year.2='<7',Geocluster='TAM',Varroa='YES',count=0,total=0,percent=0) %>%
  add_row(year.2='0', Geocluster='WEL',Varroa='YES',count=0,total=0,percent=0) %>%
  add_row(year.2='0', Geocluster='BOR',Varroa='YES',count=0,total=0,percent=0) %>%
  add_row(year.2='0', Geocluster='TAM',Varroa='YES',count=0,total=0,percent=0)
  

data_g$year.2 <- ordered(data_g$year.2, levels=c('0','<1','<2','<3','<4','<5','<6','<7'))
data_g$Geocluster <- ordered(data_g$Geocluster, level=c("WEL","BOR","TAM","VER"))

# Varroa: YES

a<-ggplot(data_g[which(data_g$Varroa=="YES"), ], aes(x=year.2, y=percent))+
  geom_bar(stat="identity",fill='black',color='black')+
  ylim(c(0,100))+
  facet_grid(rows=vars(Geocluster))+
  labs(x="",y="")+
  theme_classic()+
  theme(plot.margin=unit(c(0,0.5,-0.3,0),"lines"),
        axis.text.y=element_text(size=8))

b<-ggplot(data_all[which(data_all$Varroa=="YES"), ], aes(x=year.2, y=percent))+
  geom_bar(stat="identity",fill='black',color='black')+
  ylim(c(0,100))+
  labs(x="",y="")+
  theme_classic()+
  theme(plot.margin=unit(c(1,1.75,-0.3,-0.1),"lines"),
        axis.text.y=element_text(size=8))

c <- ggarrange(b,a,nrow=2,ncol=1,heights=c(0.3, 0.7))
c <- annotate_figure(c,
                     top=text_grob(paste("Prevalence of Varroa destructor"),size=12),
                     bottom=text_grob("Varroa exposure time (years)"),
                     left = text_grob("Percentage (%)",rot = 90))
print(c)



ggplot(data_all[which(data_all$Varroa=="YES"&data_all$year.2!='<6'&data_all$year.2!='<7'),], aes(x=year.2, y=percent))+
  geom_bar(stat="identity",fill='grey',color='black')+
  ylim(c(0,100))+
  geom_text(aes(label=percent),position = position_dodge(0.8),
            vjust = -0.3, size = 3.5)+
  labs(title="Prevalence of Varroa desturctor",x="Varroa exposure time (years)",y="Percentage (%)")+
  theme_classic()+
  theme(legend.position='none',
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=10))+
  scale_x_discrete(labels=c('0','<1','<2','<3','<4','<5','<6','<7'))

#-------------------------------------------------------------------------------------------------------
data_all <- data %>%
  group_by(year.3,Varroa) %>%
  summarise(count=length(ID)) %>%
  group_by(year.3) %>%
  mutate(total=sum(count),percent=round(100*count/total,2))

data_all <- ungroup(data_all)
data_all  <- data_all %>% add_row(year.3='0',Varroa='YES',count=0,total=0,percent=0)
data_all$year.3 <- factor(data_all$year.3,levels=c('0','<1','<2','>2'))

ggplot(data_all[which(data_all$Varroa=="YES"),], aes(x=year.3, y=percent))+
  geom_bar(stat="identity",fill='grey30',color='black')+
  ylim(c(0,100))+
  geom_text(aes(label=percent),position = position_dodge(0.8),
            vjust = -0.3, size = 3.5)+
  labs(title="Prevalence of Varroa desturctor",x="Varroa exposure time (years)",y="Percentage (%)")+
  theme_classic()+
  theme(legend.position='none',
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=10))+
  scale_x_discrete(labels=c('0','<1','<2','>2'))



result <- glm(as.factor(Varroa)~year.2, data=data_all, family=binomial(link='logit'))

data$year.3 <- factor(data$year.3, levels=c('0','<1','<2','>2'))
dt <- table(data$Varroa,data$year.3); dt
a <- chisq.test(dt,simulate.p.value = T)


 