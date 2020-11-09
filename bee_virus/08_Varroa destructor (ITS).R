library(data.table) #fread(): read the data
library(tools) #file_path_sans_ext(): without extension
library(dplyr)

rm(list = ls())

# getting files' name without extension (the files' names are the sample's ID) 
setwd("D:/Sherry/文件/ResearchProject/BeeVirus")
unzip.filename <- list.files("00_bigbee")                                     
unzip.filename.new <- file_path_sans_ext(list.files("00_bigbee")) 

data <- read.csv('new_data.csv')

target <- read.csv('00_taxa2transcript.tsv',sep='',header=F)
target_Varroa <- target[target$V1=='Varroa_destructor', ]

setwd("00_bigbee")
c <- data.frame()
for (i in 1:1043){                                               
  d <- fread(unzip.filename[i])
  d <- filter(d,target_id %in% target_Varroa$V2)
  d$ID <- unzip.filename.new[i]
  c <- rbind(c,d)
}   

negative <- c[which(grepl("Negative",c$ID)&c$tpm!=0), ]
neg_id <- negative[-which(duplicated(negative$target_id)),]

tpm <- as.data.frame(dcast(c,ID~target_id,value.var='tpm'))
detected <- which(apply(tpm[ ,2:35642],2,sum,na.rm=T)!=0)+1
tpm <- tpm[,c(1,detected)]

#---------------------------------------------------------------------
setwd("D:/Sherry/文件/ResearchProject/BeeVirus")
data_merge <- merge(data, tpm, by='ID')


# ITS in tobit model
# setup data: add variables -> time,level,trend
data_ITS <- data_merge[-which(data_merge$Geocluster=='VER'),]

##time
a <- c(min(data_ITS$year.1):max(data_ITS$year.1))
data_ITS$time <- NA
for (i in 1:length(data_ITS$year.1)) {
  data_ITS$time[i] <- which(a==data_ITS$year.1[i])
}

## trend
b <- c(0:max(data_ITS$year.1))
data_ITS$trend <- 0
for (i in which(data_ITS$year.1>=0)) {
  data_ITS$trend[i] <- which(b==data_ITS$year.1[i])-1
}

## level
data_ITS$level <- 0
data_ITS$level[which(data_ITS$year.1>=0)] <- 1

library(VGAM)
library(ggplot2)

data_mean <- data_ITS[,c(11,44:3761)] %>% group_by(year.1) %>% summarise_all(mean)
data_median <- data_ITS[,c(11,44:3761)] %>% group_by(year.1) %>% summarise_all(median)
for (i in 44:3761){
  tobit.y <- try(vglm(data=data_ITS,data_ITS[,i]~time+level+trend,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored')),TRUE)     
  if(isTRUE(class(tobit.y)=='try-error')) {next}
  data_ITS$fitted_tobit_y <- fitted(tobit.y)
  if (sum(data_ITS$fitted_tobit_y)==0) {next}
  c <- ggplot()+
    geom_point(data=data_mean,aes_string(x='year.1',y=names(data_ITS)[i]),shape=16,size=2)+
    geom_point(data=data_median,aes_string(x='year.1',y=names(data_ITS)[i]),shape=2,size=3)+
    geom_vline(xintercept = 0, linetype='dashed',color='darkgrey')+
    geom_line(data=data_ITS[which(data_ITS$level==0),],aes(x=year.1,y=fitted_tobit_y),size=1,color='blue')+
    geom_line(data=data_ITS[which(data_ITS$level==1),],aes(x=year.1,y=fitted_tobit_y),size=1,color='#f8766d')+
    labs(title=names(data_ITS)[i],
         x='Varroa exposure time (years)',y='Transcript per Million')+
    #ylim(c(0,10^6))+
    theme_classic()+
    scale_x_continuous(breaks = c(seq(from = -8, to = 7, by = 1)))
  ggsave(paste(names(data_ITS)[i],'.png'))
  
}


tobit.y <- vglm(data=data_ITS,~time+level+trend,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(tobit.y)
data_ITS$fitted_tobit_y <- fitted(tobit.y)

