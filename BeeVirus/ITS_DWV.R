library(ggplot2)
library(dplyr)
library(zoo)
library(TTR)
library(tseries)
library(forecast)
library(tframePlus)
library(imputeTS)
library(seasonal)


rm(list = ls())

data <- read.csv("00_meta.csv")

colnames(data)[1] <- "ID"
data <- data[-c(which(duplicated(data$ID))), ]

summary(data)
summary(as.factor(data$Varroa)) # NO(178), YES(204), NA(2424)
summary(as.factor(data$Country)) # BELIZE(12), MEXICO(1457), USA(1337)
summary(as.factor(data$State)) # BELIZE(12), CANPECHE(9), N.VERACRUZ(84), QUINTANA ROO(3), TAMAULIPAS(1359), 
# Texas(3), TEXAS(1334), VERACRUZ(2)
summary(as.factor(data$Geocluster)) # BEL,BOR,CAM,CAN,HOU,SAT,TAM,VER,WEL

data$Varroa <-fct_explicit_na(data$Varroa, na_level = "Unknown")
data$State <- as.factor(data$State)
levels(data$State)[3] <- "VERACRUZ"
levels(data$State)[6] <- "TEXAS"

data$Date.1 <- as.yearmon(paste(data$Year, data$Mo, sep="-"), "%Y-%m") #Sys.setlocale("LC_TIME", "English")
data$Date.2 <- as.Date(paste0(data$Year,"-12-31"), "%Y-%m-%d")

data_g <- data %>% 
  group_by(Date.1, Geocluster, Varroa) %>% 
  summarise(count=length(ID)) %>% 
  group_by(Date.1, Geocluster) %>% 
  mutate(total=sum(count), percent=round(100*count/total,2))

data_g <- data_g[which(data_g$Varroa=="YES"), ]

## first time detected Varroa: TAM-1995/Apr, BOR-1992/Apr, WEL-1994/Apr

#-----------------------------------------------------------------------------------

rm(list = ls())

data <- read.csv("new_data.csv")
data <- data[ ,c(3,4,5,10,18:34)]

data$Date.1 <- as.yearmon(paste(data$Year, data$Mo, sep="-"), "%Y-%m") #Sys.setlocale("LC_TIME", "English")
data$Date.2 <- as.Date(paste0(data$Year,"-12-31"), "%Y-%m-%d")

data.Date.1 <- data[ ,c(22,4:21)]
data.Date.1<- data.Date.1 %>% group_by(Geocluster, Date.1) %>% summarise_all(mean, na.rm=T)

data.T <- data.Date.1[which(data.Date.1$Geocluster=="TAM"), ]
data.T$Varroa <- 0
data.T$Varroa[which(data.T$Date.1>=as.yearmon("1995-4", "%Y-%m"))] <- 1

data.B <- data.Date.1[which(data.Date.1$Geocluster=="BOR"), ]
data.B$Varroa <- 0
data.B$Varroa[which(data.B$Date.1>=as.yearmon("1992-4", "%Y-%m"))] <- 1

data.W <- data.Date.1[which(data.Date.1$Geocluster=="WEL"), ]
data.W$Varroa <- 0
data.W$Varroa[which(data.W$Date.1>=as.yearmon("1994-4", "%Y-%m"))] <- 1

#---------------------------------------------------------------------------------------
## TAM
## Sample for post-intervention period (test)
test.T <- data.T[which(data.T$Varroa==1), 2:19]
## Sample for pre-intervention period (train)
train.T <- data.T[which(data.T$Varroa==0), 2:19]

train.T <- train.T %>%
  add_row(Date.1=as.yearmon("1988-11", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1989-1", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1989-2", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1989-12", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1990-2","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1991-12", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1992-4", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1992-5", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1992-6", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1992-11", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1992-12", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1993-2","%Y-%m"))

train.T <- train.T[order(train.T$Date.1), ]

test.T <- test.T %>%
  add_row(Date.1=as.yearmon("1995-6", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-7", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-8", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-9", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-10","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-11", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-12", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-1", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-2", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-3", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-4", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-6", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-7", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-8", "%Y-%m"))

test.T <- test.T[order(test.T$Date.1), ]

# identify the moel for the process over the pre-intervention period
TStrain.T <- ts(train.T$DWV,frequency=12,start=c(1988,3),end=c(1993,4))
TStest.T <- ts(test.T$DWV,frequency=12,start=c(1995,4),end=c(1996,9))

## imputation of missing data
m1.1<-na_kalman(TStrain.T)  #same as na.StrucTS
m1.1[which(m1.1<0)] <- 0
m1.2<-na_kalman(TStest.T)
m1.2[which(m1.2<0)] <- 0

plotNA.imputations(TStrain.T,m1.1)
plotNA.imputations(TStest.T,m1.2)

## monthly
adf.test(m1.1)  #whether it is stationary? H0=non-stationary, H1=stationary, p-value:0.03459
ndiffs(m1.1)  # =0
plot(decompose(m1.1))
tsdisplay(m1.1)
tsdisplay(m1.2)
### remove seasonal effect (time series has not or less than 2 periods)
#m1.2_seasonaladjusted <- seasadj(decompose(m1.2))
#m1.2_seasonaladjusted[which(m1.2_seasonaladjusted<0)] <- 0


# Arima models for pre-intervention period using auto.arima for DWV
arma1 <- auto.arima(m1.1, trace=T)  

# forecast based on the model for the pre-intervention period for DWV
arma1.f <- forecast(arma1, 41)
plot(arma1.f)

# plot
y <- as.yearmon(time(arma1.f$mean))
df1.1 <- data.frame(virus_pre=arma1.f$mean,Date=y)
df1.1[which(df1.1$virus_pre<0),1] <- 0
x <- as.yearmon(time(m1.2))
df1.2 <- data.frame(virus_obs=m1.2,Date=x)

df1 <- merge(df1.1,df1.2, by="Date",all=T)
df1 <- df1 %>% mutate(diff=virus_obs-virus_pre)
mean(df1$diff,na.rm=T)

z <- as.yearmon(time(m1.1))
df1.3 <- data.frame(Date=z,virus=m1.1)

ggplot(data=df1)+
  geom_smooth(method=lm,data=df1.3,se=F,aes(x=Date,y=virus))+
  geom_point(data=df1.3,aes(x=Date,y=virus))+
  #geom_line(size=1,aes(x=Date,y=virus_pre,color="Prediction"))+
  geom_smooth(method=lm,se=F,aes(x=Date,y=virus_pre,color='Prediction'))+
  geom_point(shape=21,fill='#2297e6',color='black',size=1.5,aes(x=Date,y=virus_pre))+
  #geom_line(size=1,aes(x=Date,y=virus_obs,color='Observation'))+
  geom_smooth(method=lm,se=F,aes(x=Date,y=virus_obs,color='Observation'))+
  geom_point(shape=21,fill='#d1600f',color="black",size=1.5,aes(x=Date,y=virus_obs))+
  geom_vline(xintercept = as.yearmon("1995-4","%Y-%m"),linetype="dashed")+
  labs(title="Intervention plot_TAM",y='viral TPH',subtitle="Deformed wing virus")+
  theme(legend.position = "bottom")+
  scale_color_manual(values=c('#d1600f','#2297e6'))

ggplot(na.omit(df1),aes(x=Date,y=diff))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept = 0, color="red", linetype="dashed")+
  labs(title="Time Series Plot of Difference",subtitle="Deformed wing virus ")

#--------------------------------------------------------------------------------------------------------
## BOR
## Sample for post-intervention period (test)
test.B <- data.B[which(data.B$Varroa==1), 2:19]
## Sample for pre-intervention period (train)
train.B <- data.B[which(data.B$Varroa==0), 2:19]

train.B <- train.B %>%
  add_row(Date.1=as.yearmon("1988-6", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1988-9", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1988-10", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1988-12", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1989-1","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1989-2", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1989-3", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1989-4", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1989-6", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1989-8", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1989-9", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1989-10","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1989-11","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1989-12","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1990-1","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1990-2","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1990-5","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1990-8","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1990-10","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1990-11","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1990-12","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1991-1","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1991-2","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1991-3","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1991-9","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1991-11","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1991-12","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1992-1","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1992-2","%Y-%m"))

train.B <- train.B[order(train.B$Date.1), ]


test.B <- test.B %>%
  add_row(Date.1=as.yearmon("1992-6", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1992-7", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1992-8", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1992-9", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1992-10","%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1992-11", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1992-12", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1993-1", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1993-2", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1993-3", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1993-11", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1993-12", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1994-1", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1994-2", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1994-5", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1994-8", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1994-9", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1994-10", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1994-11", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1994-12", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-1", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-2", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-3", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-4", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-5", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-6", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-7", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-8", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-9", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-10", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-11", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-12", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-1", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-2", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-3", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-4", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-5", "%Y-%m"))

test.B <- test.B[order(test.B$Date.1), ]

# identify the moel for the process over the pre-intervention period
TStrain.B <- ts(train.B$DWV,frequency=12,start=c(1988,5),end=c(1992,3))
TStest.B <- ts(test.B$DWV,frequency=12,start=c(1992,4),end=c(1996,6))

## imputation of missing data
m1.1<-na_kalman(TStrain.B)  #same as na.StrucTS
m1.1[which(m1.1<0)] <- 0
m1.2<-na_kalman(TStest.B)
m1.2[which(m1.2<0)] <- 0

plotNA.imputations(TStrain.B,m1.1)
plotNA.imputations(TStest.B,m1.2)

## monthly
adf.test(m1.1)  #whether it is stationary? H0=non-stationary, H1=stationary, p-value:0.3075
ndiffs(m1.1)
plot(decompose(m1.1))
tsdisplay(m1.1)
tsdisplay(m1.2)

## remove seasonal effect
m1.2_seasonaladjusted <- seasadj(mstl(m1.2))
m1.2_seasonaladjusted[which(m1.2_seasonaladjusted<0)] <- 0

# Arima models for pre-intervention period using auto.arima for DWV
arma1 <- auto.arima(m1.1, trace=T)

# forecast based on the model for the pre-intervention period for DWV
arma1.f <- forecast(arma1, 51)
plot(arma1.f)

# plot
y <- as.yearmon(time(arma1.f$mean))
df1.1 <- data.frame(virus_pre=arma1.f$mean,Date=y)
df1.1[which(df1.1$virus_pre<0),1] <- 0
x <- as.yearmon(time(m1.2_seasonaladjusted))
df1.2 <- data.frame(virus_obs=m1.2_seasonaladjusted,Date=x)

df1 <- merge(df1.1,df1.2, by="Date",all=T)
df1 <- df1 %>% mutate(diff=virus_obs-virus_pre)
mean(df1$diff,na.rm=T)

z <- as.yearmon(time(m1.1))
df1.3 <- data.frame(Date=z,virus=m1.1)

ggplot(data=df1)+
  geom_smooth(data=df1.3,se=F,aes(x=Date,y=virus))+
  geom_point(data=df1.3,aes(x=Date,y=virus))+
  #geom_line(size=1,aes(x=Date,y=virus_pre,color="Prediction"))+
  geom_smooth(se=F,aes(x=Date,y=virus_pre,color='Prediction'))+
  geom_point(shape=21,fill='#2297e6',color='black',size=1.5,aes(x=Date,y=virus_pre))+
  #geom_line(size=1,aes(x=Date,y=virus_obs,color='Observation'))+
  geom_smooth(se=F,aes(x=Date,y=virus_obs,color='Observation'))+
  geom_point(shape=21,fill='#d1600f',color="black",size=1.5,aes(x=Date,y=virus_obs))+
  geom_vline(xintercept = as.yearmon("1992-4","%Y-%m"),linetype="dashed")+
  labs(title="Intervention plot_BOR",y='viral TPH',subtitle="Deformed wing virus")+
  theme(legend.position = "bottom")+
  scale_color_manual(values=c('#d1600f','#2297e6'))

ggplot(na.omit(df1),aes(x=Date,y=diff))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept = 0, color="red", linetype="dashed")+
  labs(title="Time Series Plot of Difference",subtitle="Deformed wing virus ")

#--------------------------------------------------------------------------------------------------------
## WEL
## Sample for post-intervention period (test)
test.W <- data.W[which(data.W$Varroa==1), 2:19]
## Sample for pre-intervention period (train)
train.W<- data.W[which(data.W$Varroa==0), 2:19]

train.W <- train.W %>%
  add_row(Date.1=as.yearmon("1992-7", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1992-8", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1992-9", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1992-10", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1992-11", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1992-12", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1993-2", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1993-3", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1993-4", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1993-5", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1993-6", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1993-7", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1993-9", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1993-10", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1993-11", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1993-12", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1994-1", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1994-2", "%Y-%m"))

train.W <- train.W[order(train.W$Date.1), ]

test.W <- test.W %>%
  add_row(Date.1=as.yearmon("1994-7", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1994-9", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1994-11", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1994-12", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-1", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-2", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-4", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-5", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-6", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-7", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-8", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-9", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-10", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-11", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1995-12", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-1", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-2", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-3", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-4", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-5", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-7", "%Y-%m")) %>%
  add_row(Date.1=as.yearmon("1996-8", "%Y-%m"))
  
test.W <- test.W[order(test.W$Date.1), ]  
  
# identify the moel for the process over the pre-intervention period
TStrain.W <- ts(train.W$DWV,frequency=12,start=c(1992,3),end=c(1994,3))
TStest.W <- ts(test.W$DWV,frequency=12,start=c(1994,4),end=c(1996,9))

## imputation of missing data
m1.1<-na_kalman(TStrain.W)  #same as na.StrucTS
m1.1[which(m1.1<0)] <- 0
m1.2<-na_kalman(TStest.W)
m1.2[which(m1.2<0)] <- 0

plotNA.imputations(TStrain.W,m1.1)
plotNA.imputations(TStest.W,m1.2)

## monthly
adf.test(m1.1)  #whether it is stationary? H0=non-stationary, H1=stationary, p-value:0.3075
ndiffs(m1.1)
plot(decompose(m1.1))
tsdisplay(m1.1)
tsdisplay(m1.2)

## remove reasonal effect
m1.2_seasonaladjusted <- seasadj(mstl(m1.2))
m1.2_seasonaladjusted[which(m1.2_seasonaladjusted<0)] <- 0

# Arima models for pre-intervention period using auto.arima for DWV
arma1 <- auto.arima(m1.1, trace=T)

# forecast based on the model for the pre-intervention period for DWV
arma1.f <- forecast(arma1, 29)
plot(arma1.f)

# plot
y <- as.yearmon(time(arma1.f$mean))
df1.1 <- data.frame(virus_pre=arma1.f$mean,Date=y)
df1.1[which(df1.1$virus_pre<0),1] <- 0
x <- as.yearmon(time(m1.2))
df1.2 <- data.frame(virus_obs=m1.2_seasonaladjusted,Date=x)

df1 <- merge(df1.1,df1.2, by="Date",all=T)
df1 <- df1 %>% mutate(diff=virus_obs-virus_pre)
mean(df1$diff,na.rm=T)

z <- as.yearmon(time(m1.1))
df1.3 <- data.frame(Date=z,virus=m1.1)

ggplot(data=df1)+
  geom_smooth(method=lm,data=df1.3,se=F,aes(x=Date,y=virus))+
  geom_point(data=df1.3,aes(x=Date,y=virus))+
  #geom_line(size=1,aes(x=Date,y=virus_pre,color="Prediction"))+
  geom_smooth(method=lm,se=F,aes(x=Date,y=virus_pre,color='Prediction'))+
  geom_point(shape=21,fill='#2297e6',color='black',size=1.5,aes(x=Date,y=virus_pre))+
  #geom_line(size=1,aes(x=Date,y=virus_obs,color='Observation'))+
  geom_smooth(method=lm,se=F,aes(x=Date,y=virus_obs,color='Observation'))+
  geom_point(shape=21,fill='#d1600f',color="black",size=1.5,aes(x=Date,y=virus_obs))+
  geom_vline(xintercept = as.yearmon("1994-4","%Y-%m"),linetype="dashed")+
  labs(title="Intervention plot_WEL",y='viral TPH',subtitle="Deformed wing virus")+
  theme(legend.position = "bottom")+
  scale_color_manual(values=c('#d1600f','#2297e6'))

ggplot(na.omit(df1),aes(x=Date,y=diff))+
  geom_line()+
  geom_point()+
  geom_hline(yintercept = 0, color="red", linetype="dashed")+
  labs(title="Time Series Plot of Difference",subtitle="Deformed wing virus ")

  
  
  
  






