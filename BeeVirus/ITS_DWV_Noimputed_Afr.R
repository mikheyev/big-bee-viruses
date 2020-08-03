library(dplyr)
library(car) #dwt()
library(forecast) #tsdisplay()
library(ggplot2)
library(nlme)
library(lmtest)
library(imputeTS) #na_kalman()
library(reshape2) #melt()
library(ggsci)

rm(list = ls())

data <- read.csv("new_data.csv")
data <- data[ ,c(3,4,5,10,18:34)]
data <- data[-which(data$Geocluster=="VER"), ]


Sys.setlocale("LC_TIME", "English")
data$Date <- as.Date(paste(data$Year,data$Mo,data$Day,sep='-'), "%Y-%m-%d")
data$Date.V <- '1995-4-11'
data$Date.V[which(data$Geocluster=="BOR")] <- "1992-4-8"
data$Date.V[which(data$Geocluster=="WEL")] <- "1994-4-27"
data$Date.V <- as.Date(data$Date.V)

data$Year.V <- 1995
data$Year.V[which(data$Geocluster=="BOR")] <- 1992
data$Year.V[which(data$Geocluster=="WEL")] <- 1994

data$month <- round(difftime(data$Date,data$Date.V,units="days")/30)
data$quarter <- round(difftime(data$Date,data$Date.V,units='days')/90)
data$year.1 <- round(difftime(data$Date,data$Date.V,units='days')/365)
data$year.2 <- data$Year-data$Year.V

data.m <- data[ ,c(25,5:21)] %>% group_by(month) %>% summarise_all(mean,na.rm=T)
data.m <- data.m[-84, ]
data.q <- data[ ,c(26,5:21)] %>% group_by(quarter) %>% summarise_all(mean,na.rm=T)
data.q <- data.q[-40, ]
data.y1 <- data[ ,c(27,5:21)] %>% group_by(year.1) %>% summarise_all(mean,na.rm=T)
data.y1 <- data.y1[-12, ]
data.y2 <- data[ ,c(28,5:21)] %>% group_by(year.2) %>% summarise_all(mean,na.rm=T)
#=======================================================================================

# setup data: add variables -> time,level,trend

# monthly
a <- c(min(data.m$month):max(data.m$month))
data.m$time <- 0
for (i in 1:length(data.m$month)) {
  data.m$time[i] <- which(a==data.m$month[i])
}
data.m$level <- 0
data.m$level[which(data.m$month>=0)] <- 1

b <- c(0:max(data.m$month))
data.m$trend <- 0
for (i in which(data.m$month>=0)) {
  data.m$trend[i] <- which(b==data.m$month[i])
}

# quarterly
a <- c(min(data.q$quarter):max(data.q$quarter))
data.q$time <- 0
for (i in 1:length(data.q$quarter)) {
  data.q$time[i] <- which(a==data.q$quarter[i])
}
data.q$level <- 0
data.q$level[which(data.q$quarter>=0)] <- 1

b <- c(0:max(data.q$quarter))
data.q$trend <- 0
for (i in which(data.q$quarter>=0)) {
  data.q$trend[i] <- which(b==data.q$quarter[i])
}

# yearly.1
a <- c(min(data.y1$year.1):max(data.y1$year.1))
data.y1$time <- 0
for (i in 1:length(data.y1$year.1)) {
  data.y1$time[i] <- which(a==data.y1$year.1[i])
}
data.y1$level <- 0
data.y1$level[which(data.y1$year.1>=0)] <- 1

b <- c(0:max(data.y1$year.1))
data.y1$trend <- 0
for (i in which(data.y1$year.1>=0)) {
  data.y1$trend[i] <- which(b==data.y1$year.1[i])
}

# yearly.2
a <- c(min(data.y2$year.2):max(data.y2$year.2))
data.y2$time <- 0
for (i in 1:length(data.y2$year.2)) {
  data.y2$time[i] <- which(a==data.y2$year.2[i])
}
data.y2$level <- 0
data.y2$level[which(data.y2$year.2>=0)] <- 1

b <- c(0:max(data.y2$year.2))
data.y2$trend <- 0
for (i in which(data.y2$year.2>=0)) {
  data.y2$trend[i] <- which(b==data.y2$year.2[i])
}

#===================================================================================================
# monthly,DWV
DWV.m <- data.m[ ,c(1,7,19:21,18)]

############
# modelling
############

# Fit the GLS regression model with p=1, q=1 as the previous analysis 
model_p1q1_m <- gls(data=DWV.m,DWV~time+level+trend+Afr.Extent,
                    correlation=corARMA(p=1,q=1,form=~time),
                    method="ML")
summary(model_p1q1_m) #AIC=660.77
DWV.m$fitted_gls_p1q1 <- c(fitted(model_p1q1_m))

###############
# Plot results
###############
offset <- mean(DWV.m$Afr.Extent)*model_p1q1_m$coef[5]

ggplot()+
  geom_point(data=DWV.m,aes(x=month,y=DWV,shape='observed'))+
  geom_vline(xintercept = 0, linetype='dashed',color='darkgrey')+
  # original fiited value
  #geom_line(data=DWV.m[which(DWV.m$level==0), ],aes(x=month, y=fitted_gls_p1q1))+
  #geom_line(data=DWV.m[which(DWV.m$level==1), ],aes(x=month, y=fitted_gls_p1q1))+
  
  # adjusted fitted value
  geom_segment(x=-86,   y   =model_p1q1_m$coef[1]+model_p1q1_m$coef[2]+offset,
               xend=-1, yend=model_p1q1_m$coef[1]+model_p1q1_m$coef[2]*86+offset,
               color='#7cae00',size=1,aes(linetype='fitted'))+
  geom_segment(x=0,     y   =model_p1q1_m$coef[1]+model_p1q1_m$coef[2]*87 +model_p1q1_m$coef[3]+model_p1q1_m$coef[4]+offset,
               xend=51, yend=model_p1q1_m$coef[1]+model_p1q1_m$coef[2]*138+model_p1q1_m$coef[3]+model_p1q1_m$coef[4]*52+offset,
               color='#7cae00',size=1,aes(linetype='fitted'))+
  geom_segment(x=-1,   y   =model_p1q1_m$coef[1]+model_p1q1_m$coef[2]*86+offset,
               xend=51,yend=model_p1q1_m$coef[1]+model_p1q1_m$coef[2]*138+offset,
               color='#7cae00',size=1,aes(linetype='counterfactul'))+
  scale_linetype_manual(values=c('dotted','solid'))+
  labs(title="Monthly ITS analysis: GLS_ARMA(1,1) model with Africanized effect",subtitle='Deformed wing virus',
       x='Months after Varroa arrived',y='Viral Transcript per Hundred')+
  guides(shape = guide_legend(order=1),linetype=guide_legend(order=2))+
  scale_x_continuous(breaks = c(seq(from = -86, to = 51, by = 4)))

#------------------------------------------------------------------------------------
# quarterly,DWV
DWV.q <- data.q[ ,c(1,7,19:21,18)]

############
# modelling
############

# Fit the GLS regression model with p=2, q=2 as the previous analysis 
model_p2q2_q <- gls(data=DWV.q,DWV~time+level+trend+Afr.Extent,
                    correlation=corARMA(p=2,q=2,form=~time),
                    method="ML")
summary(model_p2q2_q)#AIC=299.04
DWV.q$fitted_gls_p2q2 <- c(fitted(model_p2q2_q))


###############
# Plot results
###############
offset <- mean(DWV.q$Afr.Extent)*model_p2q2_q$coef[5]

ggplot()+
  geom_point(data=DWV.q,aes(x=quarter,y=DWV,shape='observed'))+
  geom_vline(xintercept = 0, linetype='dashed',color='darkgrey')+
  # original fiited value
  #geom_line(data=DWV.q[which(DWV.q$level==0), ],aes(x=quarter, y=fitted_gls_p2q2))+
  #geom_line(data=DWV.q[which(DWV.q$level==1), ],aes(x=quarter, y=fitted_gls_p2q2))+
  
  # adjusted fitted value
  geom_segment(x=-29,   y   =model_p2q2_q$coef[1]+model_p2q2_q$coef[2]+offset,
               xend=-1, yend=model_p2q2_q$coef[1]+model_p2q2_q$coef[2]*29+offset,
               color='#7cae00',size=1,aes(linetype='fitted'))+
  geom_segment(x=0,     y   =model_p2q2_q$coef[1]+model_p2q2_q$coef[2]*30+model_p2q2_q$coef[3]+model_p2q2_q$coef[4]+offset,
               xend=17, yend=model_p2q2_q$coef[1]+model_p2q2_q$coef[2]*47+model_p2q2_q$coef[3]+model_p2q2_q$coef[4]*18+offset,
               color='#7cae00',size=1,aes(linetype='fitted'))+
  geom_segment(x=-1,   y   =model_p2q2_q$coef[1]+model_p2q2_q$coef[2]*29+offset,
               xend=17,yend=model_p2q2_q$coef[1]+model_p2q2_q$coef[2]*47+offset,
               color='#7cae00',size=1,aes(linetype='counterfactul'))+
  scale_linetype_manual(values=c('dotted','solid'))+
  labs(title="Quarterly ITS analysis: GLS_ARMA(2,2) model with Africanized effect",subtitle='Deformed wing virus',
       x='Quarters after Varroa arrived',y='Viral Transcript per Hundred')+
  guides(shape = guide_legend(order=1),linetype=guide_legend(order=2))+
  scale_x_continuous(breaks = c(seq(from = -29, to = 17, by = 1)))

#------------------------------------------------------------------------------------
# yearly.1,DWV
DWV.y1 <- data.y1[ ,c(1,7,19:21,18)]

############
# modelling
############

# Fit the GLS regression model with p=2, as the previous analysis 
model_p2_y1 <- gls(data=DWV.y1,DWV~time+level+trend+Afr.Extent,
                    correlation=corARMA(p=2,form=~time),
                    method="ML")
summary(model_p2_y1)#AIC=69.47
DWV.y1$fitted_gls_p2 <- c(fitted(model_p2_y1))

###############
# Plot results
###############
offset <- mean(DWV.y1$Afr.Extent)*model_p2_y1$coef[5]

ggplot()+
  geom_point(data=DWV.y1,aes(x=year.1,y=DWV,shape='observed'))+
  geom_vline(xintercept = 0, linetype='dashed',color='darkgrey')+
  # original fiited value
  #geom_line(data=DWV.y1[which(DWV.y1$level==0), ],aes(x=year.1, y=fitted_gls_p2))+
  #geom_line(data=DWV.y1[which(DWV.y1$level==1), ],aes(x=year.1, y=fitted_gls_p2))+
  
  # adjusted fitted value
  geom_segment(x=-7,   y   =model_p2_y1$coef[1]+model_p2_y1$coef[2]+offset,
               xend=-1,yend=model_p2_y1$coef[1]+model_p2_y1$coef[2]*7+offset,
               color='#7cae00',size=1,aes(linetype='fitted'))+
  geom_segment(x=0,    y   =model_p2_y1$coef[1]+model_p2_y1$coef[2]*8 +model_p2_y1$coef[3]+model_p2_y1$coef[4]+offset,
               xend=4, yend=model_p2_y1$coef[1]+model_p2_y1$coef[2]*12+model_p2_y1$coef[3]+model_p2_y1$coef[4]*5+offset,
               color='#7cae00',size=1,aes(linetype='fitted'))+
  geom_segment(x=-1,   y   =model_p2_y1$coef[1]+model_p2_y1$coef[2]*7+offset,
               xend=4, yend=model_p2_y1$coef[1]+model_p2_y1$coef[2]*12+offset,
               color='#7cae00',size=1,aes(linetype='counterfactul'))+
  scale_linetype_manual(values=c('dotted','solid'))+
  labs(title="Yearly ITS analysis: GLS_AR(2) model with Africanized effect",subtitle='Deformed wing virus',
       x='Years after Varroa arrived',y='Viral Transcript per Hundred')+
  guides(shape = guide_legend(order=1),linetype=guide_legend(order=2))+
  scale_x_continuous(breaks = c(seq(from = -7, to = 4, by = 1)))

#------------------------------------------------------------------------------------
# yearly.1,DWV
DWV.y2 <- data.y2[ ,c(1,7,19:21,18)]

############
# modelling
############

# Fit the OLS regression model
model_ols_y2 <- lm(data=DWV.y2,DWV~time+level+trend+Afr.Extent)
summary(model_ols_y2)
AIC(model_ols_y2)  #122.20
DWV.y2$fitted_ols <- model_ols_y2$fitted.values

###############
# Plot results
###############
offset <- mean(DWV.y2$Afr.Extent)*model_ols_y2$coef[5]

ggplot()+
  geom_point(data=DWV.y2,aes(x=year.2,y=DWV,shape='observed'))+
  geom_vline(xintercept = 0, linetype='dashed',color='darkgrey')+
  # original fiited value
  #geom_line(data=DWV.y2[which(DWV.y2$level==0), ],aes(x=year.2, y=fitted_ols))+
  #geom_line(data=DWV.y2[which(DWV.y2$level==1), ],aes(x=year.2, y=fitted_ols))+
  
  # adjusted fitted value
  geom_segment(x=-7,   y   =model_ols_y2$coef[1]+model_ols_y2$coef[2]+offset,
               xend=-1,yend=model_ols_y2$coef[1]+model_ols_y2$coef[2]*7+offset,
               color='#7cae00',size=1,aes(linetype='fitted'))+
  geom_segment(x=0,    y   =model_ols_y2$coef[1]+model_ols_y2$coef[2]*8 +model_ols_y2$coef[3]+model_ols_y2$coef[4]+offset,
               xend=7, yend=model_ols_y2$coef[1]+model_ols_y2$coef[2]*15+model_ols_y2$coef[3]+model_ols_y2$coef[4]*8+offset,
               color='#7cae00',size=1,aes(linetype='fitted'))+
  geom_segment(x=-1,   y   =model_ols_y2$coef[1]+model_ols_y2$coef[2]*7+offset,
               xend=7, yend=model_ols_y2$coef[1]+model_ols_y2$coef[2]*15+offset,
               color='#7cae00',size=1,aes(linetype='counterfactul'))+
  scale_linetype_manual(values=c('dotted','solid'))+
  labs(title="Yearly ITS analysis: OLS model with Africanized effect",subtitle='Deformed wing virus',
       x='Years after Varroa arrived',y='Viral Transcript per Hundred')+
  guides(shape = guide_legend(order=1),linetype=guide_legend(order=2))+
  scale_x_continuous(breaks = c(seq(from = -7, to = 7, by = 1)))






