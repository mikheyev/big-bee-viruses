library(dplyr)
library(car) #dwt()
library(forecast) #tsdisplay()
library(ggplot2)
library(nlme)
library(lmtest)
library(imputeTS) #na_kalman()
library(reshape2) #melt()
library(ggsci)
library(ResourceSelection)  #hoslem.test()

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
DWV.m <- data.m[ ,c(1,7,19:21)]

# A preliminary OLS regression
model_ols_m <- lm(data=DWV.m,DWV~time+level+trend)
summary(model_ols_m)
AIC(model_ols_m)  #657.98
DWV.m$fitted_ols <- model_ols_m$fitted.values

###########################
# Assessing Autocorrelation
###########################

# Durbin-Watson test: test for correlated residuals
dwt(model_ols_m,alternative='two.sided')  #p-value=0.454

# Graph the residuals from the OLS regression to check for serially correlated error
tsdisplay(model_ols_m$residuals)  # no significant spikes in ACF or PACF
auto.arima(DWV.m$DWV, xreg=as.matrix(DWV.m[ ,3:5]),trace=T,max.p = 20,max.q = 20)  #ARIMA(0,0,0) with constana is same as OLS
                                                                                   #ARIMA(2,0,2),(1,0,1):inf

######################
# Run the final model
#####################

# fit the GLS regression model
## ARMA(1,1)
model_p1q1_m <- gls(data=DWV.m,DWV~time+level+trend,
                   correlation=corARMA(p=1,q=1,form=~time),
                   method="ML")
summary(model_p1q1_m)  #AIC=658.82
## ARMA(2,2)
model_p2q2_m <- gls(data=DWV.m,DWV~time+level+trend,
                   correlation=corARMA(p=2,q=2,form=~time),
                   method="ML")
summary(model_p2q2_m)  #AIC=661.22

##################
# Diagnostic tests
##################

# model checking: likelihood-ration tests to check ARMA process
anova(model_p1q1_m,model_p2q2_m)  #p-value=0.45
model_p3q3_m <- update(model_p2q2_m,correlation=corARMA(p=3,q=3,form=~time))
anova(model_p2q2_m,model_p3q3_m)  #p-value=0.037
# -> ARMA(1,1) is better

DWV.m$fitted_gls_p1q1 <- c(fitted(model_p1q1_m))

# plot
DWV.m_melt <- melt(DWV.m, id=c('month','DWV','time','level','trend'))

ggplot()+
  geom_point(data=DWV.m,aes(x=month,y=DWV,shape='observed'))+
  geom_vline(xintercept = 0, linetype='dashed',color='darkgrey')+
  geom_line(data=DWV.m_melt[which(DWV.m_melt$level==0), ],aes(x=month,y=value,color=variable,linetype='fitted'),size=1)+
  geom_line(data=DWV.m_melt[which(DWV.m_melt$level==1), ],aes(x=month,y=value,color=variable,linetype='fitted'),size=1)+
  geom_segment(x=-1,   y   =model_ols_m$coef[1]+model_ols_m$coef[2]*86,
               xend=51,yend=model_ols_m$coef[1]+model_ols_m$coef[2]*138,
               color='#f8766d',size=1,aes(linetype='counterfactul'))+
  geom_segment(x=-1,   y   =model_p1q1_m$coef[1]+model_p1q1_m$coef[2]*86,
               xend=51,yend=model_p1q1_m$coef[1]+model_p1q1_m$coef[2]*138,
               color='#00bfc4',size=1,aes(linetype='counterfactul'))+
  scale_color_hue(labels=c('OLS',"GLS_ARMA(1,1)"))+
  scale_linetype_manual(values=c('dotted','solid'))+
  labs(title="Monthly ITS analysis",subtitle='Deformed wing virus',
       x='Months after Varroa arrived',y='Viral Transcript per Hundred',
       color='models:')+
  guides(shape = guide_legend(order=1),linetype=guide_legend(order=2),color=guide_legend(order=3))+
  scale_x_continuous(breaks = c(seq(from = -86, to = 51, by = 4)))

# Hosmer-Lemeshow Goodness of Fit
hoslem.test(DWV.m$DWV, fitted(model_p1q1_m))  #p-value=1

#--------------------------------------------------------------------------------------------------------
# quarterly,DWV
DWV.q <- data.q[ ,c(1,7,19:21)]

# A preliminary OLS regression
model_ols_q <- lm(data=DWV.q,DWV~time+level+trend)
summary(model_ols_q)
AIC(model_ols_q)  #299.99
DWV.q$fitted_ols <- model_ols_q$fitted.values

###########################
# Assessing Autocorrelation
###########################

# Durbin-Watson test: test for correlated residuals
dwt(model_ols_q,alternative='two.sided')  #p-value=0.456

# Graph the residuals from the OLS regression to check for serially correlated error
tsdisplay(model_ols_q$residuals)  # AR(2)/MA(2)/ARMA(2,2)
auto.arima(DWV.q$DWV, xreg=as.matrix(DWV.q[ ,3:5]),trace=T,max.p = 15,max.q = 15)  #ARIMA(0,0,0) with constant is same as OLS
                                                                                   #ARIMA(1,0,1)/(2,0,2):inf

######################
# Run the final model
#####################

# fit the GLS regression model
## ARMA(1,1)
model_p1q1_q <- gls(data=DWV.q,DWV~time+level+trend,
                    correlation=corARMA(p=1,q=1,form=~time),
                    method="ML")
summary(model_p1q1_q)  #AIC=303.29
## ARMA(2,2)
model_p2q2_q <- gls(data=DWV.q,DWV~time+level+trend,
                    correlation=corARMA(p=2,q=2,form=~time),
                    method="ML")
summary(model_p2q2_q)  #AIC=301.72

##################
# Diagnostic tests
##################

# model checking: likelihood-ration tests to check ARMA process
anova(model_p1q1_q,model_p2q2_q)  #p-value=0.0617
model_p3q3_q <- update(model_p2q2_q,correlation=corARMA(p=3,q=3,form=~time))
anova(model_p2q2_q,model_p3q3_q)  #p-value=0.0645

DWV.q$fitted_gls_p1q1 <- c(fitted(model_p1q1_q))
DWV.q$fitted_gls_p2q2 <- c(fitted(model_p2q2_q))

#plot
DWV.q_melt <- melt(DWV.q, id=c('quarter','DWV','time','level','trend'))

ggplot()+
  geom_point(data=DWV.q,aes(x=quarter,y=DWV,shape='observed'))+
  geom_vline(xintercept = 0, linetype='dashed',color='darkgrey')+
  geom_line(data=DWV.q_melt[which(DWV.q_melt$level==0), ],aes(x=quarter,y=value,color=variable,linetype='fitted'),size=1)+
  geom_line(data=DWV.q_melt[which(DWV.q_melt$level==1), ],aes(x=quarter,y=value,color=variable,linetype='fitted'),size=1)+
  geom_segment(x=-1,  y=    model_ols_q$coef[1]+model_ols_q$coef[2]*28,
               xend=17,yend=model_ols_q$coef[1]+model_ols_q$coef[2]*46,
               color='#f8766d',size=1,aes(linetype='counterfactul'))+
  geom_segment(x=-1,  y=    model_p1q1_q$coef[1]+model_p1q1_q$coef[2]*28,
               xend=17,yend=model_p1q1_q$coef[1]+model_p1q1_q$coef[2]*46,
               color='#00ba38',size=1,aes(linetype='counterfactul'))+
  geom_segment(x=-1,  y=    model_p2q2_q$coef[1]+model_p2q2_q$coef[2]*28,
               xend=17,yend=model_p2q2_q$coef[1]+model_p2q2_q$coef[2]*46,
               color='#619cff',size=1,aes(linetype='counterfactul'))+
  scale_color_hue(labels=c('OLS',"GLS_ARMA(1,1)","GLS_ARMA(2,2)"))+
  scale_linetype_manual(values=c('dotted','solid'))+
  labs(title="Quarterly ITS analysis",subtitle='Deformed wing virus',
       x='Quarters after Varroa arrived',y='Viral Transcript per Hundred',
       color='models:')+
  guides(shape = guide_legend(order=1),linetype=guide_legend(order=2),color=guide_legend(order=3))+
  scale_x_continuous(breaks = c(seq(from = -28, to = 17, by = 1)))

# Hosmer-Lemeshow Goodness of Fit
hoslem.test(DWV.q$DWV, fitted(model_p2q2_q))  #p-value=1

#--------------------------------------------------------------------------------------------------------
# yearly.1,DWV
DWV.y1 <- data.y1[ ,c(1,7,19:21)]

# A preliminary OLS regression
model_ols_y1 <- lm(data=DWV.y1,DWV~time+level+trend)
summary(model_ols_y1)
AIC(model_ols_y1)  #76.78
DWV.y1$fitted_ols <- model_ols_y1$fitted.values

###########################
# Assessing Autocorrelation
###########################

# Durbin-Watson test: test for correlated residuals
dwt(model_ols_y1,alternative='two.sided')  #p-value=0.334

# Graph the residuals from the OLS regression to check for serially correlated error
tsdisplay(model_ols_y1$residuals)  # AR(2)/MA(2)/ARMA(2,2)
auto.arima(DWV.y1$DWV, xreg=as.matrix(DWV.y1[ ,3:5]),trace=T)  #ARIMA(0,0,0)
                                                               #ARIMA(1,0,1)/(2,0,2)/(0,0,1):inf

######################
# Run the final model
#####################

# fit the GLS regression model
## AR(2)
model_p2_y1 <- gls(data=DWV.y1,DWV~time+level+trend,
                    correlation=corARMA(p=2,form=~time),
                    method="ML")
summary(model_p2_y1)  #AIC=74.56
## MA(1)
model_q1_y1 <- gls(data=DWV.y1,DWV~time+level+trend,
                   correlation=corARMA(q=1,form=~time),
                   method="ML")
summary(model_q1_y1)  #AIC=75.59
## MA(2)
model_q2_y1 <- gls(data=DWV.y1,DWV~time+level+trend,
                   correlation=corARMA(q=2,form=~time),
                   method="ML")
summary(model_q2_y1)  #AIC=75.69
## ARMA(1,1)
model_p1q1_y1 <- gls(data=DWV.y1,DWV~time+level+trend,
                     correlation=corARMA(p=1,q=1,form=~time),
                     method="ML")
summary(model_p1q1_y1)  #AIC=77.29
## ARMA(2,2)
model_p2q2_y1 <- gls(data=DWV.y1,DWV~time+level+trend,
                    correlation=corARMA(p=2,q=2,form=~time),
                    method="ML")
summary(model_p2q2_y1)  #AIC=75.29

##################
# Diagnostic tests
##################

# model checking: likelihood-ration tests to check ARMA process
model_p3_y1 <- update(model_p2_y1,correlation=corARMA(p=3,form=~time))
anova(model_p2_y1,model_p3_y1)  #p-value=0.6692
anova(model_q1_y1,model_q2_y1)  #p-value=0.1685
model_q3_y1 <- update(model_q2_y1,correlation=corARMA(q=3,form=~time))
anova(model_q2_y1,model_q3_y1)  #p-value=0.0793
anova(model_p1q1_y1,model_p2q2_y1)  #p-value=0.0497
model_p3q3_y1 <- update(model_p2q2_y1,correlation=corARMA(p=3,q=3,form=~time))
anova(model_p2q2_y1,model_p3q3_y1)  #p-value=0.1397
# -> remove ARMA(1,1)

DWV.y1$fitted_gls_p2 <- c(fitted(model_p2_y1))
DWV.y1$fitted_gls_p2q2 <- c(fitted(model_p2q2_y1))

# plot
DWV.y1_melt <- melt(DWV.y1, id=c('year.1','DWV','time','level','trend'))

ggplot()+
  geom_point(data=DWV.y1,aes(x=year.1,y=DWV,shape='observed'))+
  geom_vline(xintercept = 0, linetype='dashed',color='darkgrey')+
  geom_line(data=DWV.y1_melt[which(DWV.y1_melt$level==0), ],aes(x=year.1,y=value,color=variable,linetype='fitted'),size=1)+
  geom_line(data=DWV.y1_melt[which(DWV.y1_melt$level==1), ],aes(x=year.1,y=value,color=variable,linetype='fitted'),size=1)+
  geom_segment(x=-1,  y=    model_ols_y1$coef[1]+model_ols_y1$coef[2]*7,
               xend=4,yend=model_ols_y1$coef[1]+model_ols_y1$coef[2]*12,
               color='#f8766d',size=1,aes(linetype='counterfactul'))+
  geom_segment(x=-1,  y=    model_p2_y1$coef[1]+model_p2_y1$coef[2]*7,
               xend=4,yend=model_p2_y1$coef[1]+model_p2_y1$coef[2]*12,
               color='#00ba38',size=1,aes(linetype='counterfactul'))+
  geom_segment(x=-1,  y=    model_p2q2_y1$coef[1]+model_p2q2_y1$coef[2]*7,
               xend=4,yend=model_p2q2_y1$coef[1]+model_p2q2_y1$coef[2]*12,
               color='#619cff',size=1,aes(linetype='counterfactul'))+
  scale_color_hue(labels=c('OLS',"GLS_AR(2)","GLS_ARMA(2,2)"))+
  scale_linetype_manual(values=c('dotted','solid'))+
  labs(title="Yearly ITS analysis",subtitle='Deformed wing virus',
       x='Years after Varroa arrived',y='Viral Transcript per Hundred',
       color='models:')+
  guides(shape = guide_legend(order=1),linetype=guide_legend(order=2),color=guide_legend(order=3))+
  scale_x_continuous(breaks = c(seq(from = -7, to = 4, by = 1)))

# Hosmer-Lemeshow Goodness of Fit
hoslem.test(DWV.y1$DWV, fitted(model_p2_y1))  #p-value=1

#--------------------------------------------------------------------------------------------------------
# yearly.2,DWV
DWV.y2 <- data.y2[ ,c(1,7,19:21)]

# A preliminary OLS regression
model_ols_y2 <- lm(data=DWV.y2,DWV~time+level+trend)
summary(model_ols_y2)
AIC(model_ols_y2)  #122.20
DWV.y2$fitted_ols <- model_ols_y2$fitted.values

###########################
# Assessing Autocorrelation
###########################

# Durbin-Watson test: test for correlated residuals
dwt(model_ols_y2,alternative='two.sided')  #p-value=0.122

# Graph the residuals from the OLS regression to check for serially correlated error
tsdisplay(model_ols_y2$residuals)  # AR(1)/MA(1)/ARMA(1,1)
auto.arima(DWV.y2$DWV, xreg=as.matrix(DWV.y2[ ,3:5]),trace=T)  #ARIMA(0,0,0)
                                                               #ARIMA(1,0,1)/(2,0,2)/(0,0,1):inf

######################
# Run the final model
#####################

# fit the GLS regression model
## AR(1)
model_p1_y2 <- gls(data=DWV.y2,DWV~time+level+trend,
                   correlation=corARMA(p=1,form=~time),
                   method="ML")
summary(model_p1_y2)  #AIC=117.76
## MA(1)
model_q1_y2 <- gls(data=DWV.y2,DWV~time+level+trend,
                   correlation=corARMA(q=1,form=~time),
                   method="ML")
summary(model_q1_y2)  #AIC=114.14
## ARMA(1,1)
model_p1q1_y2 <- gls(data=DWV.y2,DWV~time+level+trend,
                   correlation=corARMA(p=1,q=1,form=~time),
                   method="ML")
summary(model_p1q1_y2)  #AIC=112.82
## ARMA(2,2)
model_p2q2_y2 <- gls(data=DWV.y2,DWV~time+level+trend,
                     correlation=corARMA(p=2,q=2,form=~time),
                     method="ML")
summary(model_p2q2_y2)  #AIC=98.88

##################
# Diagnostic tests
##################

# model checking: likelihood-ration tests to check ARMA process
model_p2_y2 <- update(model_p1_y2,correlation=corARMA(p=2,form=~time))
anova(model_p1_y2,model_p2_y2)  #p-value=0.002
model_q2_y2 <- update(model_q1_y2,correlation=corARMA(q=2,form=~time))
anova(model_q1_y2,model_q2_y2)  #p-value=0.0085
anova(model_p1q1_y2,model_p2q2_y2)  #p-value=1e-04
model_p3q3_y2 <- update(model_p2q2_y2,correlation=corARMA(p=3,q=3,form=~time))
anova(model_p2q2_y2,model_p3q3_y2)  #p-value=0.0079
# -> remove all models

#plot
DWV.y2_melt <- melt(DWV.y2, id=c('year.2','DWV','time','level','trend'))

ggplot()+
  geom_point(data=DWV.y2,aes(x=year.2,y=DWV,shape='observed'))+
  geom_vline(xintercept = 0, linetype='dashed',color='darkgrey')+
  geom_line(data=DWV.y2[which(DWV.y2$level==0), ],aes(x=year.2,y=fitted_ols,linetype='fitted'),size=1,color='#f8766d')+
  geom_line(data=DWV.y2[which(DWV.y2$level==1), ],aes(x=year.2,y=fitted_ols,linetype='fitted'),size=1,color='#f8766d')+
  geom_segment(x=-1,  y   =model_ols_y2$coef[1]+model_ols_y2$coef[2]*7,
               xend=7,yend=model_ols_y2$coef[1]+model_ols_y2$coef[2]*15,
               color='#f8766d',size=1,aes(linetype='counterfactul'))+
  scale_linetype_manual(values=c('dotted','solid'))+
  labs(title="Yearly ITS analysis: OLS model",subtitle='Deformed wing virus',
       x='Years after Varroa arrived',y='Viral Transcript per Hundred')+
  guides(shape = guide_legend(order=1),linetype=guide_legend(order=2))+
  scale_x_continuous(breaks = c(seq(from = -7, to = 7, by = 1)))

# Hosmer-Lemeshow Goodness of Fit
hoslem.test(DWV.y2$DWV, fitted(model_ols_y2))  #p-value=1
