library(rsample)
library(dplyr)
library(rpart)
library(rpart.plot)

rm(list = ls())

# loading data
data <- read.csv("new_data.csv")

# data exploration: summary(), str()=glimpse(), execute in console
# data preprocessing: remove or change variables 
data$Geocluster <- as.factor(data$Geocluster)
data <- data[ ,-c(2,3,6,7,8,9,11,12,13,14,15,16,17)]
rownames(data) <- data$ID
data <- data[ ,-1]


# building regression tree with all variables
# using tree models to see interaction is my goal, 
# so I did not split the data into train & test and not tuning trees

##ABPV
set.seed(120)
m1 <- rpart(
  formula = ABPV ~ .,
  data    = data,
  method  = "anova"
)
rpart.plot(m1, tweak=1.2) # leaf node=9
plotcp(m1)
m1$cptable # nsplit=8, cp=0.01, xerror=0.6383622
m1$variable.importance # AM,BQCV,CBPV,DWV,KBV,LSV,Mo,SV,Yeae,Afr.Extent,Geocluster,V.arrival...

#---------------------------------------------------------------------

##AM
set.seed(120)
m2 <- rpart(
  formula = AM ~ .,
  data    = data,
  method  = "anova"
)
rpart.plot(m2, tweak=1.2) # leaf node=7
plotcp(m2)
m2$cptable # nsplit=6, cp=0.01, xerror=0.3238932
m2$variable.importance  # SV,DWV,BQCV,ABPV,LSV,Year,Geocluster,V.arrival...

#------------------------------------------------------------------

##BQCV
set.seed(120)
m3 <- rpart(
  formula = BQCV ~ .,
  data    = data,
  method  = "anova"
)
rpart.plot(m3, tweak=1.2) # leaf node=11
plotcp(m3)
m3$cptable # nsplit=10, cp=0.01, xerror=0.1618887
m3$variable.importance  # LSV,SV,ABPV,AM,DWV,CBPV,Year,KBV,VDV1,IAPV,Geocluster,Afr.Extent

