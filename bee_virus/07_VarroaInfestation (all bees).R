library(ggplot2)
library(vegan)
library(RVAideMemoire)
library(dplyr)
library(reshape2)
library(lme4)
rm(list=ls())

data <- read.csv('new_data.csv')
data$V.arrival <- as.factor(data$V.arrival)
data$year.2 <- factor(data$year.2, levels=c('0','<1','<2','<3','<5','<7'))
data$year.3 <- factor(data$year.3, levels=c('0','<1','<2','>2'))

data.p <- data
data.t <- data

## viral prevalence: chi-squared test
## convert continuous variable (viral tpm) into categorical variable
for(i in 25:41){
  for(l in 1:1034){
    if(data.p[l,i]>0){data.p[l,i] <- 1}
    else{data.p[l,i] <- 0}
  }
}

b <- c()
for (i in c(26,35,28,37,25,30)) {
  set.seed(123)
  dt <- table(data.p$V.arrival,data.p[ ,i])
  a <- chisq.test(dt,simulate.p.value = T)
  p <- paste(names(data.p)[i],'x-squared:',a$statistic,'p-value:',a$p.value)
  print(dt)
  print(p)
  b <- append(b,a$p.value)
}
# adjusted p-value: false discover rate -> significant: BQCV,SBV
p.adjust(b,method='fdr')

## viral prevalence: bar chart-> x=virus, y=prevalence
## step1: calculate the viral prevalence
data_prevalence <- data.p %>% group_by(V.arrival) %>% 
  summarise(ABPV=length(which(ABPV==1))/length(ID),Amfv=length(which(Amfv==1))/length(ID),
            Arv2=length(which(Arv2==1))/length(ID),BQCV=length(which(BQCV==1))/length(ID),
            CBPV=length(which(CBPV==1))/length(ID),DWV_A=length(which(DWV_A==1))/length(ID),
            DWV_B=length(which(DWV_B==1))/length(ID),DWV_C=length(which(DWV_C==1))/length(ID),
            IAPV=length(which(IAPV==1))/length(ID),KBV=length(which(KBV==1))/length(ID),
            LSV=length(which(LSV==1))/length(ID),SBPV=length(which(SBPV==1))/length(ID),
            SBV=length(which(SBV==1))/length(ID),VDV2=length(which(VDV2==1))/length(ID),
            VDV3=length(which(VDV3==1))/length(ID),Vov1=length(which(Vov1==1))/length(ID),
            VTV=length(which(VTV==1))/length(ID))

data_prevalence_all <- data.p %>%
  summarise(ABPV=length(which(ABPV==1))/length(ID),Amfv=length(which(Amfv==1))/length(ID),
            Arv2=length(which(Arv2==1))/length(ID),BQCV=length(which(BQCV==1))/length(ID),
            CBPV=length(which(CBPV==1))/length(ID),DWV_A=length(which(DWV_A==1))/length(ID),
            DWV_B=length(which(DWV_B==1))/length(ID),DWV_C=length(which(DWV_C==1))/length(ID),
            IAPV=length(which(IAPV==1))/length(ID),KBV=length(which(KBV==1))/length(ID),
            LSV=length(which(LSV==1))/length(ID),SBPV=length(which(SBPV==1))/length(ID),
            SBV=length(which(SBV==1))/length(ID),VDV2=length(which(VDV2==1))/length(ID),
            VDV3=length(which(VDV3==1))/length(ID),Vov1=length(which(Vov1==1))/length(ID),
            VTV=length(which(VTV==1))/length(ID))
data_prevalence_all <- data_prevalence_all[ ,names(sort(data_prevalence_all[1,],decreasing=T))]

## step2: ggplot plot the barcharts
data_prevalence_melt <- melt(data_prevalence[,c('V.arrival',names(sort(data_prevalence_all[1,],decreasing=T)))], id=c('V.arrival'))


ggplot(data_prevalence_melt,aes(x=variable,y=value*100,fill=as.factor(V.arrival)))+
  geom_bar(stat='identity', width=0.65, position=position_dodge(),color='black')+
  labs(x='',y='Pathogen prevalence (% positive individual bees)',title='Africanized bee (N=735)')+
  scale_fill_manual(labels=c('Varroa-free (n=647)','Varroa-infested (n=88)'),values=c('white','grey'))+
  ylim(c(0,100))+
  theme_classic()+
  theme(legend.position = c(0.905,0.99),
        legend.background = element_rect(linetype="solid",colour ="darkgrey"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.text.x = element_text(size=11))+
  scale_x_discrete(labels=c('Amfv\n91.3','LSV\n72.8','BQCV\n55.4','SBV\n34.0','ABPV\n19.3','DWV-A\n14.1','CBPV\n5.58',
                            'KBV\n1.63','Arv2\n1.50','IAPV\n1.22','DWV-B\n1.09','DWV-C\n0.41','Vov1\n0.41','SBPV\n0.27',
                            'VTV\n0.27','VDV2\n0.14','VDV3\n0.14'))


ggplot(data_prevalence_melt[1:12,],aes(x=variable,y=value*100,fill=as.factor(V.arrival)))+
  geom_bar(stat='identity',width=0.65, position=position_dodge(),color='black')+
  labs(x='',y='Pathogen prevalence (% positive individual bees)',title='Honey bees (N=1034)')+
  scale_fill_manual(labels=c('Varroa-free (n=936)','Varroa-infested (n=98)'),values=c('white','grey'))+
  ylim(c(0,100))+
  geom_text(aes(label=round(value*100,1)),position = position_dodge(0.7),
            vjust = -0.3, size = 3.5)+
  theme_classic()+
  theme(legend.position = c(0.87,0.99),
        legend.background = element_rect(linetype="solid",colour ="darkgrey"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.text.x = element_text(size=11))+
  scale_x_discrete(labels=c('Amfv\n91.5','LSV\n75.5','BQCV\n58.3','SBV\n40.1',
                            'ABPV\n17.5','DWV-A\n12.6'))

#-----------------------------------------------------------------------------------

## GLMM: viral prevalence, logistics regression: Amfv,LSV,BQCV,SBV,ABPV,DWV-A
## Varroa presence: BQCV, ABPV, DWV-A
a<- glmer(DWV_A~V.arrival+(1|Geocluster)+(1|ID), 
          family=binomial(link = 'logit'), data=data.p)
summary(a)

## BQCV
BQCV <- glmer(BQCV~V.arrival+Afr.Extent+Latitude+PRCP+TMAX+TMIN+(1|Geocluster)+(1|ID), 
              family=binomial(link = 'logit'), data=data.p)
summary(BQCV) #AIC=1357.4

BQCV_reduced <- glmer(BQCV~V.arrival+Afr.Extent+TMAX+(1|Geocluster)+(1|ID), 
                      family=binomial(link = 'logit'), data=data.p)
summary(BQCV_reduced)

anova(BQCV_reduced, BQCV, test='F') #p=0.2606

## ABPV
ABPV <- glmer(ABPV~V.arrival+Afr.Extent+Latitude+PRCP+TMAX+TMIN+(1|Geocluster)+(1|ID), 
              family=binomial(link = 'logit'), data=data.p)
summary(ABPV) #AIC=957.3

ABPV_reduced <- glmer(ABPV~V.arrival+Afr.Extent+(1|Geocluster)+(1|ID), 
              family=binomial(link = 'logit'), data=data.p)
summary(ABPV_reduced)

anova(ABPV_reduced, ABPV, test='F') #p=0.9036


## DWV-A
DWV <- glmer(DWV_A~V.arrival+Afr.Extent+Latitude+PRCP+TMAX+TMIN+(1|Geocluster)+(1|ID), 
              family=binomial(link = 'logit'), data=data.p)
summary(DWV) #AIC=749.8

DWV_reduced <- glmer(DWV_A~V.arrival+Afr.Extent+TMAX+(1|Geocluster)+(1|ID), 
             family=binomial(link = 'logit'), data=data.p)
summary(DWV_reduced)

anova(DWV_reduced, DWV, test='F') #p=0.9036

#-----------------------------------------------------------------------------------------
# Varria exposure time: 0, <1, <2, >2
# statistical significance: chi-square test
b <- c()
for (i in c(26,35,28,37,25,30)) {
  set.seed(123)
  dt <- table(data.p[ ,i],data.p$year.3)
  a <- chisq.test(dt,simulate.p.value = T)
  p <- paste(names(data.p)[i],'x-squared:',a$statistic,'p-value:',a$p.value)
  print(dt)
  print(p)
  b <- append(b,a$p.value)
}
# adjusted p-value: false discover rate -> significant: LSV,BQCV,ABPV,DWV-A
p.adjust(b,method='fdr')

# pairwise chi-square test
library(rcompanion)
for (i in c(26,35,28,37,25,30)) {
  dt <- table(data.p$year.3, data.p[ ,i])
  a <- pairwiseNominalIndependence(dt,fisher=F,gtest=F,chisq=T,method='fdr')
  p <- paste(names(data.p)[i],a$Comparison,'p-value:',a$p.adj.Chisq)
  print(dt)
  print(p)
}


name_full <- c('Acute bee paralysis virus','Apis mellifera filamentous virus',
               'Apis rhabdovirus 2','Black queen cell virus', 'Chronic bee paralysis virus',
               'Deformed wing virus-A','Deformed wing virus-B','Deformed wing virus-C',
               'Israeli acute paralysis virus','Kashmir bee virus','Lake Sinai virus',
               'Slow bee paralysis virus','Sacbrood virus','Varroa destructor virus 2',
               'Varroa destructor virus 3','Varroa orthomyxovirus-1','Varroa Tymo-like virus')

data_prevalence_year.3 <- data.p %>% group_by(year.3,V.arrival) %>% 
  summarise(ABPV=length(which(ABPV==1))/length(ID),Amfv=length(which(Amfv==1))/length(ID),
            Arv2=length(which(Arv2==1))/length(ID),BQCV=length(which(BQCV==1))/length(ID),
            CBPV=length(which(CBPV==1))/length(ID),DWV_A=length(which(DWV_A==1))/length(ID),
            DWV_B=length(which(DWV_B==1))/length(ID),DWV_C=length(which(DWV_C==1))/length(ID),
            IAPV=length(which(IAPV==1))/length(ID),KBV=length(which(KBV==1))/length(ID),
            LSV=length(which(LSV==1))/length(ID),SBPV=length(which(SBPV==1))/length(ID),
            SBV=length(which(SBV==1))/length(ID),VDV2=length(which(VDV2==1))/length(ID),
            VDV3=length(which(VDV3==1))/length(ID),Vov1=length(which(Vov1==1))/length(ID),
            VTV=length(which(VTV==1))/length(ID))

data_prevalence_year.3[,3:19] <- data_prevalence_year.3[,3:19]*100

for (i in 3:19) {
  a <-ggplot(data_prevalence_year.3,aes_string(x='year.3',y=names(data_prevalence_year.3)[i],fill='year.3'))+
    geom_bar(stat='identity',color='black')+
    scale_fill_manual(labels=c('  0 year   (n=647)','<1 year   (n=46)','<2 years (n=18)','>2 years (n=24)'),
                      values = c('white','grey90','grey60','grey30'))+
    labs(title=paste(name_full[i-2],',',names(data_prevalence_year.3)[i]),
         x='Varroa exposure time (years)',y='Pathogen prevalence (% positive individual bees)',
         fill='Africanized bee (N=735)\n\nVarroa exposure time')+
    ylim(c(0,110))+
    geom_text(aes(label=c(round(data_prevalence_year.3[1,i],1),round(data_prevalence_year.3[2,i],1),round(data_prevalence_year.3[3,i],1),round(data_prevalence_year.3[4,i],1))),
              position = position_dodge(0.7),vjust = -0.3, size = 3.5)+
    theme_classic()+
    theme(legend.position = 'none',
          legend.title = element_text(size=12),
          legend.background = element_rect(linetype="solid",colour ="darkgrey"),
          legend.text = element_text(size=12),
          axis.text.x = element_text(size=11))
  print(a)
}

#------------------------------------------------------------------------------------
## Varroa exposure time: 0, <1, <2, >2
## GLMM: viral prevalence, logistics regression -> LSV, BQCV, DWV-A
a <- glmer(DWV_A~year.3+(1|Geocluster)+(1|ID), 
           family=binomial(link = 'logit'), data=data.p)
summary(a)

## LSV 
LSV <- glmer(LSV~year.3+Afr.Extent+Latitude+PRCP+TMAX+TMIN+(1|Geocluster)+(1|ID), 
             family=binomial(link = 'logit'), data=data.p)
summary(LSV) #AIC=1148.6

LSV_reduced <- glmer(LSV~year.3+Afr.Extent+Latitude+PRCP+(1|Geocluster)+(1|ID), 
                     family=binomial(link = 'logit'), data=data.p)
summary(LSV_reduced)

## BQCV
BQCV <- glmer(BQCV~year.3+Afr.Extent+Latitude+PRCP+TMAX+TMIN+(1|Geocluster)+(1|ID), 
              family=binomial(link = 'logit'), data=data.p)
summary(BQCV) #AIC=1351.2

BQCV_reduced <- glmer(BQCV~year.3+Afr.Extent+Latitude+TMAX+(1|Geocluster)+(1|ID), 
              family=binomial(link = 'logit'), data=data.p)
summary(BQCV_reduced)

## DWV-A
DWV_A <-  glmer(DWV_A~year.3+Afr.Extent+Latitude+PRCP+TMAX+TMIN+(1|Geocluster)+(1|ID), 
                family=binomial(link = 'logit'), data=data.p)
summary(DWV_A) #AIC=737.4

DWV_A_reduced <-  glmer(DWV_A~year.3+Afr.Extent+TMAX+(1|Geocluster)+(1|ID), 
                family=binomial(link = 'logit'), data=data.p)
summary(DWV_A_reduced)

a <- glmer(LSV~year.3+Afr.Extent+(1|Geocluster)+(1|ID), 
           family=binomial(link = 'logit'), data=data.p)
summary(a)

#========================================================================================================================

## viral tpm: permutation median test (are there any tpm differences between Before & After Varroa)
a <- c()
for (i in c(26,35,28,37,25,30)) {
  set.seed(123)
  n <- length(data.t$V.arrival)
  P <- 100000
  variable <- data.t[ ,i]
  test.stat.1 <- abs(median(data.t[data.t$V.arrival=='0',i])-median(data.t[data.t$V.arrival=='1',i]))
  
  PerSamples <- matrix(0, nrow=n, ncol=P)
  for(j in 1:P){
    PerSamples[ ,j] <- sample(variable, size=n, replace=F)
  }
  
  Perm.test.stat.1 <- rep(0,P)
  for(k in 1:P){
    Perm.test.stat.1[k] <- abs(median(PerSamples[data.t$V.arrival=='0',k])-
                                 median(PerSamples[data.t$V.arrival=='1',k]))
  }
  print(paste(names(data.t)[i],"Median:",mean(Perm.test.stat.1>=test.stat.1)))
  a <- append(a,mean(Perm.test.stat.1>=test.stat.1))
}

# adjusted p-value: false discover rate -> significant: SBV
p.adjust(a,method='fdr')


## viral tpm: boxplot-> x=viruses, y=tpm
data_median <- data.t[,25:41] %>% summarise_all(median)
data_mean <- data.t[,25:41] %>% summarise_all(mean)

names(sort(data_mean[1,],decreasing=T))

data_melt <- melt(data.t[,c(26,35,28,37,25,30,24)], id=c('V.arrival'))

ggplot(data_melt,aes(x=variable,y=value,fill=as.factor(V.arrival)))+
  geom_boxplot(outlier.shape = NA)+
  labs(x='',y='Transcript per million',title='Honey bees (N=1034)')+
  scale_y_sqrt()+
  scale_fill_manual(labels=c('Varroa-free (n=936)','Varroa-infested (n=98)'),values=c('white','grey'))+
  theme_classic()+
  theme(legend.position = c(0.87,0.99),
        legend.background = element_rect(linetype="solid",colour ="darkgrey"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size=11))+
  scale_x_discrete(labels=c('Amfv\n3350.72','LSV\n265291.5','BQCV\n24747.7','SBV\n0','ABPV\n0','DWV-A\n0'))

#-----------------------------------------------------------------------------------------------------------

# tobit model: LSV,BQCV,ABPV,DWV-a
library(VGAM)
a <- vglm(data=data.t,Amfv~V.arrival,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(a)

library(MCMCpack)
## LSV
LSV_tobit <- vglm(data=data.t,LSV~V.arrival+Afr.Extent+Latitude+PRCP+TMAX+TMIN,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(LSV_tobit)
step4vglm(LSV_tobit,direction='both',steps=10000,scope=list(lower=LSV~V.arrival+Afr.Extent)) 

LSV_tobit_reduced <- vglm(data=data.t,LSV~V.arrival+Afr.Extent+PRCP,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(LSV_tobit_reduced)

## BQCV
BQCV_tobit <- vglm(data=data.t,BQCV~V.arrival+Afr.Extent+Latitude+PRCP+TMAX+TMIN,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(BQCV_tobit)
step4vglm(BQCV_tobit,direction='both',steps=10000,scope=list(lower=BQCV~V.arrival+Afr.Extent)) 

BQCV_tobit_reduced <- vglm(data=data.t,BQCV~V.arrival+Afr.Extent+Latitude+TMIN+TMAX,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(BQCV_tobit_reduced)

## ABPV
ABPV_tobit <- vglm(data=data.t,ABPV~V.arrival+Afr.Extent+Latitude+PRCP+TMAX+TMIN,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(ABPV_tobit)
step4vglm(ABPV_tobit,direction='both',steps=10000,scope=list(lower=ABPV~V.arrival+Afr.Extent)) 

ABPV_tobit_reduced <- vglm(data=data.t,ABPV~V.arrival+Afr.Extent,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(ABPV_tobit_reduced)

## DWV_A
DWV_tobit <- vglm(data=data.t,DWV_A~V.arrival+Afr.Extent+Latitude+PRCP+TMAX+TMIN,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(DWV_tobit)
step4vglm(DWV_tobit,direction='both',steps=10000,scope=list(lower=DWV_A~V.arrival+Afr.Extent)) 

DWV_tobit_reduced <- vglm(data=data.t,DWV_A~V.arrival+Afr.Extent+PRCP+TMIN+TMAX,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(DWV_tobit_reduced)

#------------------------------------------------------------------------------------------------------------------------------

## viral tpm: pairwise permutation median test (are there any tpm differences between Before Varroa vs. Varroa exposure time)
data1 <- data[which(data$year.3=='0'|data$year.3=='<1'),]
a <- c()
for (i in c(26,35,28,37,25,30)) {
  set.seed(123)
  n <- length(data1$year.3)
  P <- 100000
  variable <- data1[ ,i]
  test.stat.1 <- abs(median(data1[data1$year.3=='0',i])-median(data1[data1$year.3=='<1',i]))
  
  PerSamples <- matrix(0, nrow=n, ncol=P)
  for(j in 1:P){
    PerSamples[ ,j] <- sample(variable, size=n, replace=F)
  }
  
  Perm.test.stat.1 <- rep(0,P)
  for(k in 1:P){
    Perm.test.stat.1[k] <- abs(median(PerSamples[data1$year.3=='0',k])-
                                 median(PerSamples[data1$year.3=='<1',k]))
  }
  print(paste(names(data1)[i],"Median:",mean(Perm.test.stat.1>=test.stat.1)))
  a <- append(a,mean(Perm.test.stat.1>=test.stat.1))
}

data2 <- data[which(data$year.3=='0'|data$year.3=='<2'),]
b <- c()
for (i in c(26,35,28,37,25,30)) {
  set.seed(123)
  n <- length(data2$year.3)
  P <- 100000
  variable <- data2[ ,i]
  test.stat.1 <- abs(median(data2[data2$year.3=='0',i])-median(data2[data2$year.3=='<2',i]))
  
  PerSamples <- matrix(0, nrow=n, ncol=P)
  for(j in 1:P){
    PerSamples[ ,j] <- sample(variable, size=n, replace=F)
  }
  
  Perm.test.stat.1 <- rep(0,P)
  for(k in 1:P){
    Perm.test.stat.1[k] <- abs(median(PerSamples[data2$year.3=='0',k])-
                                 median(PerSamples[data2$year.3=='<2',k]))
  }
  print(paste(names(data2)[i],"Median:",mean(Perm.test.stat.1>=test.stat.1)))
  b <- append(b,mean(Perm.test.stat.1>=test.stat.1))
}

data3 <- data[which(data$year.3=='0'|data$year.3=='>2'),]
c <- c()
for (i in c(26,35,28,37,25,30)) {
  set.seed(123)
  n <- length(data3$year.3)
  P <- 100000
  variable <- data3[ ,i]
  test.stat.1 <- abs(median(data3[data3$year.3=='0',i])-median(data3[data3$year.3=='>2',i]))
  
  PerSamples <- matrix(0, nrow=n, ncol=P)
  for(j in 1:P){
    PerSamples[ ,j] <- sample(variable, size=n, replace=F)
  }
  
  Perm.test.stat.1 <- rep(0,P)
  for(k in 1:P){
    Perm.test.stat.1[k] <- abs(median(PerSamples[data3$year.3=='0',k])-
                                 median(PerSamples[data3$year.3=='>2',k]))
  }
  print(paste(names(data3)[i],"Median:",mean(Perm.test.stat.1>=test.stat.1)))
  c <- append(c,mean(Perm.test.stat.1>=test.stat.1))
}

data4 <- data[which(data$year.3=='<1'|data$year.3=='<2'),]
d <- c()
for (i in c(26,35,28,37,25,30)) {
  set.seed(123)
  n <- length(data4$year.3)
  P <- 100000
  variable <- data4[ ,i]
  test.stat.1 <- abs(median(data4[data4$year.3=='<1',i])-median(data4[data4$year.3=='<2',i]))
  
  PerSamples <- matrix(0, nrow=n, ncol=P)
  for(j in 1:P){
    PerSamples[ ,j] <- sample(variable, size=n, replace=F)
  }
  
  Perm.test.stat.1 <- rep(0,P)
  for(k in 1:P){
    Perm.test.stat.1[k] <- abs(median(PerSamples[data4$year.3=='<1',k])-
                                 median(PerSamples[data4$year.3=='<2',k]))
  }
  print(paste(names(data4)[i],"Median:",mean(Perm.test.stat.1>=test.stat.1)))
  d <- append(d,mean(Perm.test.stat.1>=test.stat.1))
}


data5 <- data[which(data$year.3=='<1'|data$year.3=='>2'),]
e <- c()
for (i in c(26,35,28,37,25,30)) {
  set.seed(123)
  n <- length(data5$year.3)
  P <- 100000
  variable <- data5[ ,i]
  test.stat.1 <- abs(median(data5[data5$year.3=='<1',i])-median(data5[data5$year.3=='>2',i]))
  
  PerSamples <- matrix(0, nrow=n, ncol=P)
  for(j in 1:P){
    PerSamples[ ,j] <- sample(variable, size=n, replace=F)
  }
  
  Perm.test.stat.1 <- rep(0,P)
  for(k in 1:P){
    Perm.test.stat.1[k] <- abs(median(PerSamples[data5$year.3=='<1',k])-
                                 median(PerSamples[data5$year.3=='>2',k]))
  }
  print(paste(names(data5)[i],"Median:",mean(Perm.test.stat.1>=test.stat.1)))
  e <- append(e,mean(Perm.test.stat.1>=test.stat.1))
}


data6 <- data[which(data$year.3=='<2'|data$year.3=='>2'),]
f <- c()
for (i in c(26,35,28,37,25,30)) {
  set.seed(123)
  n <- length(data6$year.3)
  P <- 100000
  variable <- data6[ ,i]
  test.stat.1 <- abs(median(data6[data6$year.3=='<2',i])-median(data6[data6$year.3=='>2',i]))
  
  PerSamples <- matrix(0, nrow=n, ncol=P)
  for(j in 1:P){
    PerSamples[ ,j] <- sample(variable, size=n, replace=F)
  }
  
  Perm.test.stat.1 <- rep(0,P)
  for(k in 1:P){
    Perm.test.stat.1[k] <- abs(median(PerSamples[data6$year.3=='<2',k])-
                                 median(PerSamples[data6$year.3=='>2',k]))
  }
  print(paste(names(data6)[i],"Median:",mean(Perm.test.stat.1>=test.stat.1)))
  f <- append(f,mean(Perm.test.stat.1>=test.stat.1))
}

pairwise_perm <- rbind(a,b,c,d,e,f)

for(i in 1:6){
  pairwise_perm[,i] <- p.adjust(pairwise_perm[,i],method='fdr')
}

pairwise_perm <- as.data.frame(pairwise_perm)
names(pairwise_perm) <- names(data6[,c(26,35,28,37,25,30)])


name_full <- c('Acute bee paralysis virus','Apis mellifera filamentous virus',
               'Apis rhabdovirus 2','Black queen cell virus', 'Chronic bee paralysis virus',
               'Deformed wing virus-A','Deformed wing virus-B','Deformed wing virus-C',
               'Israeli acute paralysis virus','Kashmir bee virus','Lake Sinai virus',
               'Slow bee paralysis virus','Sacbrood virus','Varroa destructor virus 2',
               'Varroa destructor virus 3','Varroa orthomyxovirus-1','Varroa Tymo-like virus')

for (i in c(26,35,28,37,25,30)) {
  a <-ggplot(data=data, aes_string(x='year.3',y=names(data[i]), fill='year.3'))+
    geom_boxplot(outlier.shape = NA)+
    labs(title=paste(name_full[i-24],',',names(data[i])),x='Varroa exposure time (years)',y='Transcript per million',
         fill='Africanized bee (N=735)\n\nVarroa exposure time')+
    scale_y_sqrt()+
    scale_fill_manual(labels=c('  0 year   (n=647)','<1 year   (n=46)','<2 years (n=18)','>2 years (n=24)'),
                      values=c('white','grey90','grey60','grey30'))+
    theme_classic()+
    theme(legend.position='none',
          legend.title = element_text(size=12),
          legend.background = element_rect(linetype="solid",colour ="darkgrey"),
          legend.text = element_text(size=12),
          axis.text.x = element_text(size=11))
  print(a)
}

#-------------------------------------------------------------------------------------

# Varroa exposure time: 0, <1, <2, >2
# tobit: LSV,BQCV,DWV_A
a <- vglm(data=data,DWV_A~year.3,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(a)


## LSV
LSV_tobit <- vglm(data=data,LSV~year.3+Afr.Extent+Latitude+PRCP+TMAX+TMIN,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(LSV_tobit)
step4vglm(LSV_tobit,direction='both',steps=10000,scope=list(lower=LSV~year.3+Afr.Extent)) 

LSV_tobit_reduced <- vglm(data=data,LSV~year.3+Afr.Extent,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(LSV_tobit_reduced)

## BQCV
BQCV_tobit <- vglm(data=data,BQCV~year.3+Afr.Extent+Latitude+PRCP+TMAX+TMIN,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(BQCV_tobit)
step4vglm(BQCV_tobit,direction='both',steps=10000,scope=list(lower=BQCV~year.3+Afr.Extent)) 

BQCV_tobit_reduced <- vglm(data=data,BQCV~year.3+Afr.Extent+Latitude+TMIN+TMAX,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(BQCV_tobit_reduced)

## DWV
DWV_tobit <- vglm(data=data,DWV_A~year.3+Afr.Extent+Latitude+PRCP+TMAX+TMIN,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(DWV_tobit)
step4vglm(DWV_tobit,direction='both',steps=10000,scope=list(lower=DWV_A~year.3+Afr.Extent)) 

DWV_tobit_reduced <- vglm(data=data,DWV_A~year.3+Afr.Extent+TMAX,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(DWV_tobit_reduced)

#==========================================================================================================================================

# Diversity
library(vegan)
## Alpha diversity
## step1: calculate alpha diversity index: species richness, Shannon, Simpson
data.t$Species_richness <-apply(data.t[25:41]>0, 1, sum)
data.t$Shannon_vegan <- diversity(data.t[25:41], index = "shannon")

## step2: Mann-Whitney U test: compare non-africanized with africanized
wilcox.test(data.t$Shannon_vegan[data.t$V.arrival==0],data.t$Shannon_vegan[data.t$V.arrival==1])  #p-value:0.028

data_g <- data.t[which(data.t$Geocluster=='VER'), ]  #look into Geocluster
wilcox.test(data_g$Shannon_vegan[data_g$V.arrival==0],data_g$Shannon_vegan[data_g$V.arrival==1])

## step3: plot boxtplot
### under Varroa-free condition
data.t$Geocluster <- ordered(data.t$Geocluster, levels=c('WEL','BOR','TAM','VER'))
ggplot(data.t,aes(x=V.arrival,y=Shannon_vegan,fill=V.arrival))+
  geom_boxplot()+
  #facet_grid(rows=vars(Geocluster))+
  labs(x='',y='Shannon index',title='Shannon Diversity',subtitle='Honey bees (N=1034)')+
  scale_fill_manual(values = c('white','grey'))+
  theme_classic()+
  theme(legend.position="none",axis.text.x=element_text(size=12))+
  scale_x_discrete(labels=c('Varroa-free\n(n=936)','Varroa-infested\n(n=98)'))

### Varroa exposure time: 0,<1,<2,>2
ggplot(data.t,aes(x=year.3,y=Shannon_vegan,fill=year.3))+
  geom_boxplot()+
  #facet_grid(rows=vars(Geocluster))+
  labs(x='Varroa exposure time (years)',y='Shannon index',title='Shannon Diversity',subtitle='Honey bees (N=1034)')+
  scale_fill_manual(values = c('white','grey90','grey60','grey30'))+
  theme_classic()+
  theme(legend.position="none",axis.text.x=element_text(size=12))+
  scale_x_discrete(labels=c('0\n(n=936)','<1\n(n=56)','<2\n(n=18)','>2\n(n=24)'))

kruskal.test(Shannon_vegan~year.3, data=data.t)
pairwise.wilcox.test(data.t$Shannon_vegan, data.t$year.3, p.adjust.method = "fdr")

#====================================================================================================================
## mixed virus infection
### t-test
t.test(data.t$Species_reichness[data.t$V.arrival==0],data.t$Species_richness[data.t$V.arrival==1])

### chi-sequared 
dt <- table(data.t$Species_richness,data.t$V.arrival); dt
set.seed(123); chisq.test(dt,simulate.p.value=T) #p-value:0.0041

dt_freq <- as.data.frame(dt)
dt_freq$Freq[dt_freq$Var2=='0'] <- dt_freq$Freq[dt_freq$Var2=='0']/936*100
dt_freq$Freq[dt_freq$Var2=='1'] <- dt_freq$Freq[dt_freq$Var2=='1']/98*100
ggplot(data=dt_freq, aes(x=Var1, y=Freq, fill=Var2))+
  geom_bar(stat='identity',position=position_dodge(),color='black',width=0.7)+
  labs(x='The number of viral species',y='Percentage (%)',title='Mixed virus infections in the bee colonies',subtitle='Honey bees (N=1034)')+
  scale_fill_manual(labels=c('Varroa-free        (n=936)','Varroa-infested  (n=98)'),values=c('white','grey'))+
  theme_classic()+
  theme(legend.position = c(0.84,1.03),
        legend.background = element_rect(linetype="solid",colour ="darkgrey"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.text.x = element_text(size=11))

# Varroa exposure time: 0,<1,<2,>2
dt <- table(data.t$year.3, data.t$Species_richness); dt
set.seed(123); chisq.test(dt,simulate.p.value=T) #p-value:0.007

dt_freq <- as.data.frame(dt)
dt_freq$Freq[dt_freq$Var1=='0'] <- dt_freq$Freq[dt_freq$Var1=='0']/936*100
dt_freq$Freq[dt_freq$Var1=='<1'] <- dt_freq$Freq[dt_freq$Var1=='<1']/56*100
dt_freq$Freq[dt_freq$Var1=='<2'] <- dt_freq$Freq[dt_freq$Var1=='<2']/18*100
dt_freq$Freq[dt_freq$Var1=='>2'] <- dt_freq$Freq[dt_freq$Var1=='>2']/24*100

ggplot(data=dt_freq, aes(x=Var2, y=Freq, fill=Var1, group=Var1))+
  geom_bar(stat='identity',position=position_dodge(),color='black',width=0.7)+
  labs(x='The number of viral species',y='Percentage (%)',title='Mixed virus infections in the bee colonies',subtitle='Honey bees (N=1034)')+
  scale_fill_manual(labels=c('Varroa-free (n=936)','Varroa exposure <1 year (n=56)','Varroa exposure <2 years (n=18)','Varroa exposure >2 years (n=24)'),
                    values=c('white','grey90','grey60','grey30'))+
  theme_classic()+
  theme(legend.position = c(0.84,0.99),
        legend.background = element_rect(linetype="solid",colour ="darkgrey"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.text.x = element_text(size=11))

pairwiseNominalIndependence(dt,fisher=F,gtest=F,chisq=T,method='fdr')




