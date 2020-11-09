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


data_NoV <- data[which(data$V.arrival==0), ]

data_NoV.p <- data_NoV
detected <-which(colSums(data_NoV.p[,25:41])==0)+24
data_NoV.p <- data_NoV.p[,-detected]
data_NoV.p$Africanization <- as.factor(data_NoV.p$Africanization)

data_NoV.t <- data_NoV[,-detected]
data_NoV.t$Africanization <- as.factor(data_NoV.t$Africanization)

library(ggcorrplot)
corr <- round(cor(data_NoV.p[ ,25:39],method='spearman'), 2)
ggcorrplot(corr, hc.order = TRUE, type = "lower",
           outline.col = "white",
           ggtheme = ggplot2::theme_classic,
           colors = c("#6D9EC1", "white", "#E46726"),
           lab=T)

## viral prevalence: chi-squared test
## convert continuous variable (viral tpm) into categorical variable
for(i in 25:38){
  for(l in 1:936){
    if(data_NoV.p[l,i]>0){data_NoV.p[l,i] <- 1}
    else{data_NoV.p[l,i] <- 0}
  }
}
### statistical significance
b <- c()
for (i in c(26,35,28,37,25,30)) {
  set.seed(123)
  dt <- table(data_NoV.p$Africanization,data_NoV.p[ ,i])
  a <- chisq.test(dt,simulate.p.value = T)
  p <- paste(names(data_NoV.p)[i],'x-squared:',a$statistic,'p-value:',a$p.value)
  print(dt)
  print(p)
  b <- append(b,a$p.value)
}
# adjusted p-value: false discover rate -> significant: BQCV,SBV
p.adjust(b,method='fdr')

## viral prevalence: bar chart-> x=virus, y=prevalence
## step1: calculate the viral prevalence
data_prevalence <- data_NoV.p %>% group_by(Africanization) %>% 
  summarise(ABPV=length(which(ABPV==1))/length(ID),Amfv=length(which(Amfv==1))/length(ID),
            Arv2=length(which(Arv2==1))/length(ID),BQCV=length(which(BQCV==1))/length(ID),
            CBPV=length(which(CBPV==1))/length(ID),DWV_A=length(which(DWV_A==1))/length(ID),
            DWV_B=length(which(DWV_B==1))/length(ID),DWV_C=length(which(DWV_C==1))/length(ID),
            IAPV=length(which(IAPV==1))/length(ID),KBV=length(which(KBV==1))/length(ID),
            LSV=length(which(LSV==1))/length(ID),SBPV=length(which(SBPV==1))/length(ID),
            SBV=length(which(SBV==1))/length(ID),VTV=length(which(VTV==1))/length(ID))

data_prevalence_all <- data_NoV.p %>%
  summarise(ABPV=length(which(ABPV==1))/length(ID),Amfv=length(which(Amfv==1))/length(ID),
            Arv2=length(which(Arv2==1))/length(ID),BQCV=length(which(BQCV==1))/length(ID),
            CBPV=length(which(CBPV==1))/length(ID),DWV_A=length(which(DWV_A==1))/length(ID),
            DWV_B=length(which(DWV_B==1))/length(ID),DWV_C=length(which(DWV_C==1))/length(ID),
            IAPV=length(which(IAPV==1))/length(ID),KBV=length(which(KBV==1))/length(ID),
            LSV=length(which(LSV==1))/length(ID),SBPV=length(which(SBPV==1))/length(ID),
            SBV=length(which(SBV==1))/length(ID),VTV=length(which(VTV==1))/length(ID))
data_prevalence_all <- data_prevalence_all[ ,names(sort(data_prevalence_all[1,],decreasing=T))]

## step2: ggplot plot the barcharts
data_prevalence_melt <- melt(data_prevalence[,c('Africanization',names(sort(data_prevalence_all[1,],decreasing=T)))], id=c('Africanization'))
ggplot(data_prevalence_melt,aes(x=variable,y=value*100,fill=as.factor(Africanization)))+
  geom_bar(stat='identity',width=0.65, position=position_dodge(),color='black')+
  labs(x='',y='Pathogen prevalence (% positive individual bees)',title='Varroa-free honeybees (N=936)')+
  scale_fill_manual(labels=c('Non-Africanized (n=289)','Africanized (n=647)'),values=c('white','grey'))+
  ylim(c(0,100))+
  theme_classic()+
  theme(legend.position = c(0.905,0.99),
        legend.background = element_rect(linetype="solid",colour ="darkgrey"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.text.x = element_text(size=11))+
  scale_x_discrete(labels=c('Amfv\n91.2','LSV\n76.2','BQCV\n56.0','SBV\n40.5','ABPV\n16.5',
                            'DWV-A\n10.3','CBPV\n5.88','KBV\n1.28','Arv2\n0.85','IAPV\n0.85',
                            'DWV-B\n0.53','SBPV\n0.32','VTV\n0.21','DWV-C\n0.11'))

ggplot(data_prevalence_melt[1:12,],aes(x=variable,y=value*100,fill=as.factor(Africanization)))+
  geom_bar(stat='identity',width=0.65, position=position_dodge(),color='black')+
  labs(x='',y='Pathogen prevalence (% positive individual bees)',title='Varroa-free honey bees (N=936)')+
  scale_fill_manual(labels=c('European (n=289)','Africanized (n=647)'),values=c('white','grey'))+
  ylim(c(0,100))+
  geom_text(aes(label=round(value*100,1)),position = position_dodge(0.7),
            vjust = -0.3, size = 3.5)+
  theme_classic()+
  theme(legend.position = c(0.89,0.99),
        legend.background = element_rect(linetype="solid",colour ="darkgrey"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.text.x = element_text(size=11))+
  scale_x_discrete(labels=c('Amfv\n91.2','LSV\n76.2','BQCV\n56.0','SBV\n40.5',
                            'ABPV\n16.5','DWV-A\n10.3'))


Amfv<- glmer(Amfv~Afr.Extent+(1|Geocluster)+(1|ID), 
            family=binomial(link = 'logit'), data=data_NoV.p)
summary(Amfv)

#-------------------------------------------------------------------------------------

## GLMM: viral prevalence, logistics regression: Amfv,LSV,BQCV,SBV,ABPV,DWV-A
## Africanization: LSV,BQCV,SBV
a<- glmer(DWV_A~Africanization+(1|Geocluster)+(1|ID), 
             family=binomial(link = 'logit'), data=data_NoV.p)
summary(a)

### LSV
LSV <- glmer(LSV~Africanization+Latitude+PRCP+TMAX+TMIN+(1|Geocluster)+(1|ID), 
            family=binomial(link = 'logit'), data=data_NoV.p)
summary(LSV) #AIC=1030.4

LSV_reduced <- glmer(LSV~Africanization+(1|Geocluster)+(1|ID), 
                    family=binomial(link = 'logit'), data=data_NoV.p)
summary(LSV_reduced) #AIC=1028.6

anova(LSV_reduced, LSV, test='F') #p=0.1792

### BQCV
BQCV <- glmer(BQCV~Africanization+Latitude+PRCP+TMAX+TMIN+(1|Geocluster)+(1|ID), 
            family=binomial(link = 'logit'), data=data_NoV.p)
summary(BQCV) #AIC=1261.3

BQCV_reduced <- glmer(BQCV~Africanization+Latitude+TMAX+(1|Geocluster)+(1|ID), 
              family=binomial(link = 'logit'), data=data_NoV.p)
summary(BQCV_reduced) #AIC=1258.1

anova(BQCV_reduced, BQCV, test='F') #p=0.6847

### SBV
SBV <- glmer(SBV~Africanization+Latitude+PRCP+TMAX+TMIN+(1|Geocluster)+(1|ID), 
             family=binomial(link = 'logit'), data=data_NoV.p)
summary(SBV) #AIC=1217.4

SBV_reduced <- glmer(SBV~Africanization+PRCP+TMAX+TMIN+(1|Geocluster)+(1|ID), 
             family=binomial(link = 'logit'), data=data_NoV.p)
summary(SBV_reduced) #AIC=1217.0

anova(SBV_reduced, SBV, test='F') #p=0.2129

## Africanized extent: Amfv,LSV,BQCV,SBV,ABPV
b <- glmer(ABPV~Afr.Extent+(1|Geocluster)+(1|ID), 
           family=binomial(link = 'logit'), data=data_NoV.p)
summary(b)

### Amfv
Amfv2 <- glmer(Amfv~Afr.Extent+Latitude+PRCP+TMAX+TMIN+(1|Geocluster)+(1|ID), 
             family=binomial(link = 'logit'), data=data_NoV.p)
summary(Amfv) #AIC=543.2

Amfv_reduced2 <- glmer(Amfv~Afr.Extent+Latitude+TMIN+(1|Geocluster)+(1|ID), 
                      family=binomial(link = 'logit'), data=data_NoV.p)
summary(Amfv_reduced2) #AIC=540.7

anova(Amfv_reduced2, Amfv2, test='F') #p=0.4689

### LSV
LSV2 <- glmer(LSV~Afr.Extent+Latitude+PRCP+TMAX+TMIN+(1|Geocluster)+(1|ID), 
              family=binomial(link = 'logit'), data=data_NoV.p)
summary(LSV2) #AIC=1032.1

LSV_reduced2 <- glmer(LSV~Afr.Extent+(1|Geocluster)+(1|ID), 
                      family=binomial(link = 'logit'), data=data_NoV.p)
summary(LSV_reduced2) #AIC=1030.1

anova(LSV_reduced2, LSV2, test='F') #p=0.2

## BQCV
BQCV2 <- glmer(BQCV~Afr.Extent+Latitude+PRCP+TMAX+TMIN+(1|Geocluster)+(1|ID), 
             family=binomial(link = 'logit'), data=data_NoV.p)
summary(BQCV2) #AIC=1261.2

BQCV_reduced2 <- glmer(BQCV~Afr.Extent+Latitude+TMAX+(1|Geocluster)+(1|ID), 
              family=binomial(link = 'logit'), data=data_NoV.p)
summary(BQCV_reduced2) #AIC=1257.8

anova(BQCV_reduced2, BQCV2, test='F') #p=0.7352

## SBV
SBV2 <- glmer(SBV~Afr.Extent+Latitude+PRCP+TMAX+TMIN+(1|Geocluster)+(1|ID), 
              family=binomial(link = 'logit'), data=data_NoV.p)
summary(SBV2) #AIC=1209.5

SBV_reduced2 <- glmer(SBV~Afr.Extent+PRCP+TMAX+TMIN+(1|Geocluster)+(1|ID), 
             family=binomial(link = 'logit'), data=data_NoV.p)
summary(SBV_reduced2) #AIC=1209.5

anova(SBV_reduced2, SBV2, test='F') #p=0.1677

## ABPV
ABPV2 <- glmer(ABPV~Afr.Extent+Latitude+PRCP+TMAX+TMIN+(1|Geocluster)+(1|ID), 
             family=binomial(link = 'logit'), data=data_NoV.p)
summary(ABPV2)

ABPV_reduced2 <- glmer(ABPV~Afr.Extent+(1|Geocluster)+(1|ID), 
              family=binomial(link = 'logit'), data=data_NoV.p)
summary(ABPV_reduced2)

anova(ABPV_reduced2, ABPV2, test='F') #p=<.001 reduced model is not better than full

#================================================================================================

## viral tpm: permutation median test (are there any tpm differences between non-africanized and africanized)
a <- c()
for (i in c(26,35,28,37,25,30)) {
  set.seed(123)
  n <- length(data_NoV.t$Africanization)
  P <- 100000
  variable <- data_NoV.t[ ,i]
  test.stat.1 <- abs(median(data_NoV.t[data_NoV.t$Africanization==0,i])-median(data_NoV.t[data_NoV.t$Africanization==1,i]))
  
  PerSamples <- matrix(0, nrow=n, ncol=P)
  for(j in 1:P){
    PerSamples[ ,j] <- sample(variable, size=n, replace=F)
  }
  
  Perm.test.stat.1 <- rep(0,P)
  for(k in 1:P){
    Perm.test.stat.1[k] <- abs(median(PerSamples[data_NoV.t$Africanization==0,k])-
                                 median(PerSamples[data_NoV.t$Africanization==1,k]))
  }
  print(paste(names(data_NoV.t)[i],"Median:",mean(Perm.test.stat.1>=test.stat.1)))
  a <- append(a,mean(Perm.test.stat.1>=test.stat.1))
}

# adjusted p-value: false discover rate -> significant: SBV
p.adjust(a,method='fdr')

## viral tpm: boxplot-> x=viruses, y=tpm
data_melt <- melt(data_NoV.t[,c(26,35,28,37,25,30,40)], id=c('Africanization'))

ggplot(data_melt,aes(x=variable,y=value,fill=as.factor(Africanization)))+
  geom_boxplot(outlier.shape = NA)+
  labs(x='',y='Transcript per million',title='Varroa-free honey bees (N=936)')+
  scale_y_sqrt()+
  scale_fill_manual(labels=c('European (n=289)','Africanized (n=647)'),values=c('white','grey'))+
  theme_classic()+
  theme(legend.position = c(0.89,0.99),
        legend.background = element_rect(linetype="solid",colour ="darkgrey"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size=11))+
  scale_x_discrete(labels=c('Amfv\n3914.38','LSV\n304226','BQCV\n10802.65','SBV\n0','ABPV\n0','DWV-A\n0'))

data_median <- data_NoV.t[,25:38] %>% summarise_all(median)
data_mean <- data_NoV.t[,25:38] %>% summarise_all(mean)

names(sort(data_mean[1,],decreasing=T))
      
#--------------------------------------------------------------------------------------------------------

# tobit model: LSV,BQCV,SBV,Amfv
## Africanization: Amfv,SBV, ABPV,DWV-A
a <- vglm(data=data_NoV.t,DWV_A~Africanization,family=tobit(Lowe=0,Upper=10^6,type.fitted = 'censored'))
summary(a)

## Africanized extent: Amfv,SBV,ABPV
b <- vglm(data=data_NoV.t,Amfv~Afr.Extent,family=tobit(Lowe=0,Upper=10^6,type.fitted = 'censored'))
summary(b)

library(VGAM)
## ABPV
ABPV_tobit <- vglm(data=data_NoV.t,ABPV~Africanization+Latitude+PRCP+TMAX+TMIN,family=tobit(Lowe=0,Upper=10^6,type.fitted = 'censored'))
summary(ABPV_tobit)
step4vglm(ABPV_tobit,direction='both',steps=10000,scope=list(lower=ABPV~Africanization))

ABPV_tobit_reduced <- vglm(data=data_NoV.t,ABPV~Afr.Extent,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(ABPV_tobit_reduced)

## LSV
LSV_tobit <- vglm(data=data_NoV.t,LSV~Africanization+Latitude+PRCP+TMAX+TMIN,family=tobit(Lowe=0,Upper=10^6,type.fitted = 'censored'))
summary(LSV_tobit)
step4vglm(LSV_tobit,direction='both',steps=10000,scope=list(lower=LSV~Africanization))

LSV_tobit_reduced <- vglm(data=data_NoV.t,LSV~Africanization,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(LSV_tobit_reduced)

## BQCV
BQCV_tobit <- vglm(data=data_NoV.t,BQCV~Africanization+Latitude+PRCP+TMAX+TMIN,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(BQCV_tobit)
step4vglm(BQCV_tobit,direction='both',steps=10000,scope=list(lower=LSV~Africanization))

BQCV_tobit_reduced <- vglm(data=data_NoV.t,BQCV~Africanization+Latitude+TMAX+TMIN,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(BQCV_tobit_reduced)

## SBV
SBV_tobit1 <- vglm(data=data_NoV.t,SBV~Africanization+Latitude+PRCP+TMAX+TMIN,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(SBV_tobit1)
step4vglm(SBV_tobit1,direction='both',steps=10000)

SBV_tobit_reduced1 <- vglm(data=data_NoV.t,SBV~Africanization+Latitude+PRCP+TMAX,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(SBV_tobit_reduced1)    

SBV_tobit2 <- vglm(data=data_NoV.t,SBV~Afr.Extent+Latitude+PRCP+TMAX+TMIN,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(SBV_tobit2)
step4vglm(SBV_tobit2,direction='both',steps=10000)

SBV_tobit_reduced2 <- vglm(data=data_NoV.t,SBV~Afr.Extent+Latitude+PRCP+TMAX,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(SBV_tobit_reduced2)    

## Amfv
Amfv_tobit <- vglm(data=data_NoV.t,Amfv~Africanization+Latitude+PRCP+TMAX+TMIN,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(Amfv_tobit)
step4vglm(Amfv_tobit,direction='both',steps=10000) ## did not remain any variables

Amfv_tobit_reduced1 <- vglm(data=data_NoV.t,Amfv~Africanization+Latitude+TMAX,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(Amfv_tobit_reduced1)


## DWV-A
DWV_tobit <- vglm(data=data_NoV.t,DWV_A~Africanization+Latitude+PRCP+TMAX+TMIN,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(DWV_tobit) #AIC=3290.55
step4vglm(DWV_tobit,direction='both',steps=10000,scope=list(lower=DWV_A~Africanization)) ## did not remain any variables

DWV_tobit_reduced <- vglm(data=data_NoV.t,DWV_A~Africanization+TMIN+TMAX,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(DWV_tobit_reduced) #AIC=3287.57


## CBPV
CBPV_tobit <- vglm(data=data_NoV.t,CBPV~Africanization+Latitude+PRCP+TMAX+TMIN,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(CBPV_tobit) #AIC=3290.55
step4vglm(CBPV_tobit,direction='both',steps=10000,scope=list(lower=CBPV~Africanization)) ## did not remain any variables

CBPV_tobit_reduced <- vglm(data=data_NoV.t,CBPV~Africanization+Latitude+PRCP,family=tobit(Lower=0,Upper=10^6,type.fitted = 'censored'))
summary(CBPV_tobit_reduced) 

#==========================================================================================================================================

# Diversity
library(vegan)
data_NoV.t <- data_NoV[,-detected]
data_NoV.t$Africanization <- as.factor(data_NoV.t$Africanization)
## Alpha diversity
## step1: calculate alpha diversity index: species richness, Shannon, Simpson
data_NoV.t$Species_richness <-apply(data_NoV.t[25:38]>0, 1, sum)
data_NoV.t$Shannon_vegan <- diversity(data_NoV.t[25:38], index = "shannon")

## step2: Mann-Whitney U test: compare non-africanized with africanized
wilcox.test(data_NoV.t$Shannon_vegan[data_NoV.t$Africanization==0],data_NoV.t$Shannon_vegan[data_NoV.t$Africanization==1])  #p-value:0.046

data_NoV_g <- data_NoV.t[which(data_NoV$Geocluster=='VER'), ]  #look into Geocluster
wilcox.test(data_NoV_g$Shannon_vegan[data_NoV_g$Africanized==0],data_NoV_g$Shannon_vegan[data_NoV_g$Africanized==1])

## step3: plot boxtplot
### under Varroa-free condition
data_NoV.t$Geocluster <- ordered(data_NoV.t$Geocluster, levels=c('WEL','BOR','TAM','VER'))
ggplot(data_NoV.t,aes(x=Africanization,y=Shannon_vegan,fill=Africanization))+
  geom_boxplot()+
  #facet_grid(rows=vars(Geocluster))+
  labs(x='',y='Shannon index',title='Shannon Diversity',subtitle='Varroa-free honeybees (N=936)')+
  scale_fill_manual(values = c('white','grey'))+
  theme_classic()+
  theme(legend.position="none",axis.text.x=element_text(size=12))+
  scale_x_discrete(labels=c('Non-Africanized\n(n=289)','Africanized\n(n=647)'))

#---------------------------------------------------------------------------------------------------------------------------
# Species richness

## mixed virus infections in the bee colonies
## t-test and boxplot
t.test(data_NoV.t$Species_richness[which(data_NoV.t$Africanization==0)],data_NoV.t$Species_richness[which(data_NoV.t$Africanization==1)])

ggplot(data=data_NoV.t, aes(x=Africanization,y=Species_richness))+
  geom_boxplot()

## under Varroa-free condition: chi-squared test, plot
dt <- table(data_NoV.t$Species_richness,data_NoV.t$Africanization); dt
set.seed(123); chisq.test(dt,simulate.p.value=T) #p-value:0.0005

dt_freq <- as.data.frame(dt)
dt_freq$Freq[dt_freq$Var2=='0'] <- dt_freq$Freq[dt_freq$Var2=='0']/289*100
dt_freq$Freq[dt_freq$Var2=='1'] <- dt_freq$Freq[dt_freq$Var2=='1']/647*100
ggplot(data=dt_freq,aes(x=Var1,y=Freq,fill=Var2))+
  geom_bar(stat='identity',position=position_dodge(),color='black',width=0.7)+
  labs(x='The number of viral species',y='Percentage (%)',title='Mixed virus infections in the honey bees',subtitle='Varroa-free honeybees (N=936)')+
  scale_fill_manual(labels=c('European  (n=289)','Africanized  (n=647)'),values=c('white','grey'))+
  theme_classic()+
  theme(legend.position = c(0.87,1.036),
        legend.background = element_rect(linetype="solid",colour ="darkgrey"),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.text.x = element_text(size=11))










