library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)

rm(list = ls())

data <- read.csv("new_data.csv")

data_NoV <- data[which(data$V.arrival==0), ] # n=926
data_NoV <- data_NoV[ ,c(18,19,21,22,23,26,28,34,35,1,5,10)]

data_NoV$Afr <- 0
data_NoV$Afr[which(data_NoV$Afr.Extent>10^-7)] <- 1
data_NoV$Afr <- as.factor(data_NoV$Afr)

name_list <- list("Acute bee paralysis virus","Apis mellifera filamentous virus","Black queen cell virus",
                  "Chronic bee paralysis virus","Deformed wing virus","Lake Sinai virus","Sacbrood virus")

data.1 <- data_NoV[which(data_NoV$Afr==0), ] # n=281
data.2 <- data_NoV[which(data_NoV$Afr==1), ] # n=645
for(i in 1:7){
  a <- data.1[ ,i]
  b <- data.2[ ,i]
  c <- wilcox.test(a,b, exact = F, correct = F)
  p <- paste(colnames(data.1)[i], c$p.value)
  print(p)
}

##boxplot-overview
for(i in 1:7){
  box <- ggplot(data_NoV,aes_string(x="Afr",y=colnames(data_NoV)[i],color='Afr'))+
    geom_boxplot()+
    scale_y_sqrt()+
    scale_color_jco()+
    labs(title=name_list[i],x="Africanization",y="Transcript per hundred")+
    scale_x_discrete(labels=c("No","Yes"))+
    stat_compare_means(label.x=1.325)+
    theme(legend.position = "none")
    
  print(box)
}

##boxplot-subgroup by geocluster
for(i in 1:7){
  box <- ggplot(data_NoV,aes_string(x="Afr",y=colnames(data_NoV)[6],color='Afr'))+
    geom_boxplot()+
    scale_y_sqrt()+
    scale_color_jco()+
    labs(title=name_list[6],x="Africanization",y="Transcript per hundred")+
    scale_x_discrete(labels=c("No","Yes"))+
    stat_compare_means(label = "p.format",label.x=1.4,label.y=9,size=3)+
    facet_grid(rows=vars(Geocluster))+
    theme(legend.position = "none")
  
  print(box)
}

#####################################################################

data_NoV.W <- data_NoV[which(data_NoV$Geocluster=="WEL"), ] # n=66

data.3 <- data_NoV.W[which(data_NoV.W$Afr==0), ] # n=21
data.4 <- data_NoV.W[which(data_NoV.W$Afr==1), ] # n=45
for(i in 1:7){
  a <- data.3[ ,i]
  b <- data.4[ ,i]
  c <- wilcox.test(a,b, exact = F, correct = F)
  p <- paste(colnames(data.1)[i], c$p.value)
  print(p)
}

############################################################################

data_NoV.B <- data_NoV[which(data_NoV$Geocluster=="BOR"), ] # n=53

data.5 <- data_NoV.B[which(data_NoV.B$Afr==0), ] # n=29
data.6 <- data_NoV.B[which(data_NoV.B$Afr==1), ] # n=24
for(i in 1:7){
  a <- data.5[ ,i]
  b <- data.6[ ,i]
  c <- wilcox.test(a,b, exact = F, correct = F)
  p <- paste(colnames(data.5)[i], c$p.value)
  print(p)
}

##################################################################################

data_NoV.T <- data_NoV[which(data_NoV$Geocluster=="TAM"), ] # n=724

data.7 <- data_NoV.T[which(data_NoV.T$Afr==0), ] # n=223
data.8 <- data_NoV.T[which(data_NoV.T$Afr==1), ] # n=501
for(i in 1:7){
  a <- data.7[ ,i]
  b <- data.8[ ,i]
  c <- wilcox.test(a,b, exact = F, correct = F)
  p <- paste(colnames(data.7)[i], c$p.value)
  print(p)
}

#########################################################################################

data_NoV.V <- data_NoV[which(data_NoV$Geocluster=="VER"), ] # n=83

data.9 <- data_NoV.V[which(data_NoV.V$Afr==0), ] # n=8
data.10 <- data_NoV.V[which(data_NoV.V$Afr==1), ] # n=75

for(i in 1:7){
  a <- data.9[ ,i]
  b <- data.10[ ,i]
  c <- wilcox.test(a,b, exact = F, correct = F)
  p <- paste(colnames(data.9)[i], c$p.value)
  print(p)
}

######################################################################################

# whether the distribution of data is normal
##shapiro.test()  
# F test: to determine whether 2 variance are the same
##kruskal.test(list(A, B))


