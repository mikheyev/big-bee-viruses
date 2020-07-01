library(ggplot2)
library(ggsci)
library(ggpubr)

rm(list = ls())

data <- read.csv("new_data.csv")
data$Geocluster <- ordered(data$Geocluster, level=c("WEL","BOR","TAM","VER"))
data$V.arrival <- as.factor(data$V.arrival)

data.a <- data[-c(which(data$Afr.Extent<=10^-7)), ]

name_full <- c('Acute bee paralysis virus','Apis mellifera filamentous virus',
               'Apis rhabdovirus 2','Black queen cell virus', 'Chronic bee paralysis virus',
               'Deformed wing virus','Israeli acute paralysis virus','Kashmir bee virus',
               'Lake Sinai virus','Slow bee paralysis virus','Sacbrood virus',
               'Varroa destructor virus 1','Varroa destructor virus 2','Varroa destructor virus 3',
               'Varroa orthomyxovirus-1','Varroa Tymo-like virus','Africanization')

data.1 <- data.a[which(data.a$V.arrival==0), ]
data.2 <- data.a[which(data.a$V.arrival==1), ]

for(i in 18:33){
  a <- data.1[ ,i]
  b <- data.2[ ,i]
  c <- wilcox.test(a,b, exact = F, correct = F)
  p <- paste(colnames(data.1)[i], c$p.value)
  print(p)
}

##boxplot-overview
for(i in 18:33){
  box <- ggplot(data.a,aes_string(x="V.arrival",y=colnames(data.a)[i],color='V.arrival'))+
    geom_boxplot()+
    scale_y_sqrt()+
    scale_color_jco()+
    labs(title=name_full[i-17],x="Varroa arrival",y="Transcript per hundred")+
    scale_x_discrete(labels=c("No","Yes"))+
    stat_compare_means(label.x=1.325)+
    theme(legend.position = "none")
  
  print(box)
}

##boxplot-subgroup by geocluster
for(i in 18:33){
  box <- ggplot(data.a,aes_string(x="V.arrival",y=colnames(data.a)[33],color='V.arrival'))+
    geom_boxplot()+
    scale_y_sqrt()+
    scale_color_jco()+
    labs(title=name_full[33-17],x="Varroa arrival",y="Transcript per hundred")+
    scale_x_discrete(labels=c("No","Yes"))+
    stat_compare_means(label = "p.format",label.x=1.4,label.y=4.5,size=3)+
    facet_grid(rows=vars(Geocluster))+
    theme(legend.position = "none")
  
  print(box)
}
#------------------------------------------------------------------------

data_1 <- data[which(data$V.arrival==0), ]
data_2 <- data[which(data$V.arrival==1), ]
wilcox.test(data_1$Afr.Extent,data_2$Afr.Extent, exact = F, correct = F)

#################################################################
data_WEL <- data.a[which(data.a$Geocluster=="WEL"), ]

data.3 <- data_WEL[which(data_WEL$V.arrival==0), ]
data.4 <- data_WEL[which(data_WEL$V.arrival==1), ]


for(i in 18:33){
  a <- data.3[ ,i]
  b <- data.4[ ,i]
  c <- wilcox.test(a,b, exact = F, correct = F)
  p <- paste(colnames(data.3)[i], c$p.value)
  print(p)
}

#-------------------------------------------------------------

data_W <- data[which(data$Geocluster=="WEL"), ]

data_3 <- data_W[which(data_W$V.arrival==0), ]
data_4 <- data_W[which(data_W$V.arrival==1), ]

wilcox.test(data_3$Afr.Extent,data_4$Afr.Extent, exact = F, correct = F)

#############################################################
data_BOR <- data.a[which(data.a$Geocluster=="BOR"), ]

data.5 <- data_BOR[which(data_BOR$V.arrival==0), ]
data.6 <- data_BOR[which(data_BOR$V.arrival==1), ]

for(i in 18:33){
  a <- data.5[ ,i]
  b <- data.6[ ,i]
  c <- wilcox.test(a,b, exact = F, correct = F)
  p <- paste(colnames(data.5)[i], c$p.value)
  print(p)
}
#------------------------------------------------------------

data_B <- data[which(data$Geocluster=="BOR"), ]

data_5 <- data_B[which(data_B$V.arrival==0), ]
data_6 <- data_B[which(data_B$V.arrival==1), ]

wilcox.test(data_5$Afr.Extent,data_6$Afr.Extent, exact = F, correct = F)


#############################################################
data_TAM <- data.a[which(data.a$Geocluster=="TAM"), ]

data.7 <- data_TAM[which(data_TAM$V.arrival==0), ]
data.8 <- data_TAM[which(data_TAM$V.arrival==1), ]

for(i in 18:33){
  a <- data.7[ ,i]
  b <- data.8[ ,i]
  c <- wilcox.test(a,b, exact = F, correct = F)
  p <- paste(colnames(data.7)[i], c$p.value)
  print(p)
}
#-----------------------------------------------------------

data_T <- data[which(data$Geocluster=="TAM"), ]

data_7 <- data_T[which(data_T$V.arrival==0), ]
data_8 <- data_T[which(data_T$V.arrival==1), ]
wilcox.test(data_7$Afr.Extent,data_8$Afr.Extent, exact = F, correct = F)

############################################################
data_VER <- data.a[which(data$Geocluster=="VER"), ]

data.9 <- data_VER[which(data_VER$V.arrival==0), ]
data.10 <- data_VER[which(data_VER$V.arrival==1), ] # data in VER does not have any samples with Varroa

################################################################

# whether the distribution of data is normal
##shapiro.test()  
# F test: to determine whether 2 variance are the same
##kruskal.test(list(A, B))





