########################################
##                                    ##   
## Author: Sherry Lin                 ##
## Student ID: u6890762               ##
## Project: evolution of bee viruses  ##
## Purpose: data arrangement          ##
## Version: 14                        ##
## Date: 22/06/20                     ##
##                                    ##
########################################

library(data.table) #fread(): read the data
library(tools) #file_path_sans_ext(): without extension

rm(list = ls())

# before import"seq_info.csv, change 'SAMPLE_ID':'4183-label-inside-334183' to '334183' 
# and '1206' Varroa: YES to NA (because they could be tracheal mites, for instance)
# if they were really varroa, my feeling is that we'd get more in the next year or two. instead there is a gap of several years

# import the data: "seq_info.csv" which contains the basic information about the sequenced samples
# change the first column's name
data <- read.csv("seq_info.csv")      
names(data)[1] <- "ID"                                           

# getting files' name without extension (the files' names are the sample's ID) 
unzip.filename <- list.files("00_known_viruses")                                     
unzip.filename.new <- file_path_sans_ext(list.files("00_known_viruses")) 

# change the the working directory into a file named 'unzip'
# merge each sequenced results into a dataframe
setwd("00_known_viruses")                                                   
TPM <- data.frame()                                             
for (i in 1:1043){                                               
  a <- fread(unzip.filename[i])                                  
  a <- a[ ,c(1, 5)]                                              
  a <- as.data.frame(t(a), stringsAsFactors = F)                 
  names(a) <- a[1, ]                                             
  a <- a[-1, ]   
  a$ID <- unzip.filename.new[i]                              
  TPM <- rbind(TPM, a)                                         
}                                                                

for (i in 1:28){                                                 
  TPM[[i]] <- as.numeric(TPM[[i]])                             
}                                                                

remove_list <- c()
for(i in 1:28){
  a <- sum(TPM[i],na.rm = T)
  if (a==0){
    remove_list <- append(remove_list,i)
  }
}

TPM <- TPM[ ,-remove_list]

setwd("D:/Sherry/文件/RStudio/Bee's viruses evolution")
# when merge(full join)'data' & 'TPM', I found :
# 'TPM' includes 9 negative controls and 334713, but 'data' does not
# 'data' includes 7227, but 'TPM' does not
# thus, used the metadata to get the information of 334713 and insert it into "seq_info.csv" (in excel)
# then, run the codes again
# removed the negative controls and 7227 by 'merge' (inner join)
data_merge <- merge(data,TPM,by="ID")

# used the website: "ncbi.nlm.nih.gov/search/" to get virus description by accession
# for example: the IDs beginning with 'MK0324' are the same virus, so should sum the TPMs for them
# also, change the columns' name from accession number to viruses' name
data_merge$DWV <- data_merge$NC_004830.2+data_merge$`ENA|CEND01000001|CEND01000001.1` #DWV
data_merge <- data_merge[ ,-c(18,19)]

names(data_merge)[18:32] <- c("IAPV","KBV","CBPV","ABPV","SBPV", 
                              "LSV","BQCV","SV", "VDV1","VDV2",  
                              "VDV3","AR","AM","VTV","VO1")                                            

data_merge[ ,18:33] <- data_merge[ ,18:33]/10^4 #scale the TPM -> Transcript per hundred(0~100)
data_merge <- data_merge[ ,c(1:17,21,30,29,24,20,33,18,19,23,22,25:28,32,31)]

# import data: "k2.txt" is the extent of Africanization of each bee
# then, merge this data into seq_tpm_data by ID
africanized <- read.table("00_k2.txt")                                                                                      
names(africanized)[1] <- "ID"                                                                                            
africanized$ID <- as.character(africanized$ID)                                                                           
africanized$ID[759] <- "334183"                                                                                          

data_merge <- merge(data_merge, africanized, by = "ID")  # merge by ID, africanized data does not have "7227" 
data_merge <- data_merge[ ,-34]                           
colnames(data_merge)[34] <- "Afr.Extent"
data_merge[34] <- data_merge[34]*100

data_merge$V.arrival <- 0
data_merge$V.arrival[which(data_merge$Geocluster=="WEL"&data_merge$Year>=1994)] <- 1
data_merge$V.arrival[which(data_merge$Geocluster=="BOR"&data_merge$Year>=1992)] <- 1
data_merge$V.arrival[which(data_merge$Geocluster=="TAM"&data_merge$Year>=1995)] <- 1

write.csv(data_merge, file = "new_data.csv", row.names = F)

#################################################################################################################
## average of TPM
library(dplyr)
rm(list = ls())
data <- read.csv("new_data.csv", header = TRUE, sep = ",")
a <- group_by(data, Year) %>% summarise(ABPV=mean(ABPV, na.rm=T ), AM=mean(AM, na.rm=T), 
                                        AR=mean(AR, na.rm=T), BQCV=mean(BQCV, na.rm=T),
                                        CBPV=mean(CBPV,na.rm=T),DWV=mean(DWV, na.rm=T), 
                                        IAPV=mean(IAPV, na.rm=T),KBV=mean(KBV, na.rm=T),
                                        LSV=mean(LSV, na.rm=T),SBPV=mean(SBPV, na.rm=T),
                                        SV=mean(SV, na.rm=T),VDV1=mean(VDV1, na.rm=T),  
                                        VDV2=mean(VDV2, na.rm=T),VDV3=mean(VDV3, na.rm=T),    
                                        VO1=mean(VO1, na.rm=T),VTV=mean(VTV, na.rm=T),
                                        Afr.Extent = mean(Afr.Extent, na.rm = T))

b <- group_by(data, Year, Geocluster) %>% summarise(ABPV=mean(ABPV, na.rm=T ), AM=mean(AM, na.rm=T), 
                                                    AR=mean(AR, na.rm=T), BQCV=mean(BQCV, na.rm=T),
                                                    CBPV=mean(CBPV,na.rm=T),DWV=mean(DWV, na.rm=T), 
                                                    IAPV=mean(IAPV, na.rm=T),KBV=mean(KBV, na.rm=T),
                                                    LSV=mean(LSV, na.rm=T),SBPV=mean(SBPV, na.rm=T),
                                                    SV=mean(SV, na.rm=T),VDV1=mean(VDV1, na.rm=T),  
                                                    VDV2=mean(VDV2, na.rm=T),VDV3=mean(VDV3, na.rm=T),    
                                                    VO1=mean(VO1, na.rm=T),VTV=mean(VTV, na.rm=T),
                                                    Afr.Extent = mean(Afr.Extent, na.rm = T))

write.csv(a, file ="new_data_mean_Year.csv", row.names = F)
write.csv(b, file ="new_data_mean_YearGeo.csv", row.names = F)

###############################################################################################################
## The prevalence
library(dplyr)
rm(list = ls())

data <- read.csv("new_data.csv")

for(i in 18:33){
  for(l in 1:1034){
    if(is.na(data[l,i])){data[l,i] <- NA}
    else if(data[l,i]>0){data[l,i] <- 1}
    else{data[l,i] <- 0}
  }
}

data$Afr.Extent[which(data$Afr.Extent<=10^-7)] <- 0
data$Afr.Extent[which(data$Afr.Extent>10^-7)] <- 1

data$Varroa[which(data$Varroa=="YES")] <- 1
data$Varroa[which(data$Varroa=="NO")] <- 0
data$Varroa[which(is.na(data$Varroa))] <- 0

cols <- colnames(data[ ,c(14,18:35)])
data[cols] <- lapply(data[cols], factor)

data_pre <- group_by(data, Geocluster, Year) %>% summarise(ABPV=length(which(ABPV=="1"))/length(ID),
                                                           AM=length(which(AM=="1"))/length(ID),
                                                           AR=length(which(AR=="1"))/length(ID),
                                                           BQCV=length(which(BQCV=="1"))/length(ID),
                                                           CBPV=length(which(CBPV=="1"))/length(ID),
                                                           DWV=length(which(DWV=="1"))/length(ID),
                                                           IAPV=length(which(IAPV=="1"))/length(ID),
                                                           KBV=length(which(KBV=="1"))/length(ID),
                                                           LSV=length(which(LSV=="1"))/length(ID),
                                                           SBPV=length(which(SBPV=="1"))/length(ID),
                                                           SV=length(which(SV=="1"))/length(ID),
                                                           VDV1=length(which(VDV1=="1"))/length(ID),
                                                           VDV2=length(which(VDV2=="1"))/length(ID),
                                                           VDV3=length(which(VDV3=="1"))/length(ID),
                                                           VO1=length(which(VO1=="1"))/length(ID),
                                                           VTV=length(which(VTV=="1"))/length(ID),
                                                           Varroa=length(which(Varroa=="1"))/length(ID),
                                                           African=length(which(Afr.Extent=="1"))/length(ID))
data_pre[3:20] <- data_pre[3:20]*100

data_pre_all <- group_by(data,Year) %>% summarise(ABPV=length(which(ABPV=="1"))/length(ID),
                                                  AM=length(which(AM=="1"))/length(ID),
                                                  AR=length(which(AR=="1"))/length(ID),
                                                  BQCV=length(which(BQCV=="1"))/length(ID),
                                                  CBPV=length(which(CBPV=="1"))/length(ID),
                                                  DWV=length(which(DWV=="1"))/length(ID),
                                                  IAPV=length(which(IAPV=="1"))/length(ID),
                                                  KBV=length(which(KBV=="1"))/length(ID),
                                                  LSV=length(which(LSV=="1"))/length(ID),
                                                  SBPV=length(which(SBPV=="1"))/length(ID),
                                                  SV=length(which(SV=="1"))/length(ID),
                                                  VDV1=length(which(VDV1=="1"))/length(ID),
                                                  VDV2=length(which(VDV2=="1"))/length(ID),
                                                  VDV3=length(which(VDV3=="1"))/length(ID),
                                                  VO1=length(which(VO1=="1"))/length(ID),
                                                  VTV=length(which(VTV=="1"))/length(ID),
                                                  Varroa=length(which(Varroa=="1"))/length(ID),
                                                  African=length(which(Afr.Extent=="1"))/length(ID))
data_pre_all[2:19] <- data_pre_all[2:19]*100

write.csv(data_pre, file ="new_data_prevalence_YearGeo.csv", row.names = F)
write.csv(data_pre_all, file="new_data_prevalence_Year.csv",row.names=F)

##############################################################################################################

library(forcats)

rm(list = ls())

data <- read.csv("new_data.csv", header=T, sep=",")

summary(data) #Factors: Varroa, State, Country
summary(as.factor(data$Varroa)) #Varroa: NO(31),YES(22),NA(981)
summary(as.factor(data$State)) #States: TAMAULIPAS,TEXAS,VERACRUZ,N.VERACRUZ
summary(as.factor(data$Country)) #Country: MEXICO,USA

data$State[which(data$State=="N.VERACRUZ")] <- "VERACRUZ"
data$Varroa <-fct_explicit_na(data$Varroa, na_level = "Unknown")

data <- data[ ,c(1,3:12,14)]

data_all <- data %>%
  group_by(Year, Varroa) %>%
  summarise(count=length(ID)) %>%
  group_by(Year) %>%
  mutate(total=sum(count), percent=round(100*count/total,2))

data_c <- data %>% 
  group_by(Year,Country,Varroa) %>% 
  summarise(count=length(ID)) %>% 
  group_by(Year,Country) %>% 
  mutate(total=sum(count), percent=round(100*count/total,2))

data_s <- data %>% 
  group_by(Year,State,Varroa) %>% 
  summarise(count=length(ID)) %>% 
  group_by(Year,State) %>% 
  mutate(total=sum(count), percent=round(100*count/total,2))

data_g <- data %>% 
  group_by(Year, Geocluster, Varroa) %>% 
  summarise(count=length(ID)) %>% 
  group_by(Year, Geocluster) %>% 
  mutate(total=sum(count), percent=round(100*count/total,2))

write.csv(data_all, file="new_data_Varroa_Year.csv", row.names=F)
write.csv(data_c, file ="new_data_Varroa_YearCountry.csv", row.names = F)
write.csv(data_s, file ="new_data_Varroa_YearState.csv", row.names = F)
write.csv(data_g, file ="new_data_Varroa_YearGeo.csv", row.names = F)

##################################################################################################

library(forcats)

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


data_c <- data %>% 
  group_by(Year,Country,Varroa) %>% 
  summarise(count=length(ID)) %>% 
  group_by(Year,Country) %>% 
  mutate(total=sum(count), percent=round(100*count/total,2))

data_s <- data %>% 
  group_by(Year,State,Varroa) %>% 
  summarise(count=length(ID)) %>% 
  group_by(Year,State) %>% 
  mutate(total=sum(count), percent=round(100*count/total,2))

data_g <- data %>% 
  group_by(Year, Geocluster, Varroa) %>% 
  summarise(count=length(ID)) %>% 
  group_by(Year, Geocluster) %>% 
  mutate(total=sum(count), percent=round(100*count/total,2))

write.csv(data, file="new_meta.csv", row.names = F)
write.csv(data_c, file="new_meta_Varroa_YearCountry.csv", row.names = F)
write.csv(data_s, file="new_meta_Varroa_YearState.csv", row.names = F)
write.csv(data_g, file="new_meta_Varroa_YearGeo.csv", row.names = F)

########################################################################################################

library(forcats)

rm(list = ls())

data <- read.csv("00_Big Dataset.csv")

data$Year[which(data$Year==88)] <- 1988
data$Year[which(data$Year==89)] <- 1989
data$Year[which(data$Year==90)] <- 1990
data$Year[which(data$Year==91)] <- 1991
data$Year[which(data$Year==92)] <- 1992
data$Year[which(data$Year==93)] <- 1993
data$Year[which(data$Year==94)] <- 1994
data$Year[which(data$Year==95)] <- 1995
data$Year[which(data$Year==96)] <- 1996
data$Year[which(data$Year==97)] <- 1997
data$Year[which(data$Year==98)] <- 1998

colnames(data)[7] <- "Geocluster"
data$Geocluster <- as.character(data$Geocluster)
data$Geocluster[which(data$Geocluster=="")] <- "Unknown"

data$Varroa <-fct_explicit_na(data$Varroa, na_level = "Unknown")
data$Country <-fct_explicit_na(data$Country, na_level = "Unknown")
data$State <-fct_explicit_na(data$State, na_level = "Unknown")
data$City_original <-fct_explicit_na(data$City_original, na_level = "Unknown")
data$Geocluster <-fct_explicit_na(data$Geocluster, na_level = "Unknown")


levels(data$State)[3] <- "VERACRUZ"
levels(data$State)[6] <- "TEXAS"


data$Country[which(data$Geocluster=="WEL")] <- "USA"
data$State[which(data$Geocluster=="WEL")] <- "TEXAS"

data$Geocluster <- as.character(data$Geocluster)

data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="BASELINE")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="BORDER")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="LOS INDIOS")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="LEVEE")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="GRANJENO")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="RESACA")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="HARLINGEN")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="HRLNGN")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="MCALLAN")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="MCALLEN")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="LAS PALOMAS")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="ESCONDIDOS")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="LAS CURVAS")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$Latitude==26.40306)] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="CUCHILLA")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="LACHATA")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="LACOPITA")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="RANCHO LA FERIA")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="PIRULI")] <- "TAM"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="RANCHO LA PALMA")] <- "TAM"
data$Geocluster[which(data$Geocluster=="Unknown" & data$State=="VERACRUZ")] <- "VER"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="LONE TREE")] <- "SAT"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="DOS MESQUITES")] <- "DAL"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="PUNTA")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="RANCHO BAZAN")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="MERCY MISSION")] <- "BOR"
data$Geocluster[which(data$Geocluster=="Unknown" & data$City_original=="MRCYMISS")] <- "BOR"

data$Geocluster[which(data$Geocluster=="Unknown" & data$State=="TAMAULIPAS")] <- "TAM"


data_c <- data %>% 
  group_by(Year,Country,Varroa) %>% 
  summarise(count=length(ID)) %>% 
  group_by(Year,Country) %>% 
  mutate(total=sum(count), percent=round(100*count/total,2))

data_s <- data %>% 
  group_by(Year,State,Varroa) %>% 
  summarise(count=length(ID)) %>% 
  group_by(Year,State) %>% 
  mutate(total=sum(count), percent=round(100*count/total,2))

data_g <- data %>% 
  group_by(Year, Geocluster, Varroa) %>% 
  summarise(count=length(ID)) %>% 
  group_by(Year, Geocluster) %>% 
  mutate(total=sum(count), percent=round(100*count/total,2))



write.csv(data, file="new_BigDataset.csv", row.names = F)
write.csv(data_c, file="new_BigDataset_Varroa_YearCountry.csv", row.names = F)
write.csv(data_s, file="new_BigDataset_Varroa_YearState.csv", row.names = F)
write.csv(data_g, file="new_BigDataset_Varroa_YearGeo.csv", row.names = F)

#######################################################################################################

