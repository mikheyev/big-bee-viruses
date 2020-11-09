########################################
##                                    ##   
## Author: Sherry Lin                 ##
## Student ID: u6890762               ##
## Project: evolution of bee viruses  ##
## Purpose: data arrangement          ##
## Version: 20                        ##
## Date: 22/09/20                     ##
##                                    ##
########################################

# addressing data: find the time first detect Varroa
library(forcats)

rm(list = ls())

# before import"00_seq_info.csv, change 'SAMPLE_ID':'4183-label-inside-334183' to '334183' 
# and '1206' Varroa: YES to NA (because they could be tracheal mites, for instance)
# import the data: "seq_info.csv" which contains the basic information about the sequenced samples
data <- read.csv("00_seq_info.csv", header=T)
# change the first column's name
colnames(data)[1] <- "ID"

summary(data) #Factors: Varroa, State, Country
summary(as.factor(data$Varroa)) #Varroa: NO(31),YES(21),NA(983)
summary(as.factor(data$State)) #States: TAMAULIPAS,TEXAS,VERACRUZ,N.VERACRUZ
summary(as.factor(data$Country)) #Country: MEXICO,USA

data$State[which(data$State=="N.VERACRUZ")] <- "VERACRUZ"

write.csv(data, file ="new_seq_info.csv", row.names = F)

#------------------------------------------------------------------------

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

data$State <- as.factor(data$State)
levels(data$State)[3] <- "VERACRUZ"
levels(data$State)[6] <- "TEXAS"

write.csv(data, file="new_meta.csv", row.names = F)
#----------------------------------------------------------------------

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
data$Geocluster[which(data$Geocluster=="")] <- NA
data$State[data$State=='N.VERACRUZ'] <- "VERACRUZ"
data$State[data$State=='Texas'] <- "TEXAS"
data$Country[data$Geocluster=="WEL"] <- "USA"
data$State[data$Geocluster=="WEL"] <- "TEXAS"

data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="BASELINE")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="BORDER")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="LOS INDIOS")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="LEVEE")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="GRANJENO")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="RESACA")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="HARLINGEN")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="HRLNGN")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="MCALLAN")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="MCALLEN")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="LAS PALOMAS")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="ESCONDIDOS")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="LAS CURVAS")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$Latitude==26.40306)] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="CUCHILLA")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="LACHATA")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="LACOPITA")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="RANCHO LA FERIA")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="PIRULI")] <- "TAM"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="RANCHO LA PALMA")] <- "TAM"
data$Geocluster[which(is.na(data$Geocluster) & data$State=="VERACRUZ")] <- "VER"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="LONE TREE")] <- "SAT"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="DOS MESQUITES")] <- "DAL"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="PUNTA")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="RANCHO BAZAN")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="MERCY MISSION")] <- "BOR"
data$Geocluster[which(is.na(data$Geocluster) & data$City_original=="MRCYMISS")] <- "BOR"

data$Geocluster[which(is.na(data$Geocluster) & data$State=="TAMAULIPAS")] <- "TAM"

write.csv(data, file="new_BigDataset.csv", row.names = F)

#================================================================================
library(data.table) #fread(): read the data
library(tools) #file_path_sans_ext(): without extension

rm(list = ls())

data <- read.csv("new_seq_info.csv")      

# add new column 'Date' when Varroa was first detected by observing                                          
Sys.setlocale("LC_TIME", "English")
data$Date <- as.Date(paste(data$Year,data$Mo,data$Day,sep='-'), "%Y-%m-%d")
data$Date.V <- '1995-4-11'
data$Date.V[which(data$Geocluster=="BOR")] <- "1992-4-8"
data$Date.V[which(data$Geocluster=="WEL")] <- "1994-4-27"
data$Date.V[which(data$Geocluster=="VER")] <- NA
data$Date.V <- as.Date(data$Date.V)

data$Year.V <- NA
data$Year.V[which(data$Geocluster=="BOR")] <- 1992
data$Year.V[which(data$Geocluster=="WEL")] <- 1994
data$Year.V[which(data$Geocluster=="TAM")] <- 1995

data$month <- difftime(data$Date,data$Date.V,units="days")/30
data$quarter <- difftime(data$Date,data$Date.V,units='days')/90
data$year <- difftime(data$Date,data$Date.V,units='days')/365

# If some samples only have 'Year', we can not calculate 'days' and data will show NA. So, I substitute 'Year'-'Varroa arrived year' for yearly data
for (i in 1:1035) {
  if(is.na(data$year[i])&data$Geocluster[i]!='VER'){data$year[i] <- data$Year[i]-data$Year.V[i]}
}

# rounding of numbers:months,quarters,years
# the negative numbers means the months, quarters, or years before Varroa arrived
data$month.1 <- NA
data$quarter.1 <- NA
data$year.1 <- NA
for(i in 1:1035){
  if(data$month[i]<0 & !is.na(data$month[i])){data$month.1[i] <- -ceiling(abs(data$month[i]))}
  if(data$month[i]>=0 & !is.na(data$month[i])){data$month.1[i] <- ceiling(data$month[i])}
  if(data$quarter[i]<0 & !is.na(data$quarter[i])){data$quarter.1[i] <- -ceiling(abs(data$quarter[i]))}
  if(data$quarter[i]>=0 & !is.na(data$quarter[i])){data$quarter.1[i] <- ceiling(data$quarter[i])}
  if(data$year[i]<0 & !is.na(data$year[i])){data$year.1[i] <- -ceiling(abs(data$year[i]))}
  if(data$year[i]>=0 & !is.na(data$year[i])){data$year.1[i] <- ceiling(data$year[i])}
}

# regroup sample by Varroa exposure time: <7, <5, <3, <2, <1 years
data$year.2[data$year.1==7] <-'<7'
data$year.2[data$year.1==5] <-'<5'
data$year.2[data$year.1==3] <-'<3'
data$year.2[data$year.1==2] <-'<2'
data$year.2[which(data$year.1==1|data$year.1==0)] <-'<1'
data$year.2[which(data$year.1<0|is.na(data$year.1))] <- '0'

# regroup sample by Varroa exposure time: >2, <2, <1 years
data$year.3[which(data$year.1==7|data$year.1==5|data$year.1==3)] <-'>2'
data$year.3[which(data$year.1==2)] <-'<2'
data$year.3[which(data$year.1==0|data$year.1==1)] <-'<1'
data$year.3[which(data$year.1<0|is.na(data$year.1))] <- '0'

data$V.arrival <- 1
data$V.arrival[which(data$year.3=='0')] <- 0

data <- data[ ,c(1,21,22,5,23:31,8,7,9,11,12,10,18:20,14,32)]

# getting files' name without extension (the files' names are the sample's ID) 
setwd("D:/Sherry/文件/ResearchProject/BeeVirus")
unzip.filename <- list.files("00_known_viruses")                                     
unzip.filename.new <- file_path_sans_ext(list.files("00_known_viruses")) 

# change the the working directory into a file named 'unzip'
# merge each sequenced results into a dataframe
setwd("00_known_viruses")                                                   
a <- data.frame()                                             
for (i in 1:1043){                                               
  b <- fread(unzip.filename[i])                                  
  b$ID <- unzip.filename.new[i]                              
  a <- rbind(a,b)                                         
}     

tpm <- dcast(a,ID~target_id,value.var='tpm')
tpm <- as.data.frame(tpm)

# used the website: "ncbi.nlm.nih.gov/search/" to get virus description by accession
# for example: the IDs beginning with 'MK0324' are the same virus, so should sum the TPMs for them
# also, change the columns' name from accession number to viruses' name
detected <- which(apply(tpm[ ,2:29],2,sum,na.rm=T)!=0)+1
tpm <- tpm[ ,c(1,detected)]
names(tpm)[2:18] <- c('DWV_C','Amfv','VDV3','Arv2','Vov1','SBV','ABPV',
                       'BQCV','KBV','DWV_A','DWV_B','IAPV','CBPV','SBPV',
                       'VTV','LSV','VDV2')
tpm <- tpm[ ,c('ID',sort(names(tpm[ ,2:18])))]
tpm[which(is.na(tpm$ABPV)),2:18] <- 0


setwd("D:/Sherry/文件/ResearchProject/BeeVirus")
# when merge(full join)'data' & 'TPM', I found :
# 'TPM' includes 9 negative controls and 334713, but 'data' does not
# 'data' includes 7227, but 'TPM' does not
# thus, used the metadata to get the information of 334713 and insert it into "seq_info.csv" (in excel)
# then, run the codes again
# removed the negative controls and 7227 by 'merge' (inner join)
data_merge <- merge(data,tpm,by="ID")

# import data: "k2.txt" is the extent of Africanization of each bee
# then, merge this data into seq_tpm_data by ID
africanized <- read.table("00_k2.txt")                                                                                    
names(africanized)[1] <- "ID"                                                                                            
africanized$ID <- as.character(africanized$ID)                                                                           
africanized$ID[759] <- "334183"                                                                                          

data_merge <- merge(data_merge, africanized[ ,c(1,3)], by = "ID")  # merge by ID, africanized data does not have "7227" 
colnames(data_merge)[42] <- "Afr.Extent"
data_merge[42] <- data_merge[42]*100
data_merge$Africanization <- 0
data_merge$Africanization[which(data_merge$Afr.Extent>10^-6)] <- 1


#imputation: PRCP,TMAX,TMIN
library(weathermetrics)
library(mice)
mice.data <- mice(data_merge[,c(4,17,18,20:22,24:43)],m=3,maxit=50,method='cart')
a <- complete(mice.data,1)
data_merge[,20:22] <- a[,4:6]
data_merge$TMAX <- fahrenheit.to.celsius(data_merge$TMAX, round = 2)
data_merge$TMIN <- fahrenheit.to.celsius(data_merge$TMIN, round = 2)

write.csv(data_merge, file = "new_data.csv", row.names = F)
