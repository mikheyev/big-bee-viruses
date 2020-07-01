# change'4183-label-inside-334183' to '334183' in 00_bigbee and seq_info.csv         

library(data.table) #fread(): read the data
library(tools) #file_path_sans_ext(): without extension

rm(list = ls())

# before importing data, change 'SAMPLE_ID':'4183-label-inside-334183' to '334183'
# import the data: "sequenced.csv" which contains the basic information about the sequenced samples
# change the first column's name
data <- read.csv("seq_info.csv", header = TRUE, sep = ",")      
names(data)[1] <- "ID"             
data <- data[ ,c(1,3,4,5,10)]

# getting files' name without extension (the files' names are the sample's ID) 
unzip.filename <- list.files("00_bigbee")                                     
unzip.filename.new <- file_path_sans_ext(list.files("00_bigbee")) 

# change the the working directory into a file named 'unzip'
# merge each sequenced results into a dataframe
setwd("00_bigbee")                                                   
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

TPM <- TPM[ ,c(85279,1:85278)]
TPM.1 <- TPM

for (i in 2:85279){                                                 
  TPM[[i]] <- as.numeric(TPM[[i]])
}  

remove_list <- c()
for(i in 2:85279){
  a <- sum(TPM[i],na.rm = T)
  if (a==0){
    remove_list <- append(remove_list,i)
  }
}

TPM <- TPM[ ,-remove_list]

TPM$B_apis<-rowSums(TPM[ ,c(which(grepl("BBC0122_RS",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("BBC0122_RS",names(TPM))))]

TPM$B_asteroides<-rowSums(TPM[ ,c(which(grepl("BAST_RS",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("BAST_RS",names(TPM))))] 

TPM$F_perrara<-rowSums(TPM[ ,c(which(grepl("FPB0191_RS",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("FPB0191_RS",names(TPM))))]

TPM$G_apicola<-rowSums(TPM[ ,c(which(grepl("A9G17_RS",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("A9G17_RS",names(TPM))))]

TPM$L_apis<-rowSums(TPM[ ,c(which(grepl("JF72_RS",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("JF72_RS",names(TPM))))]

TPM$L_helsingborgensis<-rowSums(TPM[ ,c(which(grepl("JF73_RS",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("JF73_RS",names(TPM))))]

TPM$L_kullabergensis<-rowSums(TPM[ ,c(which(grepl("JF76_RS",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("JF76_RS",names(TPM))))]

TPM$L_kunkeei<-rowSums(TPM[ ,c(which(grepl("APS55_RS",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("APS55_RS",names(TPM))))]

TPM$L_melliventris<-rowSums(TPM[ ,c(which(grepl("JF74_RS",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("JF74_RS",names(TPM))))]

TPM$S_alvi<-rowSums(TPM[ ,c(which(grepl("SALWKB2_RS",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("SALWKB2_RS",names(TPM))))]

TPM$Melissococcus_plutonius<-rowSums(TPM[ ,c(which(grepl("NZ_CP0066",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("NZ_CP0066",names(TPM))))]

TPM$Ascosphaera_apis<-rowSums(TPM[ ,c(which(grepl("AZGZ",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("AZGZ",names(TPM))))]

TPM$Aspergillus_flavus<-rowSums(TPM[ ,c(which(grepl("NW_00247",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("NW_00247",names(TPM))))]

TPM$Aspergillus_niger<-rowSums(TPM[ ,c(which(grepl("NT_1665",names(TPM))|grepl("NC_007445.1",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("NT_1665",names(TPM))|grepl("NC_007445.1",names(TPM))))]

TPM$Aspergillus_fumigatus<-rowSums(TPM[ ,c(which(grepl("NC_0071",names(TPM))|grepl("NC_0072",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("NC_0071",names(TPM))|grepl("NC_0072",names(TPM))))]

TPM$Nosema_apis<-rowSums(TPM[ ,c(which(grepl("KE64",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("KE64",names(TPM))))]

TPM$Nosema_ceranae<-rowSums(TPM[ ,c(which(grepl("NW_020169",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("NW_020169",names(TPM))))]

TPM$Varroa_destructor<-rowSums(TPM[ ,c(which(grepl("XM_022",names(TPM))|grepl("XR_00267",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("XM_022",names(TPM))|grepl("XR_00267",names(TPM))))]

TPM$Apis_mellifera<-rowSums(TPM[ ,c(which(grepl("NM",names(TPM))|grepl("NR",names(TPM))|grepl("XM",names(TPM))|grepl("XR",names(TPM))))])
TPM <- TPM[ ,-c(which(grepl("NM",names(TPM))|grepl("NR",names(TPM))|grepl("XM",names(TPM))|grepl("XR",names(TPM))))]

names(TPM)[2:19] <- c("DWV","DWV_C","IAPV","KBV","CBPV","ABPV","SBPV","LSV","BQCV","SV","VDV1","VDV2",
                      "VDV3","AR","AM","VTV","VO1","Paenibacillus_larvae")

TPM <- TPM[ ,c(1,38,37,7,16,15,10,6,2:5,9,8,11:14,18,17,20:28,30,19,29,31,32,34,33,35,36)]
TPM[2:38] <- TPM[2:38]/10^4  #scale the TPM -> Transcript per hundred(0~100)

setwd("D:/Sherry/文件/RStudio/Bee's viruses evolution")
data_merge <- merge(data,TPM,by="ID")

# import data: "k2.txt" is the extent of Africanization of each bee
# merge this data by ID
africanized <- read.table("00_k2.txt")                                                                                      
names(africanized)[1] <- "ID"                                                                                            
africanized$ID <- as.character(africanized$ID)                                                                           
africanized$ID[759] <- "334183"                                                                                          

data_merge <- merge(data_merge, africanized, by = "ID")  # merge by ID, africanized data does not have "7227" 
data_merge <- data_merge[ ,-43]                           
colnames(data_merge)[43] <- "Afr.Extent"
data_merge[43] <- data_merge[43]*100

data_merge$V.arrival <- 0
data_merge$V.arrival[which(data_merge$Geocluster=="WEL"&data_merge$Year>=1994)] <- 1
data_merge$V.arrival[which(data_merge$Geocluster=="BOR"&data_merge$Year>=1992)] <- 1
data_merge$V.arrival[which(data_merge$Geocluster=="TAM"&data_merge$Year>=1995)] <- 1


write.csv(data_merge, file = "new_bigbee_data.csv", row.names = F)



