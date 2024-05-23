library(CASdatasets)
library(anytime)
data("canlifins")
setwd("~/Desktop/Canlifins")
dataset <- read.delim("C:/Users/z3532149/Desktop/LFCR/AugmentedVM/UOFM.txt", header = FALSE, sep = " ", dec = ".")

dataset <- UOFM

unique(dataset$V1)
unique(dataset$V2)
unique(dataset$V3)
unique(dataset$V5)
unique(dataset$V6)
unique(dataset$V7)
unique(dataset$V8)
unique(dataset$V9)


unique(dataset$V18)
unique(dataset$V19)

# - Data in the Date format
dataset$V4 <- anydate(dataset$V4)   
dataset$V4 <- as.Date(dataset$V4, tryFormats = "%Y %m %d")
dataset$V10 <- anydate(dataset$V10)   
dataset$V10 <- as.Date(dataset$V10, tryFormats = "%Y %m %d")
dataset$V11 <- ifelse(dataset$V11==0, 15000101, dataset$V11)
dataset$V11 <- anydate(dataset$V11)   
dataset$V11 <- as.Date(dataset$V11, tryFormats = "%Y %m %d")
dataset$V12 <- anydate(dataset$V12)   
dataset$V12 <- as.Date(dataset$V12, tryFormats = "%Y %m %d")
dataset$V13 <- anydate(dataset$V13)   
dataset$V13 <- as.Date(dataset$V13, tryFormats = "%Y %m %d")
dataset$V16 <- ifelse(dataset$V16==0, 25000101, dataset$V16)
dataset$V16 <- anydate(dataset$V16)   
dataset$V16 <- as.Date(dataset$V16, tryFormats = "%Y %m %d")
dataset$V17 <- ifelse(dataset$V17==0, 25000101, dataset$V17)
dataset$V17 <- anydate(dataset$V17)   
dataset$V17 <- as.Date(dataset$V17, tryFormats = "%Y %m %d")

colnames(dataset) <- c("V1", "V2", "V3", "StartObs", "V5", "V6", "V7", "V8", "V9", "StartContract", "Date3", "BirthPri", "BirthSec", "GPri", "GSec", "DeathPri", "DeathSec", "IncomePri", "IncomeSec")

# - Eliminate same-gender couples
dataset <- dataset[dataset$GPri!=dataset$GSec,]

#dataset$code <- paste(dataset$V1, dataset$V2, dataset$V3, sep="")

# - Get info
## - Male as primary policyholder
dataset$maleprim <- ifelse(dataset$GPri=="M",1,0)
## - Male Birth and Death date and age when joining the observational plan, and death indicator
dataset$malebirth <- as.Date(ifelse(dataset$GPri=="M", dataset$BirthPri, dataset$BirthSec), origin="1970-01-01") 
dataset$maledeath <- as.Date(ifelse(dataset$GPri=="M", dataset$DeathPri, dataset$DeathSec), origin="1970-01-01")
dataset$maleage_C <- difftime(dataset$StartContract, dataset$malebirth, units="days")/365.25 # - age at the beginning of the observational period
dataset$DeathMInd <- ifelse(dataset$maledeath!="2500-01-01",1,0)

## - Female Birth and Death date and age when joining the observational plan, and death indicator
dataset$femalebirth <- as.Date(ifelse(dataset$GPri=="F", dataset$BirthPri, dataset$BirthSec), origin="1970-01-01")
dataset$femaledeath <- as.Date(ifelse(dataset$GPri=="F", dataset$DeathPri, dataset$DeathSec), origin="1970-01-01")
dataset$femaleage_C <- difftime(dataset$StartContract, dataset$femalebirth, units="days")/365.25 # - age at the beginning of the observational period
dataset$DeathFInd <- ifelse(dataset$femaledeath!="2500-01-01",1,0)

## - Age difference male-female
dataset$agediffMF <- -difftime(dataset$malebirth, dataset$femalebirth, units="days")/365.25
dataset$maleelder <- ifelse(dataset$agediffMF>0,1,0)

### - Added check on start observation and start contract (03-05-2024)
dataset$SO <- as.Date(ifelse(dataset$StartContract < "1988-12-29", as.Date("1988-12-29", origin="1970-01-01"), dataset$StartContract), origin="1970-01-01")

## - Eliminate people likely not to be couples
dataset <- dataset[!dataset$femaleage<2,] # - eliminate people with age lower than 2

## - Operations from Frees et. al.
dataset$a <- difftime(min(dataset$SO), dataset$StartContract, units="days")/365.25
dataset$a <- ifelse(dataset$a<0, 0, dataset$a)
dataset$x <- difftime(dataset$StartContract, dataset$malebirth, units="days")/365.25 # - Male age at contract inception
dataset$y <- difftime(dataset$StartContract, dataset$femalebirth, units="days")/365.25 # - Female age at contract inception
dataset$X <- difftime(dataset$maledeath, dataset$malebirth, units="days")/365.25 # - age at the beginning of the observational period
dataset$Y <- difftime(dataset$femaledeath, dataset$femalebirth, units="days")/365.25 # - age at the beginning of the observational period
dataset$X <- ifelse(dataset$X<0, 0, dataset$X)
dataset$Y <- ifelse(dataset$Y<0, 0, dataset$Y)

## - Check if insureds are alive at the start of the observational period
dataset$malecheck <- dataset$maledeath > "1988-12-29"
dataset$femalecheck <- dataset$femaledeath > "1988-12-29"

## - Define variable future lifetime
dataset$T1 = 0 # - Future lifetime for males
dataset$T2 = 0 # - Future lifetime for females

dataset$T1 <- ifelse(dataset$DeathMInd==0, -difftime(dataset$SO, "1994-01-01", units="days")/365.25, -difftime(dataset$SO, dataset$maledeath, units="days")/365.25)
dataset$T2 <- ifelse(dataset$DeathFInd==0, -difftime(dataset$SO, "1994-01-01", units="days")/365.25, -difftime(dataset$SO, dataset$femaledeath, units="days")/365.25)

dataset$b <- difftime("1994-01-01", max(as.Date("1988-12-29"), dataset$StartContract), units="days")/365.25

# - Dataset for analysis without duplicate records
dset <- as.data.frame(dataset$maleage_C)
colnames(dset) <- "ageM"
dset$ageM <- as.numeric(dataset$maleage_C)
dset$ageF <- dataset$femaleage_C
dset$agediff <- abs(dataset$agediffMF)
dset$agediff2 <- dataset$agediffMF
dset$maleelder <- dataset$maleelder
dset$dmale <- dataset$DeathMInd
dset$dfemale <- dataset$DeathFInd
dset$TObs1 <- dataset$T1
dset$TObs2 <- dataset$T2
#dset$TObs2 <- round(dataset$T,2)
#dset$trunc2 <- round(dataset$a,2)
dset$trunc <- dataset$a
dset <- unique(dset)

ds <- NULL
save.image("~/Desktop/Canlifins/Canlifins_elab.RData")
# - Give id to dataset to identify couples with same age. In this way we can analyse which differences are there
dset$id <- 0
id <- 1
row <- 1
for(row in row:nrow(unique(dset[,1:2]))){
  dset$id[as.numeric(dset$ageM)==as.numeric(unique(dset[,1:2])[row,1]) & as.numeric(dset$ageF)==as.numeric(unique(dset[,1:2])[row,2])] <- row
  print(paste(row, "% iteration"))
}

# - Retrieve similar observations
## - First sort row and cols of the dataset by id
dset <- dset[order(dset$id),]
dset$sim = 0
for(i in 2:(nrow(dset)-1)){
  dset$sim[i] <- ifelse((dset$id[i]==dset$id[i-1] | dset$id[i]==dset$id[i+1]),1,0)
}

## - Correct observations with similarity value equal to 1 (count how many records are there for the same observation)
dset0 <- dset
dset0$idcount=0 # - initialize
for(id in unique(dset0$id)){
  dset0$idcount[dset0$id==id] <- nrow(dset0[dset0$id==id,])
}


for(id in unique(dset0$id)){
  dset0$dmale[dset0$id==id] <- max(dset0$dmale[dset0$id==id]) # - if someone is dead, just report it
  dset0$dfemale[dset0$id==id] <- max(dset0$dfemale[dset0$id==id]) # - if someone is dead, just report it
  dset0$TObs1[dset0$id==id] <- min(dset0$TObs1[dset0$id==id]) # - if the event occurs earlier, then report it, otherwise it can be a problem of non-reporting
  dset0$TObs2[dset0$id==id] <- min(dset0$TObs2[dset0$id==id]) # - if the event occurs earlier, then report it, otherwise it can be a problem of non-reporting
  dset0$trunc[dset0$id==id] <- min(dset0$trunc[dset0$id==id])
}

data <- unique(dset0)


dset0$id <- dset0$sim <- NULL

#dset1 <- unique(dset0[dset0$sim==1,])

dset2 <- rbind(dset0[dset0$sim==0,],dset1)
save.image("~/Desktop/Canlifins/Canlifins_elab2.RData")
dset[dset$agediff<1,]


mean(data$agediff)
var(data$agediff)
