#### Load packages and read data ####
library(plyr)
library(dplyr)
library(ggplot2)
library(ggmap)
library(data.table)
library(geosphere)
library(maptools)
library(maps)
library(mapproj)
library(reshape2)
library(tidyr)
library(rgdal)
library(rgeos)
library(sp)
library(raster)
library(PBSmapping)
library(foreign)
library(pbapply)
library(pscl)
library(parallel)
library(pbmcapply)
library(gridExtra)
library(scales)
library(stringr)
library(mgcv)
library(e1071)
library(MASS)
library(viridis)
library(foreach)
library(doParallel)

DataFP <- "/Users/hannahwauchope/Documents/OneDrive - University Of Cambridge/PhD/Data/WaterbirdData_Tatsuya/FullDataSet_Edits/Hannah_Consolidation/"
ResultsFP <- "/Users/hannahwauchope/Documents/OneDrive - University Of Cambridge/PhD/Chapter1/Analysis/Temporal/"
ResultsGLMFP <- "/Users/hannahwauchope/Documents/OneDrive - University Of Cambridge/PhD/Chapter1/Analysis/Temporal/GLM/"
ResultsGAMFP <- "/Users/hannahwauchope/Documents/OneDrive - University Of Cambridge/PhD/Chapter1/Analysis/Temporal/GAM/"
FiguresFP <- "/Users/hannahwauchope/Documents/OneDrive - University Of Cambridge/PhD/Chapter1/Figures_Tmp/"
ncores <- 4

load(file=paste0(DataFP, "Spec6.RData")) #Original data (Containing site name, species name, year, count, and 'hours'(for effort))
birdcounts <- Spec6
load(file=paste0(ResultsFP, "TYB_list.RData")) #Created in "Prepare count data"
load(file=paste0(ResultsFP, "TYB_St.RData")) #Created in "Prepare count data"
GenLengthSpec <- read.csv(paste0(ResultsFP, "GenLengthSpec.csv")) #Created in "Prepare generation length data"

##Aesthetics
plotfontsize <- 10
plotlegendsize <- 10
resolution <- 100
ResFactor <- 0.3

#### Prepare count data ####

##Reduce data set to sites with at least 30 years of surveying
birdcounts <- subset(birdcounts, Dataset=="CBC")
SitesAndYears <- unique(birdcounts[c("SiteCode", "Species", "Year")]) #Get the years recorded at each survey site
SiteByYear <- dcast(SitesAndYears, SiteCode + Species ~., length, value.var = "Year") #Cast for number of years at each site
names(SiteByYear) <- c("SiteCode", "NumYears") #Add names
ThirYears <- subset(SiteByYear, NumYears>29) #Cut to sites that have at least 30 years of sampling

ThirYearBirds <- merge(birdcounts, ThirYears, by=c("SiteCode", "Species")) #Bring back in the bird counts
TYB <- ThirYearBirds[c("SiteCode", "Species", "PopulationCode", "Year", "Count", "Hours", "NumYears")] #Simplify dataset

##Because data comes from CBC, need to standardise for survey hours. 
##First we need to see which species show a linear trend between log(survey effort) and counts, and then reduce to only those
##We then include hours as an offset term in the models
TYB$SpecPop <- paste(TYB$Species, TYB$Population, sep="_") #Combine species name and population into one variable
birdcounts$SpecPop <- paste(birdcounts$Species, birdcounts$Population, sep="_") #Now do the same for bird counts, cos we want as many counts as possible for our hours check (not just those with 30 years)

CountvsHoursLogHours <- rbindlist(pblapply(unique(TYB$SpecPop), function(i){ #Cycle through each SpecPop
  birdie <- subset(birdcounts, SpecPop==i) #Subset dataframe to specpop
  glmmy <- glm(Count~log(Hours), data=birdie) #Run a GLM of counts vs. log of hours
  significant <- as.data.frame(ifelse(summary(glmmy)$coeff[-1,4]< 0.05, coef(glmmy, silent=TRUE)[[2]], NA)) #Return the coefficient IF the p value is less than 0.5, otherwise NA
  significant$SpecPop <- i #Add names
  names(significant) <- c("Slope", "SpecPop") #Add names
  return(significant)
}))

CountvsHoursLogHours <- CountvsHoursLogHours[!is.na(CountvsHoursLogHours$Slope),] #Remove the NAs (i.e. insignificant relations)
CountvsHoursLogHours <- subset(CountvsHoursLogHours, Slope>0) #Remove the negative slopes
TYB <- TYB[TYB$SpecPop %in% CountvsHoursLogHours$SpecPop,] #Reduce main dataset to just cleaned species

##Ok and now continue to clean up the data

TYB$SiteSpec <- paste(TYB$Site, TYB$SpecPop, sep="_") #Combine site and specpop into one field (So we now have a unique identifier for each 'population' i.e. each species at each site)
TYB <- TYB[,c("SiteSpec", "Year", "Count", "Hours")] #Reduce dataset to just what we need
TYB_list <- split(TYB, f = TYB$SiteSpec) #Divide into list based on sitespec

TYB_list <- lapply(TYB_list, function(x){ 
  x <- x[order(x$Year),]
  x$Count <- round(x$Count, 0)
  return(x[(nrow(x)-(29)):nrow(x),])
})#This loop reduces each list element to just 30 years of counts, using the most recent. And also rounds counts because there's a few decimal ones

TYB_list <- lapply(TYB_list, function(x) ifelse(sum(x$Count)<30, return(NULL), return(x))) #Remove any list elements with less than 30 counts over the 30 years
TYB_list <- TYB_list[!sapply(TYB_list, is.null)] #Remove null elements

save(TYB_list, file=paste0(ResultsFP, "TYB_list.RData")) #Save

##And now create a standardised list of species
TYB_Names <- as.data.frame(names(TYB_list))
names(TYB_Names) <- "SiteSpec"
TYB_Names$Spec <- as.data.frame(str_split_fixed(TYB_Names$SiteSpec, "[_]", 3))[,2] #Get species names
TYB_Names_Cast <- dcast(TYB_Names, Spec~., length, value.var="SiteSpec") #Find out how many sites of each species
TYB_Names_Cast_50 <- subset(TYB_Names_Cast, .>49) #Subset to species with at least 50 sites

#For each of those species, randomly extract 50 sites
TYB_Standardising <- TYB_Names[TYB_Names$Spec %in% TYB_Names_Cast_50$Spec,] #Get the sitespec values for each of the species
TYB_Standardising$Spec <- factor(TYB_Standardising$Spec, levels=unique(TYB_Standardising$Spec)[order(unique(TYB_Standardising$Spec))]) #Make the species factors (to allow split, next line)
TYB_Standardising <- split(TYB_Standardising, f=TYB_Standardising$Spec)

RandomNumSpec <- lapply(TYB_Names_Cast_50$., function(x) sample(1:x, 50, replace=FALSE)) #Get 50 random numbers sampled from the total number of sites for each species
TYB_St <- rbindlist(lapply(1:length(TYB_Standardising), function(x) TYB_Standardising[[x]][RandomNumSpec[[x]],]))$SiteSpec #Sample from all possible sitespec for each species, using the random numbers
save(TYB_St, file=paste0(ResultsFP, "TYB_St.RData")) #Save

#### Prepare generation length data ####
SiteSpec <- as.data.frame(names(TYB_list)) #Get all the sitespec
SiteSpec$Species <- str_split_fixed(SiteSpec$'names(TYB_list)', "[_]", 3)[,2] #Extract just species names
names(SiteSpec) <- c("SiteSpec", "Species") #Rename

Genl <- read.csv(paste0(ResultsFP, "GenLengths.csv")) #Pull in generation length data
Genl <- subset(Genl, GenerationLength!="NotGiven") #Remove any species where generation length wasn't given
Genl$GenerationLength <- as.numeric(as.character(Genl$GenerationLength)) #Class gen lengths as numeric
GenlAll <- merge(SiteSpec, Genl, by="Species") #Add generation length data to each sitespec

Genl1_5 <- subset(GenlAll, GenerationLength<6)
Genl1_5$GenLengthMax <- 5
Genl6_10 <- subset(GenlAll, GenerationLength>6 & GenerationLength<11)
Genl6_10$GenLengthMax <- 10
Genl11_15 <- subset(GenlAll, GenerationLength>11 & GenerationLength<16)
Genl11_15$GenLengthMax <- 15

SampleSize <- min(nrow(Genl1_5), nrow(Genl6_10), nrow(Genl11_15)) #Find the category with the smallest number of sitespec, and randomly reduce the others to this length (to make sure we're not unfairly comparing)

Genl1_5 <- Genl1_5[sample(1:nrow(Genl1_5), SampleSize),] #Reduce each of the categories to the smallest
Genl6_10 <- Genl6_10[sample(1:nrow(Genl6_10), SampleSize),] #Reduce each of the categories to the smallest
Genl11_15 <- Genl11_15[sample(1:nrow(Genl11_15), SampleSize),] #Reduce each of the categories to the smallest

GenLengthSpec <- rbind(Genl1_5, Genl6_10, Genl11_15)

write.csv(GenLengthSpec, paste0(ResultsFP, "GenLengthSpec.csv"))

#### Models 1. Write functions ####
#These functions take the input data from the run models, with a slope and significance level of each population, with each combination of consecutive and interval samplings COMPARED TO ONE value of complete trend lengths (run multiple times for each complete trend length)
#e.g. slope and significance of population x when sampled for 10 years starting at year 4 compared to slope and significance of population x when sampled for 20 years starting at year 3
#Direction function returns the number of samples that are correct, a false alarm, missed detection or opposing when compared to the complete trend
#Magnitude function returns the number of samples that are correct according to different tolerances when comparing the slope of the sample to the complete trend
#If Standardise is TRUE, subset the input data to the standardised dataset (98 species with 50 randomly selected sites each)
#If GenLength == 0, use all data. If it is 5 the generation length = 5 or less, if 10, = 6-10, if 15, = 11-15
#If EDFNum == 0, use all data. Otherwise 1 == EDF 1,2; 3 == EDF 3,4; 5 == EDF 5,6; 7== EDF 7,8.
#Significance sets the p value cutoff (0.05 for all published results)

Direction <- function(DATA, Standardise, GenLength, EDFNum, Significance){
  Input <- DATA
  SamplingType <- unique(Input$SamplingType)
  #Reduce the data according to the various subsets
  if(GenLength==0){
    if(Standardise==TRUE){Input <- Input[Input$SiteSpec %in% TYB_St,]}
  } else {
    GenLengthSub <- subset(GenLengthSpec, GenLengthMax==GenLength)
    Input <- Input[Input$SiteSpec %in% GenLengthSub$SiteSpec,]
  } 
  
  if(EDFNum!=0){
    Input <- DATA
    Input[is.na(Input$EDFCompleteLength),]$EDFCompleteLength <- 0
    Input <- subset(Input, EDFCompleteLength<EDFNum)
  }
  
  if(nrow(Input)==0){
    return(NULL)
  }
  
  #Change the slope to just positive, negative or insignificant
  Input$Slope <- sapply(Input$Slope, function (x) if(is.na(x)){"Insignificant"} else if(x > 0){ "Positive" } else{ "Negative" })
  Input$CompleteLengthSlope <- sapply(Input$CompleteLengthSlope, function (x) if(is.na(x)){"Insignificant"} else if(x > 0){ "Positive" } else{ "Negative" })
  
  #Create a binary value for significance based on p value
  Input$Significant <- sapply(Input$Significant, function (x) if(is.na(x)){0} else if(x < Significance){1} else{0})
  Input$SignificantCompleteLength <- sapply(Input$SignificantCompleteLength, function (x) if(is.na(x)){0} else if(x < Significance){1} else{0})
  
  #Change the slope to insignificant where relevant
  Input[Input$Significant==0]$Slope <- "Insignificant"
  Input[Input$SignificantCompleteLength==0]$CompleteLengthSlope <- "Insignificant"
  Input <- Input[,c("SiteSpec", "NumYears", "Slope", "CompleteLengthSlope")]
  
  #Calculate the false alarms, missed detections, correct and incorrect
  MD_Pos <- subset(Input, CompleteLengthSlope=="Positive")
  MD_Pos$CompleteLengthSlope <- NULL
  MD_Pos$Slope <- recode_factor(MD_Pos$Slope, "Positive"="Correct", "Negative" = "Incorrect", "Insignificant" = "Missed Detection")
  
  MD_Neg <- subset(Input, CompleteLengthSlope=="Negative")
  MD_Neg$CompleteLengthSlope <- NULL
  MD_Neg$Slope <- recode_factor(MD_Neg$Slope, "Positive"="Incorrect", "Negative" = "Correct", "Insignificant" = "Missed Detection")

  FA <- subset(Input, CompleteLengthSlope=="Insignificant")
  FA$CompleteLengthSlope <- NULL
  FA$Slope <- recode_factor(FA$Slope, "Positive" = "False Alarm", "Negative" = "False Alarm", "Insignificant" = "Correct")
  
  summarise <- rbind(FA, MD_Pos, MD_Neg) 

  rm(Input)
  
  #Cast to find number of populations of each category
  summarise <- dcast(summarise, NumYears~Slope, length, value.var="Slope")
  summarise <- melt(summarise, id.vars="NumYears")
  
  #Split numyears for interval
  if(SamplingType=="Interval"){
    summarise[,c("IntervalLength", "NumYears")] <- str_split_fixed(summarise$NumYears, "[.]", 2)
  }
  
  #Order by number of years, add metadata
  summarise$NumYears <- as.numeric(as.character(summarise$NumYears))
  summarise <- summarise[order(-summarise$NumYears),]
  summarise$SamplingType <- SamplingType
  summarise$Model <- modelfam
  if(SamplingType=="Consecutive"){
    summarise$CompleteLength <- unique(DATA$CompleteLengthNum)
  }
  summarise$GenLength <- GenLength
  if(!Standardise){
    summarise$CastType <- "All"
  } else{
    summarise$CastType <- "St"
  }
  summarise$EDF <- EDFNum
  return(summarise)
}
Magnitude <- function(DATA, Standardise, GenLength, EDFNum, Significance){
  Input <- DATA
  SamplingType <- unique(Input$SamplingType)
  
  #Reduce the data according to the various subsets
  
  if(GenLength==0){
    if(Standardise==TRUE){Input <- Input[Input$SiteSpec %in% TYB_St,]}
  } else {
    GenLengthSub <- subset(GenLengthSpec, GenLengthMax==GenLength)
    Input <- Input[Input$SiteSpec %in% GenLengthSub$SiteSpec,]
  } 
  
  if(EDFNum!=0){
    Input <- DATA
    Input[is.na(Input$EDFCompleteLength),]$EDFCompleteLength <- 0
    Input <- subset(Input, EDFCompleteLength<EDFNum)
  }
  

  Input <- subset(Input, Significant<Significance & SignificantCompleteLength<Significance)
  
  if(nrow(Input)==0){
    return(NULL)
  }
  
  if(SamplingType=="Consecutive"){
    Input[,c("SamplingType", "Model", "CompleteLengthNum", "EDF", "EDFCompleteLength")] <- NULL
  } else{
    Input <- Input[,c("SiteSpec","NumYears", "Slope", "CompleteLengthSlope")]
    Input <- Input[!is.na(Input$Slope),]
  }
  
  Input <- Input[complete.cases(Input)] #Remove NAs
  
  #Find correct and incorrect for the various tolerances
  summarise <- rbindlist(lapply(c(0.01, 0.025, 0.05, 0.1, 0.25, 0.5), function(x){
    ok <- Input
    ok$Trend <- ifelse(Input$Slope<(Input$CompleteLengthSlope+x) & Input$Slope>(Input$CompleteLengthSlope-x), "Correct", "Incorrect")
    ok$Tolerance <- x
    summarise <- dcast(ok, NumYears+Tolerance~Trend, length, value.var="Trend")
    return(summarise)
  }))
  rm(Input)
  
  #Split numyears for interval
  if(SamplingType=="Interval"){
    summarise[,c("IntervalLength")] <- str_split_fixed(summarise$NumYears, "[.]", 2)[,1]
    summarise[,c("NumYears")] <- str_split_fixed(summarise$NumYears, "[.]", 2)[,2]
  }
  
  #Add metadata
  summarise$NumYears <- as.numeric(as.character(summarise$NumYears))
  summarise$SamplingType<- SamplingType
  summarise$Model <- modelfam
  summarise$CompleteLength <- unique(DATA$CompleteLengthNum)
  summarise$GenLength <- GenLength
  if(!Standardise){
    summarise$CastType <- "All"
  } else {
    summarise$CastType <- "St"
  } 
  summarise$EDF <- EDFNum
  return(summarise)
}

#### Models 2. Consecutive ####
modelfam <- "nb" #this can be adjusted to "pois" or "quasi" if you want to test other model families
Consecutive <- rbindlist(lapply(3:30, function(j){ #This apply runs through the number of years of each model run (3 years to 30 years)
  StartYear <- pbmclapply(1:(31-j), function(i){ #And this one runs through the possible start years (so e.g. j=5, i=10 would mean a model running on 5 years of data, taken from years 10-14 for every SiteSpec)
    if(file.exists(paste0(ResultsGLMFP,"Consecutive/Consecutive_GLM_", j, "Years_",i, "StartYear_", modelfam, ".csv"))==TRUE){ #Because these models take a long time to run code could run on multiple clusters. This means anything that's already done is not done again
      return(NULL)
    } else {
      write.csv(NULL, paste0(ResultsGLMFP,"Consecutive/Consecutive_GLM_", j, "Years_",i, "StartYear_", modelfam, ".csv"), row.names=FALSE)
    }
    Degraded <- lapply(TYB_list, function(x) x[i:(i+j-1),]) #Cut the counts down to the correct number of years (j), starting at the right year (i)
    Degraded2 <- lapply(Degraded, function(x) ifelse(sum(x$Count)==0, return(NULL), return(x))) #Turn all entries with 0 at all times to NULL (no point running a model on them)
    Degraded2 <- Degraded2[!sapply(Degraded2, is.null)] #Remove null elements
    
    if(modelfam=="quasi"){ #Run the glm, using the model family defined (for the paper, model family is negative binomial - "nb")
      theglmlist <- lapply(Degraded2, function (x) glm(Count~Year, offset=log(x$Hours), quasipoisson(link="log"), data=x))
    } else if(modelfam=="pois"){
      theglmlist <- lapply(Degraded2, function (x) glm(Count~Year, offset=log(x$Hours), poisson(link="log"), data=x))
    } else {
      theglmlist <- lapply(Degraded2, function (x) tryCatch(glm.nb(Count~Year + offset(log(x$Hours)), link=log, data=x), error=function(e){NULL}))
    }

    thegamlist <- lapply(Degraded2, function (x) tryCatch(gam(Count~s(Year) + offset(log(Hours)), family = quasipoisson(link='log'), data=x), error=function(e){NULL})) #Also run a GAM on the model, to determine the curve shape (using estimate degrees of freedom, EDF). 
    #Nb due to issues with running negative binomial GAMs with sparse data, this uses a quasipoisson distribution. Checks indicated high convergence between the two model families for getting EDF. 
    
    theglmlist <- theglmlist[!sapply(theglmlist, is.null)] #Remove nulls 
    thegamlist <- thegamlist[!sapply(thegamlist, is.null)] #Remove nulls
    
    slopes <- sapply(theglmlist, function (x) tryCatch(coef(x, silent=TRUE)[[2]], error=function(e){NULL})) #Extract the coefficient of slope (i.e. r)
    slopes <- slopes[!sapply(slopes, is.null)] #Remove nulls, organise data
    slopesnames <- names(slopes)
    slopes <- as.data.frame(unlist(slopes))
    slopes$Data <- slopesnames
    
    significantglm <- sapply(theglmlist, function (x) tryCatch(summary(x)$coeff[-1,4], error=function(e){NULL})) #Extract the significance of the slope coefficient (i.e. p)
    significantglm <- significantglm[!sapply(significantglm, is.null)] #Remove nulls, organise data
    significantglmnames <- names(significantglm)
    significantglm <- as.data.frame(unlist(significantglm))
    significantglm$Data <- significantglmnames
    
    edfgam <- sapply(thegamlist, function (x) tryCatch(summary(x)$s.table[[1]], error=function(e){NA})) #Extract the estimated degrees of freedom from the GAM (EDF)
    edfgam <- edfgam[!sapply(edfgam, is.null)] #Remove nulls, organise data
    edfgamnames <- names(edfgam)
    edfgam <- as.data.frame(unlist(edfgam))
    edfgam$Data <- edfgamnames
    
    if(nrow(significantglm)!=0 & nrow(slopes)!=0){ #Providing not all models have returned null, create final dataset
      finaldata <- merge(slopes, significantglm, by="Data")
      if(nrow(edfgam)>0){ #Sometimes this is zero, because no GAMs were able to run on the data (often because data only 3 years in length)
        finaldata <- merge(finaldata, edfgam, by="Data", all=TRUE)
      } else {
        finaldata$edfgam <- NA
      }
      names(finaldata) <- c("SiteSpec", "Slope", "Significant", "EDF")
    } else {
      finaldata <- data.frame("SiteSpec"=NULL, "Slope"=NULL, "Significant"=NULL, "EDF"=NULL)
    }
    
    #And now add in the sitespec where models didn't run and/or there were zeros in all years (and therefore conclusion is insignificant)
    if(nrow(finaldata) != length(Degraded)){ #I.e. if there are some records without data
      InsigSiteSpec <- as.data.frame(names(Degraded)[!names(Degraded) %in% finaldata$SiteSpec]) #Extract these records
      names(InsigSiteSpec) <- "SiteSpec" #Define records as NA on all counts
      InsigSiteSpec$Slope <- NA
      InsigSiteSpec$Significant <- NA
      InsigSiteSpec$EDF <- NA
      finaldata <- rbind(finaldata, InsigSiteSpec) #Add them to final dataset
    }
    
    finaldata$NumYears <- as.numeric(as.character(j)) #Add metadata to final dataset 
    finaldata$StartYear <- as.numeric(as.character(i))
    finaldata$SamplingType<- "Consecutive"
    finaldata$Model <- modelfam
    write.csv(finaldata, paste0(ResultsGLMFP,"Consecutive/Consecutive_GLM_", j, "Years_",i, "StartYear_", modelfam, ".csv"), row.names=FALSE) #Write out results
    return(finaldata)
  }, mc.cores=ncores)
  return(rbindlist(StartYear))
}))

Consecutive <- rbindlist(pblapply(list.files(path=paste0(ResultsGLMFP, "Consecutive/"), pattern=paste0("*", modelfam, ".csv"), full.names=TRUE), fread), fill=TRUE) #Pull the data from the above back together (assuming it has taken multiple runs to complete, otherwise the result will just be taken from the function above)
#The following function organises the data so that samples can be compared to multiple complete trend lengths. E.g. A complete trend 10 years long taken from years 10-19 could have 3 year samples taken from it starting from years 10,11,12...17 etc)
#It then compares sample trend to complete trend according to the two comparison techniques (Direction and Magnitude)
CastConsecutive <- pblapply(c(4:30), function(CompleteLength){
  if((file.exists(paste0(ResultsGLMFP, "Summaries/Direction_Consecutive_", modelfam, "_", CompleteLength,".csv"))&
      file.exists(paste0(ResultsGLMFP, "Summaries/Magnitude_Consecutive_", modelfam, "_", CompleteLength,".csv")))==TRUE){
    return(NULL)
  } else {
    write.csv(NULL, paste0(ResultsGLMFP, "Summaries/Direction_Consecutive_", modelfam, "_", CompleteLength,".csv"))
    write.csv(NULL, paste0(ResultsGLMFP, "Summaries/Magnitude_Consecutive_", modelfam, "_", CompleteLength,".csv"))

    #Work down through the years, removing the start years that don't fit (e.g. start yr at yr21 when we're only going up to 20now)
    CompleteLengthSub <- Consecutive
    CompleteLengthSub[,c("SamplingType", "Model")] <- NULL
    CompleteLengthSlopes <- subset(CompleteLengthSub, NumYears==CompleteLength)
    names(CompleteLengthSlopes) <- c("SiteSpec","CompleteLengthSlope", "SignificantCompleteLength", "EDFCompleteLength", "CompleteLengthNum", "CompleteLengthStartYear")
    
    SamplesSub <- subset(CompleteLengthSub, NumYears<CompleteLength)
    rm(CompleteLengthSub)
    AllTheStartYears <- rbindlist(pblapply(unique(CompleteLengthSlopes$CompleteLengthStartYear), function(x){
      CompleteLengthSubStartYear <- subset(CompleteLengthSlopes, CompleteLengthStartYear==x)
      StartYearSub <- SamplesSub
      StartYearSub <- subset(StartYearSub, StartYear > (x-1))
      StartYearSub$Keep <- ifelse((StartYearSub$StartYear+(StartYearSub$NumYears)-1) < (x+CompleteLength), 1, 0)
      StartYearSub <- subset(StartYearSub, Keep==1)
      StartYearMerge <- merge(StartYearSub, CompleteLengthSubStartYear, All=T)
      return(StartYearMerge)
    }))
    rm(CompleteLengthSlopes, SamplesSub)
    AllTheStartYears$Keep <- NULL
    AllTheStartYears$Model <- modelfam
    AllTheStartYears$SamplingType <- "Consecutive"
    
    Direction_All <- rbindlist(pblapply(c(0,5,10,15), function(x) Direction(AllTheStartYears, Standardise=FALSE, GenLength=x, EDFNum=0, Significance = 0.05)))
    Direction_St <- rbindlist(lapply(c(0,5,10,15), function(x) Direction(AllTheStartYears, Standardise=TRUE, GenLength=x, EDFNum=0, Significance = 0.05)))
    Direction_EDF <- rbindlist(pblapply(c(3,5,7,9), function(x) Direction(AllTheStartYears, Standardise=FALSE, GenLength=0, EDFNum=x, Significance = 0.05)))
    write.csv(rbind(Direction_All, Direction_St, Direction_EDF), file=paste0(ResultsGLMFP, "Summaries/Direction_Consecutive_", modelfam, "_", CompleteLength,".csv"), row.names=FALSE)
    
    Magnitude_All <- rbindlist(pblapply(c(0,5,10,15), function(x) Magnitude(AllTheStartYears, Standardise=FALSE, GenLength=x, EDFNum=0, Significance = 0.05)))
    Magnitude_St <- rbindlist(pblapply(c(0,5,10,15), function(x) Magnitude(AllTheStartYears, Standardise=TRUE, GenLength=x, EDFNum=0, Significance = 0.05)))
    Magnitude_EDF <- rbindlist(pblapply(c(3,5,7,9), function(x) Magnitude(AllTheStartYears, Standardise=FALSE, GenLength=0, EDFNum=x, Significance = 0.05)))
    write.csv(rbind(Magnitude_All, Magnitude_St, Magnitude_EDF), file=paste0(ResultsGLMFP, "Summaries/Magnitude_Consecutive_", modelfam, "_", CompleteLength,".csv"), row.names=FALSE)
  }
})

#### Models 3. Intervals ####  
#Create a list of all the ways 30 years can be sampled in intervals (all possible number of years and space between years)
modelfam <- "nb" #this can be adjusted to "pois" or "quasi" if you want to test other model families

Intervals <- unlist(unlist(lapply(1:28, function(x){
  lapply(2:16, function(y){
    ok <- seq(from = x, to = 30, by = y)
    if(length(ok)>2){
      ok2 <- lapply(3:(length(ok)), function(z){
        ok[1:z]
      })
      return(ok2)
    } else{return(NULL)}
  })
}), recursive = FALSE), recursive = FALSE)

if(!file.exists(file=paste0(ResultsGLMFP, "Intervals/Intervals_", modelfam, ".csv"))){
  Interval <- rbindlist(pbmclapply(Intervals, function(IntValues){
    Degraded <- lapply(TYB_list, function(x) x[IntValues,]) #Cut the counts down to the right intervals
    Degraded2 <- lapply(Degraded, function(x) ifelse(sum(x$Count)==0, return(NULL), return(x))) #Remove entries with 0 at all times
    Degraded2 <- Degraded2[!sapply(Degraded2, is.null)]
    
    if(modelfam=="quasi"){ #Run the glm, using the model family defined (for the paper, model family is negative binomial - "nb")
      theglmlist <- lapply(Degraded2, function (x) glm(Count~Year, offset=log(x$Hours), quasipoisson(link="log"), data=x))
    } else if(modelfam=="pois"){
      theglmlist <- lapply(Degraded2, function (x) glm(Count~Year, offset=log(x$Hours), poisson(link="log"), data=x))
    } else {
      theglmlist <- lapply(Degraded2, function (x) tryCatch(glm.nb(Count~Year + offset(log(x$Hours)), link=log, data=x), error=function(e){NULL}))
    }
    
    theglmlist <- theglmlist[!sapply(theglmlist, is.null)] #Remove nulls
    
    slopes <- sapply(theglmlist, function (x) tryCatch(coef(x, silent=TRUE)[[2]], error=function(e){NULL}))  #Extract the coefficient of slope (i.e. r)
    slopes <- slopes[!sapply(slopes, is.null)] #Remove nulls and organise data
    slopes <- as.data.frame(unlist(slopes))
    slopes$Data <- rownames(slopes)
    
    significantglm <- sapply(theglmlist, function (x) tryCatch(summary(x)$coeff[-1,4], error=function(e){NULL})) #Extract the significance of the slope coefficient (i.e. p)
    significantglm <- significantglm[!sapply(significantglm, is.null)] #Remove nulls and organise data
    significantglm <- as.data.frame(unlist(significantglm))
    significantglm$Data <- rownames(significantglm)
    
    if(nrow(significantglm)!=0 & nrow(slopes)!=0){ #Providing not all models have returned null, create final dataset
      finaldata <- merge(slopes, significantglm, by="Data")
      names(finaldata) <- c("SiteSpec", "Slope", "Significant")
    } else {
      finaldata <- data.frame("SiteSpec"=NULL, "Slope"=NULL, "Significant"=NULL)
    }
    
    #And now add in the sitespec that are insignificant
    if(nrow(finaldata) != length(Degraded)){ 
      InsigSiteSpec <- as.data.frame(names(Degraded)[!names(Degraded) %in% finaldata$SiteSpec])
      names(InsigSiteSpec) <- "SiteSpec"
      InsigSiteSpec$Slope <- NA
      InsigSiteSpec$Significant <- NA
      finaldata <- rbind(finaldata, InsigSiteSpec)
    }
    
    #Add metadata
    IntervalLength <- unique(sapply(1:(length(IntValues)-1), function(x) IntValues[x+1]-IntValues[x]))  
    finaldata$IntervalLength <- IntervalLength
    finaldata$StartYear <- IntValues[1]
    finaldata$NumYears <- length(IntValues)
    finaldata$Model <- modelfam
    finaldata$SamplingType <- "Interval"
    return(finaldata)
  }, mc.cores=ncores))
  write.csv(Interval, file=paste0(ResultsGLMFP, "Intervals/Intervals_", modelfam, ".csv"), row.names=FALSE)
}

if((file.exists(paste0(ResultsGLMFP, "Summaries/Direction_Interval_", modelfam, ".csv")) &
    file.exists(paste0(ResultsGLMFP, "Summaries/Magnitude_Interval_", modelfam, ".csv")))==FALSE){
  write.csv(NULL, paste0(ResultsGLMFP, "Summaries/Direction_Interval_", modelfam, ".csv"))
  write.csv(NULL, paste0(ResultsGLMFP, "Summaries/Magnitude_Interval_", modelfam, ".csv"))

  Interval <- fread(paste0(ResultsGLMFP, "Intervals/Intervals_", modelfam, ".csv"))
  CompleteTrend <- fread(paste0(ResultsGLMFP, "Consecutive/Consecutive_GLM_30Years_1StartYear_nb.csv"))[,c("SiteSpec", "Significant", "Slope", "EDF", "NumYears")] #Read in the models run on full 30 years of data (for comparison to complete)
  names(CompleteTrend) <- c("SiteSpec", "SignificantCompleteLength", "CompleteLengthSlope" , "EDFCompleteLength", "CompleteLengthNum")
  Interval <- merge(Interval, CompleteTrend, by="SiteSpec") #Merge complete 30 years with the interval runs
  Interval$NumYears <- paste0(Interval$IntervalLength, ".", Interval$NumYears) #Pull these two together to enable the intervals to run using the functions
  
  Direction_All <- rbindlist(pblapply(c(0,5,10,15), function(x) Direction(Interval, Standardise=FALSE, GenLength=x, EDFNum=0, Significance = 0.05)))
  Direction_St <- rbindlist(pblapply(c(0,5,10,15), function(x) Direction(Interval, Standardise=TRUE, GenLength=x, EDFNum=0, Significance = 0.05)))
  Direction_EDF <- rbindlist(pblapply(c(3,5,7,9), function(x) Direction(Interval, Standardise=FALSE, GenLength=0, EDFNum=x, Significance = 0.05)))

  write.csv(rbind(Direction_All, Direction_St, Direction_EDF), file=paste0(ResultsGLMFP, "Summaries/Direction_Interval_", modelfam, ".csv"), row.names=FALSE)
  
  Magnitude_All <- rbindlist(pblapply(c(0,5,10,15), function(x) Magnitude(Interval, Standardise=FALSE, GenLength=x, EDFNum=0, Significance = 0.05)))
  Magnitude_St <- rbindlist(pblapply(c(0,5,10,15), function(x) Magnitude(Interval, Standardise=TRUE, GenLength=x, EDFNum=0, Significance = 0.05)))
  Magnitude_EDF <- rbindlist(pblapply(c(3,5,7,9), function(x) Magnitude(Interval, Standardise=FALSE, GenLength=0, EDFNum=x, Significance = 0.05)))

  write.csv(rbind(Magnitude_All, Magnitude_St, Magnitude_EDF), file=paste0(ResultsGLMFP, "Summaries/Magnitude_Interval_", modelfam, ".csv"), row.names=FALSE)
}

#### Plot #### (this code still to be annotated, stay tuned)
PropOverall <- function(Direction, Category){
  Direction2 <- if(SamplingType=="Consecutive"){split(Direction, Direction$CompleteLength)}else{split(Direction, Direction$NumYears)}
  DataFrame <- rbindlist(lapply(Direction2, function(y){
    okcast <- if(SamplingType=="Consecutive"){
      dcast(y, NumYears~variable, value.var="value")
    } else {
      dcast(y, IntervalLength~variable, value.var="value")
    }
    if(Category=="Correct"){
      okcast$Prop <- sapply(1:nrow(okcast), function(x) okcast$Correct[x]/sum(okcast[x,2:ncol(okcast)]))
      okcast$Category <- "a) Correct"
    } else if(Category=="FA"){
      okcast$Prop <- sapply(1:nrow(okcast), function(x) okcast$'False Alarm'[x]/sum(okcast[x,2:ncol(okcast)]))
      okcast$Category <- "c) False Alarm"
    } else if(Category=="MD"){
      okcast$Prop <- sapply(1:nrow(okcast), function(x) okcast$'Missed Detection'[x]/sum(okcast[x,2:ncol(okcast)]))
      okcast$Category <- "d) Missed Detection"
    } else {
      okcast$Prop <- sapply(1:nrow(okcast), function(x) okcast$Incorrect[x]/sum(okcast[x,2:ncol(okcast)]))
      okcast$Category <- "b) Opposing"
    }
    if(SamplingType=="Consecutive"){
      okcast <- okcast[,c("NumYears", "Prop", "Category")]
      okcast$CompleteLength <- unique(y$CompleteLength)
      okcast$NumYears <- as.numeric(as.character(okcast$NumYears))
    } else {
      okcast <- okcast[,c("IntervalLength", "Prop", "Category")]
      okcast$NumYears <- unique(y$NumYears)
      okcast$IntervalLength <- as.numeric(as.character(okcast$IntervalLength))
    }
    okcast$Prop <- as.numeric(okcast$Prop)
    okcast[is.na(okcast)] <- 0
    return(okcast)
  }))
  return(DataFrame)
} #This function gives us the proportion of corrects for each CompleteLength, without subsetting to a threshold
DirectionPlot <- function(Direction){
  Direction$Prop <- round_any(Direction$Prop*100, 10, floor)
  Direction$Prop <- paste0(Direction$Prop, "-", Direction$Prop+10)
  Direction$Prop <- factor(Direction$Prop, levels=sapply(seq(0,90,10), function(x) paste0(x, "-", x+10)))
  
  if(SamplingType=="Consecutive"){
    xlabel <- "Complete trend length (years)"
    ylabel <- "Sample trend length (years)"
    xbreaks <- (c(seq(4,30,4)))
    ybreaks <- (c(seq(3,30,3)))
    ex <- Direction$CompleteLength
    why <- Direction$NumYears
    aspectrat <- 1
  } else {
    xlabel <- "Number of years sampled"
    ylabel <- "Years between each sample"
    xbreaks <- (c(seq(3,29,5)))
    ybreaks <- (c(seq(1,15,3)))
    aspectrat <- 0.52
    ex <- Direction$NumYears
    why <- Direction$IntervalLength
  }
  
  ggplot(Direction, aes(x=ex,y=why))+
    geom_tile(aes(fill=Prop), colour = "grey40",size = 0.01)+
    facet_wrap(~Category, ncol=2)+
    ylab(ylabel)+
    xlab(xlabel)+
    scale_fill_manual(values = viridis(n=10), name=paste0("Percentage"), guide=guide_legend(reverse = TRUE, keywidth=unit(3, "mm"), keyheight = unit(3, "mm")), drop=FALSE)+
    scale_x_continuous(breaks=xbreaks, expand = c(0,0.15))+
    scale_y_continuous(breaks=ybreaks, expand = c(0,0.15))+
    #guides(fill=guide_legend(nrow=6, byrow=TRUE))+
    #guide_legend()+
    theme(aspect.ratio=aspectrat, 
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          text = element_text(size=plotfontsize),
          panel.border = element_rect(size = 1, fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0),
          legend.title=element_text(size=plotlegendsize),
          legend.justification = "top")
}
DirectionPlotGen <- function(Direction){
  Direction$Prop <- round_any(Direction$Prop*100, 10, floor)
  Direction$Prop <- paste0(Direction$Prop, "-", Direction$Prop+10)
  Direction$Prop <- factor(Direction$Prop, levels=sapply(seq(0,90,10), function(x) paste0(x, "-", x+10)))
  
  if(SamplingType=="Consecutive"){
    xlabel <- "Complete trend length (years)"
    ylabel <- "Sample trend length (years)"
    xbreaks <- (c(seq(4,30,6)))
    ybreaks <- (c(seq(3,30,6)))
    aspectrat <- 1
    ex <- Direction$CompleteLength
    why <- Direction$NumYears
  } else {
    xlabel <- "Number of years sampled"
    ylabel <- "Years between each sample"
    xbreaks <- (c(seq(3,29,6)))
    ybreaks <- (c(seq(1,15,3)))
    aspectrat <- 0.52
    ex <- Direction$NumYears
    why <- Direction$IntervalLength
  }
  Direction$GenLength <- as.character(as.numeric(Direction$GenLength))
  Direction[Direction$GenLength==5,]$GenLength <- "a) Short"
  Direction[Direction$GenLength==10,]$GenLength <- "b) Medium"
  Direction[Direction$GenLength==15,]$GenLength <- "c) Long"
  Direction[Direction$Category=="a) Correct",]$Category <- "Corr."
  Direction[Direction$Category=="c) False Alarm",]$Category <- "FA"
  Direction[Direction$Category=="d) Missed Detection",]$Category <- "MD"
  Direction[Direction$Category=="b) Opposing",]$Category <- "Opp."
  Direction$GenLength <- factor(Direction$GenLength)
  Direction$GenLength <- factor(Direction$GenLength,levels(Direction$GenLength)[c(1,2,3)])
  ggplot(Direction, aes(x=ex,y=why))+
    geom_tile(aes(fill=Prop), colour = "grey40",size = 0.01)+
    facet_grid(Category~GenLength)+
    ylab(ylabel)+
    xlab(xlabel)+
    scale_fill_manual(values = viridis(n=10), name=paste0("Percentage"), guide=guide_legend(reverse = TRUE, keywidth=unit(3, "mm"), keyheight = unit(3, "mm")), drop=FALSE)+
    scale_x_continuous(breaks=xbreaks, expand = c(0,0.15))+
    scale_y_continuous(breaks=ybreaks, expand = c(0,0.15))+
    theme(aspect.ratio=aspectrat, 
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          text = element_text(size=plotfontsize),
          panel.border = element_rect(size = 1, fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0),
          legend.title=element_text(size=plotlegendsize),
          legend.justification = "top")
}
DirectionPlotEDF <- function(Direction){
  Direction$Prop <- round_any(Direction$Prop*100, 10, floor)
  Direction$Prop <- paste0(Direction$Prop, "-", Direction$Prop+10)
  Direction$Prop <- factor(Direction$Prop, levels=sapply(seq(0,90,10), function(x) paste0(x, "-", x+10)))
  
  if(SamplingType=="Consecutive"){
    xlabel <- "Complete trend length (years)"
    ylabel <- "Sample trend length (years)"
    xbreaks <- (c(seq(4,30,6)))
    ybreaks <- (c(seq(3,30,6)))
    aspectrat <- 1
    ex <- Direction$CompleteLength
    why <- Direction$NumYears
  } else {
    xlabel <- "Number of years sampled"
    ylabel <- "Years between each sample"
    xbreaks <- (c(seq(3,29,6)))
    ybreaks <- (c(seq(1,15,3)))
    aspectrat <- 0.52
    ex <- Direction$NumYears
    why <- Direction$IntervalLength
  }
  Direction[Direction$Category=="a) Correct",]$Category <- "Corr."
  Direction[Direction$Category=="c) False Alarm",]$Category <- "FA"
  Direction[Direction$Category=="d) Missed Detection",]$Category <- "MD"
  Direction[Direction$Category=="b) Opposing",]$Category <- "Opp."
  Direction$EDF <- as.character(as.numeric(Direction$EDF))
  Direction[Direction$EDF==3,]$EDF <- "a) 1-3"
  Direction[Direction$EDF==5,]$EDF <- "b) 3-5"
  Direction[Direction$EDF==7,]$EDF <- "c) 5-7"
  Direction[Direction$EDF==9,]$EDF <- "d) 7-9"
  Direction$EDF <- factor(Direction$EDF)
  Direction$EDF <- factor(Direction$EDF,levels(Direction$EDF)[c(1,2,3,4)])
  ggplot(Direction, aes(x=ex,y=why))+
    geom_tile(aes(fill=Prop), colour = "grey40",size = 0.01)+
    facet_grid(Category~EDF)+
    ylab(ylabel)+
    xlab(xlabel)+
    scale_fill_manual(values = viridis(n=10), name=paste0("Percentage"), guide=guide_legend(reverse = TRUE, keywidth=unit(3, "mm"), keyheight = unit(3, "mm")), drop=FALSE)+
    scale_x_continuous(breaks=xbreaks, expand = c(0,0.15))+
    scale_y_continuous(breaks=ybreaks, expand = c(0,0.15))+
    theme(aspect.ratio=aspectrat, 
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          text = element_text(size=plotfontsize),
          panel.border = element_rect(size = 1, fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0),
          legend.title=element_text(size=plotlegendsize),
          legend.justification = "top")
}

MagnitudePlot <- function(Direction){
  MagnitudeSub <- c(0.01,0.025, 0.05, 0.1, 0.25, 0.5)
  Direction <- Direction[Direction$Tolerance %in% MagnitudeSub,]
  Direction$Tolerance <- paste0("± ", Direction$Tolerance)
  if(SamplingType=="Consecutive"){
    Direction$NumYearsProp <- (Direction$NumYears/Direction$CompleteLength)*100
    Direction$CompleteLength <- as.character(as.numeric(Direction$CompleteLength))
  }
  if(SamplingType=="Interval"){
    Direction$CompleteLength <- Direction$IntervalLength
  }
  Direction$Tolerance <- as.factor(Direction$Tolerance)
  if(SamplingType=="Consecutive"){
    Direction$CompleteLength <- sapply(Direction$CompleteLength, function(x){
      if(x==5){
        paste0("a) ", x, " years")
      } else if(x==10){
        paste0("b) ", x, " years")
      } else if(x==20){
        paste0("c) ", x, " years")
      } else if(x==30){
        paste0("d) ", x, " years")
      }
    })
    Direction$CompleteLength <- factor(Direction$CompleteLength, levels=c("a) 5 years","b) 10 years","c) 20 years","d) 30 years"))
  } else {
    Direction$CompleteLength <- sapply(Direction$CompleteLength, function(x){
      if(x==1){
        paste0("a) ", x, " year between samples")
      } else if(x==3){
        paste0("b) ", x, " years between samples")
      } else if(x==5){
        paste0("c) ", x, " years between samples")
      } else if(x==7){
        paste0("d) ", x, " years between samples")
      } else if(x==9){
        paste0("e) ", x, " years between samples")
      } else if(x==11){
        paste0("f) ", x, " years between samples")
      }
    })
    Direction$CompleteLength <- factor(Direction$CompleteLength, levels=c("a) 1 year between samples","b) 3 years between samples","c) 5 years between samples","d) 7 years between samples","e) 9 years between samples","f) 11 years between samples"))
    
  }
  if(SamplingType=="Consecutive"){
    breaksx <- seq(0,30,5)
    breaksy <- seq(0,100,10)
  } else {
    breaksx <- seq(3,15,2)
    breaksy <- seq(0,100,10) 
  }
  Direction$PropCorrect <- Direction$PropCorrect*100
  ggplot(Direction, aes(x=NumYears,y=PropCorrect, colour=Tolerance))+
    geom_line()+
    geom_point(shape=18, size=2)+
    facet_wrap(~CompleteLength, nrow=2)+
    ylab("Percentage within tolerance")+
    xlab("Number of years sampled")+
    scale_colour_manual(values = viridis(n=7), name="Tolerance", guide=guide_legend(reverse = TRUE, keywidth=unit(3, "mm"), keyheight = unit(3, "mm")), drop=FALSE)+
    scale_x_continuous(breaks=breaksx)+ #, expand = c(0,0)
    scale_y_continuous(breaks=breaksy)+ #, expand = c(0,0)
    theme(aspect.ratio=1, 
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          text = element_text(size=plotfontsize),
          panel.border = element_rect(size = 1, fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0),
          legend.title=element_text(size=plotlegendsize),
          legend.key=element_blank(),
          legend.justification = "top")
}
MagnitudePlotGen <- function(Direction){
  MagnitudeSub <- c(0.01,0.05, 0.1)
  Direction <- Direction[Direction$Tolerance %in% MagnitudeSub,]
  Direction$Tolerance <- paste0("± ", Direction$Tolerance)
  Direction$GenLength <- as.character(as.numeric(Direction$GenLength))
  Direction[Direction$GenLength==5,]$GenLength <- "a) Short"
  Direction[Direction$GenLength==10,]$GenLength <- "b) Medium"
  Direction[Direction$GenLength==15,]$GenLength <- "c) Long"
  Direction$GenLength <- factor(Direction$GenLength)
  Direction$GenLength <- factor(Direction$GenLength,levels(Direction$GenLength)[c(1,2,3)])
  if(SamplingType=="Consecutive"){
    Direction$NumYearsProp <- (Direction$NumYears/Direction$CompleteLength)*100
    Direction$CompleteLength <- as.character(as.numeric(Direction$CompleteLength))
  }
  if(SamplingType=="Interval"){
    Direction$CompleteLength <- Direction$IntervalLength
  }
  Direction$Tolerance <- as.factor(Direction$Tolerance)
  if(SamplingType=="Consecutive"){
    Direction$CompleteLength <- paste0(Direction$CompleteLength, " Years")
    Direction$CompleteLength <- factor(Direction$CompleteLength, levels=c("10 Years","20 Years","30 Years"))
  } else {
    Direction$CompleteLength <- sapply(Direction$CompleteLength, function(x){
      if(x==1){
        paste0(x, " yr btwn smpls")
      } else {
        paste0(x, " yrs btwn smpls")
      }
    })
  }
  breaksy <- seq(0,100,20)
  breaksx <- seq(0,30,5)
  Direction$PropCorrect <- Direction$PropCorrect*100
  ggplot(Direction, aes(x=NumYears,y=PropCorrect, colour=GenLength))+
    geom_line()+
    geom_point(shape=18, size=1.5)+
    facet_grid(Tolerance~CompleteLength)+
    ylab("Percentage within tolerance")+
    xlab("Number of years sampled")+
    scale_colour_manual(values = viridis(n=4), name=str_wrap("Generation Length", width=10), guide=guide_legend(reverse = FALSE, keywidth=unit(3, "mm"), keyheight = unit(3, "mm")), drop=FALSE)+
    scale_x_continuous(breaks=breaksx)+ #, expand = c(0,0)
    scale_y_continuous(breaks=breaksy)+ #, expand = c(0,0)
    theme(aspect.ratio=1, 
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          text = element_text(size=plotfontsize),
          panel.border = element_rect(size = 1, fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0),
          legend.title=element_text(size=plotlegendsize),
          legend.key=element_blank(),
          legend.justification = "top")
}
MagnitudePlotEDF <- function(Direction){
  MagnitudeSub <- c(0.01,0.05, 0.1)
  Direction <- Direction[Direction$Tolerance %in% MagnitudeSub,]
  Direction$Tolerance <- paste0("± ", Direction$Tolerance)
  Direction$EDF <- as.character(as.numeric(Direction$EDF))
  Direction[Direction$EDF==3,]$EDF <- "1-3"
  Direction[Direction$EDF==5,]$EDF <- "3-5"
  Direction[Direction$EDF==7,]$EDF <- "5-7"
  Direction[Direction$EDF==9,]$EDF <- "7-9"
  Direction$EDF <- factor(Direction$EDF)
  Direction$EDF <- factor(Direction$EDF,levels(Direction$EDF)[c(1,2,3,4)])
  if(SamplingType=="Consecutive"){
    Direction$NumYearsProp <- (Direction$NumYears/Direction$CompleteLength)*100
    Direction$CompleteLength <- as.character(as.numeric(Direction$CompleteLength))
  }
  if(SamplingType=="Interval"){
    Direction$CompleteLength <- Direction$IntervalLength
  }
  Direction$Tolerance <- as.factor(Direction$Tolerance)
  if(SamplingType=="Consecutive"){
    Direction$CompleteLength <- paste0(Direction$CompleteLength, " years")
    Direction$CompleteLength <- factor(Direction$CompleteLength, levels=c("10 years","20 years","30 years"))
  } else {
    Direction$CompleteLength <- sapply(Direction$CompleteLength, function(x){
      if(x==1){
        paste0(x, " yr btwn smpls")
      } else {
        paste0(x, " yrs btwn smpls")
      }
    })
  }
  breaksy <- seq(0,100,20)
  breaksx <- seq(0,30,5)
  Direction$PropCorrect <- Direction$PropCorrect*100
  ggplot(Direction, aes(x=NumYears,y=PropCorrect, colour=EDF))+
    geom_line()+
    geom_point(shape=18, size=1.5)+
    facet_grid(Tolerance~CompleteLength)+
    ylab("Percentage within tolerance")+
    xlab("Number of years sampled")+
    scale_colour_manual(values = viridis(n=5), name=str_wrap("EDF", width=10), guide=guide_legend(reverse = FALSE, keywidth=unit(3, "mm"), keyheight = unit(3, "mm")), drop=FALSE)+
    scale_x_continuous(breaks=breaksx)+ #, expand = c(0,0)
    scale_y_continuous(breaks=breaksy)+ #, expand = c(0,0)
    theme(aspect.ratio=1, 
          panel.background = element_blank(),
          panel.grid = element_blank(), 
          text = element_text(size=plotfontsize),
          panel.border = element_rect(size = 1, fill = NA),
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0),
          legend.title=element_text(size=plotlegendsize),
          legend.key=element_blank(),
          legend.justification = "top")
}

modelfam <- "nb"

DirectionSummary <- rbindlist(lapply(list.files(path=paste0(ResultsGLMFP, "Summaries/"), pattern=paste0("^Direction_Consecutive_",modelfam, "..*.csv$"), full.names=TRUE), fread), fill=TRUE)
MagnitudeSummary <- rbindlist(lapply(list.files(path=paste0(ResultsGLMFP, "Summaries/"), pattern=paste0("^Magnitude_Consecutive_",modelfam, "..*.csv$"), full.names=TRUE), fread), fill=TRUE)

DirectionSummaryInterval <- read.csv(file=paste0(ResultsGLMFP, "Summaries/Direction_Interval_", modelfam, ".csv"))
DirectionSummary30 <- subset(DirectionSummary, CompleteLength==30)
DirectionSummary30$IntervalLength <- 1
DirectionSummary30$CompleteLength <- NULL
DirectionSummaryInterval <- rbind(DirectionSummary30, DirectionSummaryInterval)

MagnitudeSummaryInterval <- read.csv(file=paste0(ResultsGLMFP, "Summaries/Magnitude_Interval_", modelfam, ".csv"))
MagnitudeSummary30 <- subset(MagnitudeSummary, CompleteLength==30)
MagnitudeSummary30$IntervalLength <- 1
MagnitudeSummaryInterval <- rbind(MagnitudeSummary30, MagnitudeSummaryInterval)

##Direction
for(SamplingType in c("Consecutive", "Interval")){
  for(Casttype in c("All", "St")){
    if(SamplingType=="Consecutive"){
      Direction <- subset(DirectionSummary, CastType==Casttype)
    } else {
      Direction <- subset(DirectionSummaryInterval, CastType==Casttype)
    }
    
    #All
    DirectionAll <- subset(Direction, GenLength==0 & EDF==0)
    DirectionPlots <- rbindlist(pblapply(c("Correct","Incorrect", "MD","FA"), function(Category) PropOverall(DirectionAll, Category)))
    DirectionPlots$Category <- factor(DirectionPlots$Category, levels = c("a) Correct", "b) Opposing", "c) False Alarm", "d) Missed Detection"))
    if(SamplingType=="Consecutive"){
      W <- 400
      H <- 350
    } else {
      W <- 450
      H <- 240
    }
    ggsave(paste0(FiguresFP,"Direction/Direction_", SamplingType, "_", Casttype, "_", modelfam, ".pdf"), DirectionPlot(DirectionPlots), device="pdf", width = W*ResFactor, height = H*ResFactor, units = "mm") #Save as an image

  }
  
  if(SamplingType=="Consecutive"){
    Direction <- subset(DirectionSummary, CastType=="All")
  } else {
    Direction <- subset(DirectionSummaryInterval, CastType=="All")
  }
  
  #Generation length
  DirectionGen <- subset(Direction, GenLength!=0 & EDF==0)
  DirectionPlotsGen <- rbindlist(lapply(c(5,10,15), function(Gen){
    gendata <- subset(DirectionGen, GenLength==Gen)
    PropGen <- rbindlist(lapply(c("Correct","Incorrect", "MD","FA"), function(Category) PropOverall(gendata, Category)))
    PropGen$GenLength <- Gen
    PropGen$Prop <- as.numeric(PropGen$Prop)
    return(PropGen)
  }))
  if(SamplingType=="Consecutive"){
    W <- 500
    H <- 460
  } else {
    W <- 500
    H <- 280
  }
  ggsave(paste0(FiguresFP,"Gen/DirectionGen_",SamplingType, "_", Casttype, "_", modelfam, ".pdf"), DirectionPlotGen(DirectionPlotsGen), device="pdf", width = W*ResFactor, height = H*ResFactor, units = "mm") #Save as an image
  
  #EDF
  DirectionEDF <- subset(Direction, GenLength==0 & EDF!=0)
  DirectionPlotsEDF <- rbindlist(lapply(unique(DirectionEDF$EDF), function(edf){
    edfdata <- subset(DirectionEDF, EDF==edf)
    PropEDF <- rbindlist(lapply(c("Correct","Incorrect", "MD","FA"), function(Category) PropOverall(edfdata, Category)))
    PropEDF$EDF <- edf
    PropEDF$Prop <- as.numeric(PropEDF$Prop)
    return(PropEDF)
  }))
  if(SamplingType=="Consecutive"){
    W <- 500
    H <- 460
  } else {
    W <- 500
    H <- 280
  }
  ggsave(paste0(FiguresFP,"EDF/Direction_EDF_", SamplingType, "_", modelfam, ".pdf"), DirectionPlotEDF(DirectionPlotsEDF), device="pdf", width = W*ResFactor, height = H*ResFactor, units = "mm") #Save as an image
}

##Magnitude
for(SamplingType in c("Consecutive", "Interval")){
  for(Casttype in c("All", "St")){
    if(SamplingType=="Consecutive"){
      Magnitude <- subset(MagnitudeSummary, CastType==Casttype)
    } else {
      Magnitude <- subset(MagnitudeSummaryInterval, CastType==Casttype)
    }
    Magnitude$TotalSamples <- Magnitude$Correct + Magnitude$Incorrect
    Magnitude$PropCorrect <- factor(round(Magnitude$Correct/(Magnitude$Correct+Magnitude$Incorrect), 1)*100)
    MagnitudeAll <- subset(Magnitude, GenLength==0 & EDF==0)
    if(Casttype=="All"){MagnitudeGen <- subset(Magnitude, GenLength!=0 & EDF==0)}
    if(Casttype=="All"){MagnitudeEDF <- subset(Magnitude, GenLength==0 & EDF!=0)}

    if(SamplingType=="Consecutive"){
      MagnitudeAll <- MagnitudeAll[MagnitudeAll$CompleteLength %in% c(5,10,20,30),]
      if(Casttype=="All"){MagnitudeGen <- MagnitudeGen[MagnitudeGen$CompleteLength %in% c(10,20,30),]}
      if(Casttype=="All"){MagnitudeEDF <- MagnitudeEDF[MagnitudeEDF$CompleteLength %in% c(10,20,30),]}
    } else {
      MagnitudeAll <- MagnitudeAll[MagnitudeAll$IntervalLength %in% c(1,3,5,7,9,11),]
      MagnitudeAll <- subset(MagnitudeAll, NumYears<16)
      if(Casttype=="All"){MagnitudeGen <- MagnitudeGen[MagnitudeGen$IntervalLength %in% c(1,3,5),]}
      if(Casttype=="All"){MagnitudeEDF <- MagnitudeEDF[MagnitudeEDF$IntervalLength %in% c(1,3,5,7),]}
    }
    MagnitudeAll$PropCorrect <- MagnitudeAll$Correct/(MagnitudeAll$Correct+MagnitudeAll$Incorrect)
    if(Casttype=="All"){MagnitudeGen$PropCorrect <- MagnitudeGen$Correct/(MagnitudeGen$Correct+MagnitudeGen$Incorrect)}
    if(Casttype=="All"){MagnitudeEDF$PropCorrect <- MagnitudeEDF$Correct/(MagnitudeEDF$Correct+MagnitudeEDF$Incorrect)}

    if(SamplingType=="Consecutive"){
      H <- 400
      W <- 450
      H_GenEDF <- 400
      W_GenEDF <- 500
    } else {
      H <- 350
      W <- 500
      H_GenEDF <- 400
      W_GenEDF <- 500
    }
    
    ggsave(paste0(FiguresFP,"Magnitude/Magnitude_", SamplingType, "_", Casttype, "_", modelfam, ".pdf"), MagnitudePlot(MagnitudeAll), device="pdf", width = W*ResFactor, height = H*ResFactor, units = "mm")
    if(Casttype=="All"){ggsave(paste0(FiguresFP,"Gen/MagnitudeGen_", SamplingType, "_", Casttype, "_", modelfam, ".pdf"), MagnitudePlotGen(MagnitudeGen), device="pdf", width = W_GenEDF*ResFactor, height = H_GenEDF*ResFactor, units = "mm")}
    if(Casttype=="All"){ggsave(paste0(FiguresFP,"EDF/MagnitudeEDF_", SamplingType, "_", Casttype, "_", modelfam, ".pdf"), MagnitudePlotEDF(MagnitudeEDF), device="pdf", width = W_GenEDF*ResFactor, height = H_GenEDF*ResFactor, units = "mm")}
  }
}

#### Paper stats ####
modelfam <- "nb"
SamplingType<- "Consecutive"
Casttype <- "All"

AllModels <- read.csv(paste0(ResultsFP, "GLM/Consecutive/Consecutive_GLM_30Years_1StartYear_nb.csv"))
AllModels <- AllModels[!is.na(AllModels$Slope),]
AllModelsSub <- subset(AllModels, Slope<1 & Slope>-1)
nrow(AllModelsSub)/nrow(AllModels)
nrow(subset(AllModels, EDF<1.0001))/nrow(AllModels)
nrow(subset(AllModels, Significant<=0.05))/nrow(AllModels)
mean(AllModels$Slope)
max(AllModels$Slope)
min(AllModels$Slope)

p1 <- ggplot(data=AllModelsSub)+
  geom_freqpoly(aes(x=Slope), binwidth=0.01)+
  scale_x_continuous(breaks=c(seq(-1,1,0.5)), expand = c(0,0))+
  scale_y_continuous(limits=c(0, 2000), expand = c(0,0))+
  ggtitle("a)")+
  xlab("Population Growth Rate (r)")+
  ylab("Frequency")+
  theme(aspect.ratio=1, 
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size=plotfontsize),
        panel.border = element_rect(size = 1, fill = NA),
        plot.title = element_text(size=plotfontsize))
p2 <- ggplot(data=AllModels)+
  geom_freqpoly(aes(x=EDF), binwidth=0.5)+
  scale_x_continuous(breaks=c(0,2,4,6,8,10),expand = c(0,0))+
  scale_y_continuous(limits=c(0, 6000), expand = c(0,0))+
  xlab("Estimated Degrees of Freedom")+
  ylab("Frequency")+
  ggtitle("b)")+
  theme(aspect.ratio=1, 
        panel.background = element_blank(),
        panel.grid = element_blank(), 
        text = element_text(size=plotfontsize),
        panel.border = element_rect(size = 1, fill = NA),
        plot.title = element_text(size=plotfontsize))


library(gridExtra)
plotty <- grid.arrange(p1,p2, ncol=2)

W <- 600
H <- 375

ggsave(paste0(FiguresFP,"SummaryStats.pdf"), plotty, device="pdf", width = W*ResFactor, height = H*ResFactor, units = "mm") #Save as an image

PercentSignif <- nrow(subset(AllModels, Significant<0.05))/nrow(AllModels)

TYB_listStats <- rbindlist(TYB_list)
TYB_listStats$Site <- str_split_fixed(TYB_listStats$SiteSpec, "[_]",2)[,1]
TYB_listStats$SpecPop <- str_split_fixed(TYB_listStats$SiteSpec, "[_]",2)[,2]
TYB_listStats$Spec <- str_split_fixed(TYB_listStats$SiteSpec, "[_]",3)[,2]

length(unique(TYB_listStats$Site))
length(unique(TYB_listStats$SpecPop))
length(unique(TYB_listStats$Spec))

#### Build supplementary data
PropOverallSuppData <- function(Direction, Category){
  Direction2 <- if(SamplingType=="Consecutive"){split(Direction, Direction$CompleteLength)}else{split(Direction, Direction$NumYears)}
  DataFrame <- rbindlist(lapply(Direction2, function(y){
    okcast <- if(SamplingType=="Consecutive"){
      dcast(y, NumYears~variable, value.var="value")
    } else {
      dcast(y, IntervalLength~variable, value.var="value")
    }
    if(Category=="Correct"){
      okcast$Prop <- sapply(1:nrow(okcast), function(x) okcast$Correct[x]/sum(okcast[x,2:ncol(okcast)]))
      okcast$Category <- "Correct"
    } else if(Category=="FA"){
      okcast$Prop <- sapply(1:nrow(okcast), function(x) okcast$'False Alarm'[x]/sum(okcast[x,2:ncol(okcast)]))
      okcast$Category <- "FA"
    } else if(Category=="MD"){
      okcast$Prop <- sapply(1:nrow(okcast), function(x) okcast$'Missed Detection'[x]/sum(okcast[x,2:ncol(okcast)]))
      okcast$Category <- "MD"
    } else {
      okcast$Prop <- sapply(1:nrow(okcast), function(x) okcast$Incorrect[x]/sum(okcast[x,2:ncol(okcast)]))
      okcast$Category <- "Opposing"
    }
    if(SamplingType=="Consecutive"){
      okcast <- okcast[,c("NumYears", "Prop", "Category")]
      okcast$CompleteLength <- unique(y$CompleteLength)
      okcast$NumYears <- as.numeric(as.character(okcast$NumYears))
    } else {
      okcast <- okcast[,c("IntervalLength", "Prop", "Category")]
      okcast$NumYears <- unique(y$NumYears)
      okcast$IntervalLength <- as.numeric(as.character(okcast$IntervalLength))
    }
    okcast$Prop <- as.numeric(okcast$Prop)
    okcast[is.na(okcast)] <- 0
    return(okcast)
  }))
  return(DataFrame)
} #This function gives us the proportion of corrects for each CompleteLength, without subsetting to a threshold

modelfam <- "nb"

DirectionSummary <- rbindlist(lapply(list.files(path=paste0(ResultsGLMFP, "Summaries/"), pattern=paste0("^Direction_Consecutive_",modelfam, "..*.csv$"), full.names=TRUE), fread), fill=TRUE)
MagnitudeSummary <- rbindlist(lapply(list.files(path=paste0(ResultsGLMFP, "Summaries/"), pattern=paste0("^Magnitude_Consecutive_",modelfam, "..*.csv$"), full.names=TRUE), fread), fill=TRUE)

DirectionSummaryInterval <- read.csv(file=paste0(ResultsGLMFP, "Summaries/Direction_Interval_", modelfam, ".csv"))
DirectionSummary30 <- subset(DirectionSummary, CompleteLength==30)
DirectionSummary30$IntervalLength <- 1
DirectionSummary30$CompleteLength <- NULL
DirectionSummaryInterval <- rbind(DirectionSummary30, DirectionSummaryInterval)

MagnitudeSummaryInterval <- read.csv(file=paste0(ResultsGLMFP, "Summaries/Magnitude_Interval_", modelfam, ".csv"))
MagnitudeSummary30 <- subset(MagnitudeSummary, CompleteLength==30)
MagnitudeSummary30$IntervalLength <- 1
MagnitudeSummaryInterval <- rbind(MagnitudeSummary30, MagnitudeSummaryInterval)

for (SamplingType in c("Consecutive", "Interval")){
  if(SamplingType=="Consecutive"){
    Direction <- subset(DirectionSummary, CastType=="All")
  } else {
    Direction <- subset(DirectionSummaryInterval, CastType=="All")
  }
  
  DirectionAll <- subset(Direction, GenLength==0 & EDF==0)
  DirectionPlots <- rbindlist(pblapply(c("Correct","Incorrect", "MD","FA"), function(Category) PropOverallSuppData(DirectionAll, Category)))
  
  DirectionGen <- subset(Direction, GenLength!=0 & EDF==0)
  DirectionPlotsGen <- rbindlist(lapply(c(5,10,15), function(Gen){
    gendata <- subset(DirectionGen, GenLength==Gen)
    PropGen <- rbindlist(lapply(c("Correct","Incorrect", "MD","FA"), function(Category) PropOverallSuppData(gendata, Category)))
    PropGen$GenLength <- Gen
    PropGen$Prop <- as.numeric(PropGen$Prop)
    return(PropGen)
  }))
  
  DirectionEDF <- subset(Direction, GenLength==0 & EDF!=0)
  DirectionPlotsEDF <- rbindlist(lapply(unique(DirectionEDF$EDF), function(edf){
    edfdata <- subset(DirectionEDF, EDF==edf)
    PropEDF <- rbindlist(lapply(c("Correct","Incorrect", "MD","FA"), function(Category) PropOverallSuppData(edfdata, Category)))
    PropEDF$EDF <- edf
    PropEDF$Prop <- as.numeric(PropEDF$Prop)
    return(PropEDF)
  }))
  
  DirectionPlotsGen$EDF <- NA
  DirectionPlotsEDF$GenLength <- NA
  DirectionPlots[,c("EDF", "GenLength")] <- NA
  
  DirectionData <- rbind(DirectionPlots, DirectionPlotsEDF, DirectionPlotsGen)
  if(SamplingType=="Consecutive"){
    names(DirectionData)[1] <- "SampleLength"
  } else {
    names(MagnitudeData)[1] <- "SampleLength"
    names(DirectionData)[4] <- "NumYearsInSample"
    DirectionData$CompleteLength <- 30
  }
  
  write.csv(DirectionData, paste0(FiguresFP, "Supplementary/DirectionData",SamplingType, ".csv"), row.names=FALSE)
  
  if(SamplingType=="Consecutive"){
    Magnitude <- subset(MagnitudeSummary, CastType=="All")
  } else {
    Magnitude <- subset(MagnitudeSummaryInterval, CastType=="All")
  }
  Magnitude$TotalSamples <- Magnitude$Correct + Magnitude$Incorrect
  Magnitude$PropCorrect <- factor(round(Magnitude$Correct/(Magnitude$Correct+Magnitude$Incorrect), 1)*100)
  MagnitudeAll <- subset(Magnitude, GenLength==0 & EDF==0)
  MagnitudeGen <- subset(Magnitude, GenLength!=0 & EDF==0)
  MagnitudeEDF <- subset(Magnitude, GenLength==0 & EDF!=0)
  
  MagnitudeAll$PropCorrect <- MagnitudeAll$Correct/(MagnitudeAll$Correct+MagnitudeAll$Incorrect)
  MagnitudeGen$PropCorrect <- MagnitudeGen$Correct/(MagnitudeGen$Correct+MagnitudeGen$Incorrect)
  MagnitudeEDF$PropCorrect <- MagnitudeEDF$Correct/(MagnitudeEDF$Correct+MagnitudeEDF$Incorrect)
  
  MagnitudeAll[,c("Correct", "Incorrect","SamplingType", "Model", "CastType", "TotalSamples")] <- NULL
  MagnitudeAll[,c("GenLength", "EDF")] <- NA
  MagnitudeGen[,c("Correct", "Incorrect", "SamplingType","Model", "CastType", "TotalSamples")] <- NULL
  MagnitudeGen[,c("EDF")] <- NA
  MagnitudeEDF[,c("Correct", "Incorrect", "SamplingType","Model", "CastType", "TotalSamples")] <- NULL
  MagnitudeEDF[,c("GenLength")] <- NA
  
  MagnitudeData <- rbind(MagnitudeAll, MagnitudeGen, MagnitudeEDF)
  if(SamplingType=="Consecutive"){
    names(MagnitudeData)[1] <- "SampleLength"
  } else {
    names(MagnitudeData)[1] <- "SampleLength"
    names(MagnitudeData)[7] <- "NumYearsInSample"
  }
  names(MagnitudeData)[ncol(MagnitudeData)] <- "PropWithinTolerance"
  write.csv(MagnitudeData, paste0(FiguresFP, "Supplementary/MagnitudeData",SamplingType, ".csv"), row.names=FALSE)
  
}

