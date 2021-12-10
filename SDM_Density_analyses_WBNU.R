#################################################################################################
# Chapter 4: SDM and density modeling of current and 1850s veg
#################################################################################################

## Tyler Hallman
## 10/10/18

## Summary: Take prepped detection and count files for each species to model 
## suitable habitat and then density then use the habitat suitability and 
## density models to predict current and past distributions and populations



#########################################################
# 1. Get user input
#########################################################


Working.directory <-"D:\\Chapter4_Data"

###### SDM Parameters ######
# occurrence data file
occurrence.file <- "wbnu_Chapter4.csv"

# Minimum values for splitting data
minimum.occur.training <- 14
minimum.occur.test <- 6
minimum.split <- 25

# gbm parameters
bf.zi = 0.75 
tc.zi = 3 
cores = 2 # Number of cores to use for this analysis
use.calibration = 1 # This is either 1 or 0. 0 means don't calibrate.
weights.zi = 0 # Weights should be 0 for equal or 1 for sqrt count +1 (Johnston et al. 2015)
# truncation.zi = 1200 # This needs to happen with original file and not here, because otherwise counts and occurrence will be inaccurate
radius.zi = "multi"

# Calibration parameters
manual.threshold = 0 # 0 indicates use prevalence and 1 indicates use manual threshold. 
threshold.input = 0.2 # The manual threshold you would like to use to indicate suitable habitat.

# Density parameters
observations.file = "Willamette_Observations.csv" # This is the file with the observations (each row is an individual bird). This is from the oregon 2020 database

# This should be the same as the occurrence file since I'm already using all the covariates for that
observations.covariates.file = "wbnu_Chapter4.csv" 

TruncationD = 125 # This would be the truncation distance decided upon if there is one.

# Raster folders
current.raster.folder <- "D:\\Chapter4_GLO_Crosswalking\\CurrentProcessed_clipped"
historic.raster.folder <- "D:\\Chapter4_GLO_Crosswalking\\HistoricProcessed"

# Extent Shapefile
willamette.extent.shp <- "D:\\Chapter4_Data\\ArcMap\\WillametteExtent_UTM10.shp"


#########################################################
# 2. Load packages
#########################################################

###### Load the necessary packages ######
if (!require(gbm)){
  install.packages("gbm", repos="http://cran.rstudio.com/")
  library(gbm)
}
if (!require(caret)){
  install.packages("caret", repos="http://cran.rstudio.com/")
  library(caret)
}
if (!require(scam)){
  install.packages("scam", repos="http://cran.rstudio.com/")
  library(scam)
}
if (!require(dismo)){
  install.packages("dismo", repos="http://cran.rstudio.com/")
  library(dismo)
}
if (!require(devtools)){
  install.packages("devtools", repos="http://cran.rstudio.com/")
  library(devtools)
}
if (!require(detect)){
  install.packages("detect", repos="http://cran.rstudio.com/")
  library(detect)
}
if (!require(MuMIn)){
  install.packages("MuMIn", repos="http://cran.rstudio.com/")
  library(MuMIn)
}
if (!require(reshape)){
  install.packages("reshape", repos="http://cran.rstudio.com/")
  library(reshape)
}
if (!require(ape)){
  install.packages("ape", repos="http://cran.rstudio.com/")
  library(ape)
}
if (!require(raster)){
  install.packages("raster", repos="http://cran.rstudio.com/")
  library(raster)
}
if (!require(unmarked)){
  install.packages("unmarked", repos="http://cran.rstudio.com/")
  library(unmarked)
}
if (!require(AHMbook)){
  install.packages("AHMbook", repos="http://cran.rstudio.com/")
  library(AHMbook)
}


#########################################################
# 3. Set working directory and create output directory
#########################################################

###### Set working directory ###### 
setwd(Working.directory) 

###### Set other miscellaneous inputs ###### 
file.elements <- strsplit(occurrence.file, "_")[[1]]
species <- file.elements[1]

###### Create a unique model name based on the above inputs for output file ###### 
model.name <- paste("Chapter4_", species, "_UnmarkedDensity", sep="")

output.directory <- file.path(Working.directory, model.name)

# Create output directory
dir.create(output.directory)


#########################################################
# 4. Load and prep data
#########################################################

###### Load data ###### 
occurrence.data <- read.csv(file=file.path(Working.directory, occurrence.file),header=TRUE,sep=",",dec=".", fill=TRUE)
observation.data <- read.csv(file=file.path(Working.directory, observations.file),header=TRUE,sep=",",dec=".", fill=TRUE)
observation.covariates.data <- read.csv(file=file.path(Working.directory, observations.covariates.file),header=TRUE,sep=",",dec=".", fill=TRUE)

####### Set default population and suitable habitat ###### 
#population <- "NA"
#suitable.habitat.area.hectares <- "NA"

###### Remove variables with no variation or all NAs ######
ncol(occurrence.data)
occurrence.data <- occurrence.data[,colSums(is.na(occurrence.data))<nrow(occurrence.data)-20]
ncol(occurrence.data)

occurrence.data <- occurrence.data[,sapply(occurrence.data, function(x) length(unique(na.omit(x)))) > 1]
ncol(occurrence.data)

print(paste("Total Individuals Detected in Dataset Before Split: ", sum(occurrence.data$Count), sep=""))
print(paste("Total Occupied Sites in Dataset Before Split: ", sum(occurrence.data$Occur), sep=""))

###### Change pertinent columns to factors and match site ID levels between loaded files ###### 
## I may not have to do this. I have to look ats the files that I used in unmarked

observation.covariates.data$Year <- as.factor(observation.covariates.data$Year)
#observation.covariates.data$Observer <- as.factor(observation.covariates.data$Observer)
observation.covariates.data$NewSiteID <- as.factor(observation.covariates.data$NewSiteID)
observation.covariates.data <- observation.covariates.data[order(observation.covariates.data$NewSiteID),]

lvls <- levels(observation.covariates.data$NewSiteID)
observation.data$NewSiteID <- as.factor(observation.data$NewSiteID)
observation.data$NewSiteID <- factor(observation.data$NewSiteID, levels = lvls)
observation.data$Distance <- as.numeric(observation.data$Distance)
nrow(observation.data)
observation.data <- na.omit(observation.data)
nrow(observation.data)
observation.data <- subset(observation.data, Species == species)
nrow(observation.data)

###### Truncate the observation data ###### 
plot(hist(observation.data$Distance))
distances <- observation.data[order(observation.data$Distance),]
head(distances$Distance)
truncation.proportion <- 0.1*nrow(observation.data); truncation.proportion
max(distances$Distance)

Observation.data.trunc.full <-  observation.data[which(observation.data$Distance <= TruncationD), ]
nrow(Observation.data.trunc.full)
plot(hist(Observation.data.trunc.full$Distance))
length(unique(Observation.data.trunc.full$NewSiteID))

PercentTuncated <- ((nrow(observation.data)-nrow(Observation.data.trunc.full))/nrow(observation.data))*100; PercentTuncated

TruncationD.avail <- 60
Observation.data.trunc.avail <-  observation.data[which(observation.data$Distance <= TruncationD.avail), ]
nrow(Observation.data.trunc.avail)
plot(hist(Observation.data.trunc.avail$Distance))
length(unique(Observation.data.trunc.avail$NewSiteID))



#########################################################
# 5. Run ZI gbm 
#########################################################

###### Create a new data frame of just the columns for this analysis ######
model.data <- occurrence.data

## Remove the columns that won't be helpful for this species
will_cols <- grep("will", names(model.data), value = FALSE)
urban_cols <- grep("urban", names(model.data), value = FALSE)
disturbed_cols <- grep("disturbed", names(model.data), value = FALSE)
burned_cols <- grep("burned", names(model.data), value = FALSE)
vlmr_cols <- grep("vlmr", names(model.data), value = FALSE)
oak_cols <- grep("oak", names(model.data), value = FALSE)




model.data <- model.data[,-c(will_cols, urban_cols, disturbed_cols, burned_cols, vlmr_cols, oak_cols)]
names(model.data)


###### Split data ###### 
if (sum(model.data$Occur) >= minimum.split){
  counts_check_train <- 0
  counts_check_test <- 0
  sufficient.occupancies <- "na"
  
  while(counts_check_train == 0 || counts_check_test == 0){
    inTraining <- createDataPartition(model.data$Occur, p = .8, list = FALSE)
    training.data <- model.data[ inTraining,]
    test.data  <- model.data[-inTraining,]
    
    # How many sites with greater than 0 abundance should there be? 10?
    if(sum(training.data$Occur) > minimum.occur.training){
      counts_check_train <- counts_check_train + 1
    } 
    if(sum(test.data$Occur) > minimum.occur.test){
      counts_check_test <- counts_check_test + 1
    } 
  }
} else {
  training.data <- model.data
  test.data <- model.data
  sufficient.occupancies <- "Insufficient detected occupancies to split dataset"
  print("Insufficient detected occupancies to split dataset")
}

total.sites.occupied.brt <- sum(model.data$Occur)
total.prevalence <- (sum(model.data$Occur)/nrow(model.data))*100
training.sites.occupied.brt <- sum(training.data$Occur)
test.sites.occupied.brt <- sum(test.data$Occur)

###### Create groups of variables to be used in models ######
rad165_cols <- grep("165", names(training.data), value = FALSE)
rad315_cols <- grep("315", names(training.data), value = FALSE)
rad615_cols <- grep("615", names(training.data), value = FALSE)
rad1215_cols <- grep("1215", names(training.data), value = FALSE)
riv_cols <- grep("Riv", names(training.data), value = FALSE)
zi_col <- grep("Occur", names(training.data), value = FALSE) 
count_col <- which(colnames(training.data)=="Count")
survey_cols <- c(which(colnames(training.data)=="Julian_Dat"), which(colnames(training.data)=="Year"),
                 which(colnames(training.data)=="Minutes_si"))

###### Prep radius columns ######
# In this case, the args document must also include which radius is the best model to be used here. 
if(radius.zi=="165"){
  rad_cols <- rad165_cols
} else if(radius.zi=="315"){
  rad_cols <- rad315_cols
} else if(radius.zi=="615"){
  rad_cols <- rad615_cols
} else if(radius.zi=="1215"){
  rad_cols <- rad1215_cols
} else if(radius.zi=="multi"){
  rad_cols <- c(rad165_cols, rad315_cols, rad615_cols, rad1215_cols)
}

###### Create a new data frame of just the columns for this analysis ######
training.data.prepped <- training.data[,c(rad_cols, zi_col, count_col)]

## Remove and variables that don't have enough nonzero data
training.data.prepped <- training.data.prepped[,colSums(is.na(training.data.prepped))<nrow(training.data.prepped)-20]
training.data.prepped <- training.data.prepped[,sapply(training.data.prepped, function(x) length(unique(na.omit(x)))) > 1]

## Make sure the rest.data have the same column names
test.data.prepped <- subset(test.data, select = names(training.data.prepped))


###### Prep weights ######
if(weights.zi==0){
  weights.zi.prepped <- rep(1,nrow(training.data.prepped))
} else if(weights.zi==1){
  weights.zi.prepped <- sqrt(training.data.prepped$Count)+1
}

###### Remove the count column so that the dot can be used ######
count_col_zi <- grep("Count", names(training.data.prepped), value = FALSE)
training.data.prepped <- training.data.prepped[,-count_col_zi]

###### Subset the observation covariates and the independent test data to allow for calibration cutoff ######
# This is a subset of the observation covariate dataset so it can exclude counts from unsuitable habitats
observation.covariates.data.cutoff <- subset(observation.covariates.data, select= names(training.data.prepped))


###### Run the gbm ######
lr = 0.002 # This is what I'm going to have to optimize as I go. Will take a couple tries.

gbmModel.zi = gbm(formula = Occur ~ .,
                  distribution = "bernoulli",
                  data = training.data.prepped,
                  n.trees = 5000,
                  shrinkage = lr,
                  cv.folds = 10,
                  n.minobsinnode = 20,
                  bag.fraction = bf.zi,
                  weights = weights.zi.prepped,
                  interaction.depth = tc.zi,
                  n.cores = cores)


###### Add probability of occurrence predictions to the full dataset ######
trees.number <- gbm.perf(gbmModel.zi); trees.number

gbmPredictions.zi = predict(object = gbmModel.zi,
                            newdata = training.data.prepped, 
                            n.trees = gbm.perf(gbmModel.zi),
                            type = "response")

training.data.prepped$Pred.Occur <- gbmPredictions.zi

###### Calculate AUC on withheld test data ######
gbmPredictions.test = predict(object = gbmModel.zi,
                              newdata = test.data.prepped, 
                              n.trees = gbm.perf(gbmModel.zi),
                              type = "response")

AUC <- gbm.roc.area(test.data.prepped$Occur, gbmPredictions.test); AUC

###### Export locations test file for partial ROC calculation ######
test.occur <- subset(test.data, Occur == 1)
species.test <- rep(species, nrow(test.occur))
test.coords <- as.data.frame(cbind(species.test, test.occur$Longitude, test.occur$Latitude))
names(test.coords) <- c("Species", "long", "lat")
write.csv(test.coords, file.path(output.directory, paste(model.name, "_TestOccur.csv", sep="")))

###### Export covariate influence information and figures ######
sum.model.zi <- summary(gbmModel.zi); sum.model.zi
write.csv(sum.model.zi, file.path(output.directory, paste(model.name,"_Variables_zi.csv", sep="")))

for(i in 1:nrow(sum.model.zi)) { 
  png(filename=file.path(output.directory, paste( model.name, "zi_variables_", i, ".png", sep="")), 
      type="cairo",
      units="in", 
      width=10, 
      height=14, 
      pointsize=12, 
      res=600)
  
  plot(gbmModel.zi, i.var = i)
  
  dev.off()
}



#########################################################
# 6. Run calibration of threshold
#########################################################

###### Calculate the prevalence to be used as a threshold ###### 
prevalence <- sum(training.data.prepped$Occur)/nrow(training.data.prepped)

###### Make the threshold the species prevalence or user input ######
if (manual.threshold==1){
  Threshold <- threshold.input
} else {
  Threshold <- prevalence
}

if(use.calibration==0){
  calibration_cutoff <- Threshold
} else if(use.calibration==1){
  
  ###### Create new dataframe just for calibration ######
  pred.zi <- gbmPredictions.zi
  obs.zi <- training.data.prepped$Occur
  pred.zi.seq <- seq(0, 1, by=0.01)
  calib.mod.data <- data.frame(obs.zi, pred.zi, weights.zi.prepped)
  
  ###### Run a gaussian gam on dataframe ######
  calib.mod <- tryCatch(scam(obs.zi ~ s(pred.zi, k=6, bs="mpi"), weights=weights.zi.prepped, gamma=2, data=calib.mod.data),
                        error=function(e) e
  )
  
  if(!inherits(calib.mod, "error")){
    calib.pred<-predict(calib.mod,data.frame(pred.zi=pred.zi.seq))
  }
  # If the model didn't run, try with reduced degrees of freedom:
  if(inherits(calib.mod, "error")) {
    calib.mod <- tryCatch(scam(obs.zi ~ s(pred.zi, k=4, bs="mpi"), weights=weights.calib, gamma=2, data=calib.mod.data),
                          error=function(e) e
    )
    # If the reduced model ran, force a monotonically increasing spline:
    if(!inherits(calib.mod, "error")){
      calib.pred<-predict(calib.mod,data.frame(pred.zi=pred.zi.seq))
    }
    if(inherits(calib.mod, "error")) {
      g <- ceiling(runif(1, 1, 1000))
      print(paste("random no:", g))
      write.table(calib.mod.data, paste("calib.failed.data.",g, ".txt", sep=""), row.names=F)
      print(paste("Calibration model failed. Dataset written to file ", getwd(), "/calib.failed.data.txt", sep=""))
    }
  }
  
  ###### Cap predictions at 0 and 1, required for gausssian model ######
  calib.pred <- apply(cbind(calib.pred, rep(1, length(calib.pred))), 1, min)
  calib.pred <- apply(cbind(calib.pred, rep(0, length(calib.pred))), 1, max)
  
  ###### Interpolate between two values nearest to threshhold ######
  
  # First check for values exactly equal to threshhold
  # As a few strange fits might cross threshhold twice, this code first finds the maximum and then
  # truncates the fit at the maximum value and therefore finds the value nearest threshhold which
  # occurs before the maximum
  max.pred <- which(calib.pred==max(calib.pred))[1]
  calib.pred.trunc <- calib.pred[1:max.pred]
  pred.zi.seq.trunc <- pred.zi.seq[1:max.pred]
  
  if(any(calib.pred.trunc==Threshold)) {
    calibration_cutoff <- pred.zi.seq.trunc[which(calib.pred.trunc==Threshold)]
  } else {
    if(all(calib.pred.trunc<Threshold)) {
      calibration_cutoff <- 1
    } else {
      if(all(calib.pred.trunc>Threshold)) {
        calibration_cutoff <- 0
      } else {
        pre_Threshold <- max(which(calib.pred.trunc < Threshold))
        calibration_cutoff <- pred.zi.seq.trunc[pre_Threshold] + (pred.zi.seq.trunc[pre_Threshold+1] - pred.zi.seq.trunc[pre_Threshold])*(Threshold-calib.pred.trunc[pre_Threshold])/(calib.pred.trunc[pre_Threshold+1] - calib.pred.trunc[pre_Threshold])
      }
    }
  }
}

###### Select only data that fit within the calibration threshhold ######
print(paste("Desired Threshold: ", Threshold))
print(paste("Calibration Cutoff: ", calibration_cutoff))



#########################################################
# 7. Use suitability cutoff for detection and abundance data
#########################################################

###### Now calculate the calibration cutoff for the Oregon 2020 covariates dataset ######
# This is necessary so that the detection analyses are only run on counts from areas of
# suitable habitat.

# Get gbm predictions for oregon 2020 covariate dataset
gbmPredictions.observation.covariates.data.cutoff = predict(object = gbmModel.zi,
                                                            newdata = observation.covariates.data.cutoff, 
                                                            n.trees = gbm.perf(gbmModel.zi),
                                                            type = "response")

# Add the predictions to the oregon 2020 dataset
observation.covariates.data$Pred.Occur <- gbmPredictions.observation.covariates.data.cutoff

# Use calibration cutoff to get rid of counts from unsuitable habitat
rows.cutoff <- which(observation.covariates.data$Pred.Occur>calibration_cutoff)
observation.covariates.data <- observation.covariates.data[rows.cutoff,]

# This allows you to investigate changes in density dataset due to cutoff
nrow(subset(observation.covariates.data, Pred.Occur > 0.5))
sum(subset(observation.covariates.data, Pred.Occur > 0.5)$Occur)


###### Check how the dataset in suitable habitat differs ######
print(paste("Original Count: ", sum(occurrence.data$Count)))
print(paste("Original Occupied Sites: ", sum(occurrence.data$Occur))) 
print(paste("Origal number of rows: ", nrow(occurrence.data)))

print(paste("New Total Count: ", sum(observation.covariates.data$Count)))
print(paste("New Occupied Sites: ", sum(observation.covariates.data$Occur)))
print(paste("New number of rows: ", nrow(observation.covariates.data)))


###### Remove the corresponding observation/detection data ######
nrow(Observation.data.trunc.full)
Observation.data.trunc.full <- subset(Observation.data.trunc.full, NewSiteID %in% observation.covariates.data$NewSiteID)
nrow(Observation.data.trunc.full) # This will differ slightly from the numbers above due to truncation

###### Need to change the levels of the site ID to match new dataset ######
str(observation.covariates.data$NewSiteID)
observation.covariates.data$NewSiteID <- as.numeric(observation.covariates.data$NewSiteID)
str(observation.covariates.data$NewSiteID)
observation.covariates.data$NewSiteID <- as.factor(observation.covariates.data$NewSiteID)
str(observation.covariates.data$NewSiteID)

lvls <- levels(observation.covariates.data$NewSiteID)

str(Observation.data.trunc.full$NewSiteID)
Observation.data.trunc.full$NewSiteID <- factor(Observation.data.trunc.full$NewSiteID, levels = lvls)
str(Observation.data.trunc.full$NewSiteID)




#########################################################
# 8. Prep detection data
#########################################################

###### Observation Data with Distances ###### 
## Get max distance and break data into bins
B <- max(Observation.data.trunc.full$Distance); B # Could set B to 300 to base off of distances instead of specific data
nD <- 6 # May want to use fewer bins
delta <- B/nD; delta

distBreaks <- seq(0,B,delta); distBreaks

## Create a new simpler dataframe to convert to the correct format
Obs.dist <- data.frame(Observation.data.trunc.full$NewSiteID, Observation.data.trunc.full$Distance)
names(Obs.dist) <- c("NewSiteID", "Distance")
names(Obs.dist)

## Convert distance data into the right format for analyses
yDat <- formatDistData(Obs.dist, distCol="Distance",
                       transectNameCol="NewSiteID", dist.breaks=distBreaks)
nrow(yDat)
str(yDat)
head(yDat)
colSums(yDat)

###### Covariates Data ###### 
## Choose Covariates
names(observation.covariates.data)

sum(observation.covariates.data$Occur)

## Check if there should still be all levels for observer
table(observation.covariates.data$Observer)
observation.covariates.data$Observer <- as.character(observation.covariates.data$Observer)
observation.covariates.data$Observer <- as.factor(observation.covariates.data$Observer)
table(observation.covariates.data$Observer)


# Scale specific variables to model more easily
RivDist_min <- min(observation.covariates.data$will_RivDist)
RivDist <- observation.covariates.data$will_RivDist - RivDist_min
RivDist_sd <- sd(RivDist)
observation.covariates.data$RivDist <- RivDist / RivDist_sd

observation.covariates.data$RivDistLog <- log1p(observation.covariates.data$will_RivDist)

HwyDist_min <- min(observation.covariates.data$will_HwyDist)
HwyDist <- observation.covariates.data$will_HwyDist - HwyDist_min
HwyDist_sd <- sd(HwyDist)
observation.covariates.data$HwyDist <- HwyDist / HwyDist_sd

observation.covariates.data$HwyDistLog <- log1p(observation.covariates.data$will_HwyDist)

covs <- observation.covariates.data[,c("Julian_Dat", "Year", "Observer", "Minutes_si", "urban_highmed_utm10_165", # First add detection covariates
                                       "urban_highmed_utm10_315", "urban_tot_utm10_165", "urban_tot_utm10_315", "will_CANCOV__UTM10_165", 
                                       "will_TPH_GE_3__UTM10_165", "RivDist", "RivDistLog", "HwyDist", "HwyDistLog",
                                       "riparian2531_utm10_165", "riparian2531_utm10_315", "riparian2531_utm10_615", "riparian2531_utm10_1215", 
                                       "npow16_utm10_165", "npow16_utm10_315", "npow16_utm10_615", "npow16_utm10_1215", 
                                       "wvups22_utm10_165", "wvups22_utm10_315", "wvups22_utm10_615", "wvups22_utm10_1215", 
                                       "wtf_utm10_165", "wtf_utm10_315", "wtf_utm10_615", "wtf_utm10_1215", 
                                       "devos_utm10_165", "devos_utm10_315", "devos_utm10_615", "devos_utm10_1215",
                                       "Pred.Occur")] 

JulDat_min <- min(covs$Julian_Dat)
JulDat <- covs$Julian_Dat - JulDat_min
JulDat_sd <- sd(JulDat)
covs$JulDat <- JulDat / JulDat_sd

Time_min <- min(covs$Minutes_si)
Time <- covs$Minutes_si - Time_min
Time_sd <- sd(Time)
covs$Time <- Time / Time_sd


#########################################################
# 9. Run Hierarchical Distance Sampling Models (distsamp: Poisson)
#########################################################

## Convert to an unmarked dataframe for analyses
umf.dist <- unmarkedFrameDS(y=as.matrix(yDat), siteCovs=data.frame(covs), dist.breaks=distBreaks, 
                       unitsIn="m", survey="point")

## If necessary, we could scale the covariates here to increase fit and decrease model
## run time. Is this necessary though? 
#sc <- siteCovs(umf.dist)
#sc.s <- scale(sc)
#siteCovs(umf.dist) <- sc.s

names(covs)

summary(umf.dist)

## Make a list to store the models
dist.models <- list()

## Run Distance Models and add them to list. 
###### Detection portion ###### 
dist.models$null <- distsamp(~1 ~1, umf.dist, keyfun="halfnorm", output="density")
summary(dist.models$null)

dist.models$observer <- distsamp(~Observer ~1, umf.dist, keyfun="halfnorm", output="density")
summary(dist.models$observer)

dist.models$date <- distsamp(~JulDat ~1, umf.dist, keyfun="halfnorm", output="density")
summary(dist.models$date)

dist.models$time <- distsamp(~Time ~1, umf.dist, keyfun="halfnorm", output="density")
summary(dist.models$time)

dist.models$observer_date <- distsamp(~Observer + JulDat ~1, umf.dist, keyfun="halfnorm", output="density")
summary(dist.models$observer_date)

dist.models$observer_time <- distsamp(~Observer + Time ~1, umf.dist, keyfun="halfnorm", output="density")
summary(dist.models$observer_time)

dist.models$date_time <- distsamp(~JulDat + Time ~1, umf.dist, keyfun="halfnorm", output="density")
summary(dist.models$date_time)

dist.models$observer_date_time <- distsamp(~Observer + JulDat + Time ~1, umf.dist, keyfun="halfnorm", output="density")
summary(dist.models$observer_date_time)


## Check which of these make any difference
dist.MSmodels <- modSel(fitList(fits=dist.models));dist.MSmodels


## Unused variables because they are confounded with their effects on density
#dist.models$hwy <- distsamp(~HwyDist ~1, umf.dist, keyfun="halfnorm", output="density")
#summary(dist.models$hwy)

#dist.models$hwy2 <- distsamp(~HwyDist + I(HwyDist^2) ~1, umf.dist, keyfun="halfnorm", output="density")
#summary(dist.models$hwy2)

#dist.models$hwylog <- distsamp(~HwyDistLog ~1, umf.dist, keyfun="halfnorm", output="density")
#summary(dist.models$hwylog)

#dist.models$cancov <- distsamp(~will_CANCOV__UTM10_165 ~1, umf.dist, keyfun="halfnorm", output="density") # confounded as below
#summary(dist.models$cancov)

#dist.models$treedensity <- distsamp(~will_TPH_GE_3__UTM10_165 ~1, umf.dist, keyfun="halfnorm", output="density") # confounded as below
#summary(dist.models$treedensity)

#dist.models$urban165 <- distsamp(~urban_tot_utm10_165 ~1, umf.dist, keyfun="halfnorm", output="density")
#summary(dist.models$urban165)

#dist.models$urban315 <- distsamp(~urban_tot_utm10_315 ~1, umf.dist, keyfun="halfnorm", output="density")
#summary(dist.models$urban315)

#dist.models$riv <- distsamp(~RivDist ~1, umf.dist, keyfun="halfnorm", output="density") # Looks like this is confounded with them being in these habitats. As is it says that detection probability decreases as you get further from rivers. This is untrue.
#summary(dist.models$riv)

#dist.models$riv2 <- distsamp(~RivDist + I(RivDist^2) ~1, umf.dist, keyfun="halfnorm", output="density")
#summary(dist.models$riv2)

#dist.models$rivlog <- distsamp(~RivDistLog ~1, umf.dist, keyfun="halfnorm", output="density")
#summary(dist.models$rivlog)


###### Density portion with best detection structure ###### 
names(covs)

## Starting with riparian because it was influential above
dist.models$observerXriparian165 <- update(dist.models$observer, ~. ~riparian2531_utm10_165) # Update the best detection model with density covariates
summary(dist.models$observerXriparian165)

dist.models$observerXriparian315 <- update(dist.models$observer, ~. ~riparian2531_utm10_315) # Update the best detection model with density covariates
summary(dist.models$observerXriparian315)

dist.models$observerXriparian615 <- update(dist.models$observer, ~. ~riparian2531_utm10_615, starts = c(-2.52, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXriparian615)

dist.models$observerXriparian1215 <- update(dist.models$observer, ~. ~riparian2531_utm10_1215) # Update the best detection model with density covariates
summary(dist.models$observerXriparian1215)


## Now going to try npow
dist.models$observerXnpow165 <- update(dist.models$observer, ~. ~npow16_utm10_165, starts = c(-2.52, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXnpow165)

dist.models$observerXnpow315 <- update(dist.models$observer, ~. ~npow16_utm10_315, starts = c(-2.52, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXnpow315)

dist.models$observerXnpow615 <- update(dist.models$observer, ~. ~npow16_utm10_615, starts = c(-2.52, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXnpow615)

dist.models$observerXnpow1215 <- update(dist.models$observer, ~. ~npow16_utm10_1215, starts = c(-2.52, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXnpow1215)


## Now going to try wvups
dist.models$observerXwvups165 <- update(dist.models$observer, ~. ~wvups22_utm10_165, starts = c(-2.52, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXwvups165)

dist.models$observerXwvups315 <- update(dist.models$observer, ~. ~wvups22_utm10_315, starts = c(-2.52, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXwvups315)

dist.models$observerXwvups615 <- update(dist.models$observer, ~. ~wvups22_utm10_615, starts = c(-2.52, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXwvups615)

dist.models$observerXwvups1215 <- update(dist.models$observer, ~. ~wvups22_utm10_1215, starts = c(-2.52, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXwvups1215)


## wtf
dist.models$observerXwtf165 <- update(dist.models$observer, ~. ~wtf_utm10_165, starts = c(-2.52, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXwtf165)

dist.models$observerXwtf315 <- update(dist.models$observer, ~. ~wtf_utm10_315, starts = c(-2.52, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXwtf315)

dist.models$observerXwtf615 <- update(dist.models$observer, ~. ~wtf_utm10_615, starts = c(-2.52, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXwtf615)

dist.models$observerXwtf1215 <- update(dist.models$observer, ~. ~wtf_utm10_1215, starts = c(-2.52, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXwtf1215)


## devos
dist.models$observerXdevos165 <- update(dist.models$observer, ~. ~devos_utm10_165, starts = c(-2.52, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXdevos165)

dist.models$observerXdevos315 <- update(dist.models$observer, ~. ~devos_utm10_315, starts = c(-2.52, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXdevos315)

dist.models$observerXdevos615 <- update(dist.models$observer, ~. ~devos_utm10_615, starts = c(-2.52, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXdevos615)

dist.models$observerXdevos1215 <- update(dist.models$observer, ~. ~devos_utm10_1215, starts = c(-2.52, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXdevos1215)


## Check which radius is best for this variable when also accounting for detection
dist.MSmodels <- modSel(fitList(fits=dist.models));dist.MSmodels # Have a top radius for this variable as well


## Now try with the predicted habitat suitability
dist.models$observerXpred.occur <- update(dist.models$observer, ~. ~Pred.Occur, starts = c(-2.52, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXpred.occur)

dist.models$observerXpred.occur2 <- update(dist.models$observer, ~. ~Pred.Occur + I(Pred.Occur^2), starts = c(-2.52, 0, 0, -0.0318, 3.9909, 3.9843, 3.7751)) # Update the best detection model with density covariates
summary(dist.models$observerXpred.occur2)


# Now check if the other detection structures perform better with this top abundance structure
dist.models$nullXpred.occur2 <- update(dist.models$null, ~. ~Pred.Occur + I(Pred.Occur^2)) # Update the best detection model with density covariates
summary(dist.models$nullXpred.occur2)

dist.models$observerXpred.occur2 <- update(dist.models$observer, ~. ~Pred.Occur + I(Pred.Occur^2)) # Update the best detection model with density covariates
summary(dist.models$observerXpred.occur2)

dist.models$dateXpred.occur2 <- update(dist.models$date, ~. ~Pred.Occur + I(Pred.Occur^2)) # Update the best detection model with density covariates
summary(dist.models$dateXpred.occur2)

dist.models$timeXpred.occur2 <- update(dist.models$time, ~. ~Pred.Occur + I(Pred.Occur^2)) # Update the best detection model with density covariates
summary(dist.models$timeXpred.occur2)

dist.models$observer_dateXpred.occur2 <- update(dist.models$observer_date, ~. ~Pred.Occur + I(Pred.Occur^2), starts = c(-4.28, 18.41, -19.09, -0.275, 4.226, 4.4698, 4.2805, 0.0591)) # Update the best detection model with density covariates
summary(dist.models$observer_dateXpred.occur2) 

dist.models$observer_timeXpred.occur2 <- update(dist.models$observer_time, ~. ~Pred.Occur + I(Pred.Occur^2), starts = c(-2.94, 0.0, 0.0, -0.5626, 4.4124, 4.4698, 4.2805, 0.0591)) # Update the best detection model with density covariates
summary(dist.models$observer_timeXpred.occur2)

dist.models$date_timeXpred.occur2 <- update(dist.models$date_time, ~. ~Pred.Occur + I(Pred.Occur^2), starts = c(-4.52, 24.81, -34.02, 3.7546, 0.00462, 0.02000)) # Update the best detection model with density covariates
summary(dist.models$date_timeXpred.occur2)

dist.models$observer_date_timeXpred.occur2 <- update(dist.models$observer_date_time, ~. ~Pred.Occur + I(Pred.Occur^2), starts = c(-2.94, 0.0, 0.0, -0.5626, 4.4124, 4.4698, 4.2805, 0.0, 0.0591)) # Update the best detection model with density covariates
summary(dist.models$observer_date_timeXpred.occur2)


## Check if predicted suitability is a better metric
dist.MSmodels <- modSel(fitList(fits=dist.models));dist.MSmodels


###### Check to see if model predictions make sense ######
## I am particularly concerned that the range of covariate values in the historic data may not match the current data on which
## my models are built. Unrealistic predictions could be made based on this. Need to either use a lower radius model so that the 
## range of percent land cover values is more similar between historic and current models, or check to see if predictions are
## even very different and then use the top model if they aren't.
summary(umf.dist)

## Make sequence of the variable affecting density
Pred.Occur.seq <- seq(0,1,0.1)

## Now predict with top model
newdat <- data.frame(Pred.Occur = Pred.Occur.seq)
density.pred.pred.occur2 <- predict(dist.models$observerXpred.occur2, type = "state", newdata = newdat, appendData=TRUE); density.pred.pred.occur2


jpeg(file.path(output.directory, paste( model.name, '_distsamp_predoccur2.jpg', sep="")))

plot(Pred.Occur.seq, density.pred.pred.occur2[, "Predicted"], xlab="Predicted Occurrence Squared (SDM)", ylab="Predicted density", type='l', 
     ylim = c(0,1), frame=F, lwd=2)
matlines(Pred.Occur.seq, density.pred.pred.occur2[,3:4], lty=1, col="grey", lwd=1)

dev.off()



###### Model Table ######
## Make a fistlist for selection table
dist.MSmodels <- modSel(fitList(fits=dist.models));dist.MSmodels

## Export Model Table
toExport <- as(dist.MSmodels, "data.frame")
write.csv(toExport, file.path(output.directory, paste(model.name, "_DistTable.csv", sep="")))

###### Best Model ######
## Check the goodness of fit of the best model and export summary
# Goodness of fit
pb.dist <- parboot(dist.models$observerXpred.occur2, fitstats, nsim=100, report=5); pb.dist
c.hat.dist <- pb.dist@t0[2]/mean(pb.dist@t.star[,2]); c.hat.dist

filename <- "distsamp_parboot.txt"
capture.output(pb.dist, file = file.path(output.directory,filename))

# Export best model summary
filename <- "distsamp_model_summary.txt"
capture.output(summary(dist.models$observerXpred.occur2), file = file.path(output.directory,filename))

# Detection probability for all bins within the 300m radius
p.bins <- getP(dist.models$observerXpred.occur2)
head(p.bins) # This shows that the detection probability increases then decreases. This isn't right. I may have to change the bins to fix this. 

p.distsamp <- mean(rowSums(getP(dist.models$observerXpred.occur2))); p.distsamp

test <- predict(dist.models$observerXpred.occur2, type = "det", newdata = umf.dist, appendData=FALSE)
head(test) # These are the sigmas that affect the slope of the detection function

## Plot the detection function and probability density function based on the sigma from the null
sigma.distsamp <- backTransform(dist.models$null, type="det")
sigma.distsamp <- coef(sigma.distsamp)

jpeg(file.path(output.directory, paste( model.name, '_distsamp_detfun.jpg', sep="")))

plot(function(x) gxhn(x, sigma=sigma.distsamp), 0, TruncationD, xlab="Distance (m)",
     ylab="Detection prob", cex.lab=0.7,
     cex.axis=0.7, las=1)

dev.off()

jpeg(file.path(output.directory, paste( model.name, '_distsamp_probdensfun.jpg', sep="")))

plot(function(x) drhn(x, sigma.distsamp), 0, TruncationD, xlab="distance",
     ylab="Probability density", main="Point-transect")

dev.off()

p.distsamp.integration <- integrate(grhn, 0, TruncationD, sigma=sigma.distsamp)$value * 2*pi / (pi*TruncationD^2); p.distsamp.integration


#########################################################
# 10. Run Hierarchical Distance Sampling Models (gdistsamp: Negative Binomial)
#########################################################

## Convert to an unmarked dataframe for analyses
umf.gdist <- unmarkedFrameGDS(y=as.matrix(yDat), siteCovs=data.frame(covs), dist.breaks=distBreaks, 
                        unitsIn="m", survey="point", numPrimary=1)

## If necessary, we could scale the covariates here to increase fit and decrease model
## run time. Is this necessary though? 
#sc <- siteCovs(umf.gdist)
#sc.s <- scale(sc)
#siteCovs(umf.gdist) <- sc.s


## Fit models again with negative binomial distribution
## Make a list to store the models
gdist.models <- list()

## Run Distance Models and add them to list. 
###### Detection portion ###### 
gdist.models$null <- gdistsamp(~1, ~1, ~1, data=umf.gdist, keyfun="halfnorm", output="density", mixture="NB", K=10)
summary(gdist.models$null)

gdist.models$observer <- gdistsamp(~1, ~1, ~Observer, data=umf.gdist, keyfun="halfnorm", output="density", mixture="NB", K=10, starts = c(-2.94, 3.92, 5, 5, 5, -1.14))
summary(gdist.models$observer)

gdist.models$date <- gdistsamp(~1, ~1, ~JulDat, data=umf.gdist, keyfun="halfnorm", output="density", mixture="NB", K=10, starts = c(-1.7, 3.7634, 0.0313, 0.384))
summary(gdist.models$date)

gdist.models$time <- gdistsamp(~1, ~1, ~Time, data=umf.gdist, keyfun="halfnorm", output="density", mixture="NB", K=10)
summary(gdist.models$time)

#gdist.models$observer_date <- gdistsamp(~1, ~1, ~Observer + JulDat, data=umf.gdist, keyfun="halfnorm", output="density", mixture="NB", K=10)
#summary(gdist.models$observer_date)

#gdist.models$observer_time <- gdistsamp(~1, ~1, ~Observer + Time, data=umf.gdist, keyfun="halfnorm", output="density", mixture="NB", K=10)
#summary(gdist.models$observer_time)

gdist.models$date_time <- gdistsamp(~1, ~1, ~JulDat + Time, data=umf.gdist, keyfun="halfnorm", output="density", mixture="NB", K=10)
summary(gdist.models$date_time)

#gdist.models$observer_date_time <- gdistsamp(~1, ~1, ~Observer + JulDat + Time, data=umf.gdist, keyfun="halfnorm", output="density", mixture="NB", K=10)
#summary(gdist.models$observer_date_time)



## Check which of these make any difference
gdist.MSmodels <- modSel(fitList(fits=gdist.models));gdist.MSmodels


## Unused variables because they are confounded with their effects on density
#gdist.models$hwy <- gdistsamp(~1, ~1, ~HwyDist, data=umf.gdist, keyfun="halfnorm", output="density", mixture="NB", K=10)
#summary(gdist.models$hwy)

#gdist.models$hwy2 <- gdistsamp(~1, ~1, ~HwyDist + I(HwyDist^2), data=umf.gdist, keyfun="halfnorm", output="density", mixture="NB", K=10)
#summary(gdist.models$hwy2)

#gdist.models$hwylog <- gdistsamp(~1, ~1, ~HwyDistLog, data=umf.gdist, keyfun="halfnorm", output="density", mixture="NB", K=10, starts= c(-2.47, 4.9868, -0.0959, 0.323))
#summary(gdist.models$hwylog)

#gdist.models$riv <- gdistsamp(~1, ~1, ~RivDist, data=umf.gdist, keyfun="halfnorm", output="density", mixture="NB", K=10, starts= c(-2.49, 4.82692, 0.00984, 0.279))
#summary(gdist.models$riv)

#gdist.models$riv2 <- gdistsamp(~1, ~1, ~RivDist + I(RivDist^2), data=umf.gdist, keyfun="halfnorm", output="density", mixture="NB", K=10, starts= c(-2.47, 4.9868, 0.00984, 0, 0.279))
#summary(gdist.models$riv2)

#gdist.models$rivlog <- gdistsamp(~1, ~1, ~RivDistLog, data=umf.gdist, keyfun="halfnorm", output="density", mixture="NB", K=10, starts = c(-2.47, 4.9868, -0.0959, 0.323))
#summary(gdist.models$rivlog)

#gdist.models$urban165 <- gdistsamp(~1, ~1, ~urban_tot_utm10_165, data=umf.gdist, keyfun="halfnorm", output="density", mixture="NB", K=10)
#summary(gdist.models$urban165)

#gdist.models$urban315 <- gdistsamp(~1, ~1, ~urban_tot_utm10_315, data=umf.gdist, keyfun="halfnorm", output="density", mixture="NB", K=10)
#summary(gdist.models$urban315)

#gdist.models$cancov <- gdistsamp(~1, ~1, ~will_CANCOV__UTM10_165, data=umf.gdist, keyfun="halfnorm", output="density", mixture="NB", K=10)
#summary(gdist.models$cancov)

#gdist.models$treedensity <- gdistsamp(~1, ~1, ~will_TPH_GE_3__UTM10_165, data=umf.gdist, keyfun="halfnorm", output="density", mixture="NB", K=10)
#summary(gdist.models$treedensity)


###### Density portion with best detection structure ###### 
names(covs)

## For wbnu I'm starting with Riparian
gdist.models$nullXriparian165 <- update(gdist.models$null, ~riparian2531_utm10_165, ~., ~.) # Update the best detection model with density covariates
summary(gdist.models$nullXriparian165)

gdist.models$nullXriparian315 <- update(gdist.models$null, ~riparian2531_utm10_315, ~., ~.) # Update the best detection model with density covariates
summary(gdist.models$nullXriparian315)

gdist.models$nullXriparian615 <- update(gdist.models$null, ~riparian2531_utm10_615, ~., ~.) # Update the best detection model with density covariates
summary(gdist.models$nullXriparian615)

gdist.models$nullXriparian1215 <- update(gdist.models$null, ~riparian2531_utm10_1215, ~., ~.) # Update the best detection model with density covariates
summary(gdist.models$nullXriparian1215)

## Now going to try oak
gdist.models$nullXoak165 <- update(gdist.models$null, ~oak_utm10_165, ~., ~.) # Update the best detection model with density covariates
summary(gdist.models$nullXoak165)

gdist.models$nullXoak315 <- update(gdist.models$null, ~oak_utm10_315, ~., ~.) # Update the best detection model with density covariates
summary(gdist.models$nullXoak315)

gdist.models$nullXoak615 <- update(gdist.models$null, ~oak_utm10_615, ~., ~.) # Update the best detection model with density covariates
summary(gdist.models$nullXoak615)

gdist.models$nullXoak1215 <- update(gdist.models$null, ~oak_utm10_1215, ~., ~.) # Update the best detection model with density covariates
summary(gdist.models$nullXoak1215)

## Now going to try npow
gdist.models$nullXnpow165 <- update(gdist.models$null, ~npow16_utm10_165, ~., ~.) # Update the best detection model with density covariates
summary(gdist.models$nullXnpow165)

gdist.models$nullXnpow315 <- update(gdist.models$null, ~npow16_utm10_315, ~., ~.) # Update the best detection model with density covariates
summary(gdist.models$nullXnpow315)

gdist.models$nullXnpow615 <- update(gdist.models$null, ~npow16_utm10_615, ~., ~.) # Update the best detection model with density covariates
summary(gdist.models$nullXnpow615)

gdist.models$nullXnpow1215 <- update(gdist.models$null, ~npow16_utm10_1215, ~., ~.) # Update the best detection model with density covariates
summary(gdist.models$nullXnpow1215)

## Now going to try wvups
gdist.models$nullXwvups165 <- update(gdist.models$null, ~wvups22_utm10_165, ~., ~.) # Update the best detection model with density covariates
summary(gdist.models$nullXwvups165)

gdist.models$nullXwvups315 <- update(gdist.models$null, ~wvups22_utm10_315, ~., ~.) # Update the best detection model with density covariates
summary(gdist.models$nullXwvups315)

gdist.models$nullXwvups615 <- update(gdist.models$null, ~wvups22_utm10_615, ~., ~.) # Update the best detection model with density covariates
summary(gdist.models$nullXwvups615)

gdist.models$nullXwvups1215 <- update(gdist.models$null, ~wvups22_utm10_1215, ~., ~.) # Update the best detection model with density covariates
summary(gdist.models$nullXwvups1215)



## Check which radius is best for this variable when also accounting for detection
gdist.MSmodels <- modSel(fitList(fits=gdist.models));gdist.MSmodels


## Now try with the predicted habitat suitability
gdist.models$nullXpred.occur <- update(gdist.models$null, ~Pred.Occur, ~., ~.) # Update the best detection model with density covariates
summary(gdist.models$nullXpred.occur)

gdist.models$nullXpred.occur2 <- update(gdist.models$null, ~Pred.Occur + I(Pred.Occur^2), ~., ~.) # Update the best detection model with density covariates
summary(gdist.models$nullXpred.occur2)



# Now check if the other detection structures perform better with this top abundance structure
gdist.models$nullXpred.occur2 <- update(gdist.models$null, ~Pred.Occur + I(Pred.Occur^2), ~., ~.) 
summary(gdist.models$nullXpred.occur2)

#gdist.models$observerXpred.occur2 <- update(gdist.models$observer, ~Pred.Occur + I(Pred.Occur^2), ~., ~., starts = c(-4.52, 22.92, -29.81, 3.8096, 0.0372, 6.45)) 
#summary(gdist.models$observerXpred.occur2)

gdist.models$dateXpred.occur2 <- update(gdist.models$date, ~Pred.Occur + I(Pred.Occur^2), ~., ~., starts = c(-4.52, 22.92, -29.81, 3.8096, 0.0372, 6.45)) 
summary(gdist.models$dateXpred.occur2)

gdist.models$timeXpred.occur2 <- update(gdist.models$time, ~Pred.Occur + I(Pred.Occur^2), ~., ~., starts = c(-4.52, 22.92, -29.81, 3.8096, -0.0119, 6.45)) 
summary(gdist.models$timeXpred.occur2)

#gdist.models$observer_dateXpred.occur2 <- update(gdist.models$observer_date, ~Pred.Occur + I(Pred.Occur^2), ~., ~., starts = c(-3.65, 10.87, -8.33, 3.63124, 0.08225, 0.13754, 0.24862, 0.0288, 0.569)) 
#summary(gdist.models$observer_dateXpred.occur2)

#gdist.models$observer_timeXpred.occur2 <- update(gdist.models$observer_time, ~Pred.Occur + I(Pred.Occur^2), ~., ~., starts = c(-3.65, 10.87, -8.33, 3.63124, 0.08225, 0.13754, 0.24862, 0.0485, 0.569)) 
#summary(gdist.models$observer_timeXpred.occur2)

gdist.models$date_timeXpred.occur2 <- update(gdist.models$date_time, ~Pred.Occur + I(Pred.Occur^2), ~., ~., starts = c(-4.52, 22.92, -29.81, 3.8096, 0.0372, -0.012, 7.03)) # Update the best detection model with density covariates
summary(gdist.models$date_timeXpred.occur2)

#gdist.models$observer_date_timeXpred.occur2 <- update(gdist.models$observer_date_time, ~Pred.Occur + I(Pred.Occur^2), ~., ~., starts = c(-3.65, 10.87, -8.33, 3.63124, 0.08225, 0.13754, 0.24862, 0.0288, 0.569)) 
#summary(gdist.models$observer_date_timeXpred.occur2)

## Check if predicted suitability is a better metric
gdist.MSmodels <- modSel(fitList(fits=gdist.models));gdist.MSmodels



###### Check to see if model predictions make sense ######
## I am particularly concerned that the range of covariate values in the historic data may not match the current data on which
## my models are built. Unrealistic predictions could be made based on this. Need to either use a lower radius model so that the 
## range of percent land cover values is more similar between historic and current models, or check to see if predictions are
## even very different and then use the top model if they aren't.
summary(umf.gdist)

## Make sequence of the variable affecting density
Pred.Occur.seq <- seq(0,1,0.1)

## Now predict with top model
newdat <- data.frame(Pred.Occur = Pred.Occur.seq)
density.pred.pred.occur2.gdist <- predict(gdist.models$nullXpred.occur2, type = "lambda", newdata = newdat, appendData=TRUE); density.pred.pred.occur2.gdist


jpeg(file.path(output.directory, paste( model.name, '_gdistsamp_predoccur2.jpg', sep="")))

plot(Pred.Occur.seq, density.pred.pred.occur2.gdist[, "Predicted"], xlab="Predicted Occurrence Squared (SDM)", ylab="Predicted density", type='l', 
     ylim = c(0,1), frame=F, lwd=2)
matlines(Pred.Occur.seq, density.pred.pred.occur2.gdist[,3:4], lty=1, col="grey", lwd=1)

dev.off()



###### Model Table ######
## Make a fistlist for selection table
gdist.MSmodels <- modSel(fitList(fits=gdist.models));gdist.MSmodels

## Export Model Table
toExport <- as(gdist.MSmodels, "data.frame")
write.csv(toExport, file.path(output.directory, paste(model.name, "_GDistTable.csv", sep="")))

###### Best Model ######
## Check the goodness of fit of the best model and export summary
# Goodness of fit
pb.gdist <- parboot(gdist.models$nullXpred.occur2, fitstats, nsim=100, report=5); pb.gdist # Can't do this apparently because of the issues in the model
c.hat.gdist <- pb.gdist@t0[2]/mean(pb.gdist@t.star[,2]); c.hat.gdist

filename <- "gdistsamp_parboot.txt"
capture.output(pb.gdist, file = file.path(output.directory,filename))

# Export best model summary
filename <- "gdistsamp_model_summary.txt"
capture.output(summary(gdist.models$nullXpred.occur2), file = file.path(output.directory,filename))

# Detection probability for all bins within the 300m radius
p.bins <- getP(gdist.models$nullXpred.occur2)
head(p.bins) # This shows that the detection probability increases then decreases. This isn't right. I may have to change the bins to fix this. 

p.gdistsamp <- mean(rowSums(getP(gdist.models$nullXpred.occur2))); p.gdistsamp

test <- predict(gdist.models$nullXpred.occur2, type = "det", newdata = umf, appendData=FALSE)
head(test) # These are the sigmas that affect the slope of the detection function

## Plot the detection function and probability density function based on the sigma from the null
sigma.gdistsamp <- backTransform(gdist.models$null, type="det")
sigma.gdistsamp <- coef(sigma.gdistsamp)

jpeg(file.path(output.directory, paste( model.name, '_gdistsamp_detfun.jpg', sep="")))

plot(function(x) gxhn(x, sigma=sigma.gdistsamp), 0, TruncationD, xlab="Distance (m)",
     ylab="Detection prob", cex.lab=0.7,
     cex.axis=0.7, las=1)

dev.off()

jpeg(file.path(output.directory, paste( model.name, '_gdistsamp_probdensfun.jpg', sep="")))

plot(function(x) drhn(x, sigma.gdistsamp), 0, TruncationD, xlab="distance",
     ylab="Probability density", main="Point-transect")

dev.off()

p.gdistsamp.integration <- integrate(grhn, 0, TruncationD, sigma=sigma.gdistsamp)$value * 2*pi / (pi*TruncationD^2); p.gdistsamp.integration







#########################################################
# 11. Removal Models
#########################################################

###### Observation Data with Removal Interval ###### 
obs.rem <- data.frame(Observation.data.trunc.full$NewSiteID, Observation.data.trunc.full$Removal_In)
names(obs.rem) <- c("NewSiteID", "RemovalInterval")
str(obs.rem)


rem1 <- obs.rem[which(obs.rem$RemovalInterval == 1), ]
rem1.tab <- table(rem1$NewSiteID)

rem2 <- obs.rem[which(obs.rem$RemovalInterval == 2), ]
rem2.tab <- table(rem2$NewSiteID)

rem3 <- obs.rem[which(obs.rem$RemovalInterval == 3), ]
rem3.tab <- table(rem3$NewSiteID)

rem4 <- obs.rem[which(obs.rem$RemovalInterval == 4), ]
rem4.tab <- table(rem4$NewSiteID)

rem5 <- obs.rem[which(obs.rem$RemovalInterval == 5), ]
rem5.tab <- table(rem5$NewSiteID)

removal.matrix <- cbind(rem1.tab, rem2.tab, rem3.tab, rem4.tab, rem5.tab)

max(rem1.tab)
max(rem2.tab)
max(rem3.tab)
max(rem4.tab)
max(rem5.tab)


###### Covariates Data ###### 
## Choose Covariates
names(covs)


###### Analyses ###### 
## Create an unmarked frame for removal sampling
umf.rem <- unmarkedFrameMPois(y=removal.matrix, siteCovs=covs, type="removal")

## Make a list to store the models
Rem.models <- list()

## Run Removal Models and add them to list. 
## Start out with variables affecting detection.
Rem.models$null <- multinomPois(~ 1 ~ 1, umf.rem)
summary(Rem.models$null)

Rem.models$observer <- multinomPois(~ Observer ~ 1, umf.rem)
summary(Rem.models$observer)

Rem.models$date <- multinomPois(~ JulDat ~ 1, umf.rem)
summary(Rem.models$date)

Rem.models$time <- multinomPois(~ Time ~ 1, umf.rem)
summary(Rem.models$time)

Rem.models$observer_date <- multinomPois(~ Observer + JulDat ~ 1, umf.rem)
summary(Rem.models$observer_date)

Rem.models$observer_time <- multinomPois(~ Observer + Time ~ 1, umf.rem)
summary(Rem.models$observer_time)

Rem.models$date_time <- multinomPois(~ JulDat + Time ~ 1, umf.rem)
summary(Rem.models$date_time)

Rem.models$observer_date_time <- multinomPois(~ Observer + JulDat + Time ~ 1, umf.rem)
summary(Rem.models$observer_date_time)

## Check which of these make any difference
MS.Rem.models <- modSel(fitList(fits=Rem.models));MS.Rem.models


## Unused variables because they are confounded with their effects on density
#Rem.models$hwy <- multinomPois(~ HwyDist ~ 1, umf.rem)
#summary(Rem.models$hwy)

#Rem.models$hwy2 <- multinomPois(~ HwyDist + I(HwyDist^2) ~ 1, umf.rem)
#summary(Rem.models$hwy2)

#Rem.models$hwylog <- multinomPois(~ HwyDistLog ~ 1, umf.rem)
#summary(Rem.models$hwylog)

#Rem.models$riv <- multinomPois(~ RivDist ~ 1, umf.rem)
#summary(Rem.models$riv)

#Rem.models$riv2 <- multinomPois(~ RivDist + I(RivDist^2) ~ 1, umf.rem)
#summary(Rem.models$riv)

#Rem.models$rivlog <- multinomPois(~ RivDistLog ~ 1, umf.rem)
#summary(Rem.models$rivlog)

#Rem.models$cancov <- multinomPois(~ will_CANCOV__UTM10_165 ~ 1, umf.rem)
#summary(Rem.models$cancov)

#Rem.models$treedensity <- multinomPois(~ will_TPH_GE_3__UTM10_165 ~ 1, umf.rem, starts = c(-3.91e-16, -1.0933, 0.0))
#summary(Rem.models$treedensity)

#Rem.models$urban165 <- multinomPois(~ urban_tot_utm10_165 ~ 1, umf.rem)
#summary(Rem.models$urban165)

#Rem.models$urban315 <- multinomPois(~ urban_tot_utm10_315 ~ 1, umf.rem)
#summary(Rem.models$urban315)

## Now model abundance within the truncation distance
names(covs)

Rem.models$observer_dateXriparian165 <- multinomPois(~ Observer + JulDat ~ riparian2531_utm10_165, umf.rem)
summary(Rem.models$observer_dateXriparian165)

Rem.models$observer_dateXriparian315 <- multinomPois(~ Observer + JulDat ~ riparian2531_utm10_315, umf.rem)
summary(Rem.models$observer_dateXriparian315)

Rem.models$observer_dateXriparian615 <- multinomPois(~ Observer + JulDat ~ riparian2531_utm10_615, umf.rem)
summary(Rem.models$observer_dateXriparian615)

Rem.models$observer_dateXriparian1215 <- multinomPois(~ Observer + JulDat  ~ riparian2531_utm10_1215, umf.rem)
summary(Rem.models$observer_dateXriparian1215)

## Check if these are any good and what radius is best
MS.Rem.models <- modSel(fitList(fits=Rem.models));MS.Rem.models


# npow
Rem.models$observer_dateXnpow165 <- multinomPois(~ Observer + JulDat ~ npow16_utm10_165, umf.rem)
summary(Rem.models$observer_dateXnpow165)

Rem.models$observer_dateXnpow315 <- multinomPois(~ Observer + JulDat ~ npow16_utm10_315, umf.rem)
summary(Rem.models$observer_dateXnpow315)

Rem.models$observer_dateXnpow615 <- multinomPois(~ Observer + JulDat ~ npow16_utm10_615, umf.rem)
summary(Rem.models$observer_dateXnpow615)

Rem.models$observer_dateXnpow1215 <- multinomPois(~ Observer + JulDat  ~ npow16_utm10_1215, umf.rem)
summary(Rem.models$observer_dateXnpow1215)


# Now try wvups
Rem.models$observer_dateXwvups165 <- multinomPois(~ Observer + JulDat ~ wvups22_utm10_165, umf.rem)
summary(Rem.models$observer_dateXwvups165)

Rem.models$observer_dateXwvups315 <- multinomPois(~ Observer + JulDat ~ wvups22_utm10_315, umf.rem)
summary(Rem.models$observer_dateXwvups315)

Rem.models$observer_dateXwvups615 <- multinomPois(~ Observer + JulDat ~ wvups22_utm10_615, umf.rem)
summary(Rem.models$observer_dateXwvups615)

Rem.models$observer_dateXwvups1215 <- multinomPois(~ Observer + JulDat  ~ wvups22_utm10_1215, umf.rem)
summary(Rem.models$observer_dateXwvups1215)

## Check how these did
MS.Rem.models <- modSel(fitList(fits=Rem.models));MS.Rem.models


# wtf
Rem.models$observer_dateXwtf165 <- multinomPois(~ Observer + JulDat ~ wtf_utm10_165, umf.rem)
summary(Rem.models$observer_dateXwtf165)

Rem.models$observer_dateXwtf315 <- multinomPois(~ Observer + JulDat ~ wtf_utm10_315, umf.rem)
summary(Rem.models$observer_dateXwtf315)

Rem.models$observer_dateXwtf615 <- multinomPois(~ Observer + JulDat ~ wtf_utm10_615, umf.rem)
summary(Rem.models$observer_dateXwtf615)

Rem.models$observer_dateXwtf1215 <- multinomPois(~ Observer + JulDat  ~ wtf_utm10_1215, umf.rem)
summary(Rem.models$observer_dateXwtf1215)

## Check if these are any good and what radius is best
MS.Rem.models <- modSel(fitList(fits=Rem.models));MS.Rem.models


# devos
Rem.models$observer_dateXdevos165 <- multinomPois(~ Observer + JulDat ~ devos_utm10_165, umf.rem)
summary(Rem.models$observer_dateXdevos165)

Rem.models$observer_dateXdevos315 <- multinomPois(~ Observer + JulDat ~ devos_utm10_315, umf.rem)
summary(Rem.models$observer_dateXdevos315)

Rem.models$observer_dateXdevos615 <- multinomPois(~ Observer + JulDat ~ devos_utm10_615, umf.rem)
summary(Rem.models$observer_dateXdevos615)

Rem.models$observer_dateXdevos1215 <- multinomPois(~ Observer + JulDat  ~ devos_utm10_1215, umf.rem)
summary(Rem.models$observer_dateXdevos1215)

## Check if these are any good and what radius is best
MS.Rem.models <- modSel(fitList(fits=Rem.models));MS.Rem.models

## Now try habitat suitability as a predictor of abundance
Rem.models$observer_dateXpred.occur <- multinomPois(~ Observer + JulDat ~ Pred.Occur, umf.rem)
summary(Rem.models$observer_dateXpred.occur)

Rem.models$observer_dateXpred.occur2 <- multinomPois(~ Observer + JulDat ~ Pred.Occur + I(Pred.Occur^2), umf.rem)
summary(Rem.models$observer_dateXpred.occur2)


## Try each of detection structure again with the top model of density
Rem.models$nullXpred.occur2 <- multinomPois(~ 1 ~ Pred.Occur + I(Pred.Occur^2), umf.rem)
summary(Rem.models$nullXpred.occur2)

Rem.models$observerXpred.occur2 <- multinomPois(~ Observer ~ Pred.Occur + I(Pred.Occur^2), umf.rem)
summary(Rem.models$observerXpred.occur2)

Rem.models$dateXpred.occur2 <- multinomPois(~ JulDat ~ Pred.Occur + I(Pred.Occur^2), umf.rem)
summary(Rem.models$dateXpred.occur2)

Rem.models$timeXpred.occur2 <- multinomPois(~ Time ~ Pred.Occur + I(Pred.Occur^2), umf.rem)
summary(Rem.models$timeXpred.occur2)

Rem.models$observer_dateXpred.occur2 <- multinomPois(~ Observer + Date ~ Pred.Occur + I(Pred.Occur^2), umf.rem)
summary(Rem.models$observer_dateXpred.occur2)

Rem.models$observer_timeXpred.occur2 <- multinomPois(~ Observer + Time ~ Pred.Occur + I(Pred.Occur^2), umf.rem)
summary(Rem.models$observer_timeXpred.occur2)

Rem.models$observer_date_timeXpred.occur2 <- multinomPois(~ Observer + JulDat + Time ~ Pred.Occur + I(Pred.Occur^2), umf.rem)
summary(Rem.models$observer_date_timeXpred.occur2)


## How does habitat suitability measure up
MS.Rem.models <- modSel(fitList(fits=Rem.models));MS.Rem.models


## I decided against adding more to models because 
# 1. these habitat variables were already included in suitability
# 2. predicting their effects can be large even though for wvups, the maximum is only 14%. I might well get numbers that are not within the model
# in the historic data. The error at higher percentages is huge. 



###### Check to see if model predictions make sense ######
## I am particularly concerned that the range of covariate values in the historic data may not match the current data on which
## my models are built. Unrealistic predictions could be made based on this. Need to either use a lower radius model so that the 
## range of percent land cover values is more similar between historic and current models, or check to see if predictions are
## even very different and then use the top model if they aren't.
summary(umf.rem)

## Make sequence of the variable affecting density
Pred.Occur.seq <- seq(0,1,0.1)

## Now predict with top model
newdat <- data.frame(Pred.Occur = Pred.Occur.seq)
#newdat$wvups22_utm10_165 <- 0.295
head(newdat)

abundance.pred.pred.occur2.rem <- predict(Rem.models$observerXpred.occur2, type = "state", newdata = newdat, appendData=TRUE); abundance.pred.pred.occur2.rem

## Make a figure
jpeg(file.path(output.directory, paste( model.name, '_rem_predoccur2.jpg', sep="")))

plot(Pred.Occur.seq, abundance.pred.pred.occur2.rem[, "Predicted"], xlab="Predicted Occurrence Squared (SDM)", ylab="Predicted abundance (300m radius)", type='l', 
     ylim = c(0,3.1), frame=F, lwd=2)
matlines(Pred.Occur.seq, abundance.pred.pred.occur2.rem[,3:4], lty=1, col="grey", lwd=1)

dev.off()


###### Model Table ######
## Make a fistlist for selection table
MS.Rem.models <- modSel(fitList(fits=Rem.models));MS.Rem.models

## Export Model Table
toExport <- as(MS.Rem.models, "data.frame")
write.csv(toExport, file.path(output.directory, paste(model.name, "_RemTable.csv", sep="")))

###### Best Model ######
## Check the goodness of fit of the best model and export summary
# Goodness of fit
pb.rem <- parboot(Rem.models$observerXpred.occur2, fitstats, nsim=100, report=5); pb.rem # Can't do this apparently because of the issues in the model
c.hat.rem <- pb.rem@t0[2]/mean(pb.rem@t.star[,2]); c.hat.rem

filename <- "rem_parboot.txt"
capture.output(pb.rem, file = file.path(output.directory,filename))

# Export best model summary
filename <- "rem_model_summary.txt"
capture.output(summary(Rem.models$observerXpred.occur2), file = file.path(output.directory,filename))


# Detection probability after all passes for each site. Could take an average across all sites
p.rem <- mean(rowSums(getP(Rem.models$observerXpred.occur2))); p.rem




#########################################################
# 13. Negative Binomial Removal Models
#########################################################


## Try with negative binomial to see if this even works
umf.grem <- unmarkedFrameGMM(y=removal.matrix, siteCovs=covs, numPrimary=1, type="removal")

## Make a list to store the models
gRem.models <- list()

## Run Removal Models and add them to list. 
## Start out with variables affecting detection.
gRem.models$null <- gmultmix(~ 1, ~ 1, ~ 1, umf.grem, mixture="NB", K=10)
summary(gRem.models$null)

gRem.models$observer <- gmultmix(~ 1, ~ 1, ~ Observer, umf.grem, mixture="NB", K=10)
summary(gRem.models$observer)

gRem.models$date <- gmultmix(~ 1, ~ 1, ~ JulDat, umf.grem, mixture="NB", K=10)
summary(gRem.models$date)

gRem.models$time <- gmultmix(~ 1, ~ 1, ~ Time, umf.grem, mixture="NB", K=10)
summary(gRem.models$time)

gRem.models$observer_date <- gmultmix(~ 1, ~ 1, ~ Observer + JulDat, umf.grem, mixture="NB", K=10)
summary(gRem.models$observer_date)

gRem.models$observer_time <- gmultmix(~ 1, ~ 1, ~ Observer + Time, umf.grem, mixture="NB", K=10)
summary(gRem.models$observer_time)

gRem.models$date_time <- gmultmix(~ 1, ~ 1, ~ JulDat + Time, umf.grem, mixture="NB", K=10)
summary(gRem.models$date_time)

gRem.models$observer_date_time <- gmultmix(~ 1, ~ 1, ~ Observer + JulDat + Time, umf.grem, mixture="NB", K=10)
summary(gRem.models$observer_date_time)


## Check which of these make any difference
MS.gRem.models <- modSel(fitList(fits=gRem.models));MS.gRem.models


## Unused variables because they are confounded with their effects on density
#gRem.models$hwy <- gmultmix(~ 1, ~ 1, ~ HwyDist, umf.grem, mixture="NB", K=10)
#summary(gRem.models$hwy)

#gRem.models$hwy2 <- gmultmix(~ 1, ~ 1, ~ HwyDist + I(HwyDist^2), umf.grem, mixture="NB", K=10)
#summary(gRem.models$hwy2)

#gRem.models$hwylog <- gmultmix(~ 1, ~ 1, ~ HwyDistLog, umf.grem, mixture="NB", K=10)
#summary(gRem.models$hwylog)

#gRem.models$riv <- gmultmix(~ 1, ~ 1, ~ RivDist, umf.grem, mixture="NB", K=10)
#summary(gRem.models$riv)

#gRem.models$riv2 <- gmultmix(~ 1, ~ 1, ~ RivDist + I(RivDist^2), umf.grem, mixture="NB", K=10)
#summary(gRem.models$riv)

#gRem.models$rivlog <- gmultmix(~ 1, ~ 1, ~ RivDistLog, umf.grem, mixture="NB", K=10)
#summary(gRem.models$rivlog)

#gRem.models$cancov <- gmultmix(~ 1, ~ 1, ~ will_CANCOV__UTM10_165, umf.grem, mixture="NB", K=10)
#summary(gRem.models$cancov)

#gRem.models$treedensity <- gmultmix(~ 1, ~ 1, ~ will_TPH_GE_3__UTM10_165, umf.grem, mixture="NB", K=10, starts = c(-1.28, -0.7860, 0.0, 0.223))
#summary(gRem.models$treedensity)

#gRem.models$urban165 <- gmultmix(~ 1, ~ 1, ~ urban_tot_utm10_165, umf.grem, mixture="NB", K=10)
#summary(gRem.models$urban165)

#gRem.models$urban315 <- gmultmix(~ 1, ~ 1, ~ urban_tot_utm10_315, umf.grem, mixture="NB", K=10)
#summary(gRem.models$urban315)


## Now model abundance within the truncation distance
names(covs)

gRem.models$observer_dateXriparian165 <- gmultmix(~ riparian2531_utm10_165, ~ 1, ~ Observer + JulDat , umf.grem, mixture="NB", K=10)
summary(gRem.models$observer_dateXriparian165)

gRem.models$observer_dateXriparian315 <- gmultmix(~ riparian2531_utm10_315, ~ 1, ~ Observer + JulDat , umf.grem, mixture="NB", K=10)
summary(gRem.models$observer_dateXriparian315)

gRem.models$observer_dateXriparian615 <- gmultmix(~ riparian2531_utm10_615, ~ 1, ~ Observer + JulDat , umf.grem, mixture="NB", K=10)
summary(gRem.models$observer_dateXriparian615)

gRem.models$dobserver_dateXriparian1215 <- gmultmix(~ riparian2531_utm10_1215, ~ 1, ~ Observer + JulDat , umf.grem, mixture="NB", K=10)
summary(gRem.models$dobserver_dateXriparian1215)

## Check if these are any good and what radius is best
MS.gRem.models <- modSel(fitList(fits=gRem.models));MS.gRem.models

## oak
gRem.models$observer_dateXoak165 <- gmultmix(~ oak_utm10_165, ~ 1, ~ Observer + JulDat , umf.grem, mixture="NB", K=10)
summary(gRem.models$observer_dateXoak165)

gRem.models$observer_dateXoak315 <- gmultmix(~ oak_utm10_315, ~ 1, ~ Observer + JulDat , umf.grem, mixture="NB", K=10)
summary(gRem.models$observer_dateXoak315)

gRem.models$observer_dateXoak615 <- gmultmix(~ oak_utm10_615, ~ 1, ~ Observer + JulDat , umf.grem, mixture="NB", K=10)
summary(gRem.models$observer_dateXoak615)

gRem.models$dobserver_dateXoak1215 <- gmultmix(~ oak_utm10_1215, ~ 1, ~ Observer + JulDat , umf.grem, mixture="NB", K=10)
summary(gRem.models$dobserver_dateXoak1215)

## NPOW
gRem.models$observer_dateXnpow165 <- gmultmix(~ npow16_utm10_165, ~ 1, ~ Observer + JulDat , umf.grem, mixture="NB", K=10)
summary(gRem.models$observer_dateXnpow165)

gRem.models$observer_dateXnpow315 <- gmultmix(~ npow16_utm10_315, ~ 1, ~ Observer + JulDat , umf.grem, mixture="NB", K=10)
summary(gRem.models$observer_dateXnpow315)

gRem.models$observer_dateXnpow615 <- gmultmix(~ npow16_utm10_615, ~ 1, ~ Observer + JulDat , umf.grem, mixture="NB", K=10)
summary(gRem.models$observer_dateXnpow615)

gRem.models$dobserver_dateXnpow1215 <- gmultmix(~ npow16_utm10_1215, ~ 1, ~ Observer + JulDat , umf.grem, mixture="NB", K=10)
summary(gRem.models$dobserver_dateXnpow1215)

## WVUPS
gRem.models$observer_dateXwvups165 <- gmultmix(~ wvups22_utm10_165, ~ 1, ~ Observer + JulDat , umf.grem, mixture="NB", K=10)
summary(gRem.models$observer_dateXwvups165)

gRem.models$observer_dateXwvups315 <- gmultmix(~ wvups22_utm10_315, ~ 1, ~ Observer + JulDat , umf.grem, mixture="NB", K=10)
summary(gRem.models$observer_dateXwvups315)

gRem.models$observer_dateXwvups615 <- gmultmix(~ wvups22_utm10_615, ~ 1, ~ Observer + JulDat , umf.grem, mixture="NB", K=10)
summary(gRem.models$observer_dateXwvups615)

gRem.models$dobserver_dateXwvups1215 <- gmultmix(~ wvups22_utm10_1215, ~ 1, ~ Observer + JulDat , umf.grem, mixture="NB", K=10)
summary(gRem.models$dobserver_dateXwvups1215)

## Check if these are any good and what radius is best
MS.gRem.models <- modSel(fitList(fits=gRem.models));MS.gRem.models


## Now try habitat suitability as a predictor of abundance
gRem.models$observer_dateXpred.occur <- gmultmix(~ Pred.Occur, ~ 1, ~ Observer + JulDat, umf.grem, mixture="NB", K=10)
summary(gRem.models$observer_dateXpred.occur)

gRem.models$observer_dateXpred.occur2 <- gmultmix(~ Pred.Occur + I(Pred.Occur^2), ~ 1, ~ Observer + JulDat, umf.grem, mixture="NB", K=10)
summary(gRem.models$observer_dateXpred.occur2)


## Now check if the other detection structures perform better with this top abundance structure
gRem.models$nullXpred.occur2 <- gmultmix(~ Pred.Occur + I(Pred.Occur^2), ~ 1, ~ 1, umf.grem, mixture="NB")
summary(gRem.models$nullXpred.occur2)

gRem.models$observerXpred.occur2 <- gmultmix(~ Pred.Occur + I(Pred.Occur^2), ~ 1, ~ Observer, umf.grem, mixture="NB")
summary(gRem.models$observerXpred.occur2)

gRem.models$dateXpred.occur2 <- gmultmix(~ Pred.Occur + I(Pred.Occur^2), ~ 1, ~ JulDat + Observer, umf.grem, mixture="NB", K=10)
summary(gRem.models$date_observerXpred.occur2)

gRem.models$timeXpred.occur2 <- gmultmix(~ Pred.Occur + I(Pred.Occur^2), ~ 1, ~ Time, umf.grem, mixture="NB", K=10)
summary(gRem.models$timeXpred.occur2)

gRem.models$observer_dateXpred.occur2 <- gmultmix(~ Pred.Occur + I(Pred.Occur^2), ~ 1, ~ Observer + JulDat, umf.grem, mixture="NB")
summary(gRem.models$observer_dateXpred.occur2)

gRem.models$observer_timeXpred.occur2 <- gmultmix(~ Pred.Occur + I(Pred.Occur^2), ~ 1, ~ Observer + Time, umf.grem, mixture="NB")
summary(gRem.models$observer_timeXpred.occur2)

gRem.models$date_timeXpred.occur2 <- gmultmix(~ Pred.Occur + I(Pred.Occur^2), ~ 1, ~ JulDat + Time, umf.grem, mixture="NB", K=10)
summary(gRem.models$date_timeXpred.occur2)

gRem.models$observer_date_timeXpred.occur2 <- gmultmix(~ Pred.Occur + I(Pred.Occur^2), ~ 1, ~ Observer + JulDat + Time, umf.grem, mixture="NB")
summary(gRem.models$observer_date_timeXpred.occur2)

## How does habitat suitability measure up
MS.gRem.models <- modSel(fitList(fits=gRem.models));MS.gRem.models


###### Check to see if model predictions make sense ######
## I am particularly concerned that the range of covariate values in the historic data may not match the current data on which
## my models are built. Unrealistic predictions could be made based on this. Need to either use a lower radius model so that the 
## range of percent land cover values is more similar between historic and current models, or check to see if predictions are
## even very different and then use the top model if they aren't.
summary(umf)

## Make sequence of the variable affecting density
Pred.Occur.seq <- seq(0,1,0.1)

## Now predict with top model
newdat <- data.frame(Pred.Occur = Pred.Occur.seq)
head(newdat)

abundance.pred.pred.occur2.grem <- predict(gRem.models$observer_dateXpred.occur2, type = "lambda", newdata = newdat, appendData=TRUE); abundance.pred.pred.occur2.grem

## Make a figure
jpeg(file.path(output.directory, paste( model.name, '_grem_predoccur2.jpg', sep="")))

plot(Pred.Occur.seq, abundance.pred.pred.occur2.grem[, "Predicted"], xlab="Predicted Occurrence Squared (SDM)", ylab="Predicted abundance (300m radius)", type='l', 
     ylim = c(0,3.1), frame=F, lwd=2)
matlines(Pred.Occur.seq, abundance.pred.pred.occur2.grem[,3:4], lty=1, col="grey", lwd=1)

dev.off()



###### Model Table ######
## Make a fistlist for selection table
MS.gRem.models <- modSel(fitList(fits=gRem.models));MS.gRem.models

## Export Model Table
toExport <- as(MS.gRem.models, "data.frame")
write.csv(toExport, file.path(output.directory, paste(model.name, "_gRemTable.csv", sep="")))

###### Best Model ######
## Check the goodness of fit of the best model and export summary
# Goodness of fit
pb.grem <- parboot(gRem.models$observer_dateXpred.occur2, fitstats, nsim=10, report=5); pb.grem 
c.hat.grem <- pb.grem@t0[2]/mean(pb.grem@t.star[,2]); c.hat.grem

filename <- "grem_parboot.txt"
capture.output(pb.grem, file = file.path(output.directory,filename))

# Export best model summary
filename <- "grem_model_summary.txt"
capture.output(summary(gRem.models$observer_dateXpred.occur2), file = file.path(output.directory,filename))


# Detection probability after all passes for each site. Could take an average across all sites
p.grem <- mean(rowSums(getP(gRem.models$observer_dateXpred.occur2))); p.grem







#########################################################
# 14. Time of Detection Models
#########################################################

###### Observation Data with Detection History ###### 
## Convert intervals to detection history
names(Observation.data.trunc.full)
Observation.data.trunc.full$CaptureHistory <- paste(Observation.data.trunc.full$Interval_1, Observation.data.trunc.full$Interval_2, 
                                                    Observation.data.trunc.full$Interval_3, Observation.data.trunc.full$Interval_4, 
                                                    Observation.data.trunc.full$Interval_5, sep="")
names(Observation.data.trunc.full)
head(Observation.data.trunc.full$CaptureHistory)

## Convert detection hstory to a factor
str(Observation.data.trunc.full)
Observation.data.trunc.full$CaptureHistory <- as.factor(Observation.data.trunc.full$CaptureHistory)
str(Observation.data.trunc.full)
levels(Observation.data.trunc.full$CaptureHistory)

## Need to use this for 5 intervals so that R can caluclate frequencies that each 
## detection history occured even when it did not occur at all
Observation.data.trunc.full$CaptureHistory <- factor(Observation.data.trunc.full$CaptureHistory, 
                                                     levels=c("00001", "00010", "00011", "00100", 
                                                              "00101", "00110", "00111", "01000", 
                                                              "01001", "01010", "01011","01100", 
                                                              "01101", "01110", "01111", "10000", 
                                                              "10001", "10010", "10011", "10100", 
                                                              "10101", "10110", "10111", "11000", 
                                                              "11001", "11010", "11011", "11100", 
                                                              "11101", "11110", "11111"))
str(Observation.data.trunc.full)

## To make sure that unmarked knows there are more sites than are represented by the 
## observation data you must do the following. This makes it so that sites with zero
## detections are still included
#Observation.data.trunc.full$NewSiteID <- factor(Observation.data.trunc.full$NewSiteID, 
#                                       levels=levels(observation.covariates.data$NewSiteID))


## Now tabulate so that it counts the number of each detection history at each site
Obs.hist <- table(Observation.data.trunc.full$NewSiteID, Observation.data.trunc.full$CaptureHistory)
head(Obs.hist)
class(Obs.hist) <- "matrix"

###### Analyses ###### 
## Now I have to create the multinomial cell probabilities on my own since this 
## function does not have this as a data type
crPiFun <- function(p) {
  p1 <- p[,1]
  p2 <- p[,2]
  p3 <- p[,3]
  p4 <- p[,4]
  p5 <- p[,5]
  cbind("00001" = (1-p1) * (1-p2) * (1-p3) * (1-p4) * p5,
        "00010" = (1-p1) * (1-p2) * (1-p3) * p4 * (1-p5),
        "00011" = (1-p1) * (1-p2) * (1-p3) * p4 * p5,
        "00100" = (1-p1) * (1-p2) * p3 * (1-p4) * (1-p5),
        "00101" = (1-p1) * (1-p2) * p3 * (1-p4) * p5,
        "00110" = (1-p1) * (1-p2) * p3 * p4 * (1-p5),
        "00111" = (1-p1) * (1-p2) * p3 * p4 * p5,
        "01000" = (1-p1) * p2 * (1-p3) * (1-p4) * (1-p5),
        "01001" = (1-p1) * p2 * (1-p3) * (1-p4) * p5,
        "01010" = (1-p1) * p2 * (1-p3) * p4 * (1-p5),
        "01011" = (1-p1) * p2 * (1-p3) * p4 * p5,
        "01100" = (1-p1) * p2 * p3 * (1-p4) * (1-p5),
        "01101" = (1-p1) * p2 * p3 * (1-p4) * p5,
        "01110" = (1-p1) * p2 * p3 * p4 * (1-p5),
        "01111" = (1-p1) * p2 * p3 * p4 * p5,
        "10000" = p1 * (1-p2) * (1-p3) * (1-p4) * (1-p5),
        "10001" = p1 * (1-p2) * (1-p3) * (1-p4) * p5,
        "10010" = p1 * (1-p2) * (1-p3) * p4 * (1-p5),
        "10011" = p1 * (1-p2) * (1-p3) * p4 * p5,
        "10100" = p1 * (1-p2) * p3 * (1-p4) * (1-p5),
        "10101" = p1 * (1-p2) * p3 * (1-p4) * p5,
        "10110" = p1 * (1-p2) * p3 * p4 * (1-p5),
        "10111" = p1 * (1-p2) * p3 * p4 * p5,
        "11000" = p1 * p2 * (1-p3) * (1-p4) * (1-p5),
        "11001" = p1 * p2 * (1-p3) * (1-p4) * p5,
        "11010" = p1 * p2 * (1-p3) * p4 * (1-p5),
        "11011" = p1 * p2 * (1-p3) * p4 * p5,
        "11100" = p1 * p2 * p3 * (1-p4) * (1-p5),
        "11101" = p1 * p2 * p3 * (1-p4) * p5,
        "11110" = p1 * p2 * p3 * p4 * (1-p5),
        "11111" = p1 * p2 * p3 * p4 * p5)
}

## Need to tell unmarked that if there are missing observation covariates
## the y data should be excluded.
o2y <- matrix(1, 5, 31)

## Creating an observation covariate here for each interval here? I'm not sure that this
## actually makes much sense here because they were all consecutive. On the other hand, 
## perhaps birds are stimulated by the comotion of the researcher arriving and sing more
## in the first minute.
intervalMat <- matrix(c("1","2","3","4","5"),
                      length(covs$Julian_Dat), 5, 
                      byrow=TRUE) # This allows for the calculation of Mt, M0 etc. You do this by specifying
# interval while running the model. 
str(intervalMat)



## Create an unmarked frame with all of the data above
umf.tod <- unmarkedFrameMPois(y=Obs.hist, siteCovs=covs, obsCovs=list(interval=intervalMat),
                              obsToY=o2y, piFun="crPiFun")

## Make a list to store the models
TOD.models <- list()

## Run Time of Detection Models and add them to list. 
## Start out with variables affecting detection.
TOD.models$null <- multinomPois(~ 1 ~ 1, umf.tod)
summary(TOD.models$null)

TOD.models$observer <- multinomPois(~ Observer ~ 1, umf.tod)
summary(TOD.models$observer)

TOD.models$date <- multinomPois(~ JulDat ~ 1, umf.tod)
summary(TOD.models$date)

TOD.models$time <- multinomPois(~ Time ~ 1, umf.tod)
summary(TOD.models$time)

TOD.models$observer_date <- multinomPois(~ Observer + JulDat ~ 1, umf.tod)
summary(TOD.models$observer_date)

TOD.models$observer_time <- multinomPois(~ Observer + Time ~ 1, umf.tod)
summary(TOD.models$observer_time)

TOD.models$date_time <- multinomPois(~ JulDat + Time ~ 1, umf.tod)
summary(TOD.models$date_time)

TOD.models$observer_date_time <- multinomPois(~ Observer + JulDat + Time ~ 1, umf.tod)
summary(TOD.models$observer_date_time)


## Check which of these make any difference
MS.TOD.models <- modSel(fitList(fits=TOD.models));MS.TOD.models


## Unused variables because they are confounded with their effects on density
#TOD.models$hwy <- multinomPois(~ HwyDist ~ 1, umf.tod)
#summary(TOD.models$hwy)

#TOD.models$hwy2 <- multinomPois(~ HwyDist + I(HwyDist^2) ~ 1, umf.tod)
#summary(TOD.models$hwy2)

#TOD.models$hwylog <- multinomPois(~ HwyDistLog ~ 1, umf.tod)
#summary(TOD.models$hwylog)

#TOD.models$riv <- multinomPois(~ RivDist ~ 1, umf.tod)
#summary(TOD.models$riv)

#TOD.models$riv2 <- multinomPois(~ RivDist + I(RivDist^2) ~ 1, umf.tod)
#summary(TOD.models$riv)

#TOD.models$rivlog <- multinomPois(~ RivDistLog ~ 1, umf.tod)
#summary(TOD.models$rivlog)

#TOD.models$urban165 <- multinomPois(~ urban_tot_utm10_165 ~ 1, umf.tod)
#summary(TOD.models$urban165)

#TOD.models$urban315 <- multinomPois(~ urban_tot_utm10_315 ~ 1, umf.tod)
#summary(TOD.models$urban315)


## Now model abundance within the truncation distance
names(covs)

## riparian
TOD.models$observerXriparian165 <- multinomPois(~ Observer ~ riparian2531_utm10_165, umf.tod)
summary(TOD.models$observerXriparian165)

TOD.models$observerXriparian315 <- multinomPois(~ Observer ~ riparian2531_utm10_315, umf.tod)
summary(TOD.models$observerXriparian315)

TOD.models$observerXriparian615 <- multinomPois(~ Observer ~ riparian2531_utm10_615, umf.tod)
summary(TOD.models$observerXriparian615)

TOD.models$observerXriparian1215 <- multinomPois(~ Observer ~ riparian2531_utm10_1215, umf.tod)
summary(TOD.models$observerXriparian1215)

## Check how these did
MS.TOD.models <- modSel(fitList(fits=TOD.models));MS.TOD.models


## npow
TOD.models$observerXnpow165 <- multinomPois(~ Observer ~ npow16_utm10_165, umf.tod)
summary(TOD.models$observerXnpow165)

TOD.models$observerXnpow315 <- multinomPois(~ Observer ~ npow16_utm10_315, umf.tod)
summary(TOD.models$observerXnpow315)

TOD.models$observerXnpow615 <- multinomPois(~ Observer ~ npow16_utm10_615, umf.tod)
summary(TOD.models$observerXnpow615)

TOD.models$observerXnpow1215 <- multinomPois(~ Observer ~ npow16_utm10_1215, umf.tod)
summary(TOD.models$observerXnpow1215)

## Check how these did
MS.TOD.models <- modSel(fitList(fits=TOD.models));MS.TOD.models


## wvups
TOD.models$observerXwvups165 <- multinomPois(~ Observer ~ wvups22_utm10_165, umf.tod)
summary(TOD.models$observerXwvups165)

TOD.models$observerXwvups315 <- multinomPois(~ Observer ~ wvups22_utm10_315, umf.tod)
summary(TOD.models$observerXwvups315)

TOD.models$observerXwvups615 <- multinomPois(~ Observer ~ wvups22_utm10_615, umf.tod)
summary(TOD.models$observerXwvups615)

TOD.models$observerXwvups1215 <- multinomPois(~ Observer ~ wvups22_utm10_1215, umf.tod)
summary(TOD.models$observerXwvups1215)

## Check how these did
MS.TOD.models <- modSel(fitList(fits=TOD.models));MS.TOD.models


## wtf
TOD.models$observerXwtf165 <- multinomPois(~ Observer ~ wtf_utm10_165, umf.tod)
summary(TOD.models$observerXwtf165)

TOD.models$observerXwtf315 <- multinomPois(~ Observer ~ wtf_utm10_315, umf.tod)
summary(TOD.models$observerXwtf315)

TOD.models$observerXwtf615 <- multinomPois(~ Observer ~ wtf_utm10_615, umf.tod)
summary(TOD.models$observerXwtf615)

TOD.models$observerXwtf1215 <- multinomPois(~ Observer ~ wtf_utm10_1215, umf.tod)
summary(TOD.models$observerXwtf1215)

## Check how these did
MS.TOD.models <- modSel(fitList(fits=TOD.models));MS.TOD.models


## devos
TOD.models$observerXdevos165 <- multinomPois(~ Observer ~ devos_utm10_165, umf.tod)
summary(TOD.models$observerXdevos165)

TOD.models$observerXdevos315 <- multinomPois(~ Observer ~ devos_utm10_315, umf.tod)
summary(TOD.models$observerXdevos315)

TOD.models$observerXdevos615 <- multinomPois(~ Observer ~ devos_utm10_615, umf.tod)
summary(TOD.models$observerXdevos615)

TOD.models$observerXdevos1215 <- multinomPois(~ Observer ~ devos_utm10_1215, umf.tod)
summary(TOD.models$observerXdevos1215)

## Check how these did
MS.TOD.models <- modSel(fitList(fits=TOD.models));MS.TOD.models



## Now try habitat suitability as a predictor of abundance
TOD.models$observerXpred.occur <- multinomPois(~ Observer ~ Pred.Occur, umf.tod)
summary(TOD.models$observerXpred.occur)

TOD.models$observerXpred.occur2 <- multinomPois(~ Observer ~ Pred.Occur + I(Pred.Occur^2), umf.tod)
summary(TOD.models$observerXpred.occur2)


## Now check if the other detection structures perform better with this top abundance structure
TOD.models$nullXpred.occur2 <- multinomPois(~ 1 ~ Pred.Occur + I(Pred.Occur^2), umf.tod)
summary(TOD.models$nullXpred.occur2)

TOD.models$observerXpred.occur2 <- multinomPois(~ Observer ~ Pred.Occur + I(Pred.Occur^2), umf.tod)
summary(TOD.models$observerXpred.occur2)

TOD.models$dateXpred.occur2 <- multinomPois(~ JulDat ~ Pred.Occur + I(Pred.Occur^2), umf.tod)
summary(TOD.models$dateXpred.occur2)

TOD.models$timeXpred.occur2 <- multinomPois(~ Time ~ Pred.Occur + I(Pred.Occur^2), umf.tod)
summary(TOD.models$timeXpred.occur2)

TOD.models$observer_dateXpred.occur2 <- multinomPois(~ Observer + JulDat ~ Pred.Occur + I(Pred.Occur^2), umf.tod)
summary(TOD.models$observer_dateXpred.occur2)

TOD.models$observer_timeXpred.occur2 <- multinomPois(~ Observer + Time ~ Pred.Occur + I(Pred.Occur^2), umf.tod)
summary(TOD.models$observer_timeXpred.occur2)

TOD.models$date_timeXpred.occur2 <- multinomPois(~ JulDat + Time ~ Pred.Occur + I(Pred.Occur^2), umf.tod)
summary(TOD.models$date_timeXpred.occur2)

TOD.models$observer_date_timeXpred.occur2 <- multinomPois(~ Observer + JulDat + Time ~ Pred.Occur + I(Pred.Occur^2), umf.tod)
summary(TOD.models$observer_date_timeXpred.occur2)

## How does habitat suitability measure up
MS.TOD.models <- modSel(fitList(fits=TOD.models));MS.TOD.models



###### Check to see if model predictions make sense ######
## I am particularly concerned that the range of covariate values in the historic data may not match the current data on which
## my models are built. Unrealistic predictions could be made based on this. Need to either use a lower radius model so that the 
## range of percent land cover values is more similar between historic and current models, or check to see if predictions are
## even very different and then use the top model if they aren't.
summary(umf.tod)

## Make sequence of the variable affecting density
Pred.Occur.seq <- seq(0,1,0.1)

## Now predict with top model
newdat <- data.frame(Pred.Occur = Pred.Occur.seq)
head(newdat)

abundance.pred.pred.occur2.tod <- predict(TOD.models$observerXpred.occur2 , type = "state", newdata = newdat, appendData=TRUE); abundance.pred.pred.occur2.tod

## Make a figure
jpeg(file.path(output.directory, paste( model.name, '_tod_predoccur2.jpg', sep="")))

plot(Pred.Occur.seq, abundance.pred.pred.occur2.tod[, "Predicted"], xlab="Predicted Occurrence Squared (SDM)", ylab="Predicted abundance (300m radius)", type='l', 
     ylim = c(0,3.1), frame=F, lwd=2)
matlines(Pred.Occur.seq, abundance.pred.pred.occur2.tod[,3:4], lty=1, col="grey", lwd=1)

dev.off()



###### Model Table ######
## Make a fistlist for selection table
MS.TOD.models <- modSel(fitList(fits=TOD.models));MS.TOD.models

## Export Model Table
toExport <- as(MS.TOD.models, "data.frame")
write.csv(toExport, file.path(output.directory, paste(model.name, "_TODTable.csv", sep="")))

###### Best Model ######
## Check the goodness of fit of the best model and export summary
# Goodness of fit
pb.tod <- parboot(TOD.models$observerXpred.occur2, fitstats, nsim=100, report=5); pb.tod # Can't do this apparently because of the issues in the model
c.hat.tod <- pb.tod@t0[2]/mean(pb.tod@t.star[,2]); c.hat.tod

filename <- "tod_parboot.txt"
capture.output(pb.tod, file = file.path(output.directory,filename))

# Export best model summary
filename <- "tod_model_summary.txt"
capture.output(summary(TOD.models$observerXpred.occur2), file = file.path(output.directory,filename))


# Detection probability for a single pass
#backTransform(Rem.models$date_hwy2Xpred.occur2_wvups165, type="det")
head(getP(TOD.models$observerXpred.occur2))

# Detection probability after all passes for each site. Could take an average across all sites
p.tod <- mean(rowSums(getP(TOD.models$observerXpred.occur2))); p.tod







#########################################################
# 15. Negative Binomial Time of Detection Models
#########################################################


## Try with negative binomial to see if this even works
umf.gtod <- unmarkedFrameGMM(y=Obs.hist, siteCovs=covs, numPrimary=1, obsCovs=list(interval=intervalMat),
                            obsToY=o2y, piFun="crPiFun")

## Make a list to store the models
gTOD.models <- list()

## Run Time of Detection Models and add them to list. 
## Start out with variables affecting detection.
gTOD.models$null <- gmultmix(~ 1, ~ 1, ~ 1, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$null)

gTOD.models$observer <- gmultmix(~ 1, ~ 1, ~ Observer, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer)

gTOD.models$date <- gmultmix(~ 1, ~ 1, ~ JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$date)

gTOD.models$time <- gmultmix(~ 1, ~ 1, ~ Time, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$time)

gTOD.models$observer_date <- gmultmix(~ 1, ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_date)

gTOD.models$observer_time <- gmultmix(~ 1, ~ 1, ~ Observer + Time, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_time)

gTOD.models$date_time <- gmultmix(~ 1, ~ 1, ~ JulDat + Time, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$date_time)

gTOD.models$observer_date_time <- gmultmix(~ 1, ~ 1, ~ Observer + JulDat + Time, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_date_time)

## Check which of these make any difference
MS.gTOD.models <- modSel(fitList(fits=gTOD.models));MS.gTOD.models


## Unused variables because they are confounded with their effects on density
#gTOD.models$hwy <- gmultmix(~ 1, ~ 1, ~ HwyDist, umf.gtod, mixture="NB", K=10)
#summary(gTOD.models$hwy)

#gTOD.models$hwy2 <- gmultmix(~ 1, ~ 1, ~ HwyDist + I(HwyDist^2), umf.gtod, mixture="NB", K=10)
#summary(gTOD.models$hwy2)

#gTOD.models$hwylog <- gmultmix(~ 1, ~ 1, ~ HwyDistLog, umf.gtod, mixture="NB", K=10)
#summary(gTOD.models$hwylog)

#gTOD.models$riv <- gmultmix(~ 1, ~ 1, ~ RivDist, umf.gtod, mixture="NB", K=10)
#summary(gTOD.models$riv)

#gTOD.models$riv2 <- gmultmix(~ 1, ~ 1, ~ RivDist + I(RivDist^2), umf.gtod, mixture="NB", K=10)
#summary(gTOD.models$riv2)

#gTOD.models$rivlog <- gmultmix(~ 1, ~ 1, ~ RivDistLog, umf.gtod, mixture="NB", K=10)
#summary(gTOD.models$rivlog)

#gTOD.models$urban165 <- gmultmix(~ 1, ~ 1, ~ urban_tot_utm10_165, umf.gtod, mixture="NB", K=10)
#summary(gTOD.models$urban165)

#gTOD.models$urban315 <- gmultmix(~ 1, ~ 1, ~ urban_tot_utm10_315, umf.gtod, mixture="NB", K=10)
#summary(gTOD.models$urban315)



## Now model abundance within the truncation distance
names(covs)

## riparian
gTOD.models$observer_dateXriparian165 <- gmultmix(~ riparian2531_utm10_165, ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXriparian165)

gTOD.models$observer_dateXriparian315 <- gmultmix(~ riparian2531_utm10_315, ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXriparian315)

gTOD.models$observer_dateXriparian615 <- gmultmix(~ riparian2531_utm10_615, ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXriparian615)

gTOD.models$observer_dateXriparian1215 <- gmultmix(~ riparian2531_utm10_1215, ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXriparian1215)

## oak
gTOD.models$observer_dateXoak165 <- gmultmix(~ oak_utm10_165, ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXoak165)

gTOD.models$observer_dateXoak315 <- gmultmix(~ oak_utm10_315, ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXoak315)

gTOD.models$observer_dateXoak615 <- gmultmix(~ oak_utm10_615, ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXoak615)

gTOD.models$observer_dateXoak1215 <- gmultmix(~ oak_utm10_1215, ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXoak1215)

## Check if these are any good and what radius is best
MS.gTOD.models <- modSel(fitList(fits=gTOD.models));MS.gTOD.models

## npow
gTOD.models$observer_dateXnpow165 <- gmultmix(~ npow16_utm10_165, ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXnpow165)

gTOD.models$observer_dateXnpow315 <- gmultmix(~ npow16_utm10_315, ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXnpow315)

gTOD.models$observer_dateXnpow615 <- gmultmix(~ npow16_utm10_615, ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXnpow615)

gTOD.models$observer_dateXnpow1215 <- gmultmix(~ npow16_utm10_1215, ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXnpow1215)

## wvups
gTOD.models$observer_dateXwvups165 <- gmultmix(~ wvups22_utm10_165, ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXwvups165)

gTOD.models$observer_dateXwvups315 <- gmultmix(~ wvups22_utm10_315, ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXwvups315)

gTOD.models$observer_dateXwvups615 <- gmultmix(~ wvups22_utm10_615, ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXwvups615)

gTOD.models$observer_dateXwvups1215 <- gmultmix(~ wvups22_utm10_1215, ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXwvups1215)

## Check if these are any good and what radius is best
MS.gTOD.models <- modSel(fitList(fits=gTOD.models));MS.gTOD.models



## Now try habitat suitability as a predictor of abundance
gTOD.models$observer_dateXpred.occur <- gmultmix(~ Pred.Occur, ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXpred.occur)

gTOD.models$observer_dateXpred.occur2 <- gmultmix(~ Pred.Occur + I(Pred.Occur^2), ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXpred.occur2)


## Now check if the other detection structures perform better with this top abundance structure
gTOD.models$nullXpred.occur2 <- gmultmix(~ Pred.Occur + I(Pred.Occur^2), ~ 1, ~ 1, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$nullXpred.occur2)

gTOD.models$observerXpred.occur2 <- gmultmix(~ Pred.Occur + I(Pred.Occur^2), ~ 1, ~ Observer, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observerXpred.occur2)

gTOD.models$dateXpred.occur2 <- gmultmix(~ Pred.Occur + I(Pred.Occur^2), ~ 1, ~ JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$dateXpred.occur2)

gTOD.models$timeXpred.occur2 <- gmultmix(~ Pred.Occur + I(Pred.Occur^2), ~ 1, ~ Time, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$timeXpred.occur2)

gTOD.models$observer_dateXpred.occur2 <- gmultmix(~ Pred.Occur + I(Pred.Occur^2), ~ 1, ~ Observer + JulDat, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_dateXpred.occur2)

gTOD.models$observer_timeXpred.occur2 <- gmultmix(~ Pred.Occur + I(Pred.Occur^2), ~ 1, ~ Observer + Time, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_timeXpred.occur2)

gTOD.models$date_timeXpred.occur2 <- gmultmix(~ Pred.Occur + I(Pred.Occur^2), ~ 1, ~ JulDat + Time, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$date_timeXpred.occur2)

gTOD.models$observer_date_timeXpred.occur2 <- gmultmix(~ Pred.Occur + I(Pred.Occur^2), ~ 1, ~ Observer + JulDat + Time, umf.gtod, mixture="NB", K=10)
summary(gTOD.models$observer_date_timeXpred.occur2)

## How does habitat suitability measure up
MS.gTOD.models <- modSel(fitList(fits=gTOD.models));MS.gTOD.models



###### Check to see if model predictions make sense ######
## I am particularly concerned that the range of covariate values in the historic data may not match the current data on which
## my models are built. Unrealistic predictions could be made based on this. Need to either use a lower radius model so that the 
## range of percent land cover values is more similar between historic and current models, or check to see if predictions are
## even very different and then use the top model if they aren't.
summary(umf.cr1)

## Make sequence of the variable affecting density
Pred.Occur.seq <- seq(0,1,0.1)

## Now predict with top model
newdat <- data.frame(Pred.Occur = Pred.Occur.seq)
head(newdat)

abundance.pred.pred.occur2.gtod <- predict(gTOD.models$observer_dateXpred.occur2, type = "lambda", newdata = newdat, appendData=TRUE); abundance.pred.pred.occur2.gtod

## Make a figure
jpeg(file.path(output.directory, paste( model.name, '_gtod_predoccur2.jpg', sep="")))

plot(Pred.Occur.seq, abundance.pred.pred.occur2.gtod[, "Predicted"], xlab="Predicted Occurrence Squared (SDM)", ylab="Predicted abundance (300m radius)", type='l', 
     ylim = c(0,3.1), frame=F, lwd=2)
matlines(Pred.Occur.seq, abundance.pred.pred.occur2.gtod[,3:4], lty=1, col="grey", lwd=1)

dev.off()



###### Model Table ######
## Make a fistlist for selection table
MS.gTOD.models <- modSel(fitList(fits=gTOD.models));MS.gTOD.models

## Export Model Table
toExport <- as(MS.gTOD.models, "data.frame")
write.csv(toExport, file.path(output.directory, paste(model.name, "_gTODTable.csv", sep="")))

###### Best Model ######
## Check the goodness of fit of the best model and export summary
# Goodness of fit
pb.gtod <- parboot(gTOD.models$observer_dateXpred.occur2, fitstats, nsim=100, report=5); pb.gtod # Can't do this apparently because of the issues in the model
c.hat.gtod <- pb.gtod@t0[2]/mean(pb.gtod@t.star[,2]); c.hat.gtod

filename <- "gtod_parboot.txt"
capture.output(pb.gtod, file = file.path(output.directory,filename))

# Export best model summary
filename <- "gtod_model_summary.txt"
capture.output(summary(gTOD.models$observer_dateXpred.occur2), file = file.path(output.directory,filename))


# Detection probability after all passes for each site. Could take an average across all sites
p.gtod <-mean(rowSums(getP(TOD.models$observer_dateXpred.occur2))); p.gtod




#########################################################
# 16. Predict habitat suitability to current landscapes
#########################################################

###### Load current and historic rasters ######
## Current Rasters
current.f <- list.files(current.raster.folder,pattern =".tif$", full.names=TRUE)
current.ras <- lapply(current.f,raster) 
current.ras_Stack <- stack(current.ras)
names(current.ras_Stack)

current.ras_Stack.test <- subset(current.ras_Stack, names(training.data.prepped))
names(current.ras_Stack.test)


## Historic rasters
historic.f <- list.files(historic.raster.folder,pattern =".tif$", full.names=TRUE)
historic.ras <- lapply(historic.f,raster) 
historic.ras_Stack <- stack(historic.ras)
names(historic.ras_Stack)

historic.ras_Stack.test <- subset(historic.ras_Stack, names(training.data.prepped))
names(historic.ras_Stack.test)

###### Predict suitability to both landscapes and export ######
## Current
current.preds.zi.raster <- predict(current.ras_Stack.test, gbmModel.zi, n.trees=gbm.perf(gbmModel.zi), type="response") 

writeRaster(current.preds.zi.raster, file.path(output.directory, paste(model.name,"_current_suitability.tif", sep="")), 
            format="GTiff", overwrite=TRUE) 
writeRaster(current.preds.zi.raster, file.path(output.directory, paste(model.name,"_current_suitability.asc", sep="")), 
            format="ascii", overwrite=TRUE) 

# Convert to 1 0 for suitable habitat
current.preds.zi.raster_01 <- current.preds.zi.raster
current.preds.zi.raster_01[current.preds.zi.raster_01>calibration_cutoff] <- 1
current.preds.zi.raster_01[current.preds.zi.raster_01<=calibration_cutoff] <- 0
plot(current.preds.zi.raster_01)

writeRaster(current.preds.zi.raster_01, file.path(output.directory, paste(model.name,"_current_01suitable.tif", sep="")), 
            format="GTiff" , overwrite=TRUE) 
writeRaster(current.preds.zi.raster_01, file.path(output.directory, paste(model.name,"_current_01suitable.asc", sep="")), 
            format="ascii", overwrite=TRUE) 

## Historic
historic.preds.zi.raster <- predict(historic.ras_Stack.test, gbmModel.zi, n.trees=gbm.perf(gbmModel.zi), type="response") 

writeRaster(historic.preds.zi.raster, file.path(output.directory, paste(model.name,"_historic_suitability.tif", sep="")), 
            format="GTiff", overwrite=TRUE) 
writeRaster(historic.preds.zi.raster, file.path(output.directory, paste(model.name,"_historic_suitability.asc", sep="")), 
            format="ascii", overwrite=TRUE) # Need to change the path here

# Convert to 1 0 for suitable habitat
historic.preds.zi.raster_01 <- historic.preds.zi.raster
historic.preds.zi.raster_01[historic.preds.zi.raster_01>calibration_cutoff] <- 1
historic.preds.zi.raster_01[historic.preds.zi.raster_01<=calibration_cutoff] <- 0
plot(historic.preds.zi.raster_01)

writeRaster(historic.preds.zi.raster_01, file.path(output.directory, paste(model.name,"_historic_01suitable.tif", sep="")), 
            format="GTiff", overwrite=TRUE)
writeRaster(historic.preds.zi.raster_01, file.path(output.directory, paste(model.name,"_historic_01suitable.asc", sep="")), 
            format="ascii", overwrite=TRUE) # Need to change the path here


###### Calculate suitable area after threshold for both landscape ######
## Load shapefile
willamette_shape <- shapefile(willamette.extent.shp)

## Current
current.suitable.habitat <- mask(current.preds.zi.raster_01, willamette_shape)
current.suitable.habitat.cells <- cellStats(current.suitable.habitat, sum)
current.suitable.habitat.area.hectares <- (current.suitable.habitat.cells*900)/10000
plot(current.suitable.habitat)

## Historic
historic.suitable.habitat <- mask(historic.preds.zi.raster_01, willamette_shape)
historic.suitable.habitat.cells <- cellStats(historic.suitable.habitat, sum)
historic.suitable.habitat.area.hectares <- (historic.suitable.habitat.cells*900)/10000
plot(historic.suitable.habitat)

## Calculate estimated change in habitat area
suitable.area.change <- historic.suitable.habitat.area.hectares - current.suitable.habitat.area.hectares
suitable.area.change.percent <- (current.suitable.habitat.area.hectares/historic.suitable.habitat.area.hectares)*100; suitable.area.change.percent



#########################################################
# 14. Predict density to landscapes and estimate population 
#########################################################

###### Predict densities to landscapes with unmarked ######

## Prep necessary rasters
# Need to use Arcmap to convert cell size and make polygons of just suitable habitat to increase speed of prediction. Otherwise it was going to take
# a year or so to predict to the landscape

## Prep Current
#current.pred.occur <- raster(file.path(output.directory, paste(model.name,"_current_suitability.tif", sep=""))) 
current.pred.occur <- raster(file.path(output.directory, paste(model.name,"_current_suitability180.tif", sep=""))) # This is the aggregated raster to 180 m cell size to speed things up
current_suitable <- shapefile(file.path(output.directory, "wbnu_suitable_poly_current_clipped.shp")) # This is the polygon of just suitable areas

current.pred.occur <- mask(current.pred.occur, current_suitable)
current.stack <- stack(current.pred.occur)
names(current.stack) <- "Pred.Occur"
names(current.stack)

## Prep historic
#historic.pred.occur <- raster(file.path(output.directory, paste(model.name,"_historic_suitability.tif", sep=""))) 
historic.pred.occur <- raster(file.path(output.directory, paste(model.name,"_historic_suitability180.tif", sep=""))) # This is the aggregated raster to 180 m cell size to speed things up
historic_suitable <- shapefile(file.path(output.directory, "wbnu_suitable_poly_historic_clipped.shp")) # This is the polygon of just suitable areas

historic.pred.occur <- mask(historic.pred.occur, historic_suitable)
historic.stack <- stack(historic.pred.occur)
names(historic.stack) <- "Pred.Occur"
names(historic.stack)


## Density with distsamp top model
distsamp.density.current <- predict(dist.models$observerXpred.occur2, type = "state", newdata = current.stack) 

writeRaster(subset(distsamp.density.current,subset=1), file.path(output.directory, "wbnu_distsamp_current_density180.tif"), format="GTiff", overwrite=TRUE) 
writeRaster(subset(distsamp.density.current,subset=2), file.path(output.directory, "wbnu_distsamp_current_density180_se.tif"), format="GTiff", overwrite=TRUE) 
writeRaster(subset(distsamp.density.current,subset=3), file.path(output.directory, "wbnu_distsamp_current_density180_lower.tif"), format="GTiff", overwrite=TRUE) 
writeRaster(subset(distsamp.density.current,subset=4), file.path(output.directory, "wbnu_distsamp_current_density180_upper.tif"), format="GTiff", overwrite=TRUE) 


distsamp.density.historic <- predict(dist.models$observerXpred.occur2, type = "state", newdata = historic.stack)

writeRaster(subset(distsamp.density.historic,subset=1), file.path(output.directory, "wbnu_distsamp_historic_density180.tif"), format="GTiff", overwrite=TRUE) 
writeRaster(subset(distsamp.density.historic,subset=2), file.path(output.directory, "wbnu_distsamp_historic_density180_se.tif"), format="GTiff", overwrite=TRUE) 
writeRaster(subset(distsamp.density.historic,subset=3), file.path(output.directory, "wbnu_distsamp_historic_density180_lower.tif"), format="GTiff", overwrite=TRUE) 
writeRaster(subset(distsamp.density.historic,subset=4), file.path(output.directory, "wbnu_distsamp_historic_density180_upper.tif"), format="GTiff", overwrite=TRUE) 


### Density with gdistsamp (negative binomial) top model
#gdistsamp.density.current <- predict(gdist.models$nullXpred.occur2, type = "lambda", newdata = current.stack)

#writeRaster(subset(gdistsamp.density.current,subset=1), file.path(output.directory, "wbnu_gdistsamp_current_density180.tif"), format="GTiff", overwrite=TRUE) 
#writeRaster(subset(gdistsamp.density.current,subset=2), file.path(output.directory, "wbnu_gdistsamp_current_density180_se.tif"), format="GTiff", overwrite=TRUE) 
#writeRaster(subset(gdistsamp.density.current,subset=3), file.path(output.directory, "wbnu_gdistsamp_current_density180_lower.tif"), format="GTiff", overwrite=TRUE) 
#writeRaster(subset(gdistsamp.density.current,subset=4), file.path(output.directory, "wbnu_gdistsamp_current_density180_upper.tif"), format="GTiff", overwrite=TRUE) 

#gdistsamp.density.historic <- predict(gdist.models$nullXpred.occur2, type = "lambda", newdata = historic.stack) 

#writeRaster(subset(gdistsamp.density.historic,subset=1), file.path(output.directory, "wbnu_gdistsamp_historic_density180.tif"), format="GTiff", overwrite=TRUE) 
#writeRaster(subset(gdistsamp.density.historic,subset=2), file.path(output.directory, "wbnu_gdistsamp_historic_density180_se.tif"), format="GTiff", overwrite=TRUE) 
#writeRaster(subset(gdistsamp.density.historic,subset=3), file.path(output.directory, "wbnu_gdistsamp_historic_density180_lower.tif"), format="GTiff", overwrite=TRUE) 
#writeRaster(subset(gdistsamp.density.historic,subset=4), file.path(output.directory, "wbnu_gdistsamp_historic_density180_upper.tif"), format="GTiff", overwrite=TRUE)


## Density with removal top model
# Need to do some calculations to convert abundances to densities
point.count.survey.area <- pi*TruncationD^2 # This is the area surveyed by one point count using a 300m truncation distance

rem.abundance.current <- predict(Rem.models$observerXpred.occur2, type = "state", newdata = current.stack) 

prediction <- subset(rem.abundance.current,subset=1)
prediction <- (prediction/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

se <- subset(rem.abundance.current,subset=2)
se <- (se/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

lower <- subset(rem.abundance.current,subset=3)
lower <- (lower/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

upper <- subset(rem.abundance.current,subset=4)
upper <- (upper/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

writeRaster(prediction, file.path(output.directory, "wbnu_rem_current_density180.tif"), format="GTiff", overwrite=TRUE) 
writeRaster(se, file.path(output.directory, "wbnu_rem_current_density180_se.tif"), format="GTiff", overwrite=TRUE) 
writeRaster(lower, file.path(output.directory, "wbnu_rem_current_density180_lower.tif"), format="GTiff", overwrite=TRUE) 
writeRaster(upper, file.path(output.directory, "wbnu_rem_current_density180_upper.tif"), format="GTiff", overwrite=TRUE)

rem.abundance.historic <- predict(Rem.models$observerXpred.occur2, type = "state", newdata = historic.stack) 

prediction <- subset(rem.abundance.historic,subset=1)
prediction <- (prediction/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

se <- subset(rem.abundance.historic,subset=2)
se <- (se/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

lower <- subset(rem.abundance.historic,subset=3)
lower <- (lower/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

upper <- subset(rem.abundance.historic,subset=4)
upper <- (upper/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

writeRaster(prediction, file.path(output.directory, "wbnu_rem_historic_density180.tif"), format="GTiff", overwrite=TRUE) 
writeRaster(se, file.path(output.directory, "wbnu_rem_historic_density180_se.tif"), format="GTiff", overwrite=TRUE) 
writeRaster(lower, file.path(output.directory, "wbnu_rem_historic_density180_lower.tif"), format="GTiff", overwrite=TRUE) 
writeRaster(upper, file.path(output.directory, "wbnu_rem_historic_density180_upper.tif"), format="GTiff", overwrite=TRUE)


### Density with negative binomial removal top model
#grem.abundance.current <- predict(gRem.models$observer_dateXpred.occur2, type = "lambda", newdata = current.stack) 

#prediction <- subset(grem.abundance.current,subset=1)
#prediction <- (prediction/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

#se <- subset(grem.abundance.current,subset=2)
#se <- (se/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

#lower <- subset(grem.abundance.current,subset=3)
#lower <- (lower/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

#upper <- subset(grem.abundance.current,subset=4)
#upper <- (upper/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

#writeRaster(prediction, file.path(output.directory, "wbnu_grem_current_density180.tif"), format="GTiff", overwrite=TRUE) 
#writeRaster(se, file.path(output.directory, "wbnu_grem_current_density180_se.tif"), format="GTiff", overwrite=TRUE) 
#writeRaster(lower, file.path(output.directory, "wbnu_grem_current_density180_lower.tif"), format="GTiff", overwrite=TRUE) 
#writeRaster(upper, file.path(output.directory, "wbnu_grem_current_density180_upper.tif"), format="GTiff", overwrite=TRUE)

#grem.abundance.historic <- predict(gRem.models$observer_dateXpred.occur2, type = "lambda", newdata = historic.stack)

#prediction <- subset(grem.abundance.historic,subset=1)
#prediction <- (prediction/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

#se <- subset(grem.abundance.historic,subset=2)
#se <- (se/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

#lower <- subset(grem.abundance.historic,subset=3)
#lower <- (lower/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

#upper <- subset(grem.abundance.historic,subset=4)
#upper <- (upper/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

#writeRaster(prediction, file.path(output.directory, "wbnu_grem_historic_density180.tif"), format="GTiff", overwrite=TRUE) 
#writeRaster(se, file.path(output.directory, "wbnu_grem_historic_density180_se.tif"), format="GTiff", overwrite=TRUE) 
#writeRaster(lower, file.path(output.directory, "wbnu_grem_historic_density180_lower.tif"), format="GTiff", overwrite=TRUE) 
#writeRaster(upper, file.path(output.directory, "wbnu_grem_historic_density180_upper.tif"), format="GTiff", overwrite=TRUE)


## Density with time of detection top model
tod.density.current <- predict(TOD.models$observerXpred.occur2, type = "state", newdata = current.stack) 

prediction <- subset(tod.density.current,subset=1)
prediction <- (prediction/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

se <- subset(tod.density.current,subset=2)
se <- (se/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

lower <- subset(tod.density.current,subset=3)
lower <- (lower/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

upper <- subset(tod.density.current,subset=4)
upper <- (upper/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

writeRaster(prediction, file.path(output.directory, "wbnu_tod_current_density180.tif"), format="GTiff", overwrite=TRUE) 
writeRaster(se, file.path(output.directory, "wbnu_tod_current_density180_se.tif"), format="GTiff", overwrite=TRUE) 
writeRaster(lower, file.path(output.directory, "wbnu_tod_current_density180_lower.tif"), format="GTiff", overwrite=TRUE) 
writeRaster(upper, file.path(output.directory, "wbnu_tod_current_density180_upper.tif"), format="GTiff", overwrite=TRUE)

tod.density.historic <- predict(TOD.models$observerXpred.occur2, type = "state", newdata = historic.stack) 

prediction <- subset(tod.density.historic,subset=1)
prediction <- (prediction/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

se <- subset(tod.density.historic,subset=2)
se <- (se/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

lower <- subset(tod.density.historic,subset=3)
lower <- (lower/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

upper <- subset(tod.density.historic,subset=4)
upper <- (upper/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

writeRaster(prediction, file.path(output.directory, "wbnu_tod_historic_density180.tif"), format="GTiff", overwrite=TRUE) 
writeRaster(se, file.path(output.directory, "wbnu_tod_historic_density180_se.tif"), format="GTiff", overwrite=TRUE) 
writeRaster(lower, file.path(output.directory, "wbnu_tod_historic_density180_lower.tif"), format="GTiff", overwrite=TRUE) 
writeRaster(upper, file.path(output.directory, "wbnu_tod_historic_density180_upper.tif"), format="GTiff", overwrite=TRUE)



## Density with negative binomial time of detection top model
#gtod.density.current <- predict(gTOD.models$observer_dateXpred.occur2, type = "lambda", newdata = current.stack) 

#prediction <- subset(gtod.density.current,subset=1)
#prediction <- (prediction/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

#se <- subset(gtod.density.current,subset=2)
#se <- (se/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

#lower <- subset(gtod.density.current,subset=3)
#lower <- (lower/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

#upper <- subset(gtod.density.current,subset=4)
#upper <- (upper/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

#writeRaster(prediction, file.path(output.directory, "wbnu_gtod_current_density180.tif"), format="GTiff", overwrite=TRUE) 
#writeRaster(se, file.path(output.directory, "wbnu_gtod_current_density180_se.tif"), format="GTiff", overwrite=TRUE) 
#writeRaster(lower, file.path(output.directory, "wbnu_gtod_current_density180_lower.tif"), format="GTiff", overwrite=TRUE) 
#writeRaster(upper, file.path(output.directory, "wbnu_gtod_current_density180_upper.tif"), format="GTiff", overwrite=TRUE)

#gtod.density.historic <- predict(gTOD.models$observer_dateXpred.occur2, type = "lambda", newdata = historic.stack) 

#prediction <- subset(gtod.density.historic,subset=1)
#prediction <- (prediction/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

#se <- subset(gtod.density.historic,subset=2)
#se <- (se/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

#lower <- subset(gtod.density.historic,subset=3)
#lower <- (lower/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

#upper <- subset(gtod.density.historic,subset=4)
#upper <- (upper/point.count.survey.area) * 10000 # Converts abundance to per meter density then to per hectare density

#writeRaster(prediction, file.path(output.directory, "wbnu_gtod_historic_density180.tif"), format="GTiff", overwrite=TRUE) 
#writeRaster(se, file.path(output.directory, "wbnu_gtod_historic_density180_se.tif"), format="GTiff", overwrite=TRUE) 
#writeRaster(lower, file.path(output.directory, "wbnu_gtod_historic_density180_lower.tif"), format="GTiff", overwrite=TRUE) 
#writeRaster(upper, file.path(output.directory, "wbnu_gtod_historic_density180_upper.tif"), format="GTiff", overwrite=TRUE)

###### Estimate population on both landscapes ######

cell.area <- 180^2 # This is the area of each 180 m X 180 m cell in raster so that we can convert to numbers

## distsamp
# Current
current.dens.prediction <- raster(file.path(output.directory, "wbnu_distsamp_current_density180.tif")) 
mean.density.current.distsamp <- cellStats(current.dens.prediction, "mean"); mean.density.current.distsamp
current.dens.prediction <- (current.dens.prediction/10000)*cell.area
cellStats(current.dens.prediction, "max")
cellStats(current.dens.prediction, "mean")
population.current.distsamp <- cellStats(current.dens.prediction, sum); population.current.distsamp

current.dens.lower <- raster(file.path(output.directory, "wbnu_distsamp_current_density180_lower.tif")) 
mean.density.current.lower.distsamp <- cellStats(current.dens.lower, "mean"); mean.density.current.lower.distsamp
current.dens.lower <- (current.dens.lower/10000)*cell.area
cellStats(current.dens.lower, "max")
cellStats(current.dens.lower, "mean")
population.current.lower.distsamp <- cellStats(current.dens.lower, sum); population.current.lower.distsamp

current.dens.upper <- raster(file.path(output.directory, "wbnu_distsamp_current_density180_upper.tif")) 
mean.density.current.upper.distsamp <- cellStats(current.dens.upper, "mean"); mean.density.current.upper.distsamp
current.dens.upper <- (current.dens.upper/10000)*cell.area
cellStats(current.dens.upper, "max")
cellStats(current.dens.upper, "mean")
population.current.upper.distsamp <- cellStats(current.dens.upper, sum); population.current.upper.distsamp

# Historic
historic.dens.prediction <- raster(file.path(output.directory, "wbnu_distsamp_historic_density180.tif")) 
mean.density.historic.distsamp <- cellStats(historic.dens.prediction, "mean"); mean.density.historic.distsamp
historic.dens.prediction <- (historic.dens.prediction/10000)*cell.area
cellStats(historic.dens.prediction, "max")
cellStats(historic.dens.prediction, "mean")
population.historic.distsamp <- cellStats(historic.dens.prediction, sum); population.historic.distsamp

historic.dens.lower <- raster(file.path(output.directory, "wbnu_distsamp_historic_density180_lower.tif")) 
mean.density.historic.lower.distsamp <- cellStats(historic.dens.lower, "mean"); mean.density.historic.lower.distsamp
historic.dens.lower <- (historic.dens.lower/10000)*cell.area
cellStats(historic.dens.lower, "max")
cellStats(historic.dens.lower, "mean")
population.historic.lower.distsamp <- cellStats(historic.dens.lower, sum); population.historic.lower.distsamp

historic.dens.upper <- raster(file.path(output.directory, "wbnu_distsamp_historic_density180_upper.tif")) 
mean.density.historic.upper.distsamp <- cellStats(historic.dens.upper, "mean"); mean.density.historic.upper.distsamp
historic.dens.upper <- (historic.dens.upper/10000)*cell.area
cellStats(historic.dens.upper, "max")
cellStats(historic.dens.upper, "mean")
population.historic.upper.distsamp <- cellStats(historic.dens.upper, sum); population.historic.upper.distsamp

# Change
population.change.distsamp <- population.historic.distsamp - population.current.distsamp; population.change.distsamp
population.change.lower.distsamp <- population.historic.lower.distsamp - population.current.lower.distsamp; population.change.lower.distsamp
population.change.upper.distsamp <- population.historic.upper.distsamp - population.current.upper.distsamp; population.change.upper.distsamp
population.change.min.distsamp <- population.historic.lower.distsamp - population.current.upper.distsamp; population.change.min.distsamp
population.change.percent.distsamp <- (population.current.distsamp/population.historic.distsamp)*100; population.change.percent.distsamp


### gdistsamp
## Current
#current.dens.prediction <- raster(file.path(output.directory, "wbnu_gdistsamp_current_density180.tif")) 
#mean.density.current.gdistsamp <- cellStats(current.dens.prediction, "mean"); mean.density.current.gdistsamp
#current.dens.prediction <- (current.dens.prediction/10000)*cell.area
#cellStats(current.dens.prediction, "max")
#cellStats(current.dens.prediction, "mean")
#population.current.gdistsamp <- cellStats(current.dens.prediction, sum); population.current.gdistsamp

#current.dens.lower <- raster(file.path(output.directory, "wbnu_gdistsamp_current_density180_lower.tif")) 
#mean.density.current.lower.gdistsamp <- cellStats(current.dens.lower, "mean"); mean.density.current.lower.gdistsamp
#current.dens.lower <- (current.dens.lower/10000)*cell.area
#cellStats(current.dens.lower, "max")
#cellStats(current.dens.lower, "mean")
#population.current.lower.gdistsamp <- cellStats(current.dens.lower, sum); population.current.lower.gdistsamp

#current.dens.upper <- raster(file.path(output.directory, "wbnu_gdistsamp_current_density180_upper.tif")) 
#mean.density.current.upper.gdistsamp <- cellStats(current.dens.upper, "mean"); mean.density.current.upper.gdistsamp
#current.dens.upper <- (current.dens.upper/10000)*cell.area
#cellStats(current.dens.upper, "max")
#cellStats(current.dens.upper, "mean")
#population.current.upper.gdistsamp <- cellStats(current.dens.upper, sum); population.current.upper.gdistsamp

## Historic
#historic.dens.prediction <- raster(file.path(output.directory, "wbnu_gdistsamp_historic_density180.tif")) 
#mean.density.historic.gdistsamp <- cellStats(historic.dens.prediction, "mean"); mean.density.historic.gdistsamp
#historic.dens.prediction <- (historic.dens.prediction/10000)*cell.area
#cellStats(historic.dens.prediction, "max")
#cellStats(historic.dens.prediction, "mean")
#population.historic.gdistsamp <- cellStats(historic.dens.prediction, sum); population.historic.gdistsamp

#historic.dens.lower <- raster(file.path(output.directory, "wbnu_gdistsamp_historic_density180_lower.tif")) 
#mean.density.historic.lower.gdistsamp <- cellStats(historic.dens.lower, "mean"); mean.density.historic.lower.gdistsamp
#historic.dens.lower <- (historic.dens.lower/10000)*cell.area
#cellStats(historic.dens.lower, "max")
#cellStats(historic.dens.lower, "mean")
#population.historic.lower.gdistsamp <- cellStats(historic.dens.lower, sum); population.historic.lower.gdistsamp

#historic.dens.upper <- raster(file.path(output.directory, "wbnu_gdistsamp_historic_density180_upper.tif")) 
#mean.density.historic.upper.gdistsamp <- cellStats(historic.dens.upper, "mean"); mean.density.historic.upper.gdistsamp
#historic.dens.upper <- (historic.dens.upper/10000)*cell.area
#cellStats(historic.dens.upper, "max")
#cellStats(historic.dens.upper, "mean")
#population.historic.upper.gdistsamp <- cellStats(historic.dens.upper, sum); population.historic.upper.gdistsamp

## Change
#population.change.gdistsamp <- population.historic.gdistsamp - population.current.gdistsamp; population.change.gdistsamp
#population.change.lower.gdistsamp <- population.historic.lower.gdistsamp - population.current.lower.gdistsamp; population.change.lower.gdistsamp
#population.change.upper.gdistsamp <- population.historic.upper.gdistsamp - population.current.upper.gdistsamp; population.change.upper.gdistsamp
#population.change.min.gdistsamp <- population.historic.lower.gdistsamp - population.current.upper.gdistsamp; population.change.min.gdistsamp
#population.change.percent.gdistsamp <- (population.current.gdistsamp/population.historic.gdistsamp)*100; population.change.percent.gdistsamp



## rem
# Current
current.dens.prediction <- raster(file.path(output.directory, "wbnu_rem_current_density180.tif")) 
mean.density.current.rem <- cellStats(current.dens.prediction, "mean"); mean.density.current.rem
current.dens.prediction <- (current.dens.prediction/10000)*cell.area
cellStats(current.dens.prediction, "max")
cellStats(current.dens.prediction, "mean")
population.current.rem <- cellStats(current.dens.prediction, sum); population.current.rem

current.dens.lower <- raster(file.path(output.directory, "wbnu_rem_current_density180_lower.tif")) 
mean.density.current.lower.rem <- cellStats(current.dens.lower, "mean"); mean.density.current.lower.rem
current.dens.lower <- (current.dens.lower/10000)*cell.area
cellStats(current.dens.lower, "max")
cellStats(current.dens.lower, "mean")
population.current.lower.rem <- cellStats(current.dens.lower, sum); population.current.lower.rem

current.dens.upper <- raster(file.path(output.directory, "wbnu_rem_current_density180_upper.tif")) 
mean.density.current.upper.rem <- cellStats(current.dens.upper, "mean"); mean.density.current.upper.rem
current.dens.upper <- (current.dens.upper/10000)*cell.area
cellStats(current.dens.upper, "max")
cellStats(current.dens.upper, "mean")
population.current.upper.rem <- cellStats(current.dens.upper, sum); population.current.upper.rem

# Historic
historic.dens.prediction <- raster(file.path(output.directory, "wbnu_rem_historic_density180.tif")) 
mean.density.historic.rem <- cellStats(historic.dens.prediction, "mean"); mean.density.historic.rem
historic.dens.prediction <- (historic.dens.prediction/10000)*cell.area
cellStats(historic.dens.prediction, "max")
cellStats(historic.dens.prediction, "mean")
population.historic.rem <- cellStats(historic.dens.prediction, sum); population.historic.rem

historic.dens.lower <- raster(file.path(output.directory, "wbnu_rem_historic_density180_lower.tif")) 
mean.density.historic.lower.rem <- cellStats(historic.dens.lower, "mean"); mean.density.historic.lower.rem
historic.dens.lower <- (historic.dens.lower/10000)*cell.area
cellStats(historic.dens.lower, "max")
cellStats(historic.dens.lower, "mean")
population.historic.lower.rem <- cellStats(historic.dens.lower, sum); population.historic.lower.rem

historic.dens.upper <- raster(file.path(output.directory, "wbnu_rem_historic_density180_upper.tif")) 
mean.density.historic.upper.rem <- cellStats(historic.dens.upper, "mean"); mean.density.historic.upper.rem
historic.dens.upper <- (historic.dens.upper/10000)*cell.area
cellStats(historic.dens.upper, "max")
cellStats(historic.dens.upper, "mean")
population.historic.upper.rem <- cellStats(historic.dens.upper, sum); population.historic.upper.rem

# Change
population.change.rem <- population.historic.rem - population.current.rem; population.change.rem
population.change.lower.rem <- population.historic.lower.rem - population.current.lower.rem; population.change.lower.rem
population.change.upper.rem <- population.historic.upper.rem - population.current.upper.rem; population.change.upper.rem
population.change.min.rem <- population.historic.lower.rem - population.current.upper.rem; population.change.min.rem
population.change.percent.rem <- (population.current.rem/population.historic.rem)*100; population.change.percent.rem



### grem
## Current
#current.dens.prediction <- raster(file.path(output.directory, "wbnu_grem_current_density180.tif")) 
#mean.density.current.grem <- cellStats(current.dens.prediction, "mean"); mean.density.current.grem
#current.dens.prediction <- (current.dens.prediction/10000)*cell.area
#cellStats(current.dens.prediction, "max")
#cellStats(current.dens.prediction, "mean")
#population.current.grem <- cellStats(current.dens.prediction, sum); population.current.grem

#current.dens.lower <- raster(file.path(output.directory, "wbnu_grem_current_density180_lower.tif")) 
#mean.density.current.lower.grem <- cellStats(current.dens.lower, "mean"); mean.density.current.lower.grem
#current.dens.lower <- (current.dens.lower/10000)*cell.area
#cellStats(current.dens.lower, "max")
#cellStats(current.dens.lower, "mean")
#population.current.lower.grem <- cellStats(current.dens.lower, sum); population.current.lower.grem

#current.dens.upper <- raster(file.path(output.directory, "wbnu_grem_current_density180_upper.tif")) 
#mean.density.current.upper.grem <- cellStats(current.dens.upper, "mean"); mean.density.current.upper.grem
#current.dens.upper <- (current.dens.upper/10000)*cell.area
#cellStats(current.dens.upper, "max")
#cellStats(current.dens.upper, "mean")
#population.current.upper.grem <- cellStats(current.dens.upper, sum); population.current.upper.grem

## Historic
#historic.dens.prediction <- raster(file.path(output.directory, "wbnu_grem_historic_density180.tif")) 
#mean.density.historic.grem <- cellStats(historic.dens.prediction, "mean"); mean.density.historic.grem
#historic.dens.prediction <- (historic.dens.prediction/10000)*cell.area
#cellStats(historic.dens.prediction, "max")
#cellStats(historic.dens.prediction, "mean")
#population.historic.grem <- cellStats(historic.dens.prediction, sum); population.historic.grem

#historic.dens.lower <- raster(file.path(output.directory, "wbnu_grem_historic_density180_lower.tif")) 
#mean.density.historic.lower.grem <- cellStats(historic.dens.lower, "mean"); mean.density.historic.lower.grem
#historic.dens.lower <- (historic.dens.lower/10000)*cell.area
#cellStats(historic.dens.lower, "max")
#cellStats(historic.dens.lower, "mean")
#population.historic.lower.grem <- cellStats(historic.dens.lower, sum); population.historic.lower.grem

#historic.dens.upper <- raster(file.path(output.directory, "wbnu_grem_historic_density180_upper.tif")) 
#mean.density.historic.upper.grem <- cellStats(historic.dens.upper, "mean"); mean.density.historic.upper.grem
#historic.dens.upper <- (historic.dens.upper/10000)*cell.area
#cellStats(historic.dens.upper, "max")
#cellStats(historic.dens.upper, "mean")
#population.historic.upper.grem <- cellStats(historic.dens.upper, sum); population.historic.upper.grem

## Change
#population.change.grem <- population.historic.grem - population.current.grem; population.change.grem
#population.change.lower.grem <- population.historic.lower.grem - population.current.lower.grem; population.change.lower.grem
#population.change.upper.grem <- population.historic.upper.grem - population.current.upper.grem; population.change.upper.grem
#population.change.min.grem <- population.historic.lower.grem - population.current.upper.grem; population.change.min.grem
#population.change.percent.grem <- (population.current.grem/population.historic.grem)*100; population.change.percent.grem



## tod
# Current
current.dens.prediction <- raster(file.path(output.directory, "wbnu_tod_current_density180.tif")) 
mean.density.current.tod <- cellStats(current.dens.prediction, "mean"); mean.density.current.tod
current.dens.prediction <- (current.dens.prediction/10000)*cell.area
cellStats(current.dens.prediction, "max")
cellStats(current.dens.prediction, "mean")
population.current.tod <- cellStats(current.dens.prediction, sum); population.current.tod

current.dens.lower <- raster(file.path(output.directory, "wbnu_tod_current_density180_lower.tif")) 
mean.density.current.lower.tod <- cellStats(current.dens.lower, "mean"); mean.density.current.lower.tod
current.dens.lower <- (current.dens.lower/10000)*cell.area
cellStats(current.dens.lower, "max")
cellStats(current.dens.lower, "mean")
population.current.lower.tod <- cellStats(current.dens.lower, sum); population.current.lower.tod

current.dens.upper <- raster(file.path(output.directory, "wbnu_tod_current_density180_upper.tif")) 
mean.density.current.upper.tod <- cellStats(current.dens.upper, "mean"); mean.density.current.upper.tod
current.dens.upper <- (current.dens.upper/10000)*cell.area
cellStats(current.dens.upper, "max")
cellStats(current.dens.upper, "mean")
population.current.upper.tod <- cellStats(current.dens.upper, sum); population.current.upper.tod

# Historic
historic.dens.prediction <- raster(file.path(output.directory, "wbnu_tod_historic_density180.tif")) 
mean.density.historic.tod <- cellStats(historic.dens.prediction, "mean"); mean.density.historic.tod
historic.dens.prediction <- (historic.dens.prediction/10000)*cell.area
cellStats(historic.dens.prediction, "max")
cellStats(historic.dens.prediction, "mean")
population.historic.tod <- cellStats(historic.dens.prediction, sum); population.historic.tod

historic.dens.lower <- raster(file.path(output.directory, "wbnu_tod_historic_density180_lower.tif")) 
mean.density.historic.lower.tod <- cellStats(historic.dens.lower, "mean"); mean.density.historic.lower.tod
historic.dens.lower <- (historic.dens.lower/10000)*cell.area
cellStats(historic.dens.lower, "max")
cellStats(historic.dens.lower, "mean")
population.historic.lower.tod <- cellStats(historic.dens.lower, sum); population.historic.lower.tod

historic.dens.upper <- raster(file.path(output.directory, "wbnu_tod_historic_density180_upper.tif")) 
mean.density.historic.upper.tod <- cellStats(historic.dens.upper, "mean"); mean.density.historic.upper.tod
historic.dens.upper <- (historic.dens.upper/10000)*cell.area
cellStats(historic.dens.upper, "max")
cellStats(historic.dens.upper, "mean")
population.historic.upper.tod <- cellStats(historic.dens.upper, sum); population.historic.upper.tod

# Change
population.change.tod <- population.historic.tod - population.current.tod; population.change.tod
population.change.lower.tod <- population.historic.lower.tod - population.current.lower.tod; population.change.lower.tod
population.change.upper.tod <- population.historic.upper.tod - population.current.upper.tod; population.change.upper.tod
population.change.min.tod <- population.historic.lower.tod - population.current.upper.tod; population.change.min.tod
population.change.percent.tod <- (population.current.tod/population.historic.tod)*100; population.change.percent.tod



### gtod
## Current
#current.dens.prediction <- raster(file.path(output.directory, "wbnu_gtod_current_density180.tif")) 
#mean.density.current.gtod <- cellStats(current.dens.prediction, "mean"); mean.density.current.gtod
#current.dens.prediction <- (current.dens.prediction/10000)*cell.area
#cellStats(current.dens.prediction, "max")
#cellStats(current.dens.prediction, "mean")
#population.current.gtod <- cellStats(current.dens.prediction, sum); population.current.gtod

#current.dens.lower <- raster(file.path(output.directory, "wbnu_gtod_current_density180_lower.tif")) 
#mean.density.current.lower.gtod <- cellStats(current.dens.lower, "mean"); mean.density.current.lower.gtod
#current.dens.lower <- (current.dens.lower/10000)*cell.area
#cellStats(current.dens.lower, "max")
#cellStats(current.dens.lower, "mean")
#population.current.lower.gtod <- cellStats(current.dens.lower, sum); population.current.lower.gtod

#current.dens.upper <- raster(file.path(output.directory, "wbnu_gtod_current_density180_upper.tif")) 
#mean.density.current.upper.gtod <- cellStats(current.dens.upper, "mean"); mean.density.current.upper.gtod
#current.dens.upper <- (current.dens.upper/10000)*cell.area
#cellStats(current.dens.upper, "max")
#cellStats(current.dens.upper, "mean")
#population.current.upper.gtod <- cellStats(current.dens.upper, sum); population.current.upper.gtod

## Historic
#historic.dens.prediction <- raster(file.path(output.directory, "wbnu_gtod_historic_density180.tif")) 
#mean.density.historic.gtod <- cellStats(historic.dens.prediction, "mean"); mean.density.historic.gtod
#historic.dens.prediction <- (historic.dens.prediction/10000)*cell.area
#cellStats(historic.dens.prediction, "max")
#cellStats(historic.dens.prediction, "mean")
#population.historic.gtod <- cellStats(historic.dens.prediction, sum); population.historic.gtod

#historic.dens.lower <- raster(file.path(output.directory, "wbnu_gtod_historic_density180_lower.tif")) 
#mean.density.historic.lower.gtod <- cellStats(historic.dens.lower, "mean"); mean.density.historic.lower.gtod
#historic.dens.lower <- (historic.dens.lower/10000)*cell.area
#cellStats(historic.dens.lower, "max")
#cellStats(historic.dens.lower, "mean")
#population.historic.lower.gtod <- cellStats(historic.dens.lower, sum); population.historic.lower.gtod

#historic.dens.upper <- raster(file.path(output.directory, "wbnu_gtod_historic_density180_upper.tif")) 
#mean.density.historic.upper.gtod <- cellStats(historic.dens.upper, "mean"); mean.density.historic.upper.gtod
#historic.dens.upper <- (historic.dens.upper/10000)*cell.area
#cellStats(historic.dens.upper, "max")
#cellStats(historic.dens.upper, "mean")
#population.historic.upper.gtod <- cellStats(historic.dens.upper, sum); population.historic.upper.gtod

## Change
#population.change.gtod <- population.historic.gtod - population.current.gtod; population.change.gtod
#population.change.lower.gtod <- population.historic.lower.gtod - population.current.lower.gtod; population.change.lower.gtod
#population.change.upper.gtod <- population.historic.upper.gtod - population.current.upper.gtod; population.change.upper.gtod
#population.change.min.gtod <- population.historic.lower.gtod - population.current.upper.gtod; population.change.min.gtod
#population.change.percent.gtod <- (population.current.gtod/population.historic.gtod)*100; population.change.percent.gtod



#########################################################
# 15. Calculate and gather all output
#########################################################


###### Run Moran's I on abundance model residuals ######
xy <- cbind(training.data$coords.x1,training.data$coords.x2)

# Create a distance matrix between points, take the inverse and replace diagonals with 0
training.data.dists <- as.matrix(dist(xy))

# Technically this is an inverse matrix
training.data.dists.inv <- 1/training.data.dists

# Need to set the diagonals to 0 from Inf so that the Moran's function can work
diag(training.data.dists.inv) <- 0

# Examine to make sure it looks good
training.data.dists.inv[1:5, 1:5]

###### If there are repeat counts in the same location (happens on eBird), then distance should be 0 but will be ######
###### calculated as Inf. The code below changes all Inf to 0, so that Moran's I can be run ######
training.data.dists.inv[!is.finite(training.data.dists.inv)] <- 0

# Calculate residuals as there is apparently no easy way to do this from gbm models
residual.training <- training.data$Count - gbmModel.count$fit

# Now calculate Moran's I from the model residuals
moran <- Moran.I(residual.training, training.data.dists.inv)
moran.export <- c(moran$observed, moran$expected, moran$sd, moran$p.value)
names(moran.export) <- c("observed", "expected", "sd", "p.value")
write.csv(moran.export, file.path(output.directory, paste(model.name,"_moransI.csv", sep="")))

###### Create dataframe with coordinates and residuals to export ######
species.training <- rep(species, nrow(training.data))
resid.export <- as.data.frame(cbind(species.training, training.data$Count, gbmModel.count$fit, residual.training, training.data$coords.x1, training.data$coords.x2,
                                    training.data$latitude, training.data$longitude))
names(resid.export) <- c("species", "count", "fit", "residual", "coords.x1", "coords.x2", "Latitude", "Longitude")

###### Prep to calculate the spatial autocovariate  ######
write.csv(resid.export, file.path(output.directory, paste(model.name,"_resids.csv", sep="")))



###### Write the results ######
# Now add other important fields to be combined. 
head(yDat)
nrow(yDat)
site.counts <- rowSums(yDat)
head(site.counts)
head(table(site.counts))
occupied.site.counts <- site.counts[site.counts > 0]


number.suitable.sites <- nrow(observation.covariates.data)
suitable.sites.occupied <- nrow(occupied.site.counts) 
suitable.sites.unoccupied <- number.suitable.sites-suitable.sites.occupied
mean.suitable.abundance <- mean(site.counts)
mean.occupied.suitable.abundance <- mean(occupied.site.counts)
total.suitable.counted <- sum(occupied.site.counts)

## Give na to any missing variables
#c.hat.grem <- NA

results <- cbind(species, PercentTuncated, total.sites.occupied.brt, total.prevalence, AUC, current.suitable.habitat.area.hectares, historic.suitable.habitat.area.hectares, suitable.area.change, suitable.area.change.percent,
                 lr, bf.zi, tc.zi, radius.zi, trees.number, training.sites.occupied.brt, test.sites.occupied.brt,
                 TruncationD,
                 population.change.distsamp, population.change.lower.distsamp, population.change.upper.distsamp, population.change.min.distsamp, population.change.percent.distsamp,
                 mean.density.current.distsamp, population.current.distsamp, mean.density.current.lower.distsamp, population.current.lower.distsamp, mean.density.current.upper.distsamp, population.current.upper.distsamp,
                 mean.density.historic.distsamp, population.historic.distsamp, mean.density.historic.lower.distsamp, population.historic.lower.distsamp, mean.density.historic.upper.distsamp, population.historic.upper.distsamp,
                 p.distsamp, c.hat.dist,
                 #population.change.gdistsamp, population.change.lower.gdistsamp, population.change.upper.gdistsamp, population.change.min.gdistsamp, population.change.percent.gdistsamp,
                 #mean.density.current.gdistsamp, population.current.gdistsamp, mean.density.current.lower.gdistsamp, population.current.lower.gdistsamp, mean.density.current.upper.gdistsamp, population.current.upper.gdistsamp,
                 #mean.density.historic.gdistsamp, population.historic.gdistsamp, mean.density.historic.lower.gdistsamp, population.historic.lower.gdistsamp, mean.density.historic.upper.gdistsamp, population.historic.upper.gdistsamp,
                 #p.gdistsamp, c.hat.gdist,
                 population.change.rem, population.change.lower.rem, population.change.upper.rem, population.change.min.rem, population.change.percent.rem,
                 mean.density.current.rem, population.current.rem, mean.density.current.lower.rem, population.current.lower.rem, mean.density.current.upper.rem, population.current.upper.rem,
                 mean.density.historic.rem, population.historic.rem, mean.density.historic.lower.rem, population.historic.lower.rem, mean.density.historic.upper.rem, population.historic.upper.rem,
                 p.rem, c.hat.rem,
                 #population.change.grem, population.change.lower.grem, population.change.upper.grem, population.change.min.grem, population.change.percent.grem,
                 #mean.density.current.grem, population.current.grem, mean.density.current.lower.grem, population.current.lower.grem, mean.density.current.upper.grem, population.current.upper.grem,
                 #mean.density.historic.grem, population.historic.grem, mean.density.historic.lower.grem, population.historic.lower.grem, mean.density.historic.upper.grem, population.historic.upper.grem,
                 #p.grem, c.hat.grem,
                 population.change.tod, population.change.lower.tod, population.change.upper.tod, population.change.min.tod, population.change.percent.tod,
                 mean.density.current.tod, population.current.tod, mean.density.current.lower.tod, population.current.lower.tod, mean.density.current.upper.tod, population.current.upper.tod,
                 mean.density.historic.tod, population.historic.tod, mean.density.historic.lower.tod, population.historic.lower.tod, mean.density.historic.upper.tod, population.historic.upper.tod,
                 p.tod, c.hat.tod,
                 #population.change.gtod, population.change.lower.gtod, population.change.upper.gtod, population.change.min.gtod, population.change.percent.gtod,
                 #mean.density.current.gtod, population.current.gtod, mean.density.current.lower.gtod, population.current.lower.gtod, mean.density.current.upper.gtod, population.current.upper.gtod,
                 #mean.density.historic.gtod, population.historic.gtod, mean.density.historic.lower.gtod, population.historic.lower.gtod, mean.density.historic.upper.gtod, population.historic.upper.gtod,
                 #p.gtod, c.hat.gtod,
                 number.suitable.sites, suitable.sites.occupied, suitable.sites.unoccupied, mean.suitable.abundance, mean.occupied.suitable.abundance, total.suitable.counted)




write.csv(results, file.path(output.directory, paste(model.name, "_RESULTS.csv", sep="")))


