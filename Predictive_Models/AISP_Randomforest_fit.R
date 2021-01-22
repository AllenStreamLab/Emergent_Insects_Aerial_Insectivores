citation("raster")
citation("quantregForest")
#######################
# random forest model for Annual Insect Secondary Production (AISP)
#######################

library(randomForest)
library(quantregForest)
library(sf)
library(ade4)
#Set wd to worldclime data to avoid multiple downloads...
setwd("C:/Users/AllenLabWorkstation/Dropbox/Dissertation/Chapter_3_Distance_Deposition/QuantityEstimates/Manuscript/Data/NRSA_Data_Raw")

#place to save final model
model_directory <- "C:/Users/AllenLabWorkstation/Dropbox/Dissertation/Chapter_3_Distance_Deposition/QuantityEstimates/Manuscript/Predictive_Models"

#data
#####
acsp <- read.csv(paste0(model_directory,"/Annual_Insect_Secondary_Production.csv"))#NonInsect_P

acsp <- acsp[,c("SiteID", "X", "Y", "AISP")]
names(acsp)[4] <- "Production"
#################

#Check spatial autocorrelation in response
######
head(acsp)
stn.dist<-dist(cbind(acsp$X, acsp$Y))
prod.dist<-dist(acsp$Production)

mantel.rtest(stn.dist,prod.dist, nrepet = 9999)
#################

# extract bioclim variables
#####
# load worldclim data
raster_worldclim = raster::getData("worldclim", var="bio", res=2.5)
pts <- st_as_sf(acsp, coords = c("X", "Y"), crs = 4629)
pts <- st_coordinates(st_transform(pts, raster::crs(raster_worldclim)))

points_bioclim = data.frame(raster::extract(raster_worldclim, pts[,c("X","Y")]))
points_bioclim = data.frame(acsp, points_bioclim)
acsp <- points_bioclim

#Reduce variables only include bioclim variables
acsp_rf <- acsp[,c("Production", grep("bio", names(acsp), value = T))]
acsp_rf <- acsp_rf[complete.cases(acsp_rf),]

#highly correlated variables
apply(cor(acsp_rf),1, function(x) x[abs(x)>0.7])
#bio1 = c(bio5, bio6, bio10, bio11)
#bio2
#bio3 = c(bio4, bio7) 
#bio8
#bio9
#bio12 = c(bio13, bio14, bio16, bio17, bio18, bio19)
#bio15

#limited climate set
acsp_rf <- acsp_rf[c("Production", "bio1", "bio2", "bio3", "bio8", "bio9", "bio12", "bio15")]
################

#split test&trainaing  80/20
######
#data saved below
ind <- sample(2, nrow(acsp_rf), replace = T, prob = c(0.80, 0.20))
traindat <- acsp_rf[ind == 1, ]
testdat <- acsp_rf[ind == 2, ]  

rf <- randomForest(Production ~ ., data = traindat, ntree = 3000, importance = T)
################

#optimize mtry parameter - sensitivity analysis of mtry parameter
######
results_mtry_optimization <- matrix(data = NA , nrow = 0, ncol = 3)
for (i in c(seq(from = 10, to = 1000 , by = 10))){  # values of ntree
  #print(i)
  #c(1, "2*sqrt(p)","0.2*p", "p/3", "p") where p is number of variables. p=7 P/3 is rf default
  for (j in c(1, 2, 3, 4, 5, 6, 7)){   
    rf_ij <- randomForest(Production ~ ., data = traindat,
                          importance = TRUE, proximity = TRUE, 
                          ntree = i, mtry = j)
    results_mtry_optimization <- rbind(results_mtry_optimization, 
                                       c(i, j, tail(rf_ij$rsq, 1)))
  }
}

# Clean up the file format
results_mtry_optimization <- as.data.frame(results_mtry_optimization)
colnames(results_mtry_optimization) <- c("ntree", "mtry", "PVE")

mtry <- results_mtry_optimization[which.max(results_mtry_optimization$PVE), ]
mtry_defaut.max <- max(results_mtry_optimization$PVE[results_mtry_optimization$mtry == 2])

################

#check repeatability of random forest to set ntree parameater
######
results_ntree_optimization <- matrix(data = NA , nrow = 0, ncol = 2)
for (tree in seq(100, 3000, by = 100)){
  #p <- ntree
  rf_1 <- randomForest(Production ~ ., data = traindat, 
                       importance = TRUE, proximity = TRUE, 
                       ntree = tree, mtry = mtry$mtry)
  
  rf_2 <- randomForest(Production ~ ., data = traindat, 
                       importance = TRUE, proximity = TRUE, 
                       ntree = tree, mtry = mtry$mtry)
  
  r <- cor(importance(rf_1, scale = T, type = 1), 
           importance(rf_2, scale = T, type = 1))
  
  results_ntree_optimization <- rbind(results_ntree_optimization, c(tree, r))
  if(r >= 0.99){break()}
}

results_ntree_optimization <- as.data.frame(results_ntree_optimization)
colnames(results_ntree_optimization) <- c("ntree", "corr")

################

#reun model with optimized parameters
######
rf <- randomForest(Production ~ ., data = traindat, 
                   ntree = tree, mtry = mtry$mtry,
                   importance = T, proximity = T)

#quantile random forest to define 95% prediction interval
xprod <- traindat[,-1]
yprod <- traindat[, 1]

#names(rf_model_list$randomforest)
Qunatile_RF <- quantregForest(xprod, yprod, 
                              ntree = tree, mtry = 3, 
                              importance = T, keep.inbag = T)


rf_model_list <- list(model_form = Production ~ .,
                      traindat = traindat, 
                      testdat = testdat, 
                      optimize_mtry = results_mtry_optimization,
                      optimize_ntree = results_ntree_optimization, 
                      randomforest = rf,
                      quantile_rf = Qunatile_RF)
################

#place to save final model
#save(rf_model_list, file = paste0(model_directory, "/AISP_RF_Model"))

rf_model_list$randomforest
