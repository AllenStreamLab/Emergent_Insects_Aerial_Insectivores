
rm(list = ls())
gc()

###############
#Extrapolate wetted width and aisp models to entire NHDPlusV2
###############

library(sf)
library(foreign)
library(tidyverse)
library(ggplot2)
library(randomForest)
library(quantregForest)

#rescale estimate function
rescale2 <- function(x,s,c){(x*s) + c} 
log_bias_corr <- function(x) {1 + ((sd(x)/sqrt(length(x))^2)/2)}

#Set wd to worldclime data to avoid multiple downloads...
setwd("C:/Users/AllenLabWorkstation/Dropbox/Dissertation/Chapter_3_Distance_Deposition/QuantityEstimates/Manuscript/Data/NRSA_Data_Raw")
#Need full NHDPlus Dataset
nhd.dir <- "C:/Users/AllenLabWorkstation/Dropbox/Dissertation/Network_Modeling/Actual/NHDPlus"
#create subdirecory to store results
write.dir <- "C:/Users/AllenLabWorkstation/Dropbox/Dissertation/Chapter_3_Distance_Deposition/QuantityEstimates/Manuscript/Predictive_Models/National_Estimates"

#data
#need to specify climate scenario
######
#identify vpu dirs
dirs <- strsplit(list.dirs(nhd.dir), "/")
ind <- unlist(lapply(dirs, function(x) length(x)==10))
dirs <- dirs[ind]
dirs <- unlist(lapply(dirs, function(x) x[10]))
vpus <- dirs[-c(6, 15, 24)]

#load wetted width model 
wet_wid_dir <- "C:/Users/AllenLabWorkstation/Dropbox/Dissertation/Chapter_3_Distance_Deposition/QuantityEstimates/Manuscript/Predictive_Models"
load(paste0(wet_wid_dir,"/Wet_Width_BRMS_Model_Full")) #Wet_Width_BRMS_Model_Full
load(paste0(wet_wid_dir,"/Wet_Width_BRMS_Scale_Values"))

#random forest production model 
rf_dir <- "C:/Users/AllenLabWorkstation/Dropbox/Dissertation/Chapter_3_Distance_Deposition/QuantityEstimates/Manuscript/Predictive_Models"
load(paste0(rf_dir, "/AISP_RF_Model")) #rf_model_list[["randomforest"]]

#the resolution of world clime bioclim at 2.5 is about 4.5km
raster_worldclim <- raster::getData("worldclim", var = "bio", res = 2.5)
#raster_worldclim <- raster::getData('CMIP5', var = "bio", res = 2.5, rcp = 26, model = 'HE', year = 70)
#raster_worldclim <- raster::getData('CMIP5', var = "bio", res = 2.5, rcp = 85, model = 'HE', year = 70)
#################

#iterate through each vpu
for (vpu in vpus){
  print(vpu)
  
  #data for vpu
  #####
  #vpu <- vpus[1]
  
  
  directory <- grep(paste(vpu, "/NHDPlusAttributes", sep = ""),
                    list.dirs(nhd.dir, full.names = T),
                    value = T)
  hydro_path <- grep(paste(vpu, "/NHDSnapshot/Hydrography", sep = ""),
                     list.dirs(nhd.dir, full.names = T),
                     value = T)
  
  
  #remove waterbody flowlines
  fl <- read_sf(paste0(hydro_path, "/NHDFlowline.shp"))
  names(fl)[-length(names(fl))] <- toupper(names(fl))[-length(names(fl))]
  fl <- st_zm(fl, what = "ZM", drop = T)
  fl <- fl[fl$FTYPE %in% c("StreamRiver", "ArtificialPath"),]
  wb <- read_sf(paste0(hydro_path, "/NHDWaterbody.shp"))
  names(wb)[-length(names(wb))] <- toupper(names(wb))[-length(names(wb))]
  wb <- st_zm(wb, what = "ZM", drop = T)
  wb <- st_set_geometry(wb, NULL)
  fl <- fl[!fl$WBAREACOMI %in% wb$COMID, ]
  rm(wb) #drop waterbody file 

  
  #remove divergent flow paths
  #remove flowlines with 0 upstream catchments
  vaa <- grep("PlusFlowlineVAA.dbf",
              list.files(directory[1], full.names = T),
              value = T)
  vaa <- read.dbf(vaa)
  names(vaa) <- toupper(names(vaa))
  vaa <- vaa[,c("COMID", "STREAMORDE", "STREAMLEVE", "STREAMCALC", "TOTDASQKM")]
  vaa <- vaa[vaa$STREAMORDE == vaa$STREAMCALC, ]
  vaa <- vaa[vaa$TOTDASQKM > 0, ]
  fl <- merge(fl[,c("COMID", "FTYPE")],vaa[,c("COMID", "TOTDASQKM")], by = "COMID")
  rm(vaa) #drop vaa
    
  
  #transform flowlines into projected CRS
  #extract bioclime variables 
  fl <- st_transform(fl, crs = 5070)
  flowline_pts <- st_centroid(fl)
  pts <- st_coordinates(st_transform(flowline_pts, raster::crs(raster_worldclim)))
  points_bioclim <- data.frame(st_set_geometry(flowline_pts, NULL), 
                               raster::extract(raster_worldclim, pts[,c("X","Y")]))
  points_bioclim <- points_bioclim[complete.cases(points_bioclim), ]
  #preserve names for future productions
  names(points_bioclim) <- c("COMID","FTYPE","TOTDASQKM",
                           "bio1","bio2","bio3","bio4","bio5","bio6",     
                           "bio7", "bio8","bio9","bio10","bio11","bio12",
                           "bio13","bio14","bio15","bio16","bio17","bio18","bio19")
  rm(flowline_pts, pts) #drop extras
  ###############
  
  #stream width and area predictions
  #####
  #scale values
  PPT <- scale(points_bioclim[,"bio12"], center = Wetted_Width_Scales$PPT$`scaled:center`, scale = Wetted_Width_Scales$PPT$`scaled:scale`)
  TEMP <- scale(points_bioclim[,"bio1"], center = Wetted_Width_Scales$TEMP$`scaled:center`, scale = Wetted_Width_Scales$TEMP$`scaled:scale`)
  COEFppt<-scale(points_bioclim[,"bio15"], center = Wetted_Width_Scales$COEFppt$`scaled:center`, scale = Wetted_Width_Scales$COEFppt$`scaled:scale`)
  WSAREA <- scale(log(points_bioclim$TOTDASQKM), center = Wetted_Width_Scales$WSAREA$`scaled:center`, scale = Wetted_Width_Scales$WSAREA$`scaled:scale`)
  #need extra letter for some vpus
  VPU <- ifelse(nchar(vpu) > 9,
         substring(vpu, nchar(vpu)-2, nchar(vpu)), 
         substring(vpu, nchar(vpu)-1, nchar(vpu)))
  COMID <- points_bioclim$COMID
  d <- data.frame(COMID, VPU, WSAREA, PPT, TEMP, COEFppt)
  rm(COMID, VPU, WSAREA, PPT, TEMP, COEFppt)
  
  #need loop because predict() generates large vector
  Fr <- seq(1, nrow(d), 5000)
  To <- c(seq(0, nrow(d), 5000)[-1], nrow(d))
  iter<-data.frame(Fr,To)
  
  print(paste("processing", nrow(d), "in", nrow(iter), "iterations")) 
  Stream_Width_Estimate <- data.frame()
  for (i in 1:nrow(iter)){
    #i<-1#d$VPU
    
    lnstrwid <- predict(Wet_Width_BRMS_Model_Full, newdata = d[iter$Fr[i]:iter$To[i], ])
    lnstrwid <- rescale2(lnstrwid, 
                         s = Wetted_Width_Scales$WETWID$`scaled:scale`,
                         c = Wetted_Width_Scales$WETWID$`scaled:center`)
    
    temp <- data.frame(COMID = d[iter$Fr[i]:iter$To[i], "COMID"], lnstrwid)
    Stream_Width_Estimate <- rbind(Stream_Width_Estimate,  temp)
  }
  rm(lnstrwid)
  
  #backtransform log values
  Stream_Width_Estimate$Estimate <- exp(1)^Stream_Width_Estimate$Estimate
  Stream_Width_Estimate$Q2.5  <- exp(1)^Stream_Width_Estimate$Q2.5
  Stream_Width_Estimate$Q97.5 <- exp(1)^Stream_Width_Estimate$Q97.5
  
  fl <- merge(fl, Stream_Width_Estimate[,-3], by = "COMID")
  rm(Stream_Width_Estimate, d, iter)

  
  #calculate stream area
  buff <- fl$Estimate
  stream_widths <- st_buffer(fl, buff/2)
  fl$area_Estimate <- st_area(stream_widths)
  rm(stream_widths)
  
  buff <- fl$Q2.5
  stream_widths <- st_buffer(fl, buff/2)
  fl$area_Q2.5 <- st_area(stream_widths)
  rm(stream_widths)
  
  buff <- fl$Q97.5
  stream_widths <- st_buffer(fl, buff/2)
  fl$area_Q97.5 <- st_area(stream_widths)
  rm(stream_widths)
  
  fl <- st_set_geometry(fl, NULL)
  ###############
  
  #estimate emergence production
  #####
  acsp_est <- data.frame(COMID = points_bioclim$COMID, 
                         rf_pred = predict(rf_model_list$randomforest, points_bioclim),
                         predict(rf_model_list$quantile_rf, newdata =  points_bioclim, 
                                 what = c(0.025, 0.975)))
  names(acsp_est) <- c("COMID", "AISP_mgafdm_est", "AISP_mgafdm_LowQ", "AISP_mgafdm_HiQ")
  acsp_est[,-1] <- exp(1)^acsp_est[,-1]
  fl <- merge(fl, acsp_est)
  rm(acsp_est)
  ##############

  #Finish multiplying
  #####
  #emergence grams
  fl$Emer_gAFDM <- (fl$AISP_mgafdm_est*0.1915461)/1000
  fl$Emer_gAFDM_LowQ <- (fl$AISP_mgafdm_LowQ*0.1915461)/1000
  fl$Emer_gAFDM_HiQ <- (fl$AISP_mgafdm_HiQ*0.1915461)/1000
 
  #total emergence for stream area
  fl$Total_Emer_gAFDM <- fl$area_Estimate*fl$Emer_gAFDM
  fl$Total_Emer_gAFDM_LowQ <- fl$area_Q2.5*fl$Emer_gAFDM_LowQ
  fl$Total_Emer_gAFDM_HiQ <- fl$area_Q97.5*fl$Emer_gAFDM_HiQ
  ###############

  write.csv(fl, paste0(write.dir, "/Emergence", vpu, ".csv"))
}


