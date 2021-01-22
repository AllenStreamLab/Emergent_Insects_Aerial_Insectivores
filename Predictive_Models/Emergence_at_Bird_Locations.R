###################################################
#Estimate stream width 
#calculate stream area within Pixel 
#Estimate emergence 
#calculate total emergence (flux) within pixel
###################################################

#predictions at ebird sites 
library(dplyr)
library(randomForest)
library(units)
library(foreign)
library(sf)
library(brms)
library(raster)
library(quantregForest)


setwd("C:/Users/AllenLabWorkstation/Dropbox/Insectivore_Birds")
model_directory <- "C:/Users/AllenLabWorkstation/Dropbox/Dissertation/Chapter_3_Distance_Deposition/QuantityEstimates/Manuscript/Predictive_Models"

#rescale estimate function
rescale2 <- function(x,s,c){(x*s) + c} 
log_bias_corr <- function(x) {1 + ((sd(x)/sqrt(length(x))^2)/2)}

#data
#####
bird_sites <- read_sf(paste0(model_directory,"/insectivore_polygons_HUC4.shp"))
site_cntr <- st_centroid(bird_sites)

#load streams dataset
streams <- read_sf(paste0(model_directory,"/streams.shp"))  
#remove coastlines we from analysis 
#unique(st_set_geometry(streams[streams$TOTDASQKM == -9999, "FTYPE"],NULL))
streams <- streams[streams$TOTDASQKM != -9999, ]
#remove flowlines without drainage area 
streams <- streams[streams$TOTDASQKM > 0, ]
#stream points to extract worldclime data
streams_pts <- st_centroid(streams)

#load wetted width model 
load(paste0(model_directory, "/Wet_Width_BRMS_Model_Full")) #Wet_Width_BRMS_Model_Full
load(paste0(model_directory, "/Wet_Width_BRMS_Scale_Values"))

#random forest production model 
load(paste0(model_directory, "/AISP_RF_Model")) #rf_model_list[["randomforest"]]


#the resolution of world clime bioclim at 2.5 is about 4.5km
#HE = from HADGEM2-ES model, 8.5 and 2.6 RCP for 2070
raster_worldclim <- raster::getData("worldclim", var = "bio", res = 2.5)
raster_worldclim_85 <- raster::getData('CMIP5', var = "bio", res = 2.5, rcp = 85, model = 'HE', year = 70)
raster_worldclim_26 <- raster::getData('CMIP5', var = "bio", res = 2.5, rcp = 26, model = 'HE', year = 70)

#format stream width prodiction data
pts <- st_coordinates(st_transform(streams_pts, raster::crs(raster_worldclim)))
points_bioclim <- data.frame(raster::extract(raster_worldclim, pts[,c("X","Y")]))
stream_df <- data.frame(st_set_geometry(streams_pts, NULL), 
                        points_bioclim[, c("bio1", "bio12", "bio15")])
stream_df <- stream_df[complete.cases(stream_df), ]


# extract bioclime variables at points
pts <- st_coordinates(st_transform(streams_pts, raster::crs(raster_worldclim)))
points_bioclim_26 <- data.frame(raster::extract(raster_worldclim_26, pts[,c("X","Y")]))
stream_df_26 <- data.frame(st_set_geometry(streams_pts, NULL), 
                           points_bioclim_26[,c("he26bi701", "he26bi7012", "he26bi7015")])
stream_df_26 <- stream_df_26[complete.cases(stream_df_26), ]
names(stream_df_26)[-c(1:5)] <- c("bio1", "bio12", "bio15")

#format 8.5 scenerio
pts <- st_coordinates(st_transform(streams_pts, raster::crs(raster_worldclim)))
points_bioclim_85 <- data.frame(raster::extract(raster_worldclim_85, pts[,c("X","Y")]))
stream_df_85 <- data.frame(st_set_geometry(streams_pts, NULL), 
                           points_bioclim_85[,c("he85bi701", "he85bi7012", "he85bi7015")])
stream_df_85 <- stream_df_85[complete.cases(stream_df_85), ]
names(stream_df_85)[-c(1:5)] <- c("bio1", "bio12", "bio15")


#format insect production variables 
bioclimvars <- c("bio1", "bio2", "bio3", "bio8", "bio9", "bio12", "bio15")
pts <- st_coordinates(st_transform(site_cntr, raster::crs(raster_worldclim)))
points_bioclim <- data.frame(raster::extract(raster_worldclim, pts[,c("X","Y")]))
insect_df <- data.frame(cell = site_cntr$cell, points_bioclim[,bioclimvars])
insect_df <- insect_df[complete.cases(insect_df), ]


bioclimvars_Future_26 <- c("he26bi701", "he26bi702", "he26bi703", "he26bi708", 
                           "he26bi709", "he26bi7012","he26bi7015")
pts <- st_coordinates(st_transform(site_cntr, raster::crs(raster_worldclim_26)))
points_bioclim_26 <- data.frame(raster::extract(raster_worldclim_26, pts[,c("X","Y")]))
insect_df_26 <- data.frame(cell = site_cntr$cell, points_bioclim_26[,bioclimvars_Future_26])
insect_df_26 <- insect_df_26[complete.cases(insect_df_26), ]
names(insect_df_26)[-1] <- bioclimvars #changes names


bioclimvars_Future_85 <- c("he85bi701", "he85bi702", "he85bi703", "he85bi708", 
                           "he85bi709", "he85bi7012","he85bi7015")
points_bioclim_85 <- data.frame(raster::extract(raster_worldclim_85, pts[,c("X","Y")]))
insect_df_85 <- data.frame(cell = site_cntr$cell, points_bioclim_85[,bioclimvars_Future_85])
insect_df_85 <- insect_df_85[complete.cases(insect_df_85), ]
names(insect_df_85)[-c(1)] <- bioclimvars


##########

#calculate stream area for current and future climates
#####
clim_pred <- list(stream_df, stream_df_26, stream_df_85)
names(clim_pred) <- c("bio", "he26", "he85")

for (q in 1:length(clim_pred)){
  #reset streams
  streams <- read_sf("streams.shp")  
  print(names(clim_pred)[q])
  #q<-3
  clim_dat <- clim_pred[[q]]
  
  #scale values
  #####
  #names(Wetted_Width_Scales)
  PPT <- scale(clim_dat[,"bio12"], center = Wetted_Width_Scales$PPT$`scaled:center`, scale = Wetted_Width_Scales$PPT$`scaled:scale`)
  TEMP <- scale(clim_dat[,"bio1"], center = Wetted_Width_Scales$TEMP$`scaled:center`, scale = Wetted_Width_Scales$TEMP$`scaled:scale`)
  COEFppt<-scale(clim_dat[,"bio15"], center = Wetted_Width_Scales$COEFppt$`scaled:center`, scale = Wetted_Width_Scales$COEFppt$`scaled:scale`)
  WSAREA <- scale(log(clim_dat$TOTDASQKM), center = Wetted_Width_Scales$WSAREA$`scaled:center`, scale = Wetted_Width_Scales$WSAREA$`scaled:scale`)
  VPU <- clim_dat$vpu
  COMID <- clim_dat$COMID
  cell <- clim_dat$cell
  d <- data.frame(cell, COMID, VPU, WSAREA, PPT, TEMP, COEFppt)
  ###############
  
  #predict stream width
  #####
  #estimate the wetted with of each stream 
  #need loop because predict() generates large vector
  Stream_Width_Estimate <- data.frame()
  for (j in unique(d$VPU)){
    #j<-"08"
    #print(j)
    lnstrwid <- predict(Wet_Width_BRMS_Model_Full, newdata = d[d$VPU == j, ])
    lnstrwid <- rescale2(lnstrwid, 
                         s = Wetted_Width_Scales$WETWID$`scaled:scale`,
                         c = Wetted_Width_Scales$WETWID$`scaled:center`)
    
    temp <- data.frame(COMID = d[d$VPU == j, "COMID"],
                       cell = d[d$VPU == j, "cell"], lnstrwid)
    Stream_Width_Estimate <- rbind(Stream_Width_Estimate,  temp)
  }
  
  #backtransform log values
  #account for bias in back transformation???
  #log_bias_corr(Stream_Width_Estimate$Estimate)
  Stream_Width_Estimate$Estimate <- exp(1)^Stream_Width_Estimate$Estimate
  Stream_Width_Estimate$Q2.5  <- exp(1)^Stream_Width_Estimate$Q2.5
  Stream_Width_Estimate$Q97.5 <- exp(1)^Stream_Width_Estimate$Q97.5
  
  #which.max(Stream_Width_Estimate[976, ])
  #clim_dat[clim_dat$cell == 13875270, ]
  
  #high and low estimates of the mean
  streams <- merge(streams, Stream_Width_Estimate[,-4], by = c("COMID", "cell"))
  #############
  
  #calculate stream area
  #####
  #iterate through estimates()
  count = T
  for (p in "Estimate"){#c("Q2.5", "Estimate", "Q97.5")){
    #p<-"Estimate.y"
    print(p)
    buff <- st_set_geometry(streams[,p], NULL)[,1]
    
    #create buffers reflecting 95%CrI  
    stream_widths <- st_buffer(streams, buff/2)
    
    #calcuate area of stream within bird survey locaiton
    stream_area <- st_intersection(stream_widths, bird_sites)
    
    #remove width buffers that fall into an adjacent cell
    stream_area <- stream_area[stream_area$cell == stream_area$cell.1, ]
    
    # calculate area of each overlappling stream to used in area weighted mean
    stream_area$streamarea <- units::set_units(st_area(stream_area),"km^2")
    
    #use centriod of each area polygon to extract bioclime vars
    #stream_centroid <- st_centroid(stream_area)
    
    #dissolve stream buffers reclaculate total stream area
    temp <- stream_area %>% 
      group_by(cell, vpu, HUC4) %>% 
      summarize(area = sum(streamarea))
    
    #recalculate area
    temp$area <- st_area(temp)
    ################
    #plot(st_geometry(temp[temp$cell==15792973,"area"]))
    ###################
    temp <- st_set_geometry(temp, NULL)
    
    if(count){
      names(temp)[!names(temp)%in%c("vpu", "HUC4", "cell")]<-
        paste0(names(temp)[!names(temp)%in%c("vpu", "HUC4", "cell")],"_", p)
      Stream_Area <- temp
      count = F
    } else {
      names(temp)[!names(temp)%in%c("vpu", "HUC4", "cell")]<-
        paste0(names(temp)[!names(temp)%in%c("vpu", "HUC4", "cell")],"_", p)
      Stream_Area <- merge(Stream_Area, temp, by = c("vpu", "HUC4", "cell"))
    }
  }
  ##############
  clim_pred[[q]] <- Stream_Area
}

################

#predict
#####
insectp_est <- data.frame(insect_df$cell, 
                          predict(rf_model_list$randomforest, newdata = insect_df),
                          predict(rf_model_list$quantile_rf, newdata = insect_df, 
                                  what = c(0.025, 0.975)))
names(insectp_est) <- c("cell","AISP_mgafdm_est", "AISP_mgafdm_LowQ", "AISP_mgafdm_HiQ")
insectp_est[,-1] <- exp(1)^insectp_est[,-1]

#future climate RCP 2.6
insectp_est_26 <- data.frame(insect_df_26$cell, 
                             predict(rf_model_list$randomforest, newdata = insect_df_26),
                             predict(rf_model_list$quantile_rf, newdata = insect_df_26, 
                                     what = c(0.025, 0.975)))
names(insectp_est_26) <- c("cell", "AISP_mgafdm_est", "AISP_mgafdm_LowQ", "AISP_mgafdm_HiQ")
insectp_est_26[,-1] <- exp(1)^insectp_est_26[,-1]

#future climate RCP 8.5
insectp_est_85 <- data.frame(insect_df_85$cell, 
                             predict(rf_model_list$randomforest, newdata = insect_df_85),
                             predict(rf_model_list$quantile_rf, newdata = insect_df_85, 
                                     what = c(0.025, 0.975)))
names(insectp_est_85) <- c("cell", "AISP_mgafdm_est", "AISP_mgafdm_LowQ", "AISP_mgafdm_HiQ")
insectp_est_85[,-1] <- exp(1)^insectp_est_85[,-1]

insectp_est <- merge(clim_pred[["bio"]], insectp_est, by = c("cell"))
insectp_est_26 <- merge(clim_pred[["he26"]], insectp_est_26, by = c("cell"))
insectp_est_85 <- merge(clim_pred[["he85"]], insectp_est_85, by = c("cell"))
################

#finish multiplying
#####

#emergence grams
insectp_est$Emer_gAFDM <- ((insectp_est$AISP_mgafdm_est*0.1915461)/1000)
insectp_est_26$Emer_gAFDM <- ((insectp_est_26$AISP_mgafdm_est*0.1915461)/1000)
insectp_est_85$Emer_gAFDM <- ((insectp_est_85$AISP_mgafdm_est*0.1915461)/1000)

#total emergence for stream area
insectp_est$Total_Emerg_gafdm <- insectp_est$area_Estimate*insectp_est$Emer_gAFDM
insectp_est_26$Total_Emerg_gafdm <- insectp_est_26$area_Estimate*insectp_est_26$Emer_gAFDM
insectp_est_85$Total_Emerg_gafdm <- insectp_est_85$area_Estimate*insectp_est_85$Emer_gAFDM

#energy
insectp_est$Total_Emerg_kJ <- (insectp_est$Total_Emerg_gafdm * 0.5) * 23.013
insectp_est_26$Total_Emerg_kJ <- (insectp_est_26$Total_Emerg_gafdm * 0.5) * 23.013
insectp_est_85$Total_Emerg_kJ <- (insectp_est_85$Total_Emerg_gafdm * 0.5) * 23.013

####################

#write.csv(insectp_est, paste0(model_directory,"/emergence.csv"), row.names = F)
#write.csv(insectp_est_85, paste0(model_directory,"/emergence_85.csv"), row.names = F)
#write.csv(insectp_est_26, paste0(model_directory,"/emergence_26.csv"), row.names = F)
