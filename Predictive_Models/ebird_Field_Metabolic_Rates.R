# Extract weekly (1-52) summaries from ebird abundance 

library(raster)
library(ebirdst)
library(sf)

setwd("C:/Users/AllenLabWorkstation/Dropbox/Insectivore_Birds")
model_directory <- "C:/Users/AllenLabWorkstation/Dropbox/Dissertation/Chapter_3_Distance_Deposition/QuantityEstimates/Manuscript/Predictive_Models"

# extract bird abundance for study locations.
# iterates through each species. takes time
#####
df_a_i_na_all = read.csv(paste0(model_directory, "/Avian_Insectivores.csv"))
sites <- read_sf(paste0(model_directory, "/insectivore_polygons_HUC4.shp"))
#write individual bird csv
#birs.dir <- "C:/Users/AllenLabWorkstation/Dropbox/Insectivore_Birds/Site_abundances"

out.df <- data.frame()
for (species_code in as.vector(df_a_i_na_all$species_code)){
    T1<-Sys.time()
    #species_code <- as.vector(df_a_i_na_all$species_code[1])
    print(species_code)
    sp_path <- ebirdst_download(species = species_code, force = TRUE)
    abundances <- load_raster("abundance", path = sp_path)
  
    # extract coordinates (i.e. verticies) for buffers
    # much more efficient to extract values by points
    # this seems to tbe the fastest of extract by points, and using index
    z <- raster::extract(abundances, sites$cell)
    
    temp <- data.frame(st_set_geometry(sites, NULL), round(z, 2))
    temp$total <- apply(temp[-c(1:4)], 1, function (x) sum(x, na.rm = T))
  
    # remove sites where species is absent
    temp <- temp[which(temp$total != 0), ]
    
    if(nrow(temp)>0){
      #format output         
      temp <- data.frame(species_code, temp)         
      names(temp) <- c("species_code","cell","HUC4","vpu","NRSA", 
                      paste0("WK_", 1:52), "total")

      out.df <- rbind(out.df, temp)
      #writing individual species just in case i need to stop the loop for something
      #write.csv(temp, paste0(birs.dir,"/", species_code,"_abund.csv"), row.names = F)
    }
    T2<-Sys.time()
    print(T2-T1)
}

#############
write.csv(out.df, paste0(model_directory,"/insectivore_Weekly_abund_site.csv"), row.names = F)


# calculate annual FMR.
# sum weekly FMR for each site. 
#####
#birs.dir <- "C:/Users/AllenLabWorkstation/Dropbox/Insectivore_Birds/Site_abundances"

# The FMR in kJ/day, and the body masses are in g. 
# fmr = 10.5*body_mass^0.681.
#fmr <- read.csv("all_obligate_aerial_insectivores_tyrannidae_ebirdst_weighted_mean_body_masses_fmrs.csv")
fmr <- read.csv(paste0(model_directory,"Avian_Insectivores.csv"))
bird_week <- read.csv(paste0(model_directory,"/insectivore_Weekly_abund_site.csv"))

#merge by species code 
bw<-merge(bird_week, fmr[,c("species_code","fmr")])
#multiply by FMR (kj/day) and multiply by 7 to get weekly total FRM for all birds
bw[-c(1:5,58:59)] <- round(bw[-c(1:5,58:59)]*bw$fmr*7,3)
#sum across all species within each cell(quatrat)
z <- split(bw, bw$cell)
z <- lapply(z, function(x) data.frame(unique(x[,c("cell", "HUC4", "vpu", "NRSA")]),
                          sp_rich = nrow(x),
                          t(apply(x[-c(1:5, 58:59)], 2, sum, na.rm = T))))
z <- do.call(rbind, z)

#create total - kJ needed by all birds living in a quatrat throughout the year 
z$total_kjyr <- apply(z[,-c(1:5)], 1, sum, na.rm = T)
#############
write.csv(z, paste0(model_directory, "/insectivore_Weekly_FMR_site.csv"), row.names = F)

