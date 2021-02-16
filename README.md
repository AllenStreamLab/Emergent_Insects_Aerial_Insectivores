# Emergent_Insects_Aerial_Insectivores

Data description

Scripts and files were used to predict the flux of aquatic insects from streams and rivers of the contiguous United States. This work is currently undergoing peer-review. Predictive models were fit using data from the National Rivers and Streams Assessment and a global synthesis of aquatic community secondary production (Patrick et al 2019). Below are the contents of this directory. Raw datafile provide example input files downloaded from their respected sources. The scripts heading are source code (i.e. “.r” extensions) and descriptions of their output files (bullet points). Comments are provided within the scripts for further details. Please direct any questions to darinkopp@gmial.com. 


Raw Data files 
•	NRSA_Stream_Width_Data.csv – raw NRSA wetted width datafile 
•	Annual_Insect_Secondary_Prodcution.csv – annual insect production values derived from Patrick et al (2019)
•	Avian_Insectivores.csv – species list of common aerial avian insectivores
•	insectivore_polygons_HUC4.shp – Randomly selected pixels stratified by 4-digit hydrologic units. 
Scripts
Wetted_Width_Fit.r – Bayesian mixed-effects model for stream width
•	Wet_Width_BRMS_Data – r data object containing training and testing data for wetted width model 
•	Wet_Width_BRMS_Scale_Values – r data object containing original center and scale values 
•	Wet_Width_BRMS_Model_Comparision – r data object containing loo model comparison results 
•	Wet_Width_BRMS_Model_Full – r data object containing Bayesian mixed effects modeling object

AISP_Randomforest_fit.r – random forest model for annual insect secondary production 
•	AISP_RF_Model – r data list containing training and testing data and fitted random forest models 
Emergence_at_NHDPlus_Flowlines.R – extrapolation of insect production and wetted width model to NHDPlusV2
•	National_Estimates – subdirectory containing estimates for stream area, emergent insect production and total emergence under current and future climate for approx. 2.3 million streams and rivers. 
ebird_Field_Metabolic_Rates.R – Calculate field metabolic rates form eBird status and Trends dataset (Fink et al 2020)
•	insectivore_Weekly_abund_site.csv – Weekly relative abundance for each insectivore species at all sites. 
•	insectivore_Weekly_FMR_site.csv – Sum of species weekly FMR and annual FMR for each study site. 
stream_bird_intersection.r - Overlay insectivore polygons with NHDPlus flowlines 
•	streams.shp – shapefile of stream segments that intersect insectivorous bird pixels
Emergence_at_Bird_Locations.r –estimate total emergence under current and future climate scenarios at insectivore bird polygon. 
•	emergence.csv- emergence predictions for insectivorous bird sites
•	emergence_85.csv – emergence predictions for insectivorous bird sites under RCP 8.5
•	emergence_26.csv – emergence predictions for insectivorous bird sites under RCP 2.6

