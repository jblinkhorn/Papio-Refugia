library(tidysdm)
library(sf)
library(pastclim)
library(terra)
library(tidyterra)
library(ggplot2)
library(DALEX)
library(RColorBrewer)

download_dataset("WorldClim_2.1_10m")
download_dataset("Krapp2021")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #This assumes you have saved the code and data in a single directory

####LOAD DATA####
papio <- readxl::read_excel("Supplementary Table 1.xlsx") #load papio presence locations
papio2 <- st_as_sf(papio[c(1:1401),], coords = c("longitude", "latitude"))
st_crs(papio2) <- 4326
land_mask <-  get_land_mask(time_bp = 0, dataset = "Krapp2021")
papio_split <- split(papio2, papio2$Species) #split data by species
papio_split[[7]] <- papio2 #add full dataset for all Papio

#CREATE AFRICA ARABIA MASK#
africa <- rnaturalearth::ne_countries(continent = 'africa') #africa shape
africa <- africa[-38,] #removes madagascar
arabia <- rnaturalearth::ne_countries(continent = 'asia') #asia shape
arabia <- arabia[c(5, 7, 8, 9, 10, 11, 12, 13, 44, 45), ] #select arabia only
africa_arabia <- rbind(africa, arabia) #creates map of africa and arabia
africa_arabia <- terra::vect(africa_arabia) #to spatvector
africa_arabia <- terra::aggregate(africa_arabia, dissolve=T) #disolved internal divisions
rm(africa, arabia) #tidy up
writeVector(africa_arabia, "africa_arabia.shp", overwrite=T)#export for Figures

#CROP TO MASK#
crs(africa_arabia) <- "lonlat"
# crop the extent
land_mask <- crop(land_mask, africa_arabia)
# and mask to the polygon
land_mask <- mask(land_mask, africa_arabia)
ggplot() +
  geom_spatraster(data = land_mask, aes(fill = land_mask_0)) +
  geom_sf(data = papio2)

#NAME CONVENTIONS#
papio_popn <- c("P. anubis", "P. cynocephalus", "P. hamadryas", "P. kindae", "P. papio", "P. ursinus", "Papio")


#### DOWNSCALE DATA####
#This section uses dev functions in pastclim to produce "Af_Ar_downscaled_ds4.nc" which is provided in the GitHub repository

## Temperature variables
# extract low res monthly temperature variables 
tavg_vars <- c(paste0("temperature_0",1:9),paste0("temperature_",10:12)) 
download_dataset(dataset = "Krapp2021", bio_variables = tavg_vars)

tavg_series <- region_series(bio_variables = tavg_vars,
                             time_bp =  seq(0,-131000, -1000),
                             dataset = "Krapp2021",
                             ext = ext(land_mask))

# extract single low res stack
tavg_model_lres_rast <- tavg_series$temperature_01 

# high resolution modern data
download_dataset(dataset = "WorldClim_2.1_10m", bio_variables = tavg_vars)
tavg_obs_hres_all<- region_series(bio_variables = tavg_vars,
                                  time_ce =  1985,
                                  dataset = "WorldClim_2.1_10m",
                                  ext = ext(land_mask))


# estimate the temperature range
tavg_obs_range <- range(unlist(lapply(tavg_obs_hres_all,minmax, compute=TRUE))) 

# crop all high res modern data
tavg_obs_hres_all <- terra::crop(tavg_obs_hres_all, land_mask) 

# extract single high res raster
tavg_obs_hres_rast <- tavg_obs_hres_all[[1]] 

# download relief model
download_etopo()

relief_rast <- load_etopo()

# ensure same resolution as high res raster
relief_rast <- terra::resample(relief_rast, tavg_obs_hres_rast)
# create high res landmasks
land_mask_high_res <- make_land_mask(relief_rast = relief_rast,
                                     time_bp = seq(0,-131000, -1000))

# download ice masks
ice_mask_low_res <- get_ice_mask(time_bp= seq(0,-131000, -1000),dataset="Krapp2021") 
# increase resolution to high res landmask
ice_mask_high_res <- downscale_ice_mask(ice_mask_low_res = ice_mask_low_res,land_mask_high_res = land_mask_high_res)

# mask icesheets and coastlines
land_mask_high_res <- mask(land_mask_high_res, ice_mask_high_res, 
                           inverse=TRUE)

# download internal seas 
internal_seas <- readRDS(system.file("extdata/internal_seas.RDS", package="pastclim"))
# mask internal seas
land_mask_high <- mask(land_mask_high_res, internal_seas,inverse = TRUE)
# check                       
plot(land_mask_high)

# calculate delta raster and downscale monthly temperature 
tavg_downscaled_list<-list() 
for (i in 1:12){
  delta_rast<-delta_compute(x=tavg_series[[i]], ref_time = 0,
                            obs = tavg_obs_hres_all[[i]])
  tavg_downscaled_list[[i]] <- delta_downscale (x = tavg_series[[i]], delta_rast = delta_rast,
                                                x_landmask_high = land_mask_high_res,
                                                range_limits=tavg_obs_range)
}
tavg_downscaled <- terra::sds(tavg_downscaled_list)

### Precipitation variables 
# extract monthly precipitation variables
prec_vars <- c(paste0("precipitation_0",1:9),paste0("precipitation_",10:12)) 
download_dataset(dataset = "Krapp2021", bio_variables = prec_vars)

prec_series <- region_series(bio_variables = prec_vars, 
                             time_bp =  seq(0,-131000, -1000),
                             dataset = "Krapp2021",
                             ext = ext(land_mask))
# high resolution modern data
download_dataset(dataset = "WorldClim_2.1_10m", bio_variables = prec_vars) #> [1] TRUE
prec_obs_hres_all<- region_series(bio_variables = prec_vars,
                                  time_ce =  1985,
                                  dataset = "WorldClim_2.1_10m",
                                  ext = ext(land_mask))

# estimate the precipitation range
prec_obs_range <- range(unlist(lapply(prec_obs_hres_all,minmax, compute=TRUE))) 
prec_obs_range

# downscale precipitation 
prec_downscaled_list<-list() 

for (i in 1:12){
  delta_rast<-delta_compute(x=prec_series[[i]], ref_time = 0,
                            obs = prec_obs_hres_all[[i]])
  prec_downscaled_list[[i]] <- delta_downscale (x = prec_series[[i]], delta_rast = delta_rast,
                                                x_landmask_high = land_mask_high_res,
                                                range_limits = prec_obs_range)
}
prec_downscaled <- terra::sds(prec_downscaled_list)

### compute bioclim variables 
bioclim_downscaled<-bioclim_vars(tavg = tavg_downscaled, prec = prec_downscaled) # this takes ages! 

plot(bioclim_downscaled$bio01) #check


downscaled <- check
relief_rast_lyrs <- c(rep(relief_rast, 132))
time(relief_rast_lyrs) <- seq(0,-131000, -1000)
relief_rast_lyrs <- crop(relief_rast_lyrs, downscaled)
varnames(relief_rast_lyrs) <- 'altitude'
rugosity_rast_lyrs <- focal(relief_rast_lyrs, w=3, fun="sd")
varnames(rugosity_rast_lyrs) <- 'rugosity'

#names(relief_rast_lyrs) <- paste0("altitude_", 1:132)

ds4 <-  sds(downscaled[[1]],
  downscaled[[2]],
  downscaled[[3]],
  downscaled[[4]],
  downscaled[[5]],
  downscaled[[6]],
  downscaled[[7]],
  downscaled[[8]],
  downscaled[[9]],
  downscaled[[10]],
  downscaled[[11]],
  downscaled[[12]],
  downscaled[[13]],
  downscaled[[14]],
  downscaled[[15]],
  downscaled[[16]],
  downscaled[[17]],
  relief_rast_lyrs,
  rugosity_rast_lyrs)

terra:: writeCDF(ds4, "Af_Ar_downscaled_ds4.nc", overwrite=T)
## save 
terra:: writeCDF(Bio_Alt_Rug, "Af_Ar_downscaled2.nc", overwrite=T)
terra:: writeCDF(Bio_Alt_Rug, paste0("D:/Dropbox/Blinkhorn et al Papio", "/Af_Ar_downscaled.nc"), overwrite = TRUE)
## add in altitude and rugosity

land_mask_hi <- crop(relief_rast, africa_arabia)
land_mask_hi <- mask(land_mask_hi, africa_arabia)
plot(land_mask_hi)
res(land_mask_hi)



custom_data <- region_series(bio_variables =climate_vars,
                             dataset = "custom",
                             path_to_nc = "Af_Ar_downscaled_ds4.nc")


####30 Minute SDM####
biglist <- list()
set.seed(123)
for(j in 1:7){# for each baboon species
all_papio <- list()
papio2 <- papio_split[[j]] #select species data
papio3 <- thin_by_cell(papio2, raster = land_mask)#subset so that each cell has a maximum of a single presence
all_papio$thinnedcount <- nrow(papio3)
e <- extend(ext(papio2), 6) #set extent to 6 degrees beyond maximum extent of species dataset in each direction
land_mask2 <- crop(land_mask, e) #mask to extent
papio4 <- sample_pseudoabs(papio3,
                                 n = 6*nrow(papio3),
                                 raster = land_mask2,
                                 method = c("dist_min", 50000)) #generate background data

climate_vars <- get_vars_for_dataset("Krapp2021")
climate_vars <- climate_vars[!climate_vars=="biome"] #exclude biome
climate_vars <- climate_vars[!climate_vars=="lai"] #exclude lai
climate_vars <- climate_vars[!climate_vars=="npp"] #exclude npp
climate_vars <- climate_vars[!climate_vars=="altitude"] # exclude altitude

climate_present <- pastclim::region_slice(
  time_bp = 0,
  bio_variables = climate_vars,
  data = "Krapp2021",
  crop = as.polygons(land_mask2))

papio5 <- papio4 %>%bind_cols(terra::extract(climate_present, papio4, ID = FALSE))
papio5 %>% plot_pres_vs_bg(class)

papio5 %>% dist_pres_vs_bg(class)

vars_to_keep <- papio5 %>% dist_pres_vs_bg(class)
vars_to_keep <- names(vars_to_keep[vars_to_keep > 0.2])
papio6 <- papio5 %>% select(all_of(c(vars_to_keep, "class")))
vars_to_keep

climate_present <- climate_present[[vars_to_keep]]
vars_uncor <- filter_collinear(climate_present, cutoff = 0.7)
vars_uncor

papio7 <- papio6 %>% select(all_of(c(vars_uncor, "class")))
climate_present <- climate_present[[vars_uncor]]

lacerta_rec <- recipe(papio7, formula = class ~ .)
lacerta_rec
lacerta_thin <- papio7
papio7 %>% check_sdm_presence(class)

lacerta_models <-
  # create the workflow_set
  workflow_set(
    preproc = list(default = lacerta_rec),
    models = list(
      # the standard glm specs
      glm = sdm_spec_glm(),
      # rf specs with tuning
      rf = sdm_spec_rf(),
      # boosted tree model (gbm) specs with tuning
      gbm = sdm_spec_boost_tree(),
      # maxent specs with tuning
      maxent = sdm_spec_maxent()
    ),
    # make all combinations of preproc and models,
    cross = TRUE
  ) %>%
  # tweak controls to store information needed later to create the ensemble
  option_add(control = control_ensemble_grid())

lacerta_cv <- spatial_block_cv(lacerta_thin, v = 4, n=5)
check_splits_balance(lacerta_cv,class) 
autoplot(lacerta_cv)

lacerta_models <-
  lacerta_models %>%
  workflow_map("tune_grid",
               resamples = lacerta_cv, grid = 10,
               metrics = sdm_metric_set(), verbose = TRUE
  )

autoplot(lacerta_models)

lacerta_ensemble <- simple_ensemble() %>%
  add_member(lacerta_models, metric = "boyce_cont")

autoplot(lacerta_ensemble)
prediction_present <- predict_raster(lacerta_ensemble, climate_present)
lacerta_ensemble <- calib_class_thresh(lacerta_ensemble,
                                       class_thresh = "tss_max")

climate_present_prediction <- pastclim::region_slice(
  time_bp = 0,
  bio_variables = vars_uncor,
  data = "Krapp2021",
  crop = africa_arabia)

prediction_present_binary <- predict_raster(lacerta_ensemble,
                                            climate_present_prediction,
                                            type = "class",
                                            class_thresh = c("tss_max"))

all_papio$papio_models <- lacerta_models
all_papio$papio_ensemble <- lacerta_ensemble
all_papio$papio_cv <- lacerta_cv
all_papio$papio_thin <- lacerta_thin
all_papio$vars_uncor <- vars_uncor
biglist[[j]] <- all_papio
print(j)}

####30 Minute HINDCASTING####

for(j in 1:7){
time_steps <- rev(get_time_bp_steps("Krapp2021"))
prediction_past_list <- list()
for(i in 1:131){

climate_past <- pastclim::region_slice(
  time_bp = time_steps[[i]],
  bio_variables = biglist[[j]]$vars_uncor,
  data = "Krapp2021",
  crop = africa_arabia
)

prediction_past_list[[i]] <- predict_raster(biglist[[j]]$papio_ensemble, climate_past,type = "class",
                                  class_thresh = c("tss_max"))
}
biglist[[j]]$prediction_past_list <- prediction_past_list
print(j)}
####30 Minute EXPORT####

#this saves from present to past
for(k in 1:7){writeCDF(sds(biglist[[k]]$prediction_past_list), paste0(papio_popn[[k]], "_30_minutes.nc", sep=""))}

#EXPORT MODEL METRICS

asda_list <- list()
for(k in 1:7){
  asda <- biglist[[k]]$papio_ensemble$metrics
  fasda <- as.data.frame(asda[[1]])
  for(i in 2:length(asda)){fasda <- rbind(fasda, as.data.frame(asda[[i]]))}
  asda_list[[k]] <- fasda}
for(i in 1:length(asda_list)){xlsx::write.xlsx(asda_list[[i]], "Table 1a model_metrics_30_minutes.xlsx", sheetName = papio_popn[[i]], append = T)}

#EXPORT MEAN VARIABLE IMPORTANCE
var_impotance_list <- list()
for(i in 1:7){
  vip_ensemble <- model_parts(explainer = explain_tidysdm(biglist[[i]]$papio_ensemble))
  variables <- split(vip_ensemble, vip_ensemble$variable)
  var_impotance_list[[i]] <- lapply(variables, function(x) mean(x$dropout_loss))}
ws2 <- data.frame(t(do.call(plyr::rbind.fill, lapply(var_impotance_list, as.data.frame))))
names(ws2) <- papio_popn
write.csv(ws2, "Table 1b mean_importance_30_minutes.csv")

c(biglist[[1]]$thinnedcount,
  biglist[[2]]$thinnedcount,
  biglist[[3]]$thinnedcount,
  biglist[[4]]$thinnedcount,
  biglist[[5]]$thinnedcount,
  biglist[[6]]$thinnedcount,
  biglist[[7]]$thinnedcount)

####10 Minute SDM####

land_mask_hi <- disagg(land_mask, fact=3)
biglist <- list()
set.seed(123)
for(j in 1:7){
  all_papio <- list()
  papio2 <- papio_split[[j]]
  papio3 <- thin_by_cell(papio2, raster = land_mask_hi)
  all_papio$thinnedcount <- nrow(papio3)
  limit <- c(6,6,6,5,4,6,6)
  e <- extend(ext(papio2), limit[[j]])
  land_mask2 <- crop(land_mask_hi, e)
  plot(land_mask2)
  points(papio2)
  papio4 <- sample_pseudoabs(papio3,
                             n = 6*nrow(papio3),
                             raster = land_mask2,
                             method = c("dist_min", 50000)
  )
  
  climate_vars <- get_vars_for_dataset("Krapp2021")
  climate_vars <- climate_vars[!climate_vars=="biome"]
  climate_vars <- climate_vars[!climate_vars=="lai"]
  climate_vars <- climate_vars[!climate_vars=="npp"]
  climate_vars <- climate_vars[!climate_vars=="altitude"]
  
  climate_present <- pastclim::region_slice(
    time_bp = 0,
    bio_variables = climate_vars,
    data = "custom",
    path_to_nc = "Af_Ar_downscaled_ds4.nc",
    crop = as.polygons(land_mask2)
  )
  
  papio5 <- papio4 %>%bind_cols(terra::extract(climate_present, papio4, ID = FALSE))
  papio5 <- na.omit(papio5)
  papio5 %>% plot_pres_vs_bg(class)
  papio5 %>% dist_pres_vs_bg(class)
  
  vars_to_keep <- papio5 %>% dist_pres_vs_bg(class)
  vars_to_keep <- names(vars_to_keep[vars_to_keep > 0.25])
  papio6 <- papio5 %>% select(all_of(c(vars_to_keep, "class")))
  vars_to_keep
  
  climate_present <- climate_present[[vars_to_keep]]
  vars_uncor <- filter_collinear(climate_present, cutoff = 0.7)
  vars_uncor
  
  papio7 <- papio6 %>% select(all_of(c(vars_uncor, "class")))
  climate_present <- climate_present[[vars_uncor]]
  
  
  lacerta_rec <- recipe(papio7, formula = class ~ .)
  lacerta_rec
  lacerta_thin <- papio7
  papio7 %>% check_sdm_presence(class)
  
  lacerta_models <-
    # create the workflow_set
    workflow_set(
      preproc = list(default = lacerta_rec),
      models = list(
        # the standard glm specs
        glm = sdm_spec_glm(),
        # rf specs with tuning
        rf = sdm_spec_rf(),
        # boosted tree model (gbm) specs with tuning
        gbm = sdm_spec_boost_tree(),
        # maxent specs with tuning
        maxent = sdm_spec_maxent()
      ),
      # make all combinations of preproc and models,
      cross = TRUE
    ) %>%
    # tweak controls to store information needed later to create the ensemble
    option_add(control = control_ensemble_grid())
  
  lacerta_cv <- spatial_block_cv(lacerta_thin, v = 4, n=5)
  check_splits_balance(lacerta_cv,class) 
  autoplot(lacerta_cv)
  
  lacerta_models <-
    lacerta_models %>%
    workflow_map("tune_grid",
                 resamples = lacerta_cv, grid = 10,
                 metrics = sdm_metric_set(), verbose = TRUE
    )
  
  autoplot(lacerta_models)
  
  lacerta_ensemble <- simple_ensemble() %>%
    add_member(lacerta_models, metric = "boyce_cont")
  lacerta_ensemble
  
  autoplot(lacerta_ensemble)
  
  prediction_present <- predict_raster(lacerta_ensemble, climate_present)

  lacerta_ensemble <- calib_class_thresh(lacerta_ensemble,
                                         class_thresh = "tss_max"
  )
  
  climate_present_prediction <- pastclim::region_slice(
    time_bp = 0,
    bio_variables = vars_uncor,
    data = "custom",
    path_to_nc = "Af_Ar_downscaled_ds4.nc",
    crop = as.polygons(land_mask2)
  )
  
  prediction_present_binary <- predict_raster(lacerta_ensemble,
                                              climate_present_prediction,
                                              type = "class",
                                              class_thresh = c("tss_max")
  )

  all_papio$papio_models <- lacerta_models
  all_papio$papio_ensemble <- lacerta_ensemble
  all_papio$papio_cv <- lacerta_cv
  all_papio$papio_thin <- lacerta_thin
  all_papio$vars_uncor <- vars_uncor
  biglist[[j]] <- all_papio
  print(j)}

####10 minute HINDCASTING####
for(j in 1:7){
  time_steps <- rev(get_time_bp_steps("Krapp2021"))
  prediction_past_list <- list()
  for(i in 1:131){
    
    climate_past <- pastclim::region_slice(
      time_bp = time_steps[[i]],
      bio_variables = biglist[[j]]$vars_uncor,
      data = "custom",
      path_to_nc = "Af_Ar_downscaled_ds4.nc",
      crop = africa_arabia
    )
    
    prediction_past_list[[i]] <- predict_raster(biglist[[j]]$papio_ensemble, climate_past,type = "class",
                                                class_thresh = c("tss_max"))
  }
  biglist[[j]]$prediction_past_list <- prediction_past_list
  print(j)}

plot(biglist[[j]]$prediction_past_list[[75]])

####10 minute EXPORT####

for(k in 1:7){writeCDF(sds(biglist[[k]]$prediction_past_list), paste0(papio_popn[[k]], "_10_minute.nc", sep=""), overwrite=T)}

asda_list <- list()
for(k in 1:7){
  asda <- biglist[[k]]$papio_ensemble$metrics
  fasda <- as.data.frame(asda[[1]])
  for(i in 2:length(asda)){fasda <- rbind(fasda, as.data.frame(asda[[i]]))}
  asda_list[[k]] <- fasda}
for(i in 1:length(asda_list)){xlsx::write.xlsx(asda_list[[i]], "Table 1a dd model_metrics_10_minute.xlsx", sheetName = papio_popn[[i]], append = T)}

var_impotance_list <- list()
for(i in 1:7){
  vip_ensemble <- model_parts(explainer = explain_tidysdm(biglist[[i]]$papio_ensemble))
  variables <- split(vip_ensemble, vip_ensemble$variable)
  var_impotance_list[[i]] <- lapply(variables, function(x) mean(x$dropout_loss))}
ws2 <- data.frame(t(do.call(plyr::rbind.fill, lapply(var_impotance_list, as.data.frame))))
names(ws2) <- papio_popn
write.csv(ws2, "Table 1b dd mean_importance_10_minutes.csv")

c(biglist[[1]]$thinnedcount,
  biglist[[2]]$thinnedcount,
  biglist[[3]]$thinnedcount,
  biglist[[4]]$thinnedcount,
  biglist[[5]]$thinnedcount,
  biglist[[6]]$thinnedcount,
  biglist[[7]]$thinnedcount)
