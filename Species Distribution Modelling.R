set.seed(12345)
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
papio <- readxl::read_excel("Supplementary Table 1.xslx") #load papio presence locations
papio2 <- st_as_sf(papio, coords = c("longitude", "latitude"))
st_crs(papio2) <- 4326
land_mask <-  get_land_mask(time_bp = 0, dataset = "Krapp2021")

####CREATE AFRICA ARABIA MASK####
africa <- rnaturalearth::ne_countries(continent = 'africa') #africa shape
africa <- africa[-38,] #removes madagascar
arabia <- rnaturalearth::ne_countries(continent = 'asia') #asia shape
arabia <- arabia[c(5, 7, 8, 9, 10, 11, 12, 13, 44, 45), ] #select arabia only
africa_arabia <- rbind(africa, arabia) #creates map of africa and arabia
africa_arabia <- terra::vect(africa_arabia) #to spatvector
africa_arabia <- terra::aggregate(africa_arabia, dissolve=T) #disolved internal divisions
rm(africa, arabia)
writeVector(africa_arabia, "africa_arabia.shp", overwrite=T)#export for Figures

####CROP TO MASK####
crs(africa_arabia) <- "lonlat"
# crop the extent
land_mask <- crop(land_mask, africa_arabia)
# and mask to the polygon
land_mask <- mask(land_mask, africa_arabia)
ggplot() +
  geom_spatraster(data = land_mask, aes(fill = land_mask_0)) +
  geom_sf(data = papio2)

####NAME CONVENTIONS####
papio_popn <- c("P. anubis", "P. cynocephalus", "P. hamadryas", "P. kindae", "P. papio", "P. ursinus", "Papio")

####MODELLING LOOP####
biglist <- list()
papio_split <- split(papio2, papio2$Species)
papio_split[[7]] <- papio2
#start loop here
for(j in 1:7){

all_papio <- list()
papio2 <- papio_split[[j]]
papio3 <- thin_by_cell(papio2, raster = land_mask)
e <- ext(min(papio2$Longitude)-2, max(papio2$Longitude)+2, min(papio2$Latitude)-2, max(papio2$Latitude)+2)
land_mask2 <- crop(land_mask, e)
plot(land_mask2)
points(papio2)
papio4 <- sample_pseudoabs(papio3,
                                 n = 6*nrow(papio3),
                                 raster = land_mask2,
                                 method = c("dist_min", km2m(5))
)

climate_vars <- get_vars_for_dataset("Krapp2021")
climate_vars <- climate_vars[!climate_vars=="biome"]
climate_vars <- climate_vars[!climate_vars=="altitude"]

climate_present <- pastclim::region_slice(
  time_bp = 0,
  bio_variables = climate_vars,
  data = "Krapp2021",
  crop = as.polygons(land_mask2)
)

plot(climate_present)

papio5 <- papio4 %>%bind_cols(terra::extract(climate_present, papio4, ID = FALSE))
papio5 %>% plot_pres_vs_bg(class)

papio5 %>% dist_pres_vs_bg(class)

vars_to_keep <- papio5 %>% dist_pres_vs_bg(class)
vars_to_keep <- names(vars_to_keep[vars_to_keep > 0.25])
papio6 <- papio5 %>% select(all_of(c(vars_to_keep, "class")))
vars_to_keep

climate_present <- climate_present[[vars_to_keep]]
vars_uncor <- filter_high_cor(climate_present, cutoff = 0.7)
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

lacerta_cv <- spatial_block_cv(lacerta_thin, v = 5)
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
ggplot() +
  geom_spatraster(data = prediction_present, aes(fill = mean)) +
  scale_fill_terrain_c() +
  geom_sf(data = lacerta_thin %>% filter(class == "presence"))

lacerta_ensemble <- calib_class_thresh(lacerta_ensemble,
                                       class_thresh = "tss_max"
)

climate_present_prediction <- pastclim::region_slice(
  time_bp = 0,
  bio_variables = vars_uncor,
  data = "Krapp2021",
  crop = africa_arabia
)

prediction_present_binary <- predict_raster(lacerta_ensemble,
                                            climate_present_prediction,
                                            type = "class",
                                            class_thresh = c("tss_max")
)
ggplot() +
  geom_spatraster(data = prediction_present_binary, aes(fill = binary_mean)) +
  geom_sf(data = lacerta_thin %>% filter(class == "presence"))

all_papio$papio_models <- lacerta_models
all_papio$papio_ensemble <- lacerta_ensemble
all_papio$papio_cv <- lacerta_cv
all_papio$papio_thin <- lacerta_thin
all_papio$vars_uncor <- vars_uncor
biglist[[j]] <- all_papio
print(j)}


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
####END MODELLING LOOP####

####EXPORT NC FILES OF PREDICTIONS####

for(k in 1:7){writeCDF(sds(biglist[[k]]$prediction_past_list), paste0(papio_popn[[k]], ".nc", sep=""))}

####EXPORT MODEL METRICS####
asda_list <- list()
for(k in 1:7){
  asda <- biglist[[k]]$papio_ensemble$metrics
  fasda <- as.data.frame(asda[[1]])
  for(i in 2:length(asda)){fasda <- rbind(fasda, as.data.frame(asda[[i]]))}
  asda_list[[k]] <- fasda}
for(i in 1:length(asda_list)){xlsx::write.xlsx(asda_list[[i]], paste0("Table 1a model_metrics_", Sys.Date(), ".xlsx"), sheetName = papio_popn[[i]], append = T)}

####EXPORT MEAN VARIABLE IMPORTANCE####
var_impotance_list <- list()
for(i in 1:7){
  vip_ensemble <- model_parts(explainer = explain_tidysdm(biglist[[i]]$papio_ensemble))
  variables <- split(vip_ensemble, vip_ensemble$variable)
  var_impotance_list[[i]] <- lapply(variables, function(x) mean(x$dropout_loss))}
ws2 <- data.frame(t(do.call(plyr::rbind.fill, lapply(var_impotance_list, as.data.frame))))
names(ws2) <- papio_popn
write.csv(ws2, paste0("Table 1b mean_importance_", Sys.Date(), ".csv"))