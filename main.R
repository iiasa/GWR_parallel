## Open Shell from Rstudio 
# 1. Go to the Git tab 
# 2. Click the tab More 
# 3. Select Shell... 
## Run Rscript from Shell 
# 1. In the Shell, run the command: nohup R CMD BATCH main.R > main.log & 
# 2. Close the Shell window 

library(tidyverse)
library(sf)
library(sp)
library(raster)
library(rgdal)
library(spgwr)
library(parallel)
source("pggwr.R")

# Read training_samplesamples as tibble 
training_samples <- readr::read_csv("./data/train_data_cci20.csv")

# CLS class labels 
cls_labels <- tibble::tribble(
  ~cls,              ~tr_data,
  "crops",                  4,
  "tree",                   1,
  "shrub",                  2,
  "grassland",              3,
  "wetland (herbaceous)",   5,
  "bare",                   7,
  "urban/built-up",         8,
  "water",                 10)

# Filter raster values equal to -9999 and deplace raster values 6 with 7 
training_samples <- training_samples %>% 
  dplyr::filter(RASTERVALU != -9999) %>%
  dplyr::mutate(RASTERVALU = ifelse(RASTERVALU == 6, 7, RASTERVALU))

# Join CLS values and labels 
training_samples <- training_samples %>% 
  dplyr::left_join(cls_labels) 

# Define logit function 
alogit <- function(x){exp(x)/(1+exp(x))}

# Define the new projection 
new_proj <- "+proj=longlat +datum=WGS84"

# Create spatial points 
training_samples <- training_samples %>% 
  sf::st_as_sf(coords = c("x", "y"), crs = new_proj, agr = "constant")

# Get coordinates 
coords <- sf::st_coordinates(training_samples)

# Create vgi input 
vgi_data <- 
  training_samples %>% 
  dplyr::transmute(value = ifelse(RASTERVALU == tr_data, 1, 0), tr_data = tr_data)

# Create vgi grid output  
grid_cont <- sf::st_make_grid(vgi_data, cellsize = 0.1, what = "centers") 

## Subset for testing 
# vgi_data <-
#   training_samples %>%
#   dplyr::slice(1:50) %>%
#   dplyr::transmute(value = ifelse(RASTERVALU == tr_data, 1, 0), tr_data = tr_data)
# grid_cont <- sf::st_make_grid(vgi_data, cellsize = 0.1, what = "centers")

# Run parallel 
poc_time <- system.time(gwr_model <- pggwr(formula = value~1, data = as(vgi_data, "Spatial"), adapt = 0.01, n_cores = 24L,
                                           fit.points = as(grid_cont, "Spatial"), family = binomial, longlat = TRUE))
poc_time

# Apply logistic transformation 
gwr.ov <- alogit(data.frame(gwr_model$SDF)[,2]) 

# Create spatial grid 
gwr.res.ov = SpatialPixelsDataFrame(gwr_model$SDF, data.frame(gwr.ov))

# Create raster 
r <- raster(gwr.res.ov)
plot(r)

# Write raster to disc
writeRaster(r, filename = "gwr_res_ov.tif", format='HFA', overwrite = TRUE)




