
library(here)
library(tidyverse)
library(ncdf4)
library(sf)
library(viridis)
library(ggpubr)
library(paran)
library(corrplot)
library(tidymodels)
library(glmnet)
library(nnet)
library(mgcv)
library(ranger)
library(vip)
library(factoextra)
library(brulee)

#### Load data ####

# Load covariates
feature_df <- read.csv(here::here("ML_inputs", "features_2025-01-14.csv"))
feature_df %>% dplyr::rename(Yr = X) -> feature_df

# Load each of the recruitment dfs
sablefish_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "sablefish_recdevs.csv"))
black_rf_nca_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "black_rf_nca_recdevs.csv"))
black_rf_cca_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "black_rf_cca_recdevs.csv"))
black_rf_or_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "black_rf_or_recdevs.csv"))
black_rf_wa_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "black_rf_wa_recdevs.csv"))
copper_rf_north_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "copper_rf_north_recdevs.csv"))
copper_rf_south_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "copper_rf_south_recdevs.csv"))
canary_rf_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "canary_rf_recdevs.csv"))
dover_sole_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "dover_sole_recdevs.csv"))
lingcod_north_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "lingcod_north_recdevs.csv"))
lingcod_south_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "lingcod_south_recdevs.csv"))
petrale_sole_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "petrale_sole_recdevs.csv"))
rex_sole_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "rex_sole_recdevs.csv"))
shortspine_thornyhead_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "shortspine_thornyhead_recdevs.csv"))
vermilion_rf_sca_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "vermilion_rf_sca_recdevs.csv"))
vermilion_rf_nca_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "vermilion_rf_nca_recdevs.csv"))
vermilion_rf_or_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "vermilion_rf_or_recdevs.csv"))
vermilion_rf_wa_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "vermilion_rf_wa_recdevs.csv"))
quillback_rf_ca_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "quillback_rf_ca_recdevs.csv"))
rougheye_rf_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "rougheye_rf_recdevs.csv"))
yellowtail_rf_cat_rec <- read.csv(here::here("ML_inputs", "recruitment_data", "yellowtail_rf_recdevs.csv"))

# join features with recruitment data
join_features_recdevs <- function(recdevs){
  full_join(recdevs, feature_df, by = "Yr") %>% 
    dplyr::filter(!(is.na(bbv_PC1)) & !(is.na(Value))) -> output
  
  return(output)
}

join_features_recdevs(sablefish_cat_rec) -> sablefish
join_features_recdevs(black_rf_nca_cat_rec) -> black_rf_nca
join_features_recdevs(black_rf_cca_cat_rec) -> black_rf_cca
join_features_recdevs(black_rf_or_cat_rec) -> black_rf_or
join_features_recdevs(black_rf_wa_cat_rec) -> black_rf_wa
join_features_recdevs(copper_rf_north_cat_rec) -> copper_rf_north
join_features_recdevs(copper_rf_south_cat_rec) -> copper_rf_south
join_features_recdevs(canary_rf_cat_rec) -> canary_rf
join_features_recdevs(dover_sole_cat_rec) -> dover_sole
join_features_recdevs(lingcod_north_cat_rec) -> lingcod_north
join_features_recdevs(lingcod_south_cat_rec) -> lingcod_south
join_features_recdevs(petrale_sole_cat_rec) -> petrale_sole
join_features_recdevs(rex_sole_cat_rec) -> rex_sole
join_features_recdevs(shortspine_thornyhead_cat_rec) -> shortspine_thornyhead
join_features_recdevs(vermilion_rf_sca_cat_rec) -> vermilion_rf_sca
join_features_recdevs(vermilion_rf_nca_cat_rec) -> vermilion_rf_nca
join_features_recdevs(vermilion_rf_or_cat_rec) -> vermilion_rf_or
join_features_recdevs(vermilion_rf_wa_cat_rec) -> vermilion_rf_wa
join_features_recdevs(quillback_rf_ca_cat_rec) -> quillback_rf_ca
join_features_recdevs(rougheye_rf_cat_rec) -> rougheye_rf
join_features_recdevs(yellowtail_rf_cat_rec) -> yellowtail_rf