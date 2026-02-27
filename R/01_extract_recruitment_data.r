# 01_extract_recruitment_data

# this script takes Stock Synthesis model files and extracts the recruitment residuals from the model.

# remotes::install_github("r4ss/r4ss")

library(here)
library(tidyverse)
library(r4ss)

# This script loads old assessments for the following species, although they were assessed in 2025, because the
# 2025 assessments because are not currently available:
# - Yellowtail Rockfish (not loaded, last assessed in 2017)
# - Sablefish (current data is from 2023)
# - Rougheye/Blackspotted Rockfish (cryptic species complex) (not loaded, last assessed in 2013)
# - Chilipepper Rockfish (not loaded, last assessed in 2025)
# - CA Quillback Rockfish (not loaded, 2025 was a major change to stock)

#### Define functions ####

# Identify how informed recruitment deviations are
identify_informed_recdevs <- function(SS_obj, sd_ratio_cutoff){
  # Get recruitment time series
  # extract df of parameters
  recdev_params <- SS_obj$parameters[
    substring(SS_obj$parameters[["Label"]], 1, 12) %in% c("Main_RecrDev"),
  ]
  
  # add year labels
  recdev_params[["Yr"]] <- as.numeric(substring(recdev_params[["Label"]], 14))
  # Record value of sigmaR
  sigmaR <- SS_obj$sigma_R_in
  
  # keep only key columns
  recdev_params %>% 
    dplyr::select(Yr, Value, Parm_StDev) -> recdevs
  
  # identify where recDevs are informed
  # do this based on the first value where it drops below the ratio, and the last value before
  # it goes back above the ratio. If you do it just based on the cutoff, you get 
  # a discontinuous time series where a couple of points in the middle are not as well informed.
  recdevs %>% 
    mutate(sd_ratio = Parm_StDev/sigmaR) %>% 
    mutate(informed_devs = ifelse(sd_ratio < sd_ratio_cutoff, "informed", "uninformed")) -> recdevs
  
  # find the first informed value
  first_informed_year <- min(subset(recdevs, informed_devs == "informed")$Yr)
  last_informed_year <- max(subset(recdevs, informed_devs == "informed")$Yr)
  
  recdevs %>% 
    mutate(keep = ifelse(Yr >= first_informed_year & Yr <= last_informed_year, TRUE, FALSE)) -> recdevs
  
  return(recdevs)
}

# Identify bias-adjusted recdevs
identify_bias_adj_recdevs <- function(SS_obj){
  # identify bias adjustment
  # option 1: select the full adjustment period
  first_yr_full <- round(as.numeric(SS_obj$breakpoints_for_bias_adjustment_ramp["first_yr_full"]))
  last_yr_full <- round(as.numeric(SS_obj$breakpoints_for_bias_adjustment_ramp["last_yr_full"]))
  
  # option 2: select the ramps
  first_yr_ramp <- round(as.numeric(SS_obj$breakpoints_for_bias_adjustment_ramp["last_yr_early"]))+1
  last_yr_ramp <- round(as.numeric(SS_obj$breakpoints_for_bias_adjustment_ramp["first_yr_recent"]))-1
  
  # Get recruitment time series
  # extract df of parameters
  recdev_params <- SS_obj$parameters[
    substring(SS_obj$parameters[["Label"]], 1, 12) %in% c("Main_RecrDev"),
  ]
  
  # add year labels
  recdev_params[["Yr"]] <- as.numeric(substring(recdev_params[["Label"]], 14))
  # Record value of sigmaR
  sigmaR <- SS_obj$sigma_R_in
  
  # keep only key columns
  recdev_params %>% 
    dplyr::select(Yr, Value, Parm_StDev) -> recdevs
  
  # identify which years are full bias adjustment
  recdevs %>% 
    mutate(keep = ifelse(Yr >= first_yr_ramp & Yr <= last_yr_ramp, TRUE, FALSE)) -> recdevs
  
  return(recdevs)
}

# show which recdevs were kept
plot_informed_recdevs <- function(recdevs, stock_name, SS_obj, sd_ratio_cutoff){
  
  sigmaR <- SS_obj$sigma_R_in
  
  informed_devs_sd_plot <- ggplot(recdevs, aes(x = Yr, y = Parm_StDev, color = keep)) +
    geom_point() +
    geom_hline(yintercept = sigmaR, lty = 2) +
    annotate(geom = "text", label = "SigmaR", x = max(recdevs$Yr), y = sigmaR+0.05, hjust = 1, vjust = 0) +
    geom_hline(yintercept = sd_ratio_cutoff * sigmaR, lty = 2, color = "red") +
    annotate(geom = "text", label = "Cutoff", x = max(recdevs$Yr), 
             y = sd_ratio_cutoff * sigmaR+0.05, hjust = 1, vjust = 0,
             color = "red") +
    geom_line(data = subset(recdevs, keep == TRUE), aes(x = Yr, y = Parm_StDev)) +
    scale_y_continuous(breaks = seq(0, round(max(recdevs$Parm_StDev), 1), 0.1), labels = seq(0, round(max(recdevs$Parm_StDev), 1), 0.1)) +
    scale_color_manual(values = c("gray80", "gray10")) +
    ggtitle(stock_name) +
    ylab("RecDev Estimate SD") +
    xlab("Year") +
    theme_minimal()
  
  informed_devs_plot <- ggplot(recdevs, aes(x = Yr, y = Value, color = keep)) +
    geom_point() +
    ylab("Log Recruitment Deviation") +
    geom_line(data = subset(recdevs, keep == TRUE), aes(x = Yr, y = Value)) +
    scale_color_manual(values = c("gray80", "gray10")) +
    ggtitle(stock_name) +
    theme_minimal() +
    xlab("Year")
    
  
  return(list(informed_devs_sd_plot, informed_devs_plot))
  
}




#### Models without a report file on the Council website ####

# hake_base_2024 <- SS_output(here::here("assessment_model_files", "hake-2024-model-files-zip-file"))


#### Models with a report file on the Council website ####

## Sablefish
sablefish_2025 <- SS_output(here::here("assessment_model_files", "2025_sablefish"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
sablefish_2025$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> sablefish_2025_SSB_recruitment
sablefish_2025$recruit -> sablefish_2025_recruit_detail


## Black Rockfish 2023 - Northern CA
black_rockfish_northern_ca_2023 <- SS_output(here::here("assessment_model_files", "black-rockfish-ca-2023-model-files-zip", "North"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
black_rockfish_northern_ca_2023$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> black_rockfish_northern_ca_2023_SSB_recruitment
black_rockfish_northern_ca_2023$recruit -> black_rockfish_northern_ca_2023_recruit_detail

## Black Rockfish 2023 - Central CA
black_rockfish_central_ca_2023 <- SS_output(here::here("assessment_model_files", "black-rockfish-ca-2023-model-files-zip", "Central"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
black_rockfish_central_ca_2023$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> black_rockfish_central_ca_2023_SSB_recruitment
black_rockfish_central_ca_2023$recruit -> black_rockfish_central_ca_2023_recruit_detail

## Black Rockfish 2023 - Oregon
black_rockfish_oregon_2023 <- SS_output(here::here("assessment_model_files", "black-rockfish-OR-2023"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
black_rockfish_oregon_2023$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> black_rockfish_oregon_2023_SSB_recruitment
black_rockfish_oregon_2023$recruit -> black_rockfish_oregon_2023_recruit_detail

## Black Rockfish 2023 - Washington
black_rockfish_washington_2023 <- SS_output(here::here("assessment_model_files", "black-rockfish-WA-2023"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
black_rockfish_washington_2023$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> black_rockfish_washington_2023_SSB_recruitment
black_rockfish_washington_2023$recruit -> black_rockfish_washington_2023_recruit_detail

## California Copper Rockfish 2023 - North of Point Conception
ca_copper_north_2023 <- SS_output(here::here("assessment_model_files", "ca-copper-rockfish-model-files-2023-zip-file", "nca"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
ca_copper_north_2023$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> ca_copper_north_2023_SSB_recruitment
ca_copper_north_2023$recruit -> ca_copper_north_2023_recruit_detail

## California Copper Rockfish 2023 - South
ca_copper_south_2023 <- SS_output(here::here("assessment_model_files", "ca-copper-rockfish-model-files-2023-zip-file", "sca"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
ca_copper_south_2023$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> ca_copper_south_2023_SSB_recruitment
ca_copper_south_2023$recruit -> ca_copper_south_2023_recruit_detail

## Canary Rockfish 2023
canary_rockfish_2023 <- SS_output(here::here("assessment_model_files", "canary-rockfish-2023-model-files-zip-file"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
canary_rockfish_2023$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> canary_rockfish_2023_SSB_recruitment
canary_rockfish_2023$recruit -> canary_rockfish_2023_recruit_detail


## Dover Sole 2021
dover_sole_2021 <- SS_output(here::here("assessment_model_files", "DoverSole_2021_Model_For_DeVore"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
dover_sole_2021$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> dover_sole_2021_SSB_recruitment
dover_sole_2021$recruit -> dover_sole_2021_recruit_detail

## Lingcod 2021 - North of 40°10'N
lingcod_north_2021 <- SS_output(here::here("assessment_model_files", "Lingcod_2021_Model_For_DeVore", "north"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
lingcod_north_2021$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> lingcod_north_2021_SSB_recruitment
lingcod_north_2021$recruit -> lingcod_north_2021_recruit_detail

## Lingcod 2021 - South of 40°10'N
lingcod_south_2021 <- SS_output(here::here("assessment_model_files", "Lingcod_2021_Model_For_DeVore", "south"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
lingcod_south_2021$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> lingcod_south_2021_SSB_recruitment
lingcod_south_2021$recruit -> lingcod_south_2021_recruit_detail

## Oregon Copper Rockfish 2021
or_copper_rockfish_2021 <- SS_output(here::here("assessment_model_files", "oregon-copper-rockfish-2021"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
or_copper_rockfish_2021$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> or_copper_rockfish_2021_SSB_recruitment
or_copper_rockfish_2021$recruit -> or_copper_rockfish_2021_recruit_detail

## Petrale Sole 2023
petrale_sole_2023 <- SS_output(here::here("assessment_model_files", "petrale-sole-2023-model-files-zip-file", "outputs"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
petrale_sole_2023$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> petrale_sole_2023_SSB_recruitment
petrale_sole_2023$recruit -> petrale_sole_2023_recruit_detail

## Rex Sole 2023
rex_sole_2023 <- SS_output(here::here("assessment_model_files", "rex_sole_2023_base_model", "run"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
rex_sole_2023$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> rex_sole_2023_SSB_recruitment
rex_sole_2023$recruit -> rex_sole_2023_recruit_detail

## Shortspine Thornyhead 2023
shortspine_thornyhead_2023 <- SS_output(here::here("assessment_model_files", "shortspine-thornyhead-2023-model-files-zip-file"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
shortspine_thornyhead_2023$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> shortspine_thornyhead_2023_SSB_recruitment
shortspine_thornyhead_2023$recruit -> shortspine_thornyhead_2023_recruit_detail

## Spiny Dogfish 2021
spiny_dogfish_2021 <- SS_output(here::here("assessment_model_files", "Spiny_Dogfish_Model_For_DeVore"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
spiny_dogfish_2021$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> spiny_dogfish_2021_SSB_recruitment
spiny_dogfish_2021$recruit -> spiny_dogfish_2021_recruit_detail

## Squarespot Rockfish 2021
squarespot_rockfish_2021 <- SS_output(here::here("assessment_model_files", "Squarespot_rockfish_2021_Model_For_Devore"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
squarespot_rockfish_2021$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> squarespot_rockfish_2021_SSB_recruitment
squarespot_rockfish_2021$recruit -> squarespot_rockfish_2021_recruit_detail

## Vermilion & Sunset Rockfish 2021 - Southern California
vermilion_SCA_2021 <- SS_output(here::here("assessment_model_files", "Verm_Sunset_2021_SCA_Model_For_DeVore"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
vermilion_SCA_2021$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> vermilion_SCA_2021_SSB_recruitment
vermilion_SCA_2021$recruit -> vermilion_SCA_2021_recruit_detail

## Vermilion & Sunset Rockfish 2021 - Northern California
vermilion_NCA_2021 <- SS_output(here::here("assessment_model_files", "Vermilion_NCA_2021_Model_For_DeVore"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
vermilion_NCA_2021$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> vermilion_NCA_2021_SSB_recruitment
vermilion_NCA_2021$recruit -> vermilion_NCA_2021_recruit_detail

# Vermilion Rockfish 2021 - Oregon
vermilion_OR_2021 <- SS_output(here::here("assessment_model_files", "Vermilion_OR_2021_Model_For_DeVore"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
vermilion_OR_2021$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> vermilion_OR_2021_SSB_recruitment
vermilion_OR_2021$recruit -> vermilion_OR_2021_recruit_detail

# Vermilion Rockfish 2021 - Washington
vermilion_WA_2021 <- SS_output(here::here("assessment_model_files", "Vermilion_WA_2021_Model_For_DeVore"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
vermilion_WA_2021$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> vermilion_WA_2021_SSB_recruitment
vermilion_WA_2021$recruit -> vermilion_WA_2021_recruit_detail

# Copper Rockfish 2021 - Washington
copper_WA_2021 <- SS_output(here::here("assessment_model_files", "washington-copper-rockfish-2021"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
copper_WA_2021$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> copper_WA_2021_SSB_recruitment
copper_WA_2021$recruit -> copper_WA_2021_recruit_detail

# Quillback Rockfish 2025 - California
quillback_CA_2025 <- SS_output(here::here("assessment_model_files", "ca-quillback-2025-model-files-zip-file"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
quillback_CA_2025$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> quillback_CA_2025_SSB_recruitment
quillback_CA_2025$recruit -> quillback_CA_2025_recruit_detail

# Rougheye/Blackspotted Rockfish 2025
rougheye_2025 <- SS_output(here::here("assessment_model_files", "rougheye-blackspotted-2025-model-files-zip-file"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
rougheye_2025$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> rougheye_2025_SSB_recruitment
rougheye_2025$recruit -> rougheye_2025_recruit_detail

# Yellowtail Rockfish 2025
yellowtail_2025 <- SS_output(here::here("assessment_model_files", "yellowtail-rockfish-2025"), verbose = FALSE, printstats = FALSE, hidewarn = TRUE)
yellowtail_2025$timeseries %>% 
  dplyr::select(Yr, SpawnBio, Recruit_0) -> yellowtail_2025_SSB_recruitment
yellowtail_2025$recruit -> yellowtail_2025_recruit_detail


# Pacific Hake 2025 - special treatment
#' Create a data frame containing model output recruitment from a
#' model list which is in long format ready for [ggplot2::ggplot()]
#'
#' @param model_lst A list of models, each created by [create_rds_file()]
#' @param model_names A vector of model names,the same length as `model_lst`
#' @param devs Logical. If `TRUE` return recruitment deviations, if `FALSE`,
#' return absolute recruitment
#' @param ... Arguments passed to [extract_mcmc_quant()]
#'
#' @return A list containing a [tibble::tibble()]
#'
#' @export
create_group_df_recr <- function(model_lst = NULL,
                                 model_names = NULL,
                                 devs = FALSE,
                                 relative = FALSE,
                                 ...){
  
  vals <- paste0(ifelse(devs,
                        "dev",
                        ifelse(relative,
                               "r_rel_",
                               "r")),
                 c("lower", "med", "upper"))
  if(!devs){
    # Needed for the x's on the recruitment plot
    vals <-  c(vals, ifelse(relative,
                            "r_rel_mean",
                            "rmean"))
  }
  
  d <- bind_cols(extract_mcmc_quant(model_lst,
                                    model_names,
                                    vals[1],
                                    TRUE),
                 extract_mcmc_quant(model_lst,
                                    model_names,
                                    vals[2]),
                 extract_mcmc_quant(model_lst,
                                    model_names,
                                    vals[3]))
  if(!devs){
    d <- bind_cols(d,
                   extract_mcmc_quant(model_lst,
                                      model_names,
                                      vals[4]))
  }
  d <- d |>
    mutate(model = factor(model, levels = model_names),
           year = as.numeric(year))
  
  list(d = d)
}
#' Extract MCMC quantiles from models and return them
#' in a [ggplot2::ggplot()]-ready data frame
#'
#' @param model_lst A list of models, each created by [create_rds_file()]
#' @param model_names A vector of model names,the same length as `models_lst`
#' @param type A name as found in an `mcmccalcs` object of a model object,
#' for example extraction of: `base_model$mcmccalcs$smed` would require
#' `type` = "smed"
#' @param inc_model_year Logical. If `TRUE`, include the `model` and `year`
#' columns in the data frame
#' @param end_yrs A vector of the end years for each model.
#' If one value, it will apply to all models
#'
#' @return A [tibble::tibble()]
#' @export
extract_mcmc_quant <- function(model_lst,
                               model_names,
                               type,
                               inc_model_year = FALSE,
                               end_yrs = 2025){
  
  if(length(end_yrs) == 1){
    end_yrs <- rep(end_yrs, length(model_lst))
  }
  
  if(length(end_yrs) != length(model_lst)){
    stop("Length of `end_yrs` does not equal the length ",
         "of `model_lst`")
  }
  
  out <- map2(model_lst, end_yrs, ~{
    tmp <- .x$mcmccalcs[[type]]
    type_nms <- names(tmp)
    type_nms <- type_nms[suppressWarnings(!is.na(as.numeric(type_nms)))]
    tmp <- tmp[type_nms]
    tmp[as.numeric(type_nms) <= .y]
  }) |>
    map_dfr(~{.x}) |>
    mutate(model = model_names) |>
    dplyr::select(model, everything()) |>
    pivot_longer(-"model",
                 names_to = "year",
                 values_to = type)
  
  if(inc_model_year){
    out
  }else{
    out |>
      dplyr::select(-c(model, year))
  }
}
# hake_2025 <- readRDS(here::here("assessment_model_files", "pacific_hake_2025", "01-base.rds"))
# d_obj <- create_group_df_recr(model_lst = list(hake_2025),
#                               model_names = list(hake_2025), devs = TRUE)

# NOTE: Some of these don't seem to have any main recruitment deviations estimated? I think 
# these just aren't estimating recruitment deviations but just taking values
# from the SR relationship
# or_copper_rockfish_2021
# spiny_dogfish_2021
# squarespot_rockfish_2021
# copper_WA_2021
# we should drop these from our analysis

# sablefish_2025_recdevs <- identify_informed_recdevs(SS_obj = sablefish_2025, sd_ratio_cutoff = 0.8)
# black_rockfish_northern_ca_2023_recdevs <- identify_informed_recdevs(SS_obj = black_rockfish_northern_ca_2023, sd_ratio_cutoff = 0.8)
# black_rockfish_central_ca_2023_recdevs <- identify_informed_recdevs(SS_obj = black_rockfish_central_ca_2023, sd_ratio_cutoff = 0.8)
# black_rockfish_oregon_2023_recdevs <- identify_informed_recdevs(SS_obj = black_rockfish_oregon_2023, sd_ratio_cutoff = 0.8)
# black_rockfish_washington_2023_recdevs <- identify_informed_recdevs(SS_obj = black_rockfish_washington_2023, sd_ratio_cutoff = 0.8)
# ca_copper_north_2023_recdevs <- identify_informed_recdevs(SS_obj = ca_copper_north_2023, sd_ratio_cutoff = 0.8)
# ca_copper_south_2023_recdevs <- identify_informed_recdevs(SS_obj = ca_copper_south_2023, sd_ratio_cutoff = 0.8)
# canary_rockfish_2023_recdevs <- identify_informed_recdevs(SS_obj = canary_rockfish_2023, sd_ratio_cutoff = 0.8)
# dover_sole_2021_recdevs <- identify_informed_recdevs(SS_obj = dover_sole_2021, sd_ratio_cutoff = 0.8)
# lingcod_north_2021_recdevs <- identify_informed_recdevs(SS_obj = lingcod_north_2021, sd_ratio_cutoff = 0.8)
# lingcod_south_2021_recdevs <- identify_informed_recdevs(SS_obj = lingcod_south_2021, sd_ratio_cutoff = 0.8)
# or_copper_rockfish_2021_recdevs <- identify_informed_recdevs(SS_obj = or_copper_rockfish_2021, sd_ratio_cutoff = 0.8)
# petrale_sole_2023_recdevs <- identify_informed_recdevs(SS_obj = petrale_sole_2023, sd_ratio_cutoff = 0.8)
# rex_sole_2023_recdevs <- identify_informed_recdevs(SS_obj = rex_sole_2023, sd_ratio_cutoff = 0.8)
# shortspine_thornyhead_2023_recdevs <- identify_informed_recdevs(SS_obj = shortspine_thornyhead_2023, sd_ratio_cutoff = 0.8)
# spiny_dogfish_2021_recdevs <- identify_informed_recdevs(SS_obj = spiny_dogfish_2021, sd_ratio_cutoff = 0.8)
# squarespot_rockfish_2021_recdevs <- identify_informed_recdevs(SS_obj = squarespot_rockfish_2021, sd_ratio_cutoff = 0.8)
# vermilion_SCA_2021_recdevs <- identify_informed_recdevs(SS_obj = vermilion_SCA_2021, sd_ratio_cutoff = 0.8)
# vermilion_NCA_2021_recdevs <- identify_informed_recdevs(SS_obj = vermilion_NCA_2021, sd_ratio_cutoff = 0.8)
# vermilion_OR_2021_recdevs <- identify_informed_recdevs(SS_obj = vermilion_OR_2021, sd_ratio_cutoff = 0.8)
# vermilion_WA_2021_recdevs <- identify_informed_recdevs(SS_obj = vermilion_WA_2021, sd_ratio_cutoff = 0.8)
# copper_WA_2021_recdevs <- identify_informed_recdevs(SS_obj = copper_WA_2021, sd_ratio_cutoff = 0.8)
# quillback_CA_2025_recdevs <- identify_informed_recdevs(SS_obj = quillback_CA_2025, sd_ratio_cutoff = 0.8)
# rougheye_2025_recdevs <- identify_informed_recdevs(SS_obj = rougheye_2025, sd_ratio_cutoff = 0.8)
# yellowtail_2025_recdevs <- identify_informed_recdevs(SS_obj = yellowtail_2025, sd_ratio_cutoff = 0.8)

sablefish_2025_recdevs <- identify_bias_adj_recdevs(SS_obj = sablefish_2025)
black_rockfish_northern_ca_2023_recdevs <- identify_bias_adj_recdevs(SS_obj = black_rockfish_northern_ca_2023)
black_rockfish_central_ca_2023_recdevs <- identify_bias_adj_recdevs(SS_obj = black_rockfish_central_ca_2023)
black_rockfish_oregon_2023_recdevs <- identify_bias_adj_recdevs(SS_obj = black_rockfish_oregon_2023)
black_rockfish_washington_2023_recdevs <- identify_bias_adj_recdevs(SS_obj = black_rockfish_washington_2023)
ca_copper_north_2023_recdevs <- identify_bias_adj_recdevs(SS_obj = ca_copper_north_2023)
ca_copper_south_2023_recdevs <- identify_bias_adj_recdevs(SS_obj = ca_copper_south_2023)
canary_rockfish_2023_recdevs <- identify_bias_adj_recdevs(SS_obj = canary_rockfish_2023)
dover_sole_2021_recdevs <- identify_bias_adj_recdevs(SS_obj = dover_sole_2021)
lingcod_north_2021_recdevs <- identify_bias_adj_recdevs(SS_obj = lingcod_north_2021)
lingcod_south_2021_recdevs <- identify_bias_adj_recdevs(SS_obj = lingcod_south_2021)
or_copper_rockfish_2021_recdevs <- identify_bias_adj_recdevs(SS_obj = or_copper_rockfish_2021)
petrale_sole_2023_recdevs <- identify_bias_adj_recdevs(SS_obj = petrale_sole_2023)
rex_sole_2023_recdevs <- identify_bias_adj_recdevs(SS_obj = rex_sole_2023)
shortspine_thornyhead_2023_recdevs <- identify_bias_adj_recdevs(SS_obj = shortspine_thornyhead_2023)
spiny_dogfish_2021_recdevs <- identify_bias_adj_recdevs(SS_obj = spiny_dogfish_2021)
squarespot_rockfish_2021_recdevs <- identify_bias_adj_recdevs(SS_obj = squarespot_rockfish_2021)
vermilion_SCA_2021_recdevs <- identify_bias_adj_recdevs(SS_obj = vermilion_SCA_2021)
vermilion_NCA_2021_recdevs <- identify_bias_adj_recdevs(SS_obj = vermilion_NCA_2021)
vermilion_OR_2021_recdevs <- identify_bias_adj_recdevs(SS_obj = vermilion_OR_2021)
vermilion_WA_2021_recdevs <- identify_bias_adj_recdevs(SS_obj = vermilion_WA_2021)
copper_WA_2021_recdevs <- identify_bias_adj_recdevs(SS_obj = copper_WA_2021)
quillback_CA_2025_recdevs <- identify_bias_adj_recdevs(SS_obj = quillback_CA_2025)
rougheye_2025_recdevs <- identify_bias_adj_recdevs(SS_obj = rougheye_2025)
yellowtail_2025_recdevs <- identify_bias_adj_recdevs(SS_obj = yellowtail_2025)



# visualize informed recdevs
sablefish_2025_recdev_plots <- plot_informed_recdevs(recdevs = sablefish_2025_recdevs, stock_name = "Sablefish (2025)", SS_obj = sablefish_2025, sd_ratio_cutoff = 0.8)
black_rockfish_northern_ca_2023_recdev_plots <- plot_informed_recdevs(recdevs = black_rockfish_northern_ca_2023_recdevs, stock_name = "Black Rockfish NCA (2023)", SS_obj = black_rockfish_northern_ca_2023, sd_ratio_cutoff = 0.8)
black_rockfish_central_ca_2023_recdev_plots <- plot_informed_recdevs(recdevs = black_rockfish_central_ca_2023_recdevs, stock_name = "Black Rockfish CCA (2023)", SS_obj = black_rockfish_central_ca_2023, sd_ratio_cutoff = 0.8)
black_rockfish_oregon_2023_recdev_plots <- plot_informed_recdevs(recdevs = black_rockfish_oregon_2023_recdevs, stock_name = "Black Rockfish OR (2023)", SS_obj = black_rockfish_oregon_2023, sd_ratio_cutoff = 0.8)
black_rockfish_washington_2023_recdev_plots <- plot_informed_recdevs(recdevs = black_rockfish_washington_2023_recdevs, stock_name = "Black Rockfish WA (2023)", SS_obj = black_rockfish_washington_2023, sd_ratio_cutoff = 0.8)
ca_copper_north_2023_recdev_plots <- plot_informed_recdevs(recdevs = ca_copper_north_2023_recdevs, stock_name = "Copper Rockfish NCA (2023)", SS_obj = ca_copper_north_2023, sd_ratio_cutoff = 0.8)
ca_copper_south_2023_recdev_plots <- plot_informed_recdevs(recdevs = ca_copper_south_2023_recdevs, stock_name = "Copper Rockfish SCA (2023)", SS_obj = ca_copper_south_2023, sd_ratio_cutoff = 0.8)
canary_rockfish_2023_recdev_plots <- plot_informed_recdevs(recdevs = canary_rockfish_2023_recdevs, stock_name = "Canary Rockfish (2023)", SS_obj = canary_rockfish_2023, sd_ratio_cutoff = 0.8)
dover_sole_2021_recdev_plots <- plot_informed_recdevs(recdevs = dover_sole_2021_recdevs, stock_name = "Dover Sole (2021)", SS_obj = dover_sole_2021, sd_ratio_cutoff = 0.8)
lingcod_north_2021_recdev_plots <- plot_informed_recdevs(recdevs = lingcod_north_2021_recdevs, stock_name = "Lingcod North (2021)", SS_obj = lingcod_north_2021, sd_ratio_cutoff = 0.8)
lingcod_south_2021_recdev_plots <- plot_informed_recdevs(recdevs = lingcod_south_2021_recdevs, stock_name = "Lingcod South (2021)", SS_obj = lingcod_south_2021, sd_ratio_cutoff = 0.8)
# plot_informed_recdevs(recdevs = or_copper_rockfish_2021_recdevs, stock_name = "Copper Rockfish OR (2021)")
petrale_sole_2023_recdev_plots <- plot_informed_recdevs(recdevs = petrale_sole_2023_recdevs, stock_name = "Petrale Sole (2023)", SS_obj = petrale_sole_2023, sd_ratio_cutoff = 0.8)
rex_sole_2023_recdev_plots <- plot_informed_recdevs(recdevs = rex_sole_2023_recdevs, stock_name = "Rex Sole (2023)", SS_obj = rex_sole_2023, sd_ratio_cutoff = 0.8)
shortspine_thornyhead_2023_recdev_plots <- plot_informed_recdevs(recdevs = shortspine_thornyhead_2023_recdevs, stock_name = "Shortspine Thornyhead (2023)", SS_obj = shortspine_thornyhead_2023, sd_ratio_cutoff = 0.8)
# spiny_dogfish_2021_recdevs <- identify_informed_recdevs(SS_obj = spiny_dogfish_2021, sd_ratio_cutoff = 0.8)
# squarespot_rockfish_2021_recdevs <- identify_informed_recdevs(SS_obj = squarespot_rockfish_2021, sd_ratio_cutoff = 0.8)
vermilion_SCA_2021_recdev_plots <- plot_informed_recdevs(recdevs = vermilion_SCA_2021_recdevs, stock_name = "Vermilion Rockfish SCA (2021)", SS_obj = vermilion_SCA_2021, sd_ratio_cutoff = 0.8)
vermilion_NCA_2021_recdev_plots <- plot_informed_recdevs(recdevs = vermilion_NCA_2021_recdevs, stock_name = "Vermilion Rockfish NCA (2021)", SS_obj = vermilion_NCA_2021, sd_ratio_cutoff = 0.8)
vermilion_OR_2021_recdev_plots <- plot_informed_recdevs(recdevs = vermilion_OR_2021_recdevs, stock_name = "Vermilion Rockfish OR (2021)", SS_obj = vermilion_OR_2021, sd_ratio_cutoff = 0.8)
vermilion_WA_2021_recdev_plots <- plot_informed_recdevs(recdevs = vermilion_WA_2021_recdevs, stock_name = "Vermilion Rockfish WA (2021)", SS_obj = vermilion_WA_2021, sd_ratio_cutoff = 0.8)
# copper_WA_2021_recdevs <- identify_informed_recdevs(SS_obj = copper_WA_2021, sd_ratio_cutoff = 0.8)
quillback_CA_2025_recdev_plots <- plot_informed_recdevs(recdevs = quillback_CA_2025_recdevs, stock_name = "Quillback Rockfish CA (2025)", SS_obj = quillback_CA_2025, sd_ratio_cutoff = 0.8)
rougheye_2025_recdev_plots <- plot_informed_recdevs(recdevs = rougheye_2025_recdevs, stock_name = "Rougheye/Blackspotted Rockfish (2025)", SS_obj = rougheye_2025, sd_ratio_cutoff = 0.8)
yellowtail_2025_recdev_plots <- plot_informed_recdevs(recdevs = yellowtail_2025_recdevs, stock_name = "Yellowtail Rockfish (2025)", SS_obj = yellowtail_2025, sd_ratio_cutoff = 0.8)




#### Visualize data ####

recr_cat_colors <- c("Bust" = "#762a83",
                     "low" = "#762a83",
                     "medium" = "gray90",
                     "Average" = "gray90",
                     "Baseline Recruitment" = "gray90",
                     "Boom" = "#006d2c",
                     "high" = "#006d2c",
                     "Recruitment Event" = "#006d2c")

# plot time series of recruitment
plot_recruit_ts <- function(recdevs, stock_name, recruit_var = dev, category_name = NULL){
  
  recdev_ts_plot <- ggplot(recdevs, aes(x = Yr, y = {{recruit_var}}, color = {{category_name}})) +
    geom_line(color = "black", lty = 2) +
    geom_point(size = 2.5) +
    scale_color_manual(values = recr_cat_colors) +
    theme_minimal() +
    # ggtitle(stock_name) +
    ylab("Log Recruitment Deviation") +
    xlab("Year") +
    theme(legend.position = "bottom",
          legend.direction = "vertical")
  
  
  return(recdev_ts_plot)
}

# plot standardized recruitment deviations
plot_recdevs <- function(recruit_detail, stock_name){
  recruit_detail_main <- subset(recruit_detail, era == "Main")
  
  # standardize recruitment time series (z-score)
  recruit_detail_main$pred_recr_scaled <- scale(recruit_detail_main$pred_recr)
  
  recdev_std_plot <- ggplot(recruit_detail_main, aes(x = Yr, y = pred_recr_scaled)) +
    geom_line(color = "black", lty = 2) +
    geom_point(size = 2.5) +
    theme_minimal() +
    # ggtitle(stock_name) +
    ylab("Recruitment Deviations") +
    xlab("Year")
  
  return(recdev_std_plot)
}



#### Classify recruitment as categorical ####

# Stachura et al. 2014: weakest (<25%), average (middle 50%), and strongest (>75%)
stachura_classify <- function(recdevs){
  
  recdevs_keep <- subset(recdevs, keep == TRUE)
  
  stachura_breaks <- quantile(recdevs_keep$Value, probs = c(0.25, 0.75))
  recdevs_keep %>% 
    mutate(stachura_class = cut(Value, 
                                breaks = c(-Inf, stachura_breaks, Inf),
                                labels = c("low", "medium", "high"))) -> recdevs_keep
  
  return(recdevs_keep)
  
}


stachura_classify_pred_recr <- function(recruit_detail, recdevs){
  
  years_to_keep <- subset(recdevs, keep == TRUE)$Yr
  
  recruit_detail_main <- subset(recruit_detail, Yr %in% years_to_keep)
  
  stachura_breaks <- quantile(recruit_detail_main$pred_recr, probs = c(0.25, 0.75))
  recruit_detail_main %>% 
    mutate(stachura_class = cut(pred_recr, 
                                breaks = c(-Inf, stachura_breaks, Inf),
                                labels = c("low", "medium", "high"))) -> recruit_detail_main
  
  return(recruit_detail_main)
  
}

# Boom = 1 SD above mean, bust = 1 SD below mean, average = -1 to 1 SD
# one_sd_boom_bust_classify <- function(recdevs){
#   
#   recdevs_keep <- subset(recdevs, keep == TRUE)
#   
#   # standardize recruitment time series (z-score)
#   recdevs_keep$Value_scaled <- scale(recdevs_keep$Value)
#   
#   recdevs_keep %>% 
#     mutate(one_sd_boom_bust = cut(Value_scaled,
#                                   breaks = c(-Inf, -1, 1, Inf),
#                                   labels = c("Bust", "Average", "Boom"))) -> recdevs_keep
#   
#   return(recdevs_keep)
#   
# }

# one_sd_boom_bust_classify_pred_recr <- function(recruit_detail, recdevs){
#   
#   years_to_keep <- subset(recdevs, keep == TRUE)$Yr
#   
#   recruit_detail_main <- subset(recruit_detail, Yr %in% years_to_keep)
#   
#   # standardize recruitment time series (z-score)
#   recruit_detail_main$pred_recr_scaled <- scale(recruit_detail_main$pred_recr)
#   
#   recruit_detail_main %>% 
#     mutate(one_sd_boom_bust = cut(pred_recr_scaled,
#                                   breaks = c(-Inf, -1, 1, Inf),
#                                   labels = c("Bust", "Average", "Boom"))) -> recruit_detail_main
#   
#   return(recruit_detail_main)
#   
# }

one_sd_boom_bust_classify <- function(recdevs){
  
  years_to_keep <- subset(recdevs, keep == TRUE)$Yr
  
  recdevs <- subset(recdevs, Yr %in% years_to_keep)
  
  # standardize recruitment time series (z-score)
  recdevs$recr_scaled <- scale(recdevs$Value)
  
  recdevs %>% 
    mutate(one_sd_boom_bust = cut(recr_scaled,
                                  breaks = c(-Inf, -1, 1, Inf),
                                  labels = c("Bust", "Average", "Boom"))) %>% 
    dplyr::select(-keep) -> recdevs
  
  return(recdevs)
  
}

# Episodic: Boom = 1 SD above mean
episodic_classify <- function(recdevs){
  
  recdevs_keep <- subset(recdevs, keep == TRUE)
  
  # standardize recruitment time series (z-score)
  recdevs_keep$Value_scaled <- scale(recdevs_keep$Value)
  
  recdevs_keep %>% 
    mutate(episodic_recruitment = cut(Value_scaled,
                                      breaks = c(-Inf, 1, Inf),
                                      labels = c("Baseline Recruitment", "Recruitment Event"))) -> recdevs_keep
  
  return(recdevs_keep)
  
}

episodic_classify_pred_recr <- function(recruit_detail, recdevs){
  
  years_to_keep <- subset(recdevs, keep == TRUE)$Yr
  
  recruit_detail_main <- subset(recruit_detail, Yr %in% years_to_keep)
  
  # standardize recruitment time series (z-score)
  recruit_detail_main$pred_recr_scaled <- scale(recruit_detail_main$pred_recr)
  
  recruit_detail_main %>% 
    mutate(episodic_recruitment = cut(pred_recr_scaled,
                                      breaks = c(-Inf, 1, Inf),
                                      labels = c("Baseline Recruitment", "Recruitment Event"))) -> recruit_detail_main
  
  return(recruit_detail_main)
  
}


#### Classify each time series and export ####

one_sd_boom_bust_classify(sablefish_2025_recdevs) -> sablefish_cat_rec

one_sd_boom_bust_classify(black_rockfish_northern_ca_2023_recdevs) -> black_rf_nca_cat_rec

one_sd_boom_bust_classify(black_rockfish_central_ca_2023_recdevs) -> black_rf_cca_cat_rec

one_sd_boom_bust_classify(black_rockfish_oregon_2023_recdevs) -> black_rf_or_cat_rec

one_sd_boom_bust_classify(black_rockfish_washington_2023_recdevs) -> black_rf_wa_cat_rec

one_sd_boom_bust_classify(ca_copper_north_2023_recdevs) -> copper_rf_north_cat_rec

one_sd_boom_bust_classify(ca_copper_south_2023_recdevs) -> copper_rf_south_cat_rec

one_sd_boom_bust_classify(canary_rockfish_2023_recdevs) -> canary_rf_cat_rec

one_sd_boom_bust_classify(dover_sole_2021_recdevs) -> dover_sole_cat_rec

one_sd_boom_bust_classify(lingcod_north_2021_recdevs) -> lingcod_north_cat_rec

one_sd_boom_bust_classify(lingcod_south_2021_recdevs) -> lingcod_south_cat_rec

# one_sd_boom_bust_classify(or_copper_rockfish_2021_recdevs) -> copper_rf_or_cat_rec

one_sd_boom_bust_classify(petrale_sole_2023_recdevs) -> petrale_sole_cat_rec

one_sd_boom_bust_classify(rex_sole_2023_recdevs) -> rex_sole_cat_rec

one_sd_boom_bust_classify(shortspine_thornyhead_2023_recdevs) -> shortspine_thornyhead_cat_rec

# one_sd_boom_bust_classify(spiny_dogfish_2021_recdevs) -> spiny_dogfish_cat_rec

# one_sd_boom_bust_classify(squarespot_rockfish_2021_recdevs) -> squarespot_rf_cat_rec

one_sd_boom_bust_classify(vermilion_SCA_2021_recdevs) -> vermilion_rf_sca_cat_rec

one_sd_boom_bust_classify(vermilion_NCA_2021_recdevs) -> vermilion_rf_nca_cat_rec

one_sd_boom_bust_classify(vermilion_OR_2021_recdevs) -> vermilion_rf_or_cat_rec

one_sd_boom_bust_classify(vermilion_WA_2021_recdevs) -> vermilion_rf_wa_cat_rec

# one_sd_boom_bust_classify(copper_WA_2021_recdevs) -> copper_rf_wa_cat_rec

one_sd_boom_bust_classify(quillback_CA_2025_recdevs) -> quillback_rf_ca_cat_rec

one_sd_boom_bust_classify(rougheye_2025_recdevs) -> rougheye_rf_cat_rec

one_sd_boom_bust_classify(yellowtail_2025_recdevs) -> yellowtail_rf_cat_rec

# export each time series
write.csv(sablefish_cat_rec, here::here("ML_inputs", "recruitment_data", "sablefish_recdevs.csv"), row.names = FALSE)
write.csv(black_rf_nca_cat_rec, here::here("ML_inputs", "recruitment_data", "black_rf_nca_recdevs.csv"), row.names = FALSE)
write.csv(black_rf_or_cat_rec, here::here("ML_inputs", "recruitment_data", "black_rf_or_recdevs.csv"), row.names = FALSE)
write.csv(black_rf_wa_cat_rec, here::here("ML_inputs", "recruitment_data", "black_rf_wa_recdevs.csv"), row.names = FALSE)
write.csv(copper_rf_north_cat_rec, here::here("ML_inputs", "recruitment_data", "copper_rf_north_recdevs.csv"), row.names = FALSE)
write.csv(copper_rf_south_cat_rec, here::here("ML_inputs", "recruitment_data", "copper_rf_south_recdevs.csv"), row.names = FALSE)
write.csv(canary_rf_cat_rec, here::here("ML_inputs", "recruitment_data", "canary_rf_recdevs.csv"), row.names = FALSE)
write.csv(dover_sole_cat_rec, here::here("ML_inputs", "recruitment_data", "dover_sole_recdevs.csv"), row.names = FALSE)
write.csv(lingcod_north_cat_rec, here::here("ML_inputs", "recruitment_data", "lingcod_north_recdevs.csv"), row.names = FALSE)
write.csv(lingcod_south_cat_rec, here::here("ML_inputs", "recruitment_data", "lingcod_south_recdevs.csv"), row.names = FALSE)
# write.csv(copper_rf_or_cat_rec, here::here("ML_inputs", "recruitment_data", "copper_rf_or_recdevs.csv"), row.names = FALSE)
write.csv(petrale_sole_cat_rec, here::here("ML_inputs", "recruitment_data", "petrale_sole_recdevs.csv"), row.names = FALSE)
write.csv(rex_sole_cat_rec, here::here("ML_inputs", "recruitment_data", "rex_sole_recdevs.csv"), row.names = FALSE)
write.csv(shortspine_thornyhead_cat_rec, here::here("ML_inputs", "recruitment_data", "shortspine_thornyhead_recdevs.csv"), row.names = FALSE)
# write.csv(spiny_dogfish_cat_rec, here::here("ML_inputs", "recruitment_data", "spiny_dogfish_recdevs.csv"), row.names = FALSE)
# write.csv(squarespot_rf_cat_rec, here::here("ML_inputs", "recruitment_data", "squarespot_rf_recdevs.csv"), row.names = FALSE)
write.csv(vermilion_rf_sca_cat_rec, here::here("ML_inputs", "recruitment_data", "vermilion_rf_sca_recdevs.csv"), row.names = FALSE)
write.csv(vermilion_rf_nca_cat_rec, here::here("ML_inputs", "recruitment_data", "vermilion_rf_nca_recdevs.csv"), row.names = FALSE)
write.csv(vermilion_rf_or_cat_rec, here::here("ML_inputs", "recruitment_data", "vermilion_rf_or_recdevs.csv"), row.names = FALSE)
write.csv(vermilion_rf_wa_cat_rec, here::here("ML_inputs", "recruitment_data", "vermilion_rf_wa_recdevs.csv"), row.names = FALSE)
# write.csv(copper_rf_wa_cat_rec, here::here("ML_inputs", "recruitment_data", "copper_rf_wa_recdevs.csv"), row.names = FALSE)
write.csv(quillback_rf_ca_cat_rec, here::here("ML_inputs", "recruitment_data", "quillback_rf_ca_recdevs.csv"), row.names = FALSE)
write.csv(rougheye_rf_cat_rec, here::here("ML_inputs", "recruitment_data", "rougheye_rf_recdevs.csv"), row.names = FALSE)
write.csv(yellowtail_rf_cat_rec, here::here("ML_inputs", "recruitment_data", "yellowtail_rf_recdevs.csv"), row.names = FALSE)


#### Plot recruitment data ####
rougheye_2025_recdevs %>% 
  left_join(dplyr::select(rougheye_rf_cat_rec, Yr, one_sd_boom_bust), by = "Yr") %>% 
  mutate(one_sd_boom_bust = ifelse(is.na(one_sd_boom_bust), "Outside bias adjustment", one_sd_boom_bust)) %>% 
  mutate(one_sd_boom_bust = factor(one_sd_boom_bust, 
                                       levels = c("Outside bias adjustment", "Bust", "Average", "Boom")))-> rougheye_recdevs_for_plot


recr_cat_colors <- c("Bust" = "#762a83",
                     "Average" = "#fee08b",
                     "Boom" = "#006d2c",
                     "Outside bias adjustment" = "gray90")

ggplot(rougheye_recdevs_for_plot, aes(x = Yr, y = Value, color = one_sd_boom_bust)) +
  # geom_line(color = "black", lty = 2) +
  geom_point(size = 3.5) +
  scale_color_manual(values = recr_cat_colors, name = "Boom/Bust Recruitment") +
  scale_x_continuous(breaks = seq(1900, 2020, by = 20)) +
  theme_minimal() +
  # ggtitle(stock_name) +
  ylab("Log Recruitment Deviation") +
  xlab("Year") +
  theme(legend.position = c(0.2, 0.8))
  




