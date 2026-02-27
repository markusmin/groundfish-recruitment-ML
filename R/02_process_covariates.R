# 02_process_covariates

# this script takes the potential covariates for recruitment and processes them for modeling.

library(here)
library(tidyverse)
library(ncdf4)
library(sf)
library(viridis)
library(ggpubr)
library(paran)

#### Define functions ####

# Load ncdf data, join model outputs, truncate to 100 km offshore
# note that the near-real-time product ends in October 2022, so
# need to be careful with later processing
load_reformat_roms <- function(wcra_path, wcnrt_path,
                               var_name,
                               lon_cutoff = -127,
                               lat_min = 34.5,
                               lat_max = 49,
                               dist_shore_cutoff = 100){
  
  #### Load and join two ROMS runs
  
  ### load 1980-2010
  wcra_nc <- nc_open(wcra_path)
  wcra_nc_lon <- ncvar_get(wcra_nc, "lon_rho")
  wcra_nc_lat <- ncvar_get(wcra_nc, "lat_rho")
  wcra_nc_var <- ncvar_get(wcra_nc, var_name)
  
  # store all of this in a long format df
  # round them all to three digits - this deals with the floating-point precision mismatch.
  lat_wcra = rep(as.vector(round(wcra_nc_lat,3)),
            dim(wcra_nc_var)[3])
  lon_wcra = rep(as.vector(round(wcra_nc_lon,3)),
            dim(wcra_nc_var)[3])
  year_wcra = rep(rep(1980:2010, each = 12),
             each = length(as.vector(wcra_nc_lat)))
  month_wcra = rep(rep(1:12,length(1980:2010)),
              each = length(as.vector(wcra_nc_lat)))
  var_wcra = wcra_nc_var

  wcra_df <- as.data.frame(cbind(lat_wcra, lon_wcra, year_wcra, month_wcra, var_wcra))
  names(wcra_df) <- c("lat", "lon", "year", "month", var_name)
  
  ### load 2011-2022
  wcnrt_nc <- nc_open(wcnrt_path)
  
  wcnrt_nc_year <- ncvar_get(wcnrt_nc, "year")
  wcnrt_nc_month <- ncvar_get(wcnrt_nc, "month")
  wcnrt_nc_lon <- ncvar_get(wcnrt_nc, "lon")
  wcnrt_nc_lat <- ncvar_get(wcnrt_nc, "lat")
  wcnrt_nc_var <- ncvar_get(wcnrt_nc, var_name)
  
  # store all of this in a long format df
  lat_wcnrt <- as.vector(round(wcnrt_nc_lat,3))
  lon_wcnrt <- as.vector(round(wcnrt_nc_lon,3))
  var_wcnrt <- as.vector(wcnrt_nc_var)
  
  wcnrt_df <- data.frame(lat = rep(wcnrt_nc_lat, length(wcnrt_nc_month)),
                                     lon = rep(wcnrt_nc_lon, length(wcnrt_nc_month)),
                                     year = rep(wcnrt_nc_year, each = length(lat_wcnrt)),
                                     month = rep(wcnrt_nc_month, each = length(lat_wcnrt)),
                                     var = var_wcnrt)
  
  names(wcnrt_df) <- c("lat", "lon", "year", "month", var_name)
  
  ### join wcra and wcnrt
  rbind(wcra_df,wcnrt_df) -> roms_df
  
  #### Data processing
  
  roms_df_noNA <- roms_df[!is.na(roms_df[[var_name]]), ]
  # also drop all of the data points west of 127 W as pre-processing step to save computational time
  roms_df_noNA <- subset(roms_df_noNA, lon > lon_cutoff)
  # subset only data north of point Conception
  roms_df_noNA <- subset(roms_df_noNA, lat >= lat_min & lat <= lat_max)
  

  # Pre-processing steps: Subset only coast to 100 km offshore
  # add in distance from shore - this is now your survey domain, with covariates
  
  # Convert roms_df to sf and then to UTM zone 10
  UTM_zone_10_crs <- 32610
  st_as_sf(roms_df_noNA, coords = c("lon", "lat"), crs = 4326) -> roms_df_sf
  sf::st_transform(roms_df_sf, crs = UTM_zone_10_crs) -> roms_df_sf_proj
  # convert to km
  roms_df_sf_proj_km <- roms_df_sf_proj
  roms_df_sf_proj_km$geometry <- roms_df_sf_proj$geometry/1000
  roms_df_sf_proj_km <- st_set_crs(roms_df_sf_proj_km, UTM_zone_10_crs)
  
  # make a polygon for your ocean area and subset it
  roms_df_sf_proj_km_2000_01 <- subset(roms_df_sf_proj_km, year == 2000 & month == 1)
  # calculate distance to shore
  roms_df_sf_proj_km_2000_01$dist_shore <- as.numeric(st_distance(roms_df_sf_proj_km_2000_01, US_west_coast_proj_km))
  roms_df_sf_proj_km_2000_01 %>% 
    filter(dist_shore <= dist_shore_cutoff) -> roms_df_sf_proj_km_2000_01_trimmed
  # make a concave hull polygon
  st_concave_hull(st_union(roms_df_sf_proj_km_2000_01_trimmed), ratio = 0.01) -> polygon_trimmed
  
  # intersect this with the full data
  roms_df_sf_trimmed <- st_intersection(roms_df_sf_proj_km, polygon_trimmed)
  
  return(roms_df_sf_trimmed)
  
}

# Take trimmed ROMS output and summarize temporally

# Make sure that you don't keep any partial years
summarize_annual_roms <- function(roms_sf, year_month_end, var_name){
  
  # add a time column
  roms_sf$time <- ymd(paste0(roms_sf$year, "-", paste0(roms_sf$month), "-01"))
  
  # add a new column for year, defined by end of the year
  # to align this with recruitment, this needs to be the last six months of the
  # preceding year and the first six months of the next year (if your year ends in June)
  roms_sf$year2 <- year(roms_sf$time + months(year_month_end))
  
  # check if all months are covered in a year, drop any where the full year isn't represented
  roms_sf %>% 
    st_drop_geometry() %>% 
    group_by(year2) %>% 
    # summarise(unique_months = paste(list(unique(month)), collapse = ", ")) -> roms_months
    summarise(unique_months = list(unique(month))) %>% 
    mutate(n_months = map_int(unique_months, length)) -> roms_months
  
  subset(roms_months, n_months != 12)$year2 -> years_to_drop
  
  roms_sf <- subset(roms_sf, !(year2 %in% years_to_drop))
  
  
  # calculate annual means at each point
  roms_sf %>% 
    group_by(year2, geometry) %>% 
    summarise(mean_var = mean(.data[[var_name]])) -> roms_sf_annual_mean
  
  # convert back to df
  roms_df_annual_mean <- as.data.frame(st_drop_geometry(roms_sf_annual_mean))
  roms_df_annual_mean <- cbind(roms_df_annual_mean, st_coordinates(roms_sf_annual_mean))
  
  names(roms_df_annual_mean) <- c("year", paste0(var_name, "_mean"), "X", "Y")
  
  return(roms_df_annual_mean)
}

summarize_seasonally_roms <- function(roms_sf, start_month, end_month, var_name){
  
  roms_sf_season <- subset(roms_sf, month >= start_month & month <= end_month)
  
  # calculate seasonal means at each point
  roms_sf_season %>% 
    group_by(year, geometry) %>% 
    summarise(mean_var = mean(.data[[var_name]])) -> roms_sf_seasonal_mean
  
  # convert back to df
  roms_df_seasonal_mean <- as.data.frame(st_drop_geometry(roms_sf_seasonal_mean))
  roms_df_seasonal_mean <- cbind(roms_df_seasonal_mean, st_coordinates(roms_sf_seasonal_mean))
  
  names(roms_df_seasonal_mean) <- c("year", paste0(var_name, "_mean"), "X", "Y")
  
  return(roms_df_seasonal_mean)
  
}

summarize_annual_upwelling <- function(upwelling_df, year_month_end, var_name){
  
  # add a time column
  upwelling_df$time <- ymd(paste0(upwelling_df$year, "-", paste0(upwelling_df$month), "-01"))
  
  # add a new column for year, defined by end of the year
  # to align this with recruitment, this needs to be the last six months of the
  # preceding year and the first six months of the next year (if your year ends in June)
  upwelling_df$year2 <- year(upwelling_df$time + months(year_month_end))
  
  # calculate annual means at each point
  upwelling_df %>% 
    group_by(year2) %>% 
    summarise(mean_var = mean(.data[[var_name]])) -> upwelling_df_annual_mean
  
  names(upwelling_df_annual_mean) <- c("year", paste0(var_name, "_mean"))
  
  return(upwelling_df_annual_mean)
}

# Run EOF/PCA

run_EOF <- function(roms_summary, var_name){
  # reformat into a matrix
  roms_summary %>% 
    mutate(pos = paste0(round(X,2), "-", round(Y,2))) %>% 
    dplyr::select(-c(X,Y)) %>% 
    pivot_wider(names_from = pos, values_from = var_name) %>% 
    column_to_rownames("year") %>% 
    as.matrix() -> roms_summary_matrix
  
  # run EOF
  pca <- prcomp(roms_summary_matrix, center = TRUE, scale = TRUE)
  
  # Run Horn's parallel analysis test to determine significant PCs
  horn <- paran(roms_summary_matrix, iterations = 100)
  
  n_pcs <- horn$Retained
  
  # store the EOFs in a DF with the locations
  eofs <- cbind(subset(roms_summary, year == 2015)$X,
                Y = subset(roms_summary, year == 2015)$Y,
                pca$rotation[,1:n_pcs])
  
  colnames(eofs) <- c("X", "Y",
                          paste0("PC",1:n_pcs))
  
  # return the time series of PCs
  pc_ts_df <- cbind(as.numeric(names(pca$x[,1])),
                    pca$x[,1:n_pcs])
  
  colnames(pc_ts_df) <- c("year",
                          paste0("PC",1:n_pcs))
  
  
  return(list(
    pca = pca,
    EOFs = eofs,
    pc_ts_df = pc_ts_df))
  
}


# Plot EOF/PCA spatial patterns
plot_EOF_PCs <- function(EOFs){
  nPCs <- ncol(EOFs) - 2
  
  PC_maps <- vector(mode = "list", length = nPCs)
  
  for (i in 1:nPCs){
    PC_maps[[i]] <- west_coast_basemap_km +
      geom_point(data = EOFs, aes(x = X, y = Y, color = .data[[paste0("PC", i)]])) +
      # geom_point(data = EOFs, aes(x = X, y = Y, color = PC1)) +
      scale_color_viridis(name = paste0("PC", i, " Loadings")) +
      theme(legend.position = c(0.7, 0.7))
  }
  
  # arrange all of these
  PC_maps_combined <- ggarrange(plotlist = PC_maps,
            ncol = 2, nrow = 2)
  
  
  return(PC_maps_combined)
}


#### load the spatial data and create a basemap ####
usa_spdf <- st_read(here::here("Data", "map_files", "USA_adm0.shp"))
# load BC
CAN_spdf <- st_read(here::here("Data", "map_files", "canada", "lpr_000b16a_e.shp"))
BC_spdf <- filter(CAN_spdf, PRENAME == "British Columbia")
BC_proj <- st_transform(BC_spdf, crs = 4326)


# crop them to our desired area
US_west_coast <- sf::st_crop(usa_spdf,
                             c(xmin = -134, ymin = 30, xmax = -115.5, ymax = 48))

BC_coast <- sf::st_crop(BC_proj,
                        c(xmin = -134, ymin = 30, xmax = -115.5, ymax = 48))



# convert both shapefiles to a different projection (UTM zone 10) so that they can be plotted with the sdmTMB output
UTM_zone_10_crs <- 32610

US_west_coast_proj <- sf::st_transform(US_west_coast, crs = UTM_zone_10_crs)
BC_coast_proj <- sf::st_transform(BC_coast, crs = UTM_zone_10_crs)

# make this projection into kilometers
US_west_coast_proj_km <- st_as_sf(US_west_coast_proj$geometry/1000, crs = UTM_zone_10_crs)
BC_coast_proj_km <- st_as_sf(BC_coast_proj$geometry/1000, crs = UTM_zone_10_crs)


#### create base map for visualizing data
west_coast_basemap <- ggplot(US_west_coast) +
  geom_sf() +
  geom_sf(data = BC_coast) +
  coord_sf(ylim = c(30, 48),  xlim = c(-134, -115.5)) +
  # scale_x_continuous(breaks = c(124,125,126)) +
  ylab("Latitude")+
  xlab("Longitude")+
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.14, 0.2),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave(here::here("figures", "west_coast_basemap.png"), west_coast_basemap, height = 12, width = 10)

west_coast_basemap_km <- ggplot(US_west_coast_proj_km) +
  geom_sf() +
  geom_sf(data = BC_coast_proj_km) +
  # coord_sf(ylim = c(30, 48),  xlim = c(-134, -115.5)) +
  # scale_x_continuous(breaks = c(124,125,126)) +
  ylab("Latitude")+
  xlab("Longitude")+
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.14, 0.2),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave(here::here("figures", "west_coast_basemap_km.png"), west_coast_basemap_km, height = 12, width = 10)
  
#### Load and reformat the oceanographic data ####  

### BBV
bbv_wcra_path <- here::here("Data", "ROMS", "wcra31_monthly", "wcra31_bbv_200_monthly_1980_2010.nc")
bbv_wcnrt_path <- here::here("Data", "ROMS", "wcnrt_monthly", "wcnrt_bbv_200_monthly_201101_202210.nc")

bbv_roms_sf <- load_reformat_roms(wcra_path = bbv_wcra_path, 
                               wcnrt_path = bbv_wcnrt_path,
                               var_name = "bbv_200",
                               lon_cutoff = -127,
                               lat_min = 34.5,
                               lat_max = 49,
                               dist_shore_cutoff = 100)

bbv_roms_annual_means <- summarize_annual_roms(roms_sf = bbv_roms_sf, 
                                               year_month_end = 6, 
                                               var_name = "bbv_200")

bbv_roms_winter_means <- summarize_seasonally_roms(roms_sf = bbv_roms_sf, 
                                               start_month = 1,
                                               end_month = 3,
                                               var_name = "bbv_200")

bbv_roms_spring_means <- summarize_seasonally_roms(roms_sf = bbv_roms_sf, 
                                                    start_month = 4,
                                                    end_month = 6,
                                                    var_name = "bbv_200")

bbv_roms_summer_means <- summarize_seasonally_roms(roms_sf = bbv_roms_sf, 
                                                    start_month = 7,
                                                    end_month = 9,
                                                    var_name = "bbv_200")

bbv_roms_fall_means <- summarize_seasonally_roms(roms_sf = bbv_roms_sf, 
                                                    start_month = 10,
                                                    end_month = 12,
                                                    var_name = "bbv_200")

### ILD


ild_wcra_path <- here::here("Data", "ROMS", "wcra31_monthly", "wcra31_ild_05_monthly_1980_2010.nc")
ild_wcnrt_path <- here::here("Data", "ROMS", "wcnrt_monthly", "wcnrt_ild_05_monthly_201101_202210.nc")

ild_roms_sf <- load_reformat_roms(wcra_path = ild_wcra_path, 
                                  wcnrt_path = ild_wcnrt_path,
                                  var_name = "ild_05",
                                  lon_cutoff = -127,
                                  lat_min = 34.5,
                                  lat_max = 49,
                                  dist_shore_cutoff = 100)

ild_roms_annual_means <- summarize_annual_roms(roms_sf = ild_roms_sf, 
                                               year_month_end = 6, 
                                               var_name = "ild_05")

ild_roms_winter_means <- summarize_seasonally_roms(roms_sf = ild_roms_sf, 
                                                   start_month = 1,
                                                   end_month = 3,
                                                   var_name = "ild_05")

ild_roms_spring_means <- summarize_seasonally_roms(roms_sf = ild_roms_sf, 
                                                   start_month = 4,
                                                   end_month = 6,
                                                   var_name = "ild_05")

ild_roms_summer_means <- summarize_seasonally_roms(roms_sf = ild_roms_sf, 
                                                   start_month = 7,
                                                   end_month = 9,
                                                   var_name = "ild_05")

ild_roms_fall_means <- summarize_seasonally_roms(roms_sf = ild_roms_sf, 
                                                 start_month = 10,
                                                 end_month = 12,
                                                 var_name = "ild_05")


### SST


sst_wcra_path <- here::here("Data", "ROMS", "wcra31_monthly", "wcra31_sst_monthly_1980_2010.nc")
sst_wcnrt_path <- here::here("Data", "ROMS", "wcnrt_monthly", "wcnrt_sst_monthly_201101_202210.nc")

sst_roms_sf <- load_reformat_roms(wcra_path = sst_wcra_path, 
                                  wcnrt_path = sst_wcnrt_path,
                                  var_name = "sst",
                                  lon_cutoff = -127,
                                  lat_min = 34.5,
                                  lat_max = 49,
                                  dist_shore_cutoff = 100)

sst_roms_annual_means <- summarize_annual_roms(roms_sf = sst_roms_sf, 
                                               year_month_end = 6, 
                                               var_name = "sst")

sst_roms_winter_means <- summarize_seasonally_roms(roms_sf = sst_roms_sf, 
                                                   start_month = 1,
                                                   end_month = 3,
                                                   var_name = "sst")

sst_roms_spring_means <- summarize_seasonally_roms(roms_sf = sst_roms_sf, 
                                                   start_month = 4,
                                                   end_month = 6,
                                                   var_name = "sst")

sst_roms_summer_means <- summarize_seasonally_roms(roms_sf = sst_roms_sf, 
                                                   start_month = 7,
                                                   end_month = 9,
                                                   var_name = "sst")

sst_roms_fall_means <- summarize_seasonally_roms(roms_sf = sst_roms_sf, 
                                                 start_month = 10,
                                                 end_month = 12,
                                                 var_name = "sst")


### SSH


ssh_wcra_path <- here::here("Data", "ROMS", "wcra31_monthly", "wcra31_ssh_monthly_1980_2010.nc")
ssh_wcnrt_path <- here::here("Data", "ROMS", "wcnrt_monthly", "wcnrt_ssh_monthly_201101_202210.nc")

ssh_roms_sf <- load_reformat_roms(wcra_path = ssh_wcra_path, 
                                  wcnrt_path = ssh_wcnrt_path,
                                  var_name = "ssh",
                                  lon_cutoff = -127,
                                  lat_min = 34.5,
                                  lat_max = 49,
                                  dist_shore_cutoff = 100)

ssh_roms_annual_means <- summarize_annual_roms(roms_sf = ssh_roms_sf, 
                                               year_month_end = 6, 
                                               var_name = "ssh")

ssh_roms_winter_means <- summarize_seasonally_roms(roms_sf = ssh_roms_sf, 
                                                   start_month = 1,
                                                   end_month = 3,
                                                   var_name = "ssh")

ssh_roms_spring_means <- summarize_seasonally_roms(roms_sf = ssh_roms_sf, 
                                                   start_month = 4,
                                                   end_month = 6,
                                                   var_name = "ssh")

ssh_roms_summer_means <- summarize_seasonally_roms(roms_sf = ssh_roms_sf, 
                                                   start_month = 7,
                                                   end_month = 9,
                                                   var_name = "ssh")

ssh_roms_fall_means <- summarize_seasonally_roms(roms_sf = ssh_roms_sf, 
                                                 start_month = 10,
                                                 end_month = 12,
                                                 var_name = "ssh")



#### EOF (PCA) on ROMS outputs ####

### BBV

bbv_eof <- run_EOF(roms_summary = bbv_roms_annual_means, 
        var_name = "bbv_200_mean")

bbv_eof_loadings_plot <- plot_EOF_PCs(bbv_eof$EOFs)

ggsave(here::here("figures", "feature_extraction", "bbv_eof_loadings_plot.png"), bbv_eof_loadings_plot, height = 12, width = 10)


### ILD
ild_eof <- run_EOF(roms_summary = ild_roms_annual_means, 
                   var_name = "ild_05_mean")

ild_eof_loadings_plot <- plot_EOF_PCs(ild_eof$EOFs)

ggsave(here::here("figures", "feature_extraction", "ild_eof_loadings_plot.png"), ild_eof_loadings_plot, height = 12, width = 10)


### SSH
ssh_eof <- run_EOF(roms_summary = ssh_roms_annual_means, 
                   var_name = "ssh_mean")

ssh_eof_loadings_plot <- plot_EOF_PCs(ssh_eof$EOFs)

ggsave(here::here("figures", "feature_extraction", "ssh_eof_loadings_plot.png"), ssh_eof_loadings_plot, height = 12, width = 10)


### SST
sst_eof <- run_EOF(roms_summary = sst_roms_annual_means, 
                   var_name = "sst_mean")

sst_eof_loadings_plot <- plot_EOF_PCs(sst_eof$EOFs)

ggsave(here::here("figures", "feature_extraction", "sst_eof_loadings_plot.png"), sst_eof_loadings_plot, height = 12, width = 10)




#### Upwelling indices ####

# load upwelling data from Mary Hunsicker - summarized into regions:
# North (40.5°-47°N)
# Central (34.5–40.5°N)
# South (31–34.5°N)

roms_data_north <- read.csv(here::here("Data", "ROMS", "roms_data_north.csv"))
roms_data_central <- read.csv(here::here("Data", "ROMS", "roms_data_central.csv"))
roms_data_south <- read.csv(here::here("Data", "ROMS", "roms_data_south.csv"))

# summarize into annual indices (July-June year)
roms_data_north %>% 
  dplyr::select(year, month, CUTI) %>% 
  mutate(region = "north") -> CUTI_north

roms_data_central %>% 
  dplyr::select(year, month, CUTI) %>% 
  mutate(region = "central") -> CUTI_central

roms_data_north %>% 
  dplyr::select(year, month, BEUTI) %>% 
  mutate(region = "north") -> BEUTI_north

roms_data_central %>% 
  dplyr::select(year, month, BEUTI) %>% 
  mutate(region = "central") -> BEUTI_central


CUTI_north_annual <- summarize_annual_upwelling(upwelling_df = CUTI_north, 
                                                year_month_end = 6, 
                                                var_name = "CUTI")

CUTI_central_annual <- summarize_annual_upwelling(upwelling_df = CUTI_central, 
                                                year_month_end = 6, 
                                                var_name = "CUTI")

BEUTI_north_annual <- summarize_annual_upwelling(upwelling_df = BEUTI_north, 
                                                year_month_end = 6, 
                                                var_name = "BEUTI")

BEUTI_central_annual <- summarize_annual_upwelling(upwelling_df = BEUTI_central, 
                                                  year_month_end = 6, 
                                                  var_name = "BEUTI")

#### Combine indices and export ####
bbv_predictors <- as.data.frame(bbv_eof$pc_ts_df)
colnames(bbv_predictors) <- c("year", paste0("bbv_PC", 1:(ncol(bbv_predictors)-1)))
ild_predictors <- as.data.frame(ild_eof$pc_ts_df)
colnames(ild_predictors) <- c("year", paste0("ild_PC", 1:(ncol(ild_predictors)-1)))
ssh_predictors <- as.data.frame(ssh_eof$pc_ts_df)
colnames(ssh_predictors) <- c("year", paste0("ssh_PC", 1:(ncol(ssh_predictors)-1)))
sst_predictors <- as.data.frame(sst_eof$pc_ts_df)
colnames(sst_predictors) <- c("year", paste0("sst_PC", 1:(ncol(sst_predictors)-1)))
colnames(CUTI_north_annual) <- c("year", "CUTI_north")
colnames(CUTI_central_annual) <- c("year", "CUTI_central")
colnames(BEUTI_north_annual) <- c("year", "BEUTI_north")
colnames(BEUTI_central_annual) <- c("year", "BEUTI_central")

bbv_predictors %>% 
  left_join(ild_predictors, by = "year") %>% 
  left_join(ssh_predictors, by = "year") %>% 
  left_join(sst_predictors, by = "year") %>% 
  left_join(CUTI_north_annual, by = "year") %>% 
  left_join(CUTI_central_annual, by = "year") %>% 
  left_join(BEUTI_north_annual, by = "year") %>% 
  left_join(BEUTI_central_annual, by = "year") %>% 
  column_to_rownames("year") -> feature_df

feature_df_scaled <- as.data.frame(scale(feature_df))

# examine correlation structure
feature_df_cor_mat <- cor(feature_df_scaled)
heatmap(feature_df_cor_mat, symm = TRUE)

as.data.frame(feature_df_cor_mat) %>%
  rownames_to_column("var1") %>%
  pivot_longer(-var1, names_to = "var2", values_to = "cor") %>%
  filter(var1 < var2) %>%   # removes duplicates and diagonal
  mutate(abs_cor = abs(cor)) %>%
  arrange(desc(abs_cor)) -> feature_df_cor_table

# BEUTI_central and CUTI_central are 0.915 correlated.
# BEUTI_north   CUTI_north are 0.750 correlated.
# Let's just keep BEUTI for now.

feature_df_scaled <- dplyr::select(feature_df_scaled, -c("CUTI_north", "CUTI_central"))

# export this file
write.csv(feature_df_scaled, here::here("ML_inputs", "features_2025-01-14.csv"))







