#' ---
#' title: 'Understanding puma home range characteristics'
#' author: Bernardo Niebuhr
#' output:
#'   github_document: default
#'   html_document: default
#' ---

# Load packages
if(!require(install.load)) install.packages("install.load"); library(install.load)

install.load::install_load("ezknitr", "knitr")
install.load::install_load("tidyverse", "lubridate", "purrr", "ggpubr", "GGally")
install.load::install_load("sf", "raster", "stars", "tmap")
install.load::install_load("amt", "ctmm")
install.load::install_load("ggeffects")

# load data
library(PardasIPC)

# Print options for this document
options(width = 165)
opts_knit$set(root_dir = rprojroot::find_rstudio_root_file()) # root project folder
opts_chunk$set(eval = F, error = F, message = F, warning = F, cache = F, echo = T, results = T)

# --------------- label=setup
# Set up 

# Clean everything before beginning
rm(list = ls())

# --------------- label=load_data_and_functions

# Load data
load("data/movement_data_pumas_dispersal_behavior.rda")
movement_data_pumas_behavior <- pumas_behav
movement_data_pumas_behavior %>% print(width = Inf)
str(movement_data_pumas_behavior)

# load maps
file_maps <- list.files("spatial_data/analysis/", pattern = ".tif", full.names = T)

maps_raster <- raster::stack(file_maps)
land_use_raster <- raster::raster(file_maps[4]) # just for plotting
maps_rast <- terra::rast(file_maps)

# did not work
# maps <- stars::read_stars(file_maps, proxy = T)
# maps
# str(maps)
# plot(maps[3])
# land_use <- stars::read_stars(file_maps[3], proxy = F)

# crs to use
# crs_use <- maps_raster@crs
# crs_use_rast <- terra::crs(maps_rast)
crs_use <- sf::st_crs(maps_rast)
# crs_use <- sp::CRS("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m")

# transform data into track object
mov_track <- movement_data_pumas_behavior %>% 
  amt::make_track(X, Y, timestamp, all_cols = T, crs = 4326) %>% 
  dplyr::mutate(phase = relevel(as.factor(phase), ref = "pre-dispersal"),
                selected_model = as.factor(selected_model))
str(mov_track)

# get only residents, transform coordinates
mov_track_res <- mov_track %>% 
  dplyr::filter(phase != "dispersal") %>% 
  dplyr::mutate(name_phase = paste0(id, "_", str_sub(phase, 1, 3)),
                name_phase = ifelse(grepl("pos", name_phase), paste0(name_phase, "t"),
                                    name_phase)) %>% 
  amt::transform_coords(crs_to = crs_use)
mov_track_res$name_phase

# ---------------
# one individual, test

# select 1 individual
mov_track_1 <- mov_track_res %>% 
  dplyr::filter(id == "Jussara")

# calculate MCP and KDE
mcp1 <- amt::hr_mcp(mov_track_1, levels = c(0.5, 0.95))
kde1 <- amt::hr_kde(mov_track_1, levels = c(0.5, 0.95))

# plot
plot(kde1)
plot(mcp1, add.relocations = F, add = T, border = "red")

# area
hr_area(mcp1, units = T)
hr_area(kde1)

# isopleths
hr_isopleths(mcp1)

# overlap
hr_overlap(mcp1, kde1)

# ---------------
# several estimators, all individuals 

#----------------
# plot variograms
exportFIG <- T
if(exportFIG) {
  pdf('output/04_home_range/variograms.pdf', width = 10, height = 10)
}

# ctmm, make variogram for all individuals
inds <- unique(mov_track_res$name_phase)
for(i in 1:length(inds)) {
  mov_telemetry <- mov_track_res %>% 
    dplyr::filter(name_phase == inds[i]) %>% 
    amt::as_telemetry()
  
  # build the variogram
  SVF <- variogram(mov_telemetry)
  level <- c(0.5,0.95) # 50% and 95% CIs
  # xlim <- c(0, 5 %#% "day") # 0-5 day window
  xlim <- c(0, 12 %#% "hour") # 0-12 hour window
  
  # plot
  par(mfrow = c(2,2))
  plot(SVF, xlim = xlim, level=level)
  title(inds[i])
  plot(SVF, fraction = 0.75, level=level) # 75% of the monitoring period
  plot(mov_telemetry, type = 'o', pch = 19, col = 2)
  par(mfrow = c(1,1))
}

if(exportFIG) dev.off()

# test akde
# M_IID <- ctmm.fit(mov_telemetry)
# M_IID$AICc
# m_ouf <- ctmm.guess(mov_telemetry, interactive = FALSE)
# M_OUF <- ctmm.fit(mov_telemetry, m_ouf) # autocorrelation in positions and speed
# summary(M_OUF)

# from visual interpretation of variograms
# individuals that did not reach stability - all pre-
to_remove <- c("Nick_pre", "Mineiro_pre", "Piloto_pre", "Tupa_pre", 
               "Pepira_pre", "Kurupi_pre", "Zorro_pre")

#----------------
# estimate HR

# nest
mov_track_all <- mov_track_res %>%
  dplyr::filter(!(name_phase %in% to_remove)) |> 
  # dplyr::select(id,)
  tidyr::nest(-c(id, name_phase))

# estimate HR
hr1_all <- mov_track_all %>%
  dplyr::mutate(hr_mcp = map(data, hr_mcp),
                hr_kde = map(data, hr_kde),
                hr_locoh = map(data, ~ hr_locoh(., n = ceiling(sqrt(nrow(.))))),
                hr_akde_iid = map(data, ~ hr_akde(., fit_ctmm(., "iid"))),
                hr_akde_ou = map(data, ~ hr_akde(., fit_ctmm(., "ou"))),
                hr_akde_ouf = map(data, ~ hr_akde(., fit_ctmm(., "ouf"))),
                hr_akde_auto = map(data, ~ hr_akde(., fit_ctmm(., "auto"))))

# save
saveRDS(hr1_all, "output/04_home_range/hr_calculated.rds")
# hr1_all <- readRDS("output/04_home_range/hr_calculated.rds")

# understand the fits

# calculate area and extract AICc and model parameters for akde
hr1_all_area <- hr1_all %>% 
  tidyr::pivot_longer(hr_mcp:hr_akde_auto, names_to = "estimator",
                      values_to = "hr") %>% 
  dplyr::mutate(area = map(hr, hr_area),
                area = map(area, ~ .x %>% dplyr::mutate(area = area/1e6)),
                isopleth_95 = map(hr, hr_isopleths),
                aicc = map(hr, ~ try(.$model$AICc)),
                mov_parms = map(hr, ~ .x %>% .$model %>% summary() %>% .[[3]]))
                # tau_position_days = = map(hr, ~ try(summary(.$model)$CI[2,])),
                # tau_velocity_minutes = map(hr, ~ try(summary(.$model)$CI[3,])),
                # speed_kmday = map(hr, ~ try(summary(.$model)$CI[4,])),)
# land_use = map(hr, ~ raster::extract(maps[["land_use_FBDS_8_30m"]], hr_isopleths(.))),
# road_density = map(hr, ~ raster::extract(maps[["road_density_km_per100km2"]], hr_isopleths(.))),
# dist_road = map(hr, ~ raster::extract(maps[["dist_road_m"]], hr_isopleths(.))),
# dist_hidro = map(hr, ~ raster::extract(maps[["dist_hidro_m"]], hr_isopleths(.))))

# plot

# plot home range size vs estimator, for each individual
hr1_all_area %>% 
  tidyr::unnest(area) %>% 
  dplyr::filter(what == "estimate") |> 
  ggplot(aes(estimator, area)) +
  geom_point() + 
  facet_wrap(~id, scales = "free")

# plot home range size vs individual, for each estimator
hr1_all_area %>% 
  tidyr::unnest(area) %>% 
  dplyr::filter(what == "estimate") |> 
  ggplot(aes(fct_reorder(id, area), area, color = name_phase)) +
  geom_point() +
  facet_wrap(~estimator, scales = "free") +
  theme(axis.text.x=element_blank())

# # plot movement parameters for akde
# hr1_all_area %>% 
#   dplyr::filter(estimator == "hr_akde_auto") %>% 
#   tidyr::unnest(c(area, aicc:speed)) %>% 
#   dplyr::slice(seq(2, 36, by = 3)) %>% 
#   ggplot() +
#   geom_point(aes(fct_reorder(id, area), speed, color = id))
# change y to tau_position, tau_velocity

# polygons for akde
hr_akde_polygons <- hr1_all_area %>% 
  dplyr::filter(estimator == "hr_akde_auto") %>% 
  pull(isopleth_95) %>% 
  dplyr::bind_rows() %>% 
  dplyr::bind_cols(
    hr1_all_area %>% 
      dplyr::filter(estimator == "hr_akde_auto") %>% 
      dplyr::select(1:2) |> 
      dplyr::slice(rep(1:n(), each = 3))
  ) %>% 
  dplyr::mutate(area_km2 = area/1e6) |> 
  dplyr::filter(what == "estimate")

# original points
hr_akde_points <- mov_track_res %>% 
  dplyr::filter(!(name_phase %in% to_remove)) %>%
  sf::st_as_sf(coords = c("x_", "y_"), remove = F, crs = crs_use)

# plot all individuals and hr polygons in a single map  
hr1_all_area %>% 
  dplyr::filter(estimator == "hr_akde_auto") %>% 
  # dplyr::filter()
  tidyr::unnest(data) %>% 
  ggplot() +
  geom_point(aes(x_, y_, color = id)) +
  geom_sf(data = hr_akde_polygons, fill = NA, inherit_aes = F)

# plot each individual and polygon
hr_akde_points$name <- hr_akde_points$id
hr_akde_points <- hr_akde_points %>% 
  dplyr::filter(!(name %in% c("Kurupi", "Zorro")))
tm_shape(hr_akde_points) +
  tm_dots() +
  tm_facets(by = "id") +
  tm_shape(hr_akde_polygons) +
  tm_polygons(alpha = 0, ) +
  tm_facets(by = "id")

# plot all individuals with land use raster

# colors for raster
cols_rast <- c("#80b1d3", "#fccde5", "#fb8072", "#b3de69", "#8dd3c7", "#bebada", "#fdb462", "#ffffb3")
# classes for raster
labs <- c("water", "other anthropogenic", "urban", "forest",
          "non forest natural", "forestry", "sugarcane", "pasture")
# levels for raster
levs <- c(1, 2, 3, 4, 5, 6, 7, 8)

# individuals names
ind_names <- hr_akde_polygons$id

# list of plots
plots_hr <- list()
# for each individual, do
for(i in 1:length(ind_names)) {
  
  # polygons for this individual
  pol <- hr_akde_polygons %>% 
    dplyr::filter(id == ind_names[i])
  
  # points for this individual
  pts <- hr_akde_points %>% 
    dplyr::filter(id == ind_names[i])
  
  # cut land use raster
  ext <- pts %>% 
    sf::st_buffer(dist = 3000) %>%
    extent()
  lu <- raster::crop(land_use_raster, ext)
  lu[] <- factor(lu[], levels = levs, labels = labs)
  
  # plot and save in the list
  plots_hr[[i]] <- tm_shape(lu) +
    tm_raster(palette = cols_rast, legend.show = F) +
    tm_shape(pts) +
    tm_dots(size = 0.1) +
    tm_shape(pol) +
    tm_polygons(alpha = 0, lwd = 2) +
    tm_layout(title = ind_names[i], title.bg.color = "white", 
              title.bg.alpha = 0.5) +
    tm_scale_bar(breaks = c(0, 3, 6), text.size = 0.8, position = c("left", "bottom"))
  
}

# combine plots
combined_plots_hr <- tmap_arrange(plots_hr)
combined_plots_hr

# save plots
tmap_save(combined_plots_hr, filename = "output/04_home_range/plot_akde_raw.png", width = 25, height = 25,
          units = "cm", dpi = 300)

#----------------
# do plots for each estimator

estims <- c("hr_mcp", "hr_kde", "hr_akde_auto")
for(est in estims) {
  
  # polygons for the estimator
  hr_polygons <- hr1_all_area %>% 
    dplyr::filter(estimator == est) %>% 
    dplyr::mutate(
      isopleth_95 = purrr::map2(isopleth_95, id, function(a, b) { a$id <- b; return(a)}),
      isopleth_95_pol = purrr::map(isopleth_95, ~ st_cast(., "POLYGON"))) %>% 
    dplyr::pull(isopleth_95_pol) %>% 
    dplyr::bind_rows() |> 
    dplyr::filter(what == "estimate")
    # data.table::rbindlist(.) %>% 
    # sf::st_as_sf()
  
  # original points
  hr_points <- mov_track_res %>% 
    dplyr::filter(!(name_phase %in% to_remove)) %>%
    sf::st_as_sf(coords = c("x_", "y_"), remove = F, crs = crs_use)
  
  # lines
  hr_lines <- hr_points %>% 
    group_by(id) %>% 
    summarize(do_union=FALSE) %>% 
    sf::st_cast("LINESTRING")
  
  # plot each individual and polygon
  g_norast <- tm_shape(hr_points) +
    tm_dots(size = 0.1) +
    tm_facets(by = "id") +
    tm_shape(hr_lines) +
    tm_lines() +
    tm_facets(by = "id") +
    tm_shape(hr_polygons) +
    tm_polygons(alpha = 0, lwd = 2) +
    tm_facets(by = "id")
  g_norast
  
  # save plot
  tmap_save(g_norast, filename = paste0("output/04_home_range/plot_no_rast_", est, "_raw.png"),
            width = 15, height = 15, units = "cm", dpi = 300)
  
  # plot all individuals with land use raster
  
  # colors for raster
  cols_rast <- c("#80b1d3", "#fccde5", "#fb8072", "#b3de69", "#8dd3c7", "#bebada", "#fdb462", "#ffffb3")
  # classes for raster
  labs <- c("water", "other anthropogenic", "urban", "forest",
            "non forest natural", "forestry", "sugarcane", "pasture")
  # levels for raster
  levs <- c(1, 2, 3, 4, 5, 6, 7, 8)
  
  # individuals names
  ind_names <- unique(sort(hr_polygons$id))
  
  # list of plots
  plots_hr <- list()
  # for each individual, do
  for(i in 1:length(ind_names)) {
    
    # polygons for this individual
    pol <- hr_polygons %>% 
      dplyr::filter(id == ind_names[i])
    
    # points for this individual
    pts <- hr_points %>% 
      dplyr::filter(id == ind_names[i])
    
    # lines for the trajectory
    lins <- pts %>% 
      group_by(id) %>% 
      summarize(do_union=FALSE) %>% 
      sf::st_cast("LINESTRING")
    
    # cut land use raster
    ext <- pts %>% 
      sf::st_buffer(dist = 3000) %>%
      extent()
    lu <- raster::crop(land_use_raster, ext)
    lu[] <- factor(lu[], levels = levs, labels = labs)
    
    # plot and save in the list
    plots_hr[[i]] <- tm_shape(lu) +
      tm_raster(palette = cols_rast, legend.show = F) +
      tm_shape(pts) +
      tm_dots(size = 0.05, alpha = 0.5) +
      tm_shape(lins) +
      tm_lines(alpha = 0.5) +
      tm_shape(pol) +
      tm_polygons(alpha = 0, lwd = 3) +
      tm_layout(title = ind_names[i], title.bg.color = "white", 
                title.bg.alpha = 0.5) +
      tm_scale_bar(breaks = c(0, 3, 6), text.size = 0.8, position = c("left", "bottom"))
    
  }
  
  # combine plots
  combined_plots_hr <- tmap_arrange(plots_hr)
  combined_plots_hr
  
  # save plots
  tmap_save(combined_plots_hr, filename = paste0("output/04_home_range/plot_raster_", est, "_raw.png"), 
            width = 25, height = 25, units = "cm", dpi = 300)
}

# export as shp
outfolder <- "output/04_home_range/shp/"

shp_export <- hr1_all_area %>% 
  dplyr::select(id, name_phase, estimator, hr, isopleth_95)
for(i in 1:nrow(shp_export)) {
  out_name <- paste0(outfolder, shp_export$name_phase[i], "_",
                     shp_export$estimator[i], ".shp")
  sf::st_write(shp_export$isopleth_95[[i]], out_name)
}

# export as gpkg
outfolder <- "output/04_home_range/gpkg/"

for(i in 1:nrow(shp_export)) {
  out_name <- paste0(outfolder, shp_export$name_phase[i], "_",
                     shp_export$estimator[i], ".gpkg")
  sf::st_write(shp_export$isopleth_95[[i]], out_name)
}

#----------
# cut water from polygons when the polygons go in the middle of a reservatoir

# read polygons
water_pols <- sf::st_read("spatial_data/water_remove_HR.shp") |> 
  sf::st_union() |> 
  sf::st_transform(crs = crs_use)
plot(water_pols)

# test for one polygon
# hr1_all_area$isopleth_95[[7]] |> plot()
# hr1_all_area$isopleth_95[[7]] |> 
#   sf::st_difference(water_pols) |> 
#   plot()

# cut reservatoirs for all polygons - nowater = nw
hr1_all_area_nw <- hr1_all_area |> 
  dplyr::mutate(isopleth_95_nw = purrr::map(isopleth_95, ~ sf::st_difference(.x, water_pols)),
                area_nw_km2 = purrr::map(isopleth_95_nw, ~ sf::st_area(.x)/1e6))

hr1_all_area_nw |> 
  dplyr::filter(id == c("Nick"), estimator == "hr_akde_auto") |> 
  pull(isopleth_95_nw) |> first() |> 
  plot()

hr1_all_area_nw$area_nw_km2[[1]]

#----------------
# re-do plots for each estimator

estims <- c("hr_mcp", "hr_kde", "hr_akde_auto")
for(est in estims) {
  
  # polygons for the estimator
  hr_polygons <- hr1_all_area_nw %>% 
    dplyr::filter(estimator == est) %>% 
    dplyr::mutate(
      isopleth_95_nw = purrr::map2(isopleth_95_nw, id, function(a, b) { a$id <- b; return(a)}),
      isopleth_95_pol = isopleth_95_nw) |> #purrr::map(isopleth_95_nw, ~ st_cast(., "POLYGON"))) %>% 
    dplyr::pull(isopleth_95_pol) %>% 
    dplyr::bind_rows() |> 
    dplyr::filter(what == "estimate")
  # data.table::rbindlist(.) %>% 
  # sf::st_as_sf()
  
  # original points
  hr_points <- mov_track_res %>% 
    dplyr::filter(!(name_phase %in% to_remove)) %>%
    sf::st_as_sf(coords = c("x_", "y_"), remove = F, crs = crs_use)
  
  # lines
  hr_lines <- hr_points %>% 
    group_by(id) %>% 
    summarize(do_union=FALSE) %>% 
    sf::st_cast("LINESTRING")
  
  # plot each individual and polygon
  g_norast <- tm_shape(hr_points) +
    tm_dots(size = 0.1) +
    tm_facets(by = "id") +
    tm_shape(hr_lines) +
    tm_lines() +
    tm_facets(by = "id") +
    tm_shape(hr_polygons) +
    tm_polygons(alpha = 0, lwd = 2) +
    tm_facets(by = "id")
  g_norast
  
  # save plot
  tmap_save(g_norast, filename = paste0("output/04_home_range/plot_no_rast_", est, "_nw.png"),
            width = 15, height = 15, units = "cm", dpi = 300)
  
  # plot all individuals with land use raster
  
  # colors for raster
  cols_rast <- c("#80b1d3", "#fccde5", "#fb8072", "#b3de69", "#8dd3c7", "#bebada", "#fdb462", "#ffffb3")
  # classes for raster
  labs <- c("water", "other anthropogenic", "urban", "forest",
            "non forest natural", "forestry", "sugarcane", "pasture")
  # levels for raster
  levs <- c(1, 2, 3, 4, 5, 6, 7, 8)
  
  # individuals names
  ind_names <- unique(sort(hr_polygons$id))
  
  # list of plots
  plots_hr <- list()
  # for each individual, do
  for(i in 1:length(ind_names)) {
    
    # polygons for this individual
    pol <- hr_polygons %>% 
      dplyr::filter(id == ind_names[i])
    
    # points for this individual
    pts <- hr_points %>% 
      dplyr::filter(id == ind_names[i])
    
    # lines for the trajectory
    lins <- pts %>% 
      group_by(id) %>% 
      summarize(do_union=FALSE) %>% 
      sf::st_cast("LINESTRING")
    
    # cut land use raster
    ext <- pts %>% 
      sf::st_buffer(dist = 3000) %>%
      extent()
    lu <- raster::crop(land_use_raster, ext)
    lu[] <- factor(lu[], levels = levs, labels = labs)
    
    # plot and save in the list
    plots_hr[[i]] <- tm_shape(lu) +
      tm_raster(palette = cols_rast, legend.show = F) +
      tm_shape(pts) +
      tm_dots(size = 0.05, alpha = 0.5) +
      tm_shape(lins) +
      tm_lines(alpha = 0.5) +
      tm_shape(pol) +
      tm_polygons(alpha = 0, lwd = 3) +
      tm_layout(title = ind_names[i], title.bg.color = "white", 
                title.bg.alpha = 0.5) +
      tm_scale_bar(breaks = c(0, 3, 6), text.size = 0.8, position = c("left", "bottom"))
    
  }
  
  # combine plots
  combined_plots_hr <- tmap_arrange(plots_hr)
  combined_plots_hr
  
  # save plots
  tmap_save(combined_plots_hr, filename = paste0("output/04_home_range/plot_raster_", est, "_nw.png"), 
            width = 25, height = 25, units = "cm", dpi = 300)
}

# export as shp
outfolder <- "output/04_home_range/shp_nw/"

shp_export <- hr1_all_area_nw %>% 
  dplyr::select(id, name_phase, estimator, hr, isopleth_95_nw)
for(i in 1:nrow(shp_export)) {
  out_name <- paste0(outfolder, shp_export$name_phase[i], "_",
                     shp_export$estimator[i], "_nw.shp")
  sf::st_write(shp_export$isopleth_95_nw[[i]], out_name)
}

# export as gpkg
outfolder <- "output/04_home_range/gpkg_nw/"

for(i in 1:nrow(shp_export)) {
  out_name <- paste0(outfolder, shp_export$name_phase[i], "_",
                     shp_export$estimator[i], "_nw.gpkg")
  sf::st_write(shp_export$isopleth_95_nw[[i]], out_name, delete_dsn = TRUE)
}

# Get estimates and extract covariates
hr1_all_area_nw_env <- hr1_all_area_nw |> 
  dplyr::mutate(iso_95_nw_est = purrr::map(isopleth_95_nw, ~ . |> dplyr::filter(what == "estimate")),
                land_use = map(iso_95_nw_est, ~ terra::extract(maps_rast[["FBDS_30m_map_all_states_sugarcane7_pasture8"]], terra::vect(.))),
                road_density = map(iso_95_nw_est, ~ terra::extract(maps_rast[["road_density_km_per100km2"]], terra::vect(.))),
                road_dist = map(iso_95_nw_est, ~ terra::extract(maps_rast[["dist_roads_dnit2014"]], terra::vect(.))),
                urban_density = map(iso_95_nw_est, ~ terra::extract(maps_rast[["urban_density_n_pixels_per10km2"]], terra::vect(.))),
                urban_dist = map(iso_95_nw_est, ~ terra::extract(maps_rast[["dist_cities_FBDS_30m"]], terra::vect(.))),
                dist_hidro = map(iso_95_nw_est, ~ terra::extract(maps_rast[["dist_hidro_5m"]], terra::vect(.))))

# save
# saveRDS(hr1_all_area_nw_env, file = "output/04_home_range/hr_calculated_nw_env.rds")
hr1_all_area_nw_env <- readr::read_rds("output/04_home_range/hr_calculated_nw_env.rds")

# table with HR parms
# name, sex, akde - model selected, size, tau_pos, tau_vel, speed, kde - size, mcp
labs <- c("water", "other anthropogenic", "urban", "forest",
          "non forest natural", "forestry", "sugarcane", "pasture")

hr_df <- hr1_all_area_nw_env |> 
  dplyr::mutate(best_model = map(hr, ~ try(summary(.x$model)$name)),
                best_model = map(best_model, ~ ifelse(class(.x) == "character", .x, NA)),
                area_nw_km2_est = map(area_nw_km2, ~ median(.x)),
                land_use = map(land_use, ~ factor(.x[,2], levels = 1:8, labels = labs)),
                lu_tab = map(land_use, ~ table(.x)/sum(table(.x))),
                prop_forest = map(lu_tab, ~ 100*.x[which(names(.x) == "forest")]),
                prop_nonforest = map(lu_tab, ~ 100*.x[which(names(.x) == "non forest natural")]),
                prop_forestry = map(lu_tab, ~ 100*.x[which(names(.x) == "forestry")]),
                prop_sugarcane = map(lu_tab, ~ 100*.x[which(names(.x) == "sugarcane")]),
                prop_pasture = map(lu_tab, ~ 100*.x[which(names(.x) == "pasture")]),
                road_dens_avg = map(road_density, ~ mean(.x[,2], na.rm = T)),
                road_dist_avg_km = map(road_dist, ~ mean(.x[,2], na.rm = T)/1e3),
                urban_dens_avg = map(urban_density, ~ mean(.x[,2], na.rm = T)),
                urban_dist_avg_km = map(urban_dist, ~ mean(.x[,2], na.rm = T)/1e3),
                hidro_dist_avg_km = map(dist_hidro, ~ mean(.x[,2], na.rm = T)/1e3)) |> 
  dplyr::select(id, name_phase, estimator, best_model, area_km2 = area_nw_km2_est, starts_with("prop"), contains("avg")) |> 
  tidyr::unnest(best_model:hidro_dist_avg_km) |> 
  dplyr::mutate(area_km2 = as.numeric(area_km2)) |> 
  dplyr::left_join(pumas_metadata |> 
                     dplyr::rename(id = name), by = "id")
hr_df

# explore

# forest, interesting
hr_df |> 
  dplyr::filter(id != "Jussara") |> 
  ggplot(aes(prop_forest, as.numeric(area_km2))) +
  geom_point() +
  geom_smooth(method = "gam") +
  facet_wrap(~estimator)

hr_forest_bytype <- hr_df |>
  dplyr::filter(id != "Jussara", estimator == "hr_akde_ouf") |> 
  ggplot(aes(prop_forest, as.numeric(area_km2))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = Gamma(link = "log")), col = "black") +
  facet_wrap(~dispersal_translocated) +
  labs(x = bquote("Proportion of forest (%)"),
       y = bquote("Home range size ( km"^2~")")) +
  theme_minimal()
hr_forest_bytype

# sugarcane, interesting
hr_df |> 
  dplyr::filter(id != "Jussara") |> 
  ggplot(aes(prop_sugarcane, as.numeric(area_km2))) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_wrap(~estimator)

# pasture, no
hr_df |> 
  dplyr::filter(id != "Jussara") |> 
  ggplot(aes(prop_pasture, as.numeric(area_km2))) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_wrap(~estimator)

# forestry - interesting
hr_df |> 
  dplyr::filter(id != "Jussara") |> 
  ggplot(aes(prop_forestry, as.numeric(area_km2))) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_wrap(~estimator)

# non-forest - no
hr_df |> 
  dplyr::filter(id != "Jussara") |> 
  ggplot(aes(prop_nonforest, as.numeric(area_km2))) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_wrap(~estimator)

# urb dist - interesting (not density)
hr_df |> 
  dplyr::filter(id != "Jussara") |> 
  ggplot(aes(urban_dist_avg_km, as.numeric(area_km2))) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_wrap(~estimator)

# road dens - interesting (not dist)
hr_df |> 
  dplyr::filter(id != "Jussara") |> 
  ggplot(aes(road_dens_avg, as.numeric(area_km2))) +
  geom_point() +
  geom_smooth(method = "glm") +
  facet_wrap(~estimator)

hr_roads_bytype <- hr_df |> 
  dplyr::filter(id != "Jussara", estimator == "hr_akde_ouf") |> 
  ggplot(aes(road_dens_avg, as.numeric(area_km2))) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "glm", method.args = list(family = Gamma(link = "log")), , col = "black") +
  facet_wrap(~dispersal_translocated) +
  labs(x = bquote("Density of roads ( km/km"^2~")"),
       y = bquote("Home range size ( km"^2~")")) +
  theme_minimal()
hr_roads_bytype

# hidro - no
hr_df |> 
  dplyr::filter(id != "Jussara") |> 
  ggplot(aes(hidro_dist_avg_km, as.numeric(area_km2))) +
  geom_point() +
  geom_smooth(method = "glm", col = "black") +
  facet_wrap(~estimator) + 
  theme_minimal()

# check_correlation
colmns <- c("prop|avg")
colmns <- grep(colmns, colnames(hr_df)) 
corr <- GGally::ggpairs(hr_df, columns = colmns)
corr
ggsave("correlations_among_predictors.png", plot = corr, path = "output/04_home_range//", 
       width = 25, height = 25, units = "cm", dpi = 300)

# combine figures

#xlim(0, NA)

g_hr3_bytype <- ggpubr::ggarrange(hr_forest_bytype, hr_roads_bytype, ncol = 1, labels = c("A", "B"))
g_hr3_bytype
ggsave("home_range_forest_roads_bytype.png", plot = g_hr3_bytype, path = "output/04_home_range/", 
       width = 20, height = 20, units = "cm", dpi = 300)

# hidro and urban density correlated 0.756
# non forest correlated with hidro and (0.9) and urban density (0.76)
# forest and sugarcane slightly correlated but ok -0.52

hr_df2 <- hr_df |> 
  dplyr::filter(id != "Jussara")

# no power to test difference between sex - only three females
# m_hr0 <- glm(area_km2 ~ sex, 
#             family = Gamma(link = "log"), 
#             data = filter(hr_df2, estimator == "hr_akde_auto")) 

m_hr1 <- glm(area_km2 ~ scale(prop_forest) +
              scale(road_dens_avg) + scale(urban_dist_avg_km), 
            family = Gamma(link = "log"), 
            data = filter(hr_df2, estimator == "hr_akde_auto"))

summary(m_hr1)
car::vif(m_hr1)
ggpredict(m_hr1, terms = "prop_forest") |> plot()

# with differences between individuals
m_hr1.5 <- glm(area_km2 ~ scale(prop_forest) +
                 scale(road_dens_avg) + scale(urban_dist_avg_km) +
                 dispersal_translocated - 1, 
             family = Gamma(link = "log"), 
             data = filter(hr_df2, estimator == "hr_akde_auto"))

summary(m_hr1.5)
car::vif(m_hr1.5)
ggpredict(m_hr1.5, terms = c("prop_forest", "dispersal_translocated")) |> plot()

# with differences between individuals
m_hr1.6 <- glm(area_km2 ~ scale(prop_forest):dispersal_translocated +
                 scale(road_dens_avg):dispersal_translocated + scale(urban_dist_avg_km):dispersal_translocated +
                 dispersal_translocated - 1, 
               family = Gamma(link = "log"), 
               data = filter(hr_df2, estimator == "hr_akde_auto"))

summary(m_hr1.6)
car::vif(m_hr1.6)
ggpredict(m_hr1.6, terms = c("prop_forest", "dispersal_translocated")) |> plot()
ggpredict(m_hr1.6, terms = c("road_dens_avg", "dispersal_translocated")) |> plot()

m_hr1.6 <- glm(area_km2 ~ scale(prop_forest) +
                 scale(road_dens_avg):dispersal_translocated + scale(urban_dist_avg_km) +
                 dispersal_translocated - 1, 
               family = Gamma(link = "log"), 
               data = filter(hr_df2, estimator == "hr_akde_auto"))

summary(m_hr1.6)
car::vif(m_hr1.6)
ggpredict(m_hr1.6, terms = c("prop_forest", "dispersal_translocated")) |> plot()
ggpredict(m_hr1.6, terms = c("road_dens_avg", "dispersal_translocated")) |> plot()

m_hr2 <- glm(area_km2 ~ scale(prop_nonforest) +
              scale(road_dens_avg) + scale(urban_dist_avg_km), 
            family = Gamma(link = "log"), 
            data = filter(hr_df2, estimator == "hr_akde_auto"))

summary(m_hr2)
car::vif(m_hr2)
ggpredict(m_hr2, terms = "prop_nonforest") |> plot()

m_hr3 <- glm(area_km2 ~ scale(prop_forestry) + 
               scale(road_dens_avg) + scale(urban_dist_avg_km), 
             family = Gamma(link = "log"), 
             data = filter(hr_df2, estimator == "hr_akde_auto"))

summary(m_hr3)
car::vif(m_hr3)
ggpredict(m_hr3, terms = "prop_forestry") |> plot()

m_hr4 <- glm(area_km2 ~ scale(prop_sugarcane) +
               scale(road_dens_avg) + scale(urban_dist_avg_km), 
             family = Gamma(link = "log"), 
             data = filter(hr_df2, estimator == "hr_akde_auto"))

summary(m_hr4)
car::vif(m_hr4)
ggpredict(m_hr4, terms = "prop_sugarcane") |> plot()

bbmle::AICctab(m_hr1, m_hr1.5, m_hr1.6, m_hr2, m_hr3, m_hr4, base = T, weights = T)

# export
bbmle::AICctab(m_hr1, m_hr2, m_hr3, m_hr4, base = T, weights = T) |> 
  as.data.frame() |> 
  mutate_if(is.numeric, round, digits = 4) %>% 
  dplyr::mutate(formula = as.character(c(formula(m_hr1), formula(m_hr2),
                                         formula(m_hr3), formula(m_hr4)))) |> 
  readr::write_csv(file = "output/04_home_range/AIC_tab_HR.csv")


# best fit
best_fit <- m_hr1

# coefficients
broom::tidy(best_fit)

coefs <- broom::tidy(best_fit) %>% 
  dplyr::mutate(significant = as.factor(if_else(p.value < 0.05, 1, 0)),
                p = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)))
coefs %>% 
  mutate_if(is.numeric, round, digits = 4) %>% 
  readr::write_csv(file = "output/04_home_range/coef_table_HR_forest_road.csv")

# prediction
p_hr <- ggeffects::ggpredict(best_fit, terms = c("prop_forest"))
g_hr <- plot(p_hr, show.title = F) +
  geom_point(aes(prop_forest, area_km2), inherit.aes = F, col = "grey50",
             data = filter(hr_df2, estimator == "hr_akde_auto")) +
  labs(x = "Proportion of forest (%)",
       y = bquote("Home range size ( km"^2~")"))
  #xlim(0, NA)

p_hr2 <- ggeffects::ggpredict(best_fit, terms = c("road_dens_avg"))
g_hr2 <- plot(p_hr2, show.title = F) +
  geom_point(aes(road_dens_avg, area_km2), inherit.aes = F, col = "grey50",
             data = filter(hr_df2, estimator == "hr_akde_auto")) +
  labs(x = bquote("Density of roads ( km/km"^2~")"),
       y = bquote("Home range size ( km"^2~")"))
#xlim(0, NA)

g_hr3 <- ggpubr::ggarrange(g_hr, g_hr2, ncol = 2, labels = c("A", "B"))
g_hr3
ggsave("home_range_forest_roads.png", plot = g_hr3, path = "output/04_home_range/", 
       width = 20, height = 10, units = "cm", dpi = 300)

ggsave("home_range_forest.png", plot = g_hr, path = "output/04_home_range/", 
       width = 15, height = 15, units = "cm", dpi = 300)

# estimators akde
hr_df |> 
  dplyr::filter(estimator == "hr_akde_auto") |> 
  pull(best_model) |> 
  unique()

# table HR
hr_sum <- hr_df |> 
  dplyr::select(name = id, sex, estimator, area_km2) |> 
  dplyr::filter(estimator %in% c("hr_mcp", "hr_kde", "hr_akde_ouf")) |> 
  tidyr::pivot_wider(names_from = "estimator", values_from = area_km2) |> 
  dplyr::rename(MCP = hr_mcp, KDE = hr_kde, AKDE = hr_akde_ouf) |> 
  dplyr::arrange(sex, name) |> 
  dplyr::mutate(sex = factor(sex, levels = c("Female", "Male"), labels = c("F", "M"))) 

# write
hr_sum |>
  mutate_if(is.numeric, round, digits = 2) %>%
  readr::write_csv("output/04_home_range/table_home_range_sizes.csv")

# hr_sum <- readr::read_csv("output/04_home_range/table_home_range_sizes.csv")

# summary
range(hr_sum$AKDE) 
mean(hr_sum$AKDE)

# plot
g_hr_sum <- hr_sum |> 
  tidyr::pivot_longer(cols = MCP:AKDE, names_to = "estimator", values_to = "hr") |> 
  dplyr::mutate(estimator = factor(estimator, levels = c("MCP", "KDE", "AKDE"))) |> 
  ggplot(aes(estimator, hr)) +
  geom_boxplot() +
  labs(x = "Home range estimator",
       y = bquote("Home range size ( km"^2~")")) +
  theme_ggeffects()
g_hr_sum

ggsave("home_range_methods.png", plot = g_hr_sum, path = "output/04_home_range/", 
       width = 15, height = 15, units = "cm", dpi = 300)


# coefficients of all models
broom::tidy(m_hr1)

coefs_all <- broom::tidy(m_hr1) %>% 
  dplyr::mutate(formula = as.character(c(formula(m_hr1)))) |> 
  dplyr::bind_rows(broom::tidy(m_hr2) |> 
                     dplyr::mutate(formula = as.character(c(formula(m_hr2))))) |> 
  dplyr::bind_rows(broom::tidy(m_hr3) |> 
                     dplyr::mutate(formula = as.character(c(formula(m_hr3))))) |> 
  dplyr::bind_rows(broom::tidy(m_hr4) |> 
                     dplyr::mutate(formula = as.character(c(formula(m_hr4))))) |> 
  dplyr::mutate(significant = as.factor(if_else(p.value < 0.05, 1, 0)),
                p = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)))

coefs_all %>% 
  mutate_if(is.numeric, round, digits = 4) %>% 
  readr::write_csv(file = "output/04_home_range/coef_table_HR_all_models.csv")