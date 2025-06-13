#' ---
#' title: "Puma space use and dispersal in the tropics: Explore and visualize data"
#' author: Bernardo Niebuhr
#' date: ""
#' output: 
#'   html_document: default
#'   github_document: default
#' ---

# --------------- label=setup, warning=FALSE, message=FALSE, echo=TRUE

# Load packages
if(!require(install.load)) install.packages("install.load"); library(install.load)

install.load::install_load("knitr", "rprojroot")
install.load::install_load("dplyr", "ggplot2", "forcats", "leaflet")
install.load::install_load("move", "amt", "adehabitatLT")

# Print options for this document
options(width = 165)
opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) # root project folder
opts_chunk$set(error = F, message = F, warning = F, cache = F, echo = T, results = T)


#' This script aims at calculating basic statistics of movement to understand
#' when each animal was monitored, which trajectories they followed, 
#' and ro visualize their tracks.
#' 
#' # Loading data

library(PardasIPC)
data("pumas")

# --------------- label=transform_into_movement_data

#' # Transformation of data into movement data
#' 
#' Here we transformed the raw data into movement data in different formats: `move`, `amt`, 
#' and `adehabitatLT` packages. We chose to use different movement data formats since different analyses 
#' require data to be loaded thorugh distinct packages. They may also be useful to show
#' different movement-related variables.

#' ## move

# move object
mov_data_df <- as.data.frame(pumas) # move package require a data.frame, not a tibble
move_data <-move(x = mov_data_df$X, y = mov_data_df$Y,
                 time = mov_data_df$timestamp,
                 proj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"),
                 data = mov_data_df, animal = mov_data_df$id, sensor = "GPS")

move_data
# move::plot(move_data)

# If we want to change projection

# animals from Brazil - use Albers SIRGAS2000
new_projection <- "+proj=aea +lat_1=-2 +lat_2=-22 +lat_0=-12 +lon_0=-54 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
move_data_meters <- spTransform(move_data, CRSobj = sp::CRS(new_projection))

# move::plot(move_data_meters)
move_data_meters

# If we need to re-transform the move object into data.frame
# animals <- as(move_data, "data.frame")

#' ## amt
# use SIRGAS2000, UTM 22S - EPSG 31982
mov_track <- pumas %>% 
  amt::mk_track(.x = X, .y = Y, .t = timestamp, crs = sp::CRS("+init=epsg:4326"),
                id, tagID, id_number, shortname, species, sex, weight_kg, estimated_age_months, 
                study_name, dispersal_behavior) %>% 
  amt::transform_coords(crs_to = sp::CRS(new_projection))

# movement basic statistics
summ_track <- mov_track %>% 
  dplyr::group_by(id_number, id) %>% 
  dplyr::summarise(
    begin = min(t_),
    end = max(t_),
    range = diff(range(t_)),
    n = n()
  ) 

summ_track %>% 
  knitr::kable()

# n days
summ_track$range %>% mean
summ_track$range %>% range
# n positions
mov_track %>%
  nrow
summ_track$n %>% mean
summ_track$n %>% range


#' ## adehabitat

# Coordinates
coords <- mov_track[, c("x_", "y_")] %>% 
  dplyr::rename(x = x_, y = y_)
# Time
head(mov_track$t_)

# ltraj object
mov_traj <- as.ltraj(xy = coords, date = mov_track$t_, id = mov_track$id, 
                     burst = mov_track$id, infolocs = mov_track[, -c(1:3)])

# General information
mov_traj
head(mov_traj[[1]])
mov_traj_df <- ld(mov_traj) # Transform into common data.frame
plot(mov_traj[3]); mov_traj[3]

# maximum distance crossed by the animals, in straight line
# (a first proxy for dispersal)
mov_traj_df %>% 
  dplyr::group_by(id) %>% 
  dplyr::summarise(
    r2max = sqrt(max(R2n))/1000,
    tot_dist = sum(dist, na.rm = T)/1000
  ) %>% 
  print(n = 30)

# Visualization

#' # Where are the individuals?
#' 
#' Here we explore where are the animals monitored in a few simple plots.

# --------------- label=visualization, fig.width=15, fig.height=10
# Plot - all
move::plot(move_data, pch = 19, xlab = "Longitude", ylab = "Latitude")

# --------------- label=plot_data, fig.width=9, fig.height=7
# Geographical coordinates
mov_track_ll <- mov_track %>% 
  amt::transform_coords(crs_to = sp::CRS("+init=epsg:4326"))
mov_track_ll$study_name %>% unique

# cols
pal <- leaflet::colorFactor(palette = palette(), domain = factor(mov_track_ll$id))
# cols for dispersive/non dispersive
# pal <- colorFactor(palette = palette(), domain = factor(mov_track_ll$dispersive_behavior))

# Legado
mov_track_ll %>% 
  subset(study_name == "Legado_Aguas") %>% 
  leaflet::leaflet() %>% 
  leaflet::addTiles() %>% 
  leaflet::addCircles(lng = ~x_, lat = ~y_, color = ~pal(id), fill = TRUE, 
                      stroke = T, fillOpacity = 0.7, label = ~id, group = "Animals") %>%
  leaflet::addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>% 
  leaflet::addProviderTiles(providers$Stamen.TonerLabels, group = "Labels") %>% 
  # Layers control
  leaflet::addLayersControl(
    baseGroups = c("Satellite", "Labels"),
    overlayGroups = c("Animals"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>% 
  leaflet::addMiniMap(tiles = providers$Esri.WorldStreetMap,
                      toggleDisplay = TRUE)

# Pardas
mov_track_ll %>% 
  subset(study_name == "Pardas_Tiete") %>% 
  leaflet::leaflet() %>% 
  leaflet::addTiles() %>% 
  leaflet::addCircles(lng = ~x_, lat = ~y_, color = ~pal(id), fill = TRUE,
                      stroke = T, fillOpacity = 0.7, label = ~id, group = "Animals") %>%
  leaflet::addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>% 
  leaflet::addProviderTiles(providers$Stamen.TonerLabels, group = "Labels") %>% 
  # Layers control
  leaflet::addLayersControl(
    baseGroups = c("Satellite", "Labels"),
    overlayGroups = c("Animals"),
    options = layersControlOptions(collapsed = FALSE)
  ) %>% 
  leaflet::addMiniMap(tiles = providers$Esri.WorldStreetMap,
                      toggleDisplay = TRUE)

#' The plots below show the positions of each animal separetelly, without satellite images.
#' Animals in blue were (a priori considered as) residents and animals in orange were dispersers.

# --------------- label=plot_no_satellite, fig.width=15, fig.height=10
# Plots without satellite image

# Cougars
g2 <- ggplot(mov_track) + 
  geom_point(aes(x = x_, y = y_, col = dispersal_behavior)) +
  facet_wrap(~id, scales = "free", ncol = 3) +
  # coord_map() +
  theme(legend.position="none") +
  labs(x = "x", y = "y")
g2

# lay <- matrix(c(rep(1,16), rep(c(2,2,2,2),2)), nrow = 4)
# gridExtra::grid.arrange(g1, g2, layout_matrix = lay)

# --------------- label=monitoring_period

#' # Monitoring period
#' 
#' Here we have shown for how long and when each individual was monitored.

# Monitoring period
g_samp0 <- mov_track %>% 
  ggplot() +
  geom_point(aes(x = t_, y = fct_reorder(id, id_number, .desc = T))) +
  theme_bw(base_size = 12) +
  scale_x_datetime(labels = scales::date_format("%Y"),
                   date_minor_breaks = "6 months", date_breaks = "1 year") +
  labs(x = "",
       y = "Individual")
g_samp0
