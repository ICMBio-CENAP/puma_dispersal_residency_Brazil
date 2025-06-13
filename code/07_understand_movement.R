#' ---
#' title: 'Understanding puma movement characteristics'
#' author: Bernardo Niebuhr
#' output:
#'   github_document: default
#'   html_document: default
#' ---

# Load packages
if(!require(install.load)) install.packages("install.load"); library(install.load)

install.load::install_load("ezknitr", "knitr")
install.load::install_load("tidyverse", "lubridate", "purrr", "ggpubr", "GGally")
install.load::install_load("amt", "lme4")
install.load::install_load("broom", "ggeffects")

# load data
library(PardasIPC)

# Print options for this document
options(width = 165)
opts_knit$set(root_dir = rprojroot::find_rstudio_root_file()) # root project folder
opts_chunk$set(eval = T, error = F, message = F, warning = F, cache = F, echo = T, results = T)

# --------------- label=setup
# Set up 

# Clean everything before beginning
rm(list = ls())

# --------------- label=load_data_and_functions

# Load data
load("data/movement_data_pumas_dispersal_behavior.rda")
pumas_behav

# add translocation data
pumas_behav <- pumas_behav %>% 
  dplyr::left_join(
    pumas_metadata |> 
      dplyr::select(id = name, dispersal_translocated, translocated),
    by = "id")

movement_data_pumas_behavior <- pumas_behav
movement_data_pumas_behavior %>% print(width = Inf)
str(movement_data_pumas_behavior)

# crs to use
crs_use <- "+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m"

# transform data into track object
mov_track <- movement_data_pumas_behavior %>% 
  amt::mk_track(X, Y, timestamp, all_cols = T, crs = "+init=epsg:4326") %>% 
  dplyr::mutate(phase = relevel(as_factor(phase), ref = "pre-dispersal"),
                selected_model = as_factor(selected_model)) 
str(mov_track)

# load maps
file_maps <- list.files("spatial_data/analysis/", pattern = ".tif", full.names = T)

# maps_raster <- terra::stack(file_maps)
# land_use_raster <- raster::raster(file_maps[4]) 

maps_raster <- terra::rast(file_maps)
names(maps_raster) <- file_maps |> 
  strsplit(split = "/") |> 
  sapply(last) |> 
  gsub(pattern = ".tif", replacement = "")
land_use_raster <- terra::rast(file_maps[4])
names(land_use_raster) <- names(maps_raster)[4]

# --------------- label=annotate_data

# divide trajectories into bursts if there are gaps and extract the landscape information
# reproject as well
mov_track_annotated <- mov_track %>% 
  # reproject
  amt::transform_coords(crs_to = terra::crs(maps_raster)) %>% 
  amt::time_of_day() %>% 
  # nest
  tidyr::nest_legacy(-id) %>% 
  # divide into bursts, calculate step length and turning angle, extract covariates
  dplyr::mutate(trk_rsp = purrr::map(data, function(x) 
    x %>% 
      amt::track_resample(rate = hours(1), tolerance = minutes(45)) %>% 
      amt::filter_min_n_burst(min_n = 3) %>% 
      amt::steps_by_burst(keep_cols = "start") %>% 
      amt::extract_covariates(maps_raster))) |> 
  # remove original data
  dplyr::select(-data) %>% 
  # unnest
  tidyr::unnest_legacy(trk_rsp)
mov_track_annotated

# save(mov_track_annotated, file = "data/mov_track_annotated.RData")
load("data/mov_track_annotated.RData")

# --------------- label=explore_habitat_use

colnames(mov_track_annotated)

# land use classes
labs <- c("water", "other anthropogenic", "urban", "forest",
          "non forest natural", "forestry", "sugarcane", "pasture")
mov_track_annotated <- mov_track_annotated %>% 
  dplyr::mutate(land_use_FBDS = factor(land_use_FBDS_8_30m , 
                                             levels = 1:8,
                                             labels = labs),
                log_dist_hidro = log10(dist_hidro_m + 1),
                log_dist_road = log10(dist_road_m + 1),
                log_dist_urban = log10(dist_urban_m + 1))

mov_track_annotated <- mov_track_annotated %>% 
  dplyr::mutate(phase_simplified = ifelse(phase == "dispersal", "dispersal", "ranging"),
                id_trans = ifelse(translocated == "translocated", paste(id, "(translocated)"),
                                  id))

# plot number of points per land use class
plotdata <- mov_track_annotated %>% 
  dplyr::group_by(phase_simplified, tod_, land_use_FBDS) %>%
  # dplyr::group_by(phase, land_use_FBDS) %>%
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(
    pct = n/sum(n),
    lbl = scales::percent(pct)
  )

cols_rast <- c("#80b1d3", "#fccde5", "#fb8072", "#b3de69", "#8dd3c7", 
               "#bebada", "#fdb462", "#ffffb3")

# plot for all individuals
g1 <- plotdata %>% 
  ggplot(aes(x = phase_simplified, y = pct, fill = land_use_FBDS)) + 
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~tod_) + 
  scale_y_continuous(breaks = seq(0, 1, .2), 
                     label = scales::percent) +
  geom_text(aes(label = lbl), 
            size = 3, 
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = cols_rast) +
  labs(y = "Proportion of use", 
       fill = "Land use class",
       x = "Movement phase",
       title = "") +
  theme_minimal()
g1

ggsave(filename = 'habitat_use_vs_tod_phase.png', plot = g1, 
       path = "output/03_mov_patterns_discrete/", device = 'png', 
       width = 17, height = 15, units = 'cm', dpi = 300)

# compare between translocated and non-translocated individuals
plotdata_trans <- mov_track_annotated %>% 
  dplyr::group_by(phase_simplified, dispersal_translocated, tod_, land_use_FBDS) %>%
  # dplyr::group_by(phase, land_use_FBDS) %>%
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(
    pct = n/sum(n),
    lbl = scales::percent(pct)
  )

# plot for all individuals
g1_trans <- plotdata_trans %>% 
  ggplot(aes(x = phase_simplified, y = pct, fill = land_use_FBDS)) + 
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(tod_ ~ dispersal_translocated) + 
  scale_y_continuous(breaks = seq(0, 1, .2), 
                     label = scales::percent) +
  geom_text(aes(label = lbl), 
            size = 3, 
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = cols_rast) +
  labs(y = "Proportion of use", 
       fill = "Land use class",
       x = "Movement phase",
       title = "") +
  theme_minimal()
g1_trans

ggsave(filename = 'habitat_use_vs_tod_phase_trans_vs_nontrans.png', plot = g1_trans, 
       path = "output/03_mov_patterns_discrete/", device = 'png', 
       width = 17, height = 15, units = 'cm', dpi = 300)

# # plot for all individuals - pt
# labs.pt <- c("Água", "Outras áreas antrópicas", "Áreas urbanas", "Floresta",
#           "Áreas naturais não florestais", "Silvicultura", "Cana-de-acúcar", "Pastagem")
# g1.pt <- plotdata %>% 
#   dplyr::ungroup() %>% 
#   dplyr::mutate(phase_simplified = factor(phase_simplified, labels = c("residência", "dispersão")),
#                 tod_ = factor(tod_, labels = c("dia", "noite")),
#                 land_use_FBDS = factor(land_use_FBDS, labels = labs.pt)) %>% 
#   ggplot(aes(x = phase_simplified, y = pct, fill = land_use_FBDS)) + 
#   geom_bar(stat = "identity", position = "fill") +
#   facet_wrap(~tod_) + 
#   scale_y_continuous(breaks = seq(0, 1, .2), 
#                      label = scales::percent) +
#   geom_text(aes(label = lbl), 
#             size = 3, 
#             position = position_stack(vjust = 0.5)) +
#   scale_fill_manual(values = cols_rast) +
#   labs(y = "Proporção de uso", 
#        fill = "Uso da terra",
#        x = "Fase de movimento",
#        title = "") +
#   theme_minimal()
# g1.pt
# 
# ggsave(filename = 'habitat_use_vs_tod_phase_pt.png', plot = g1.pt, 
#        path = "output", device = 'png', 
#        width = 17, height = 15, units = 'cm', dpi = 300)

plotdata_ind <- mov_track_annotated %>% 
  dplyr::group_by(phase_simplified = ifelse(phase == "dispersal", "dispersal", "ranging") %>% 
                    as_factor()) %>% 
  dplyr::group_by(id_trans, phase_simplified, tod_, land_use_FBDS) %>%
  # dplyr::group_by(phase, land_use_FBDS) %>%
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(
    pct = n/sum(n),
    lbl = scales::percent(pct)
  )

g2 <- plotdata_ind %>% 
  ggplot(aes(x = phase_simplified, y = pct, fill = land_use_FBDS)) + 
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(id_trans~tod_) + 
  scale_y_continuous(breaks = seq(0, 1, .2), 
                     label = scales::percent) +
  # geom_text(aes(label = lbl), 
  #           size = 3, 
  #           position = position_stack(vjust = 0_5)) +
  scale_fill_manual(values = cols_rast) +
  labs(y = "Proportion of use", 
       fill = "Land use class",
       x = "Movement phase",
       title = "") +
  theme_minimal()
g2

ggsave(filename = 'habitat_use_vs_tod_phase_each_individual.png', plot = g2, 
       path = "output/03_mov_patterns_discrete/", device = 'png', 
       width = 30, height = 23, units = 'cm', dpi = 300)

# some notes
# Aracatuba: more pasture and other anthropic during dispersal
# Jussara: only forest, only residency
# Kurupi: less antropic, more sugarcane, more forest at night during dispersal
# Mineiro: less sugarcane, more forest during dispersal
# Nick: more forest during dispersal
# Pepira: less sugarcane, more forest during dispersal!
# Piloto: less forest, less sugarcane (obvious), more pasture during dispersal
# Pora: resident
# rafiki: resident
# sucuri: resident
# tupa: less water, more non forest natural, more forest during dispersal
# zeus: resident
# zorro: less forest, less pasture, more sugarcane during dispersal - opposite to others

# if the amount of sugarcane and other antropic increases or decreases during dispersal
# is related to where the animal is coming from and where it is going

# plot other variables

# road dist
g3 <- ggplot(mov_track_annotated, 
       aes(x = phase_simplified, y = dist_road_m/1000, color = tod_)) + 
  stat_summary(geom="errorbar", fun_data=mean_cl_normal, fun_args=list(conf_int=0.95), 
               position = position_dodge(width = .7), size = 1, width = 0) +
  # stat_summary(geom="line", fun_y=mean, linetype="dashed")+
  stat_summary(geom="point", fun = mean, size = 3, position = position_dodge(width = 0.7)) +
  # facet_wrap(~sp, scales = 'free') +
  theme_bw() +
  labs(x = 'Movement phase', y = 'Distance to roads (km)', color = NULL)
g3
ggsave(filename = 'avg_dist_to_roads_vs_dispersal_tod.png', plot = g3, 
       path = "output/03_mov_patterns_discrete/", device = 'png', 
       width = 12, height = 9, units = 'cm', dpi = 300)

g3_simp <- ggplot(mov_track_annotated %>% filter(phase != "pre-dispersal"), 
             aes(x = phase_simplified, y = dist_road_m/1000, color = tod_)) + 
  stat_summary(geom="errorbar", fun.data=mean_cl_normal, fun.args=list(conf_int=0.95), 
               position = position_dodge(width = .7), size = 1, width = 0) +
  # stat_summary(geom="line", fun_y=mean, linetype="dashed")+
  stat_summary(geom="point", fun = mean, size = 3, position = position_dodge(width = 0.7)) +
  # facet_wrap(~sp, scales = 'free') +
  theme_bw() +
  ylim(0, NA) +
  labs(x = 'Movement phase', y = 'Distance to roads (km)', color = NULL)
g3_simp
ggsave(filename = 'avg_dist_to_roads_vs_dispersal_tod_no_pre.png', plot = g3_simp, 
       path = "output/03_mov_patterns_discrete/", device = 'png', 
       width = 12, height = 9, units = 'cm', dpi = 300)

# log scale, road dist
ggplot(mov_track_annotated, 
       aes(x = phase_simplified, y = log_dist_road, color = tod_)) + 
  stat_summary(geom="errorbar", fun.data=mean_cl_normal, fun.args=list(conf.int=0.95), 
               position = position_dodge(width = .7), size = 1, width = 0) +
  # stat_summary(geom="line", fun_y=mean, linetype="dashed")+
  stat_summary(geom="point", fun = mean, size = 3, position = position_dodge(width = 0.7)) +
  # facet_wrap(~sp, scales = 'free') +
  theme_bw() +
  labs(x = 'Movement phase', y = 'Log-distance to roads (m)', color = NULL)

# road density
g4 <- ggplot(mov_track_annotated, 
       aes(x = phase_simplified, y = road_density_km_per100km2, color = tod_)) + 
  stat_summary(geom="errorbar", fun.data=mean_cl_normal, fun.args=list(conf.int=0.95), 
               position = position_dodge(width = .7), size = 1, width = 0) +
  # stat_summary(geom="line", fun_y=mean, linetype="dashed")+
  stat_summary(geom="point", fun = mean, size = 3, position = position_dodge(width = 0.7)) +
  # facet_wrap(~sp, scales = 'free') +
  theme_bw() +
  ylim(0, NA) +
  labs(x = 'Movement phase', y = 'Road density (km/100km2)', color = NULL)
g4
ggsave(filename = 'roads_density_vs_dispersal_tod.png', plot = g4, 
       path = "output/03_mov_patterns_discrete/", device = 'png', 
       width = 12, height = 9, units = 'cm', dpi = 300)

g4_simp <- ggplot(mov_track_annotated %>% filter(phase != "pre-dispersal"), 
             aes(x = phase_simplified, y = road_density_km_per100km2, color = tod_)) + 
  stat_summary(geom="errorbar", fun.data=mean_cl_normal, fun.args=list(conf.int=0.95), 
               position = position_dodge(width = .7), size = 1, width = 0) +
  # stat_summary(geom="line", fun.y=mean, linetype="dashed")+
  stat_summary(geom="point", fun = mean, size = 3, position = position_dodge(width = 0.7)) +
  # facet_wrap(~sp, scales = 'free') +
  theme_bw() +
  ylim(0, NA) +
  labs(x = 'Movement phase', y = 'Road density (km/100km2)', color = NULL)
g4_simp
ggsave(filename = 'roads_density_vs_dispersal_tod_no_pre.png', plot = g4_simp, 
       path = "output/03_mov_patterns_discrete/", device = 'png', 
       width = 12, height = 9, units = 'cm', dpi = 300)

# roads and comparison between translocated and non-translocated
g4_simp_trans <- ggplot(mov_track_annotated %>% filter(phase != "pre-dispersal"), 
                  aes(x = phase_simplified, y = road_density_km_per100km2, color = tod_)) + 
  stat_summary(geom="errorbar", fun.data=mean_cl_normal, fun.args=list(conf.int=0.95), 
               position = position_dodge(width = .7), size = 1, width = 0) +
  # stat_summary(geom="line", fun.y=mean, linetype="dashed")+
  stat_summary(geom="point", fun = mean, size = 3, position = position_dodge(width = 0.7)) +
  facet_wrap(~dispersal_translocated, scales = 'fixed') +
  theme_bw() +
  ylim(0, NA) +
  labs(x = 'Movement phase', y = 'Road density (km/100km2)', color = NULL)
g4_simp_trans
ggsave(filename = 'roads_density_vs_dispersal_tod_no_pre_trans_vs_nontrans.png', plot = g4_simp_trans, 
       path = "output/03_mov_patterns_discrete/", device = 'png', 
       width = 12, height = 9, units = 'cm', dpi = 300)

# by individual
g4_simp_id <- ggplot(mov_track_annotated %>% filter(phase != "pre-dispersal"), 
                        aes(x = phase_simplified, y = road_density_km_per100km2, color = tod_)) + 
  stat_summary(geom="errorbar", fun.data=mean_cl_normal, fun.args=list(conf.int=0.95), 
               position = position_dodge(width = .7), size = 1, width = 0) +
  # stat_summary(geom="line", fun.y=mean, linetype="dashed")+
  stat_summary(geom="point", fun = mean, size = 3, position = position_dodge(width = 0.7)) +
  facet_wrap(~id_trans, scales = 'fixed') +
  theme_bw() +
  ylim(0, NA) +
  labs(x = 'Movement phase', y = 'Road density (km/100km2)', color = NULL)
g4_simp_id
ggsave(filename = 'roads_density_vs_dispersal_tod_no_pre_id.png', plot = g4_simp_id, 
       path = "output/03_mov_patterns_discrete/", device = 'png', 
       width = 12, height = 9, units = 'cm', dpi = 300)

# dist water
g5 <- ggplot(mov_track_annotated %>% filter(phase != "pre-dispersal"), 
             aes(x = phase_simplified, y = dist_hidro_m, color = tod_)) + 
  stat_summary(geom="errorbar", fun.data=mean_cl_normal, fun.args=list(conf.int=0.95), 
               position = position_dodge(width = .7), size = 1, width = 0) +
  # stat_summary(geom="line", fun.y=mean, linetype="dashed")+
  stat_summary(geom="point", fun = mean, size = 3, position = position_dodge(width = 0.7)) +
  # facet_wrap(~sp, scales = 'free') +
  theme_bw() +
  lims(y = c(0, NA)) +
  labs(x = 'Movement phase', y = 'Distance to water (m)', color = NULL)
g5
ggsave(filename = 'dist_water_dispersal_tod.png', plot = g4, 
       path = "output/03_mov_patterns_discrete/", device = 'png', 
       width = 12, height = 9, units = 'cm', dpi = 300)

# translocated vs non translocated
g5_trans <- ggplot(mov_track_annotated %>% filter(phase != "pre-dispersal"), 
             aes(x = phase_simplified, y = dist_hidro_m, color = tod_)) + 
  stat_summary(geom="errorbar", fun.data=mean_cl_normal, fun.args=list(conf.int=0.95), 
               position = position_dodge(width = .7), size = 1, width = 0) +
  # stat_summary(geom="line", fun.y=mean, linetype="dashed")+
  stat_summary(geom="point", fun = mean, size = 3, position = position_dodge(width = 0.7)) +
  facet_wrap(~dispersal_translocated, scales = 'fixed') +
  theme_bw() +
  lims(y = c(0, NA)) +
  labs(x = 'Movement phase', y = 'Distance to water (m)', color = NULL)
g5_trans
ggsave(filename = 'dist_water_dispersal_tod_trans.png', plot = g5_trans, 
       path = "output/03_mov_patterns_discrete/", device = 'png', 
       width = 12, height = 9, units = 'cm', dpi = 300)

# dist to urban
g6_simp <- ggplot(mov_track_annotated %>% filter(phase != "pre-dispersal"), 
                  aes(x = phase_simplified, y = dist_road_m/1000, color = tod_)) + 
  stat_summary(geom="errorbar", fun.data=mean_cl_normal, fun.args=list(conf.int=0.95), 
               position = position_dodge(width = .7), size = 1, width = 0) +
  # stat_summary(geom="line", fun.y=mean, linetype="dashed")+
  stat_summary(geom="point", fun = mean, size = 3, position = position_dodge(width = 0.7)) +
  # facet_wrap(~sp, scales = 'free') +
  theme_bw() +
  ylim(0, NA) +
  labs(x = 'Movement phase', y = 'Distance to urban areas (km)', color = NULL)
g6_simp

ggsave(filename = 'avg_dist_to_urban_vs_dispersal_tod_no_pre.png', plot = g6_simp, 
       path = "output/03_mov_patterns_discrete/", device = 'png', 
       width = 12, height = 9, units = 'cm', dpi = 300)

# urban density
g7_simp <- ggplot(mov_track_annotated %>% filter(phase != "pre-dispersal"), 
                  aes(x = phase_simplified, y = urban_density_n_pixels_per100km2, color = tod_)) + 
  stat_summary(geom="errorbar", fun.data=mean_cl_normal, fun.args=list(conf.int=0.95), 
               position = position_dodge(width = .7), size = 1, width = 0) +
  # stat_summary(geom="line", fun_y=mean, linetype="dashed")+
  stat_summary(geom="point", fun = mean, size = 3, position = position_dodge(width = 0.7)) +
  # facet_wrap(~sp, scales = 'free') +
  theme_bw() +
  ylim(0, NA) +
  labs(x = 'Movement phase', y = 'Urban density (% in 100km2)', color = NULL)
g7_simp
ggsave(filename = 'urban_density_vs_dispersal_tod_no_pre.png', plot = g7_simp, 
       path = "output/03_mov_patterns_discrete/", device = 'png', 
       width = 12, height = 9, units = 'cm', dpi = 300)

g7_trans <- ggplot(mov_track_annotated %>% filter(phase != "pre-dispersal"), 
                  aes(x = phase_simplified, y = urban_density_n_pixels_per100km2, color = tod_)) + 
  stat_summary(geom="errorbar", fun.data=mean_cl_normal, fun.args=list(conf.int=0.95), 
               position = position_dodge(width = .7), size = 1, width = 0) +
  # stat_summary(geom="line", fun_y=mean, linetype="dashed")+
  stat_summary(geom="point", fun = mean, size = 3, position = position_dodge(width = 0.7)) +
  facet_wrap(~dispersal_translocated, scales = 'fixed') +
  theme_bw() +
  ylim(0, NA) +
  labs(x = 'Movement phase', y = 'Urban density (% in 100km2)', color = NULL)
g7_trans
ggsave(filename = 'urban_density_vs_dispersal_tod_no_pre_trans_vs_nontrans.png', plot = g7_trans, 
       path = "output/03_mov_patterns_discrete/", device = 'png', 
       width = 12, height = 9, units = 'cm', dpi = 300)

# --------------- label=correlation_variables_and_exploratory_analysis

colnames(mov_track_annotated)
colmns <- c("sex", "weight_kg", "estimated_age_months", "phase", "dist_hidro_m", "dist_road_m",
            "dist_urban_m", "land_use_FBDS", "road_density_km_per100km2", 
            "urban_density_n_pixels_per100km2")
colmns <- which(colnames(mov_track_annotated) %in% colmns) 
corr <- GGally::ggpairs(mov_track_annotated, 
                        columns = colmns)
corr
ggsave("correlations_among_predictors.png", plot = corr, path = "output/03_mov_patterns_discrete/", 
       width = 25, height = 25, units = "cm", dpi = 300)

ggplot(mov_track_annotated) +
  geom_boxplot(aes(land_use_FBDS, sl_))

mov_track_annotated %>%
  ggplot(aes(dist_hidro_m, sl_, colour = land_use_FBDS)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~land_use_FBDS)

mov_track_annotated %>%
  dplyr::filter(phase != "pre-dispersal") %>% 
  ggplot(aes(x = sl_, y = land_use_FBDS, fill = phase)) +
  ggridges::geom_density_ridges() +
  facet_wrap(~phase)

mov_track_annotated %>%
  dplyr::filter(land_use_FBDS == "pasture") %>% 
  ggplot(aes(dist_hidro_m, sl_, colour = land_use_FBDS)) +
  geom_point() +
  geom_smooth()

mov_track_annotated %>%
  ggplot(aes(road_density_km_per100km2, sl_)) +
  geom_smooth(se = F)

mov_track_annotated %>%
  ggplot(aes(road_density_km_per100km2, sl_, colour = land_use_FBDS)) +
  geom_smooth(se = F)

mov_track_annotated %>%
  ggplot(aes(urban_density_n_pixels_per100km2, sl_)) +
  geom_point() +
  geom_smooth(se = F)


# --------------- label=run_models_for_sl
# run models for sl and ta

mov_track_annotated$phase_simplified |> unique()

data_fit <- mov_track_annotated %>% 
  dplyr::filter(phase != "pre-dispersal", land_use_FBDS != "urban",
                land_use_FBDS != "water", land_use_FBDS != "forestry") %>% 
  droplevels() %>% 
  dplyr::mutate(sl_ = sl_ + 1,
                phase = factor(phase_simplified, levels = c("ranging", "dispersal")),) |> 
                # phase = forcats::fct_recode(phase, "ranging" = "residency")) %>% #,
                #land_use_FBDS = as_factor(ifelse(land_use_FBDS == "pasture", "other anthropogenic", land_use_FBDS))) %>% 
  dplyr::rename(road_density = road_density_km_per100km2,
                urban_density = urban_density_n_pixels_per100km2) %>% 
  dplyr::mutate(road_density_std = scale(road_density)/2,
                urban_density_std = scale(urban_density)/2)

# check ranges of variables
summary(data_fit)
data_fit$phase
data_fit$dispersal_translocated
data_fit$type <- factor(data_fit$dispersal_translocated, levels = c("resident", "disperser", "translocated"))

# sex and body size are correlated
bs_sex <- data_fit %>% 
  dplyr::group_by(id_number) %>% 
  dplyr::summarise(sex = first(sex),
            weight_kg = first(weight_kg)) 

cor.test(as.numeric(as_factor(bs_sex$sex)), bs_sex$weight_kg)
aov(bs_sex$weight_kg ~ as_factor(bs_sex$sex))

bs_sex %>% 
  ggplot(aes(sex, weight_kg)) +
  geom_boxplot() +
  geom_point(aes(colour = id_number))

# no effect model
m0 <- glm(sl_ ~ 1, data = data_fit, family = Gamma(link = "log"))
# only additive terms
m0_5 <- glm(sl_ ~ sex + tod_ + phase + land_use_FBDS + log(dist_road_m + 1) + log(dist_urban_m + 1) + 
              log(dist_hidro_m + 1), data = data_fit, family = Gamma(link = "log"))
# only sex and land use, with interactions
m1 <- glm(sl_ ~ sex + tod_ * phase * land_use_FBDS, data = data_fit, family = Gamma(link = "log"))
# m2 <- glmer(sl_ ~ sex + tod_ * phase * land_use_FBDS + (1|name), data = data_fit, family = Gamma(link = "log"))
# adding road density
m3 <- glm(sl_ ~ sex + tod_ + phase + tod_:phase:land_use_FBDS + 
            tod_:phase:road_density, data = data_fit, family = Gamma(link = "log"))
# adding distance to water
m4 <- glm(sl_ ~ sex + tod_ + phase + tod_:phase:land_use_FBDS + 
            tod_:phase:road_density +
            log(dist_hidro_m + 1):tod_:phase, data = data_fit, family = Gamma(link = "log"))
# adding urban density
m5 <- glm(sl_ ~ sex + tod_ + phase + 
            tod_:phase:land_use_FBDS + tod_:phase:road_density +
            tod_:phase:log_dist_hidro + tod_:phase:urban_density, 
          data = data_fit, family = Gamma(link = "log"))
# using standardized urban and road desity instead
m5_5 <- glm(sl_ ~ sex + tod_ + phase + 
            tod_:phase:land_use_FBDS + tod_:phase:scale(road_density) +
            tod_:phase:log(dist_hidro_m + 1) + tod_:phase:scale(urban_density), 
            data = data_fit, family = Gamma(link = "log"))
# m6 <- glm(sl_ ~ sex + tod_ + phase + land_use_FBDS + log_dist_road + log_dist_urban + log_dist_hidro + 
# Using log_distance to roads instead of density
m6 <- glm(sl_ ~ sex + tod_ + phase +  
            tod_:phase:land_use_FBDS + tod_:phase:log(dist_road_m + 1) +
            phase:tod_:log(dist_urban_m + 1) + tod_:phase:log(dist_hidro_m + 1), 
          data = data_fit, family = Gamma(link = "log"))
# using log in the model statement, with no water which makes little difference
m7 <- glm(sl_ ~ sex + tod_ + phase +  
            tod_:phase:land_use_FBDS + tod_:phase:log(dist_road_m + 1) +
            phase:tod_:log(dist_urban_m + 1), 
          data = data_fit, family = Gamma(link = "log"))
# adding difference between types of individuals
m8 <- glm(sl_ ~ sex + tod_ + phase + type +    
            tod_:phase:type:land_use_FBDS + tod_:phase:type:log(dist_road_m + 1) +
            phase:tod_:type:log(dist_urban_m + 1), 
          data = data_fit, family = Gamma(link = "log"))
# removing phase
m9 <- glm(sl_ ~ sex + tod_ + type +    
            tod_:type:land_use_FBDS + tod_:type:log(dist_road_m + 1) +
            tod_:type:log(dist_urban_m + 1), 
          data = data_fit, family = Gamma(link = "log"))
# removing tod
m10 <- glm(sl_ ~ sex + phase + type +    
             phase:type:land_use_FBDS + phase:type:log(dist_road_m + 1) +
             phase:type:log(dist_urban_m + 1), 
          data = data_fit, family = Gamma(link = "log"))

summary(m0)
summary(m1)
summary(m3)
summary(m4)
summary(m5)
summary(m6)
summary(m8)

bbmle::AICctab(m0, m0_5, m1, m3, m4, m5, m5_5, m6, m7, m8, m9, m10)

# m1
plot(ggpredict(m1, c("phase")))
plot(ggpredict(m1, c("sex")))
plot(ggpredict(m1, c("tod_")))
plot(ggpredict(m1, c("land_use_FBDS", "phase")))
plot(ggpredict(m1, c("land_use_FBDS", "phase")))
plot(ggpredict(m1, c("land_use_FBDS", "tod_")))
plot(ggpredict(m1, c("land_use_FBDS", "tod_", "phase")))
plot(ggpredict(m1, c("road", "tod_", "phase")))

broom::tidy(m1)

# m3
plot(ggpredict(m3, c("phase")))
plot(ggpredict(m3, c("sex")))
plot(ggpredict(m3, c("tod_")))
plot(ggpredict(m3, c("land_use_FBDS")))
plot(ggpredict(m3, c("road_density")))
plot(ggpredict(m3, c("land_use_FBDS", "phase")))
plot(ggpredict(m3, c("land_use_FBDS", "phase")))
plot(ggpredict(m3, c("land_use_FBDS", "tod_")))
plot(ggpredict(m3, c("land_use_FBDS", "phase", "tod_")))
plot(ggpredict(m3, c("road_density", "tod_")))
plot(ggpredict(m3, c("road_density", "phase")))
plot(ggpredict(m3, c("road_density", "tod_", "phase")))

broom::tidy(m3)
summary(m3)

# m4
plot(ggpredict(m4, c("phase")))
plot(ggpredict(m4, c("sex")))
plot(ggpredict(m4, c("tod_")))
plot(ggpredict(m4, c("land_use_FBDS")))
plot(ggpredict(m4, c("road_density")))
plot(ggpredict(m4, c("log_dist_hidro")))
plot(ggpredict(m4, c("land_use_FBDS", "phase")))
plot(ggpredict(m4, c("land_use_FBDS", "phase")))
plot(ggpredict(m4, c("land_use_FBDS", "tod_")))
plot(ggpredict(m4, c("land_use_FBDS", "phase", "tod_")))
plot(ggpredict(m4, c("road_density", "tod_")))
plot(ggpredict(m4, c("road_density", "phase")))
plot(ggpredict(m4, c("road_density", "tod_", "phase")))
plot(ggpredict(m4, c("log_dist_hidro", "tod_")))
plot(ggpredict(m4, c("log_dist_hidro", "phase")))
plot(ggpredict(m4, c("log_dist_hidro", "phase", "tod_")))

broom::tidy(m4) %>% 
  print(n = 100)
summary(m4)

# m5
plot(ggpredict(m5, c("phase")))
plot(ggpredict(m5, c("sex")))
plot(ggpredict(m5, c("tod_")))
plot(ggpredict(m5, c("land_use_FBDS")))
plot(ggpredict(m5, c("road_density")))
plot(ggpredict(m5, c("urban_density")))
plot(ggpredict(m5, c("log_dist_hidro")))
plot(ggpredict(m5, c("land_use_FBDS", "phase")))
plot(ggpredict(m5, c("land_use_FBDS", "tod_")))
plot(ggpredict(m5, c("land_use_FBDS", "phase", "tod_")))
plot(ggpredict(m5, c("road_density", "tod_")))
plot(ggpredict(m5, c("road_density", "phase")))
plot(ggpredict(m5, c("road_density", "tod_", "phase")))
plot(ggpredict(m5, c("log_dist_hidro", "tod_")))
plot(ggpredict(m5, c("log_dist_hidro", "phase")))
plot(ggpredict(m5, c("log_dist_hidro", "phase", "tod_")))
plot(ggpredict(m5, c("urban_density", "phase")))
plot(ggpredict(m5, c("urban_density", "tod_")))
plot(ggpredict(m5, c("urban_density", "phase", "tod_")))

broom::tidy(m5) %>% 
  print(n = 100)
summary(m5)

# m6
cols <- c("steelblue", "firebrick")

plot(ggemmeans(m6, c("sex")))
# plot(ggemmeans(m6, c("phase")))
plot(ggpredict(m6, c("phase"))) 
# plot(ggemmeans(m6, c("tod_")))
# average difference between sexes
ggemmeans(m6, c("sex")) * 24
# average difference between day and night
ggemmeans(m6, c("tod_")) * 24 # both sexes
ggemmeans(m6, c("tod_"), condition = c(sex = "Female"))
ggemmeans(m6, c("tod_"), condition = c(sex = "Male"))
plot(ggpredict(m6, c("tod_")))
plot(ggemmeans(m6, c("tod_", "phase")))
plot(ggpredict(m6, c("land_use_FBDS")))
plot(ggpredict(m6, c("dist_road_m")))
plot(ggpredict(m6, c("dist_urban_m")))
plot(ggpredict(m6, c("dist_hidro_m")))
plot(ggpredict(m6, c("land_use_FBDS", "phase")))
plot(ggpredict(m6, c("land_use_FBDS", "tod_")))

res1 <- ggpredict(m6, c("land_use_FBDS", "phase", "tod_"), 
                  condition = c(sex = "Male"))
res1$x <- res1$x %>% 
  forcats::fct_relevel("other anthropogenic", after = 2) %>% 
  forcats::fct_rev()
res1$group <- res1$group |> 
  forcats::fct_recode("resident" = "post-dispersal")
res1$predicted <- res1$predicted * 24/1000
res1$conf.low <- res1$conf.low * 24/1000
res1$conf.high <- res1$conf.high * 24/1000
g_res1 <- plot(res1,
               show.title = F, show.x.title = F, show.y.title = F) + 
  labs(y = "Movement rate (km/day)",
       x = "Land use class", 
       colour = "Movement phase") +
  ylim(0, 15) +
  scale_color_manual(values = cols)
g_res1 + coord_flip()

ggsave("results_landuse_phase_tod.png", plot = g_res1, 
       path = "output/03_mov_patterns_discrete/", 
         width =  28, height = 12, units = "cm", dpi = 300)

ggsave("results_landuse_phase_tod2.png", plot = g_res1 +
           coord_flip(), 
         path = "output/03_mov_patterns_discrete/", 
         width =  20, height = 10, units = "cm", dpi = 300)


# # pt
# classes <- c("Floresta", "Áreas naturais não florestais", 
#              "Outras áreas antrópicas", "Cana-de-açúcar",
#              "Pastagem")
# res1_pt <- ggpredict(m6, c("land_use_FBDS", "phase", "tod_"), 
#                   condition = c(sex = "Male"))
# res1_pt$x <- res1_pt$x %>% 
#   forcats::fct_relevel("other anthropogenic", after = 2) %>% 
#   forcats::fct_rev() %>% 
#   factor(labels = classes)
# res1_pt$group <- factor(res1_pt$group, labels = c("dispersão", "pós-dispersão"))
# res1_pt$facet <- factor(res1_pt$facet, labels = c("dia", "noite"))
# attr(res1_pt, "x.axis.labels") <- classes[c(3,1,2,4,5)]
# g_res1_pt <- plot(res1_pt,
#                show.title = F, show.x.title = F, show.y.title = F) + 
#   labs(y = "Taxa de movimento (m/h)",
#        x = "Uso da terra", 
#        colour = "Fase") +
#   scale_color_manual(values = cols)
# g_res1_pt + coord_flip()
# 
# ggsave("results_landuse_phase_tod_pt.png", plot = g_res1_pt, 
#        path = "output/mov_patterns_discrete/", 
#        width =  28, height = 12, units = "cm", dpi = 300)
# 
# ggsave("results_landuse_phase_tod2_pt.png", plot = g_res1_pt +
#          coord_flip(), 
#        path = "output/mov_patterns_discrete/", 
#        width =  20, height = 10, units = "cm", dpi = 300)

plot(ggpredict(m6, c("dist_road_m", "tod_")))
plot(ggpredict(m6, c("dist_road_m", "phase")), lim = c(0, NA))

res2 <- ggpredict(m6, c("dist_road_m", "phase", "tod_"), 
                  condition = c(sex = "Male"))
res2$x <- res2$x/1e3
res2$predicted <- res2$predicted * 24/1e3
res2$conf.low <- res2$conf.low * 24/1e3
res2$conf.high <- res2$conf.high * 24/1e3
res2$group <- res2$group |> 
  forcats::fct_recode("resident" = "post-dispersal")
g_res2 <- plot(res2, lim = c(0, 800), show.title = F) + 
  xlim(0, 5) + 
  labs(y = "Movement rate (km/day)",
       x = "Distance to roads (km)",
       colour = "Movement phase") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols)
g_res2
ggsave("results_roads_phase_tod.png", plot = g_res2, 
       path = "output/03_mov_patterns_discrete/", 
       width =  20, height = 15, units = "cm", dpi = 300)

plot(ggpredict(m6, c("dist_road_m", "tod_", "phase", "land_use_FBDS")), 
     lim = c(0, 600), show.title = F)

plot(ggpredict(m6, c("dist_hidro_m", "tod_")))
plot(ggpredict(m6, c("dist_hidro_m", "phase")))
plot(ggpredict(m6, c("dist_hidro_m", "phase", "tod_")))

plot(ggpredict(m6, c("dist_urban_m", "phase")))
plot(ggpredict(m6, c("dist_urban_m", "tod_")), lim = c(0, NA))
plot(ggpredict(m6, c("dist_urban_m", "phase", "tod_", "land_use_FBDS")), lim = c(0, NA))

res3 <- ggpredict(m6, c("dist_urban_m", "phase", "tod_"), condition = c(sex = "Male"))
res3$x <- res3$x/1e3
res3$predicted <- res3$predicted * 24/1e3
res3$conf.low <- res3$conf.low * 24/1e3
res3$conf.high <- res3$conf.high * 24/1e3
res3$group <- res3$group |> 
  forcats::fct_recode("resident" = "post-dispersal")
g_res3 <- plot(res3, lim = c(0, NA), show.title = F) + 
  xlim(0, 5) + 
  labs(y = "Movement rate (km/day)",
       x = "Distance to urban areas (km)",
       colour = "Movement phase") +
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols)
g_res3
ggsave("results_urban_phase_tod.png", plot = g_res3, 
       path = "output/03_mov_patterns_discrete/", 
       width =  20, height = 15, units = "cm", dpi = 300)

g_res4 <- ggpubr::ggarrange(g_res2, g_res3, nrow = 2, labels = c("A", "B"),
                            common.legend = T, legend = "right")

ggsave("results_roads_urban_phase_tod.png", plot = g_res4, 
       path = "output/03_mov_patterns_discrete/", 
       width =  20, height = 18, units = "cm", dpi = 300)

# Joint figure
g_res4_1 <- ggpubr::ggarrange(g_res1_flip, 
                              g_res2, g_res3,
                              nrow = 3, labels = c("A", "B", "C"),
                              common.legend = T, legend = "right", align = "v")

ggsave("results_lc_roads_urban_phase_tod.png", plot = g_res4_1, 
       path = "output/03_mov_patterns_discrete/", 
       width =  20, height = 22, units = "cm", dpi = 300)


# export coefficients
coefs <- broom::tidy(m6) %>% 
  dplyr::mutate(significant = as.factor(if_else(p.value < 0.05, 1, 0)),
                p = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)))
coefs %>% 
  mutate_if(is.numeric, round, digits = 4) %>% 
  readr::write_csv(file = "output/03_mov_patterns_discrete/coef_table_movement_rate.csv")

summary(m6)

car::vif(m0_5)
car::vif(m1)
car::vif(m6)
plot(m6)

# residual check
plot(density(resid(m0, type='response')))
lines(density(resid(m6, type='response')), col='red')

plot(density(resid(m0, type='pearson')))
lines(density(resid(m6, type='pearson')), col='red')

plot(density(rstandard(m0, type='pearson')))
lines(density(rstandard(m6, type='pearson')), col='red')

plot(density(resid(m0, type='deviance')))
lines(density(resid(m6, type='deviance')), col='red')

plot(density(rstandard(m0, type='deviance')))
lines(density(rstandard(m6, type='deviance')), col='red')

library(statmod)
plot(density(qresid(m0)))
lines(density(qresid(m6)), col='red')

par(mfrow=c(1,2))
scatter.smooth(1:length(rstandard(m0, type='deviance')), rstandard(m0, type='deviance'), col='gray')
scatter.smooth(1:length(rstandard(m6, type='deviance')), rstandard(m6, type='deviance'), col='gray')

# residuals vs response
par(mfrow=c(1,2))
scatter.smooth(predict(m0, type='response'), rstandard(m0, type='deviance'), col='gray')
scatter.smooth(predict(m6, type='response'), rstandard(m6, type='deviance'), col='gray')

par(mfrow=c(1,2))
scatter.smooth(sqrt(predict(m0, type='response')), qresid(m0), col='gray')
scatter.smooth(sqrt(predict(m6, type='response')), qresid(m1), col='gray')

# residuals vs predictors

# qqnorm
par(mfrow=c(1,2))
qqnorm(qresid(m0)); qqline(qresid(m0))
qqnorm(qresid(m6)); qqline(qresid(m6))

# cooks distance
par(mfrow=c(1,2))
plot(cooks.distance(m0), type='h')
plot(cooks.distance(m6), type='h')

cooksd_m0 <- cooks.distance(m0)
cooksd_m6 <- cooks.distance(m6)

length(cooksd_m0[cooksd_m0 > mean(cooksd_m0) * 2])
length(cooksd_m6[cooksd_m6 > mean(cooksd_m6) * 2])

which(influence.measures(m0)$is.inf[,'cook.d'] )

library(boot)
m6_diag <- glm.diag(m6)
glm.diag.plots(m6, m6_diag)

#------------------------
# m8
cols <- c("steelblue", "firebrick")
cols <- c("#117733", "#88CCEE", "#CC6677")
cols <- c("#000000", "#009E73", "#D55E00")

summary(m8)

plot(ggemmeans(m8, c("sex")))
# plot(ggemmeans(m6, c("phase")))
plot(ggpredict(m8, c("phase"))) 
plot(ggpredict(m8, c("type"))) 
# average difference between sexes
ggemmeans(m8, c("sex")) * 24
# average difference between day and night
ggemmeans(m8, c("tod_")) * 24 # both sexes
ggemmeans(m8, c("tod_"), condition = c(sex = "Female"))
ggemmeans(m8, c("tod_"), condition = c(sex = "Male"))
plot(ggpredict(m8, c("tod_")))
plot(ggemmeans(m8, c("tod_", "phase")))
plot(ggpredict(m8, c("land_use_FBDS")))
plot(ggpredict(m8, c("dist_road_m")))
plot(ggpredict(m8, c("dist_urban_m")))
# plot(ggpredict(m8, c("dist_hidro_m")))
plot(ggpredict(m8, c("land_use_FBDS", "phase")))
plot(ggpredict(m8, c("land_use_FBDS", "type")))
plot(ggpredict(m8, c("land_use_FBDS", "tod_")))
plot(ggpredict(m8, c("land_use_FBDS", "type", "phase")))
# plot(ggpredict(m8, c("land_use_FBDS", "phase", "type")))

res1 <- ggpredict(m8, c("land_use_FBDS", "type", "phase"), 
                  condition = c(sex = "Male"))
res1$x <- res1$x %>% 
  forcats::fct_relevel("other anthropogenic", after = 2) %>% 
  forcats::fct_rev()
res1$predicted <- res1$predicted * 24/1000
res1$conf.low <- res1$conf.low * 24/1000
res1$conf.high <- res1$conf.high * 24/1000
g_res1 <- plot(res1,
               show_title = F, show_x_title = F, show_y_title = F) +
  labs(y = "Movement rate (km/day)",
       x = "Land use class", 
       colour = "Movement phase") +
  # ylim(0, 10) +
  scale_color_manual(values = cols) +
  coord_cartesian(ylim = c(0, 20))
g_res1
g_res1_flip <- g_res1 + coord_flip(ylim = c(0, 20))
g_res1_flip

ggsave("results_landuse_phase_tod_review.png", plot = g_res1, 
       path = "output/03_mov_patterns_discrete/", 
       width =  28, height = 12, units = "cm", dpi = 300)

ggsave("results_landuse_phase_tod2_review.png", plot = g_res1_flip, 
       path = "output/03_mov_patterns_discrete/", 
       width =  20, height = 10, units = "cm", dpi = 300)

# roads
plot(ggpredict(m8, c("dist_road_m", "tod_")))
plot(ggpredict(m8, c("dist_road_m", "phase")), lim = c(0, NA))
plot(ggpredict(m8, c("dist_road_m", "type"))) + ylim(0, 500)
plot(ggpredict(m8, c("dist_road_m", "type", "phase", "tod_"))) + ylim(0, 500)
plot(ggpredict(m8, c("dist_road_m", "type", "phase"))) + ylim(0, 500)
plot(ggpredict(m8, c("dist_road_m", "phase", "type"))) + ylim(0, 500)

res2 <- ggpredict(m8, c("dist_road_m", "type", "phase"), 
                  condition = c(sex = "Male"))
res2$x <- res2$x/1e3
res2$predicted <- res2$predicted * 24/1e3
res2$conf.low <- res2$conf.low * 24/1e3
res2$conf.high <- res2$conf.high * 24/1e3
# res2$predicted[res2$group == "resident" & res2$facet == "residency"] <- NA
# res2$predicted[res2$group == "resident" & res2$facet == "residency"] <- NA
# res2$predicted[res2$group == "resident" & res2$facet == "residency"] <- NA
g_res2 <- plot(res2, show_title = F) + 
  xlim(0, 10) + 
  ylim(0, 20) +
  labs(y = "Movement rate (km/day)",
       x = "Distance to roads (km)",
       colour = "Movement phase") +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols)# +
  # coord_cartesian(ylim = c(0, 20)) 
g_res2
ggsave("results_roads_phase_tod_review.png", plot = g_res2, 
       path = "output/03_mov_patterns_discrete/", 
       width =  20, height = 15, units = "cm", dpi = 300)

# plot(ggpredict(m8, c("dist_hidro_m", "tod_")))
# plot(ggpredict(m8, c("dist_hidro_m", "phase")))
# plot(ggpredict(m8, c("dist_hidro_m", "phase", "tod_")))

plot(ggpredict(m8, c("dist_urban_m", "phase")))
plot(ggpredict(m8, c("dist_urban_m", "tod_")), lim = c(0, NA))
plot(ggpredict(m8, c("dist_urban_m", "phase", "tod_", "land_use_FBDS")), lim = c(0, NA))
plot(ggpredict(m8, c("dist_urban_m", "type")))+ ylim(0, 2000)
plot(ggpredict(m8, c("dist_urban_m", "type", "phase")))+ ylim(0, 2000)

res3 <- ggpredict(m8, c("dist_urban_m", "type", "phase"),#, "tod_"), 
                  condition = c(sex = "Male"))#, tod_ = "day"))
res3$x <- res3$x/1e3
res3$predicted <- res3$predicted * 24/1e3
res3$conf.low <- res3$conf.low * 24/1e3
res3$conf.high <- res3$conf.high * 24/1e3
# res3$group <- res3$group |> 
#   forcats::fct_recode("resident" = "post-dispersal")
g_res3 <- plot(res3, lim = c(0, 30), show_title = F) + 
  xlim(0, 5) + 
  # ylim(0, 100) +
  labs(y = "Movement rate (km/day)",
       x = "Distance to urban areas (km)",
       colour = "Movement phase") +
  scale_color_manual(values = cols) + 
  scale_fill_manual(values = cols) +
  coord_cartesian(ylim = c(0, 20))
g_res3 
ggsave("results_urban_phase_tod_review.png", plot = g_res3, 
       path = "output/03_mov_patterns_discrete/", 
       width =  20, height = 15, units = "cm", dpi = 300)

g_res4 <- ggpubr::ggarrange(g_res2, g_res3, nrow = 2, labels = c("A", "B"),
                            common.legend = T, legend = "right")

ggsave("results_roads_urban_phase_tod_review.png", plot = g_res4, 
       path = "output/03_mov_patterns_discrete/", 
       width =  20, height = 18, units = "cm", dpi = 300)

# Joint figure
g_res4_1 <- ggpubr::ggarrange(g_res1_flip, 
                              g_res2, g_res3,
                              nrow = 3, labels = c("A", "B", "C"),
                              common.legend = T, legend = "right", align = "v")
g_res4_1

ggsave("results_lc_roads_urban_phase_tod_review.png", plot = g_res4_1, 
       path = "output/03_mov_patterns_discrete/", 
       width =  20, height = 22, units = "cm", dpi = 300)


# export coefficients
coefs <- broom::tidy(m8) %>% 
  dplyr::mutate(significant = as.factor(if_else(p.value < 0.05, 1, 0)),
                p = ifelse(p.value < 0.001, "<0.001", round(p.value, 3)))
coefs %>% 
  mutate_if(is.numeric, round, digits = 4) %>% 
  readr::write_csv(file = "output/03_mov_patterns_discrete/coef_table_movement_rate_review.csv")
View(coefs |> dplyr::filter(significant == 1))

summary(m8)

car::vif(m0_5)
car::vif(m1)
car::vif(m8)
plot(m8)

#-------------
# # using MuMiN

# install.packages("MuMIn", dep = T)
# library(MuMIn)

# m6
# 
# options(na_action = "na_fail")
# # na_action = "na_omit"
# 
# models <- MuMIn::dredge(m6)
# 
# sw(models)
# sw(subset(model_sel(models), delta <= 200))
# sw(model_avg(models, subset = delta <= 200))
# 
# MuMIn::importance(models, delta <= 200)

#-----------------
# fit model with random intercepts

# library(lme4)
# 
# m8 <- glmer(sl_ ~ sex + tod_ + phase +  
#             tod_:phase:land_use_FBDS + tod_:phase:scale(dist_road_m) +
#             phase:tod_:scale(dist_urban_m) + tod_:phase:scale(dist_hidro_m) +
#             (1|id),
#             data = data_fit, family = Gamma(link = "log"))
# m8

#-----------------
# fit each individuals

# m6 <- glm(sl_ ~ sex -1 + tod_:phase:land_use_FBDS + tod_:phase:log_dist_road +
#             phase:tod_:log_dist_urban + tod_:phase:log_dist_hidro, data = data_fit, family = Gamma)
# 
# ind_fit <- data_fit %>% 
#   dplyr::select(id, sl_, phase, tod_, land_use_FBDS, 
#                 dist_road_m, dist_urban_m, dist_hidro_m) %>% 
#   tidyr::nest(data = c(sl_, tod_, land_use_FBDS, dist_road_m, dist_urban_m, 
#                        dist_hidro_m)) %>% 
#   dplyr::mutate(fit = purrr::map(data, function(x) 
#     try(glm(sl_ ~ tod_:land_use_FBDS + tod_:log(dist_road_m + 1) +
#           tod_:log(dist_urban_m + 1) + tod_:log(dist_hidro_m + 1), 
#         data = x, family = Gamma(link = "log")))),
#     parms = purrr::map(fit, ~ try(broom::tidy(.))))
# 
# ind_fit_parms <- ind_fit %>%
#   dplyr::mutate(class = purrr::map(fit, ~ class(.)[1])) |> 
#   tidyr::unnest(class) |> 
#   dplyr::filter(class != "try-error") |> 
#   dplyr::select(id, phase, parms) %>% 
#   tidyr::unnest(parms)
# 
# ind_fit_parms %>% 
#   dplyr::filter(phase == "dispersal") %>%  
#   ggplot() +
#   geom_point(aes(estimate, term, color = id)) +
#   geom_linerange(aes(xmin = estimate - 1.96*std.error,
#                      xmax = estimate + 1.96*std.error,
#                      y = term, colour = id))
#   # xlim(-0.1, 0.1)
# 
# ind_fit_parms %>% 
#   dplyr::group_by(phase, term) %>% 
#   summarise(avg = mean(estimate),
#           inf = avg - 1.96*sd(estimate),
#           upp = avg + 1.96*sd(estimate)) %>% 
#   ggplot() +
#   geom_point(aes(avg, term, colour = phase),
#              position = position_dodge(0.5)) +
#   geom_linerange(aes(xmin = inf, xmax = upp, y = term, colour = phase),
#                  position = position_dodge(0.5)) +
#   xlim(-5, 5)
