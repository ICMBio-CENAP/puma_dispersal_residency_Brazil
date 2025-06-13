#' ---
#' title: "Understanding puma dispersal characteristics"
#' author: Bernardo Niebuhr
#' output:
#'   github_document: default
#'   html_document: default
#' ---

# Load packages
if(!require(install.load)) install.packages("install.load"); library(install.load)

install.load::install_load("ezknitr", "knitr")
install.load::install_load("tidyverse", "lubridate", "purrr", "ggpubr")
install.load::install_load("amt", "circular")
install.load::install_load("broom", "broom.mixed", "ggeffects")

# load data
library(PardasIPC)

# Print options for this document
options(width = 165)
opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) # root project folder
opts_chunk$set(eval = TRUE, error = F, message = F, warning = F, cache = F, echo = T, results = T)

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

# crs to use
crs_use <- sp::CRS("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m")

# transform data into track object
mov_track <- movement_data_pumas_behavior %>% 
  amt::mk_track(X, Y, timestamp, all_cols = T, crs = "epsg:4326") %>% 
  dplyr::mutate(phase = relevel(as.factor(phase), ref = "pre-dispersal"),
                selected_model = as.factor(selected_model)) %>% 
  dplyr::left_join(
    dplyr::select(pumas_metadata, 
                  id = name, dispersal_translocated, translocated),
  by = "id")
str(mov_track)

dispersers <- mov_track %>% 
  dplyr::filter(selected_model != "resident") %>% 
  pull(id) %>% 
  unique() %>% sort()

dispersers_non_translocated <- mov_track %>% 
  dplyr::filter(dispersal_translocated == "disperser") %>% 
  pull(id) %>% 
  unique() %>% sort()

dispersers_translocated <- mov_track %>% 
  dplyr::filter(dispersal_translocated == "translocated") %>% 
  pull(id) %>% 
  unique() %>% sort()

mov_track$phase %>% levels()
mov_track$phase <- forcats::fct_recode(
  mov_track$phase,
  "ranging" = "post-dispersal"
)

# ----- label=plots_individuals

# plot phases
# cols <- c("springgreen4", "steelblue", "firebrick")
# cols <- c("#009E73", "#0072B2", "#D55E00")
cols <- c("#117733", "#88CCEE", "#CC6677")

g1 <- ggplot(mov_track) + 
  geom_point(aes(x = x_, y = y_, col = phase)) +
  facet_wrap(~id, scales = "free", ncol = 3) +
  # coord_map() +
  theme_classic() +
  scale_color_manual(values = cols) +
  labs(x = "x", y = "y", color = "Movement phase")
g1
ggsave("individuals_dispersal_classified.png", plot = g1, path = "output/02_understand_dispersal/", 
       width = 20, height = 20, units = "cm", dpi = 300)

# in portuguese
# g1pt <- mov_track %>% 
#   dplyr::mutate(fase = factor(phase, labels = c("pré-dispersão", "dispersão", "pós-dispersão"))) %>% 
#   ggplot() + 
#   geom_point(aes(x = x_, y = y_, col = fase)) +
#   facet_wrap(~id, scales = "free", ncol = 3) +
#   # coord_map() +
#   theme_classic() +
#   scale_color_manual(values = cols) +
#   labs(x = "x", y = "y", color = "Fase do movimento")
# g1pt
# ggsave("individuals_dispersal_classified_pt.png", plot = g1pt, path = "output/02_understand_dispersal/", 
#        width = 20, height = 20, units = "cm", dpi = 300)

# plot models
g2 <- ggplot(mov_track) + 
  geom_point(aes(x = x_, y = y_, col = selected_model)) +
  facet_wrap(~id, scales = "free", ncol = 3) +
  # coord_map() +
  theme_classic() +
  labs(x = "x", y = "y", color = "Fitted behavior")
g2
ggsave("individuals_dispersal_models.png", plot = g2, path = "output/02_understand_dispersal/", 
       width = 20, height = 20, units = "cm", dpi = 300)

# transform coordinates
mov_track <- mov_track %>% 
  amt::transform_coords(crs_to = crs_use)

# calculate sl and ta
mov_track_parms <- mov_track %>% 
  tidyr::nest(data = c(x_, y_, t_, hDOP)) %>% 
  dplyr::mutate(
    data_rsp = map(data, function(x) {
      x %>% 
        amt::track_resample(rate = hours(1), tolerance = hours(1)) %>% 
        amt::filter_min_n_burst(min_n = 3) %>% 
        amt::steps_by_burst(keep_cols = "start") 
      })) %>% 
  tidyr::unnest(cols = data_rsp) %>% 
  dplyr::mutate(z = x1_ + 1i * y1_)

# some individual summaries
inds_summary <- mov_track %>% 
  dplyr::group_by(id) %>% 
  dplyr::summarise(
    sex = unique(sex),
    weight_kg = unique(weight_kg),
    age_months = unique(estimated_age_months),
    begin = min(t_),
    end = max(t_),
    days_monitored = as.integer(diff(range(t_))),
    n = n(),
    behavior = unique(selected_model),
    behav_disp_trans = unique(dispersal_translocated),
    translocated = unique(translocated))
inds_summary  

ind_summary_vals <- mov_track_parms %>% 
  dplyr::group_by(id, phase) %>%
  dplyr::summarise(
    dispersal_start = min(t1_),
    dispersal_end = max(t2_),
    dispersal_duration = as.integer(difftime(dispersal_end, dispersal_start, units = "days")),
    zi = first(z),
    zf = last(z),
    euclidean_distance_km = round(Mod(zf - zi)/1000, digits = 1),
    total_distance_km = round(sum(sl_)/1000, digits = 1),
    tot_dist_sl_km = round(sum(sl_)/1000, digits = 1))
ind_summary_vals

individual_info <- dplyr::left_join(inds_summary, ind_summary_vals, by = "id") %>% 
  dplyr::mutate(
    dispersal_start_day = as.numeric(floor(difftime(dispersal_start, begin - days(1), units = "days"))),
    dispersal_end_day = as.numeric(floor(difftime(dispersal_end, begin - days(1), units = "days"))),
    dispersal_age = round(age_months + dispersal_start_day/30, digits = 1)
  ) %>% 
  dplyr::select(-starts_with("z"))

# # mean and sd
# avg_row <- individual_info %>%
#   dplyr::mutate_if(is.numeric, mean, na.rm = T) %>%
#   slice(1)
# avg_row[1, c(1,2,5,6,9)] <- NA
# avg_row$id <- "avg"
# 
# sd_row <- individual_info %>%
#   dplyr::mutate_if(is.numeric, sd, na.rm = T) %>%
#   slice(1)
# sd_row[1, c(1,2,5,6,9)] <- NA
# sd_row$id <- "sd"
# 
# individual_info <- bind_rows(individual_info, avg_row, sd_row)

# save complete table
# save(individual_info, file = "output/02_understand_dispersal/individual_info_complete.rda")
# individual_info %>% 
#   readr::write_csv("output/02_understand_dispersal/individual_info_complete.csv")

# clean table
individual_info_short <- individual_info %>% 
  dplyr::filter(!(behavior != "resident" & phase != "dispersal"))

individual_info_short[individual_info_short$behavior == "resident", 11:18] <- NA

individual_info_short <- individual_info_short %>% 
  dplyr::mutate(behavior = if_else(behavior == "resident", "resident", "dispersed")) %>% 
  dplyr::select(-phase)

# save short table
individual_info_short$total_distance_km[individual_info_short$id == "Pepira"] <- NA
individual_info_short$tot_dist_sl_km[individual_info_short$id == "Pepira"] <- NA
# save(individual_info_short, file = "output/02_understand_dispersal/individual_info_concise.rda")
# individual_info_short %>% 
#   readr::write_csv("output/02_understand_dispersal/individual_info_concise.csv")

# avg parameters for all
(indiv_all <- individual_info_short %>% 
    dplyr::select(weight_kg, age_months, days_monitored, n))

# mean
apply(indiv_all, 2, mean)
# min
apply(indiv_all, 2, min)
# max
apply(indiv_all, 2, max)
# sum
apply(indiv_all, 2, sum)

# avg parameters for dispersers
(indiv_disps <- individual_info_short %>% 
  dplyr::filter(id %in% dispersers) %>% 
  # dplyr::filter(id != "Pepira") %>% 
  dplyr::select(weight_kg, contains("age"), dispersal_duration, contains("km")))

colMeans(indiv_disps, na.rm = T)
apply(indiv_disps, 2, min, na.rm = T)
apply(indiv_disps, 2, max, na.rm = T)

# avg parameters for dispersers which were not translocated
(indiv_disps_non_trans <- individual_info_short %>% 
    dplyr::filter(id %in% dispersers_non_translocated) %>% 
    # dplyr::filter(id != "Pepira") %>% 
    dplyr::select(weight_kg, contains("age"), dispersal_duration, contains("km")))

colMeans(indiv_disps_non_trans, na.rm = T)
apply(indiv_disps_non_trans, 2, min, na.rm = T)
apply(indiv_disps_non_trans, 2, max, na.rm = T)

# avg parameters for translocated dispersers
(indiv_disps_trans <- individual_info_short %>% 
    dplyr::filter(id %in% dispersers_translocated) %>% 
    # dplyr::filter(id != "Pepira") %>% 
    dplyr::select(weight_kg, contains("age"), dispersal_duration, contains("km")))

colMeans(indiv_disps_trans, na.rm = T)
apply(indiv_disps_trans, 2, min, na.rm = T)
apply(indiv_disps_trans, 2, max, na.rm = T)

# avg parameters for residents
(indiv_resid <- individual_info_short %>% 
    dplyr::filter(!(id %in% dispersers)) %>% 
    # dplyr::filter(id != "Pepira") %>% 
    dplyr::select(weight_kg, age_months, dispersal_duration, contains("km")))

colMeans(indiv_resid, na.rm = T)

# Females vs males
100*table(individual_info_short$sex)/nrow(individual_info_short)

# calculate step length and turning angle

# check sampling
mov_track %>% 
  amt::summarize_sampling_rate_many("id")
# one individual, Pepira, has gaps, that is why the mean is higher

#---
# Analysis considering cumulative distance within each day

# calculate step length and angle
mov_track_days <- mov_track_parms %>%
  dplyr::mutate(day = lubridate::as_date(t1_)) %>% 
  dplyr::group_by(id, phase, day) %>% 
  dplyr::summarise(daily_dist = sum(sl_)/1000,
                   n = n(),
                   mean_ta = ifelse(n > 15, mean(ta_, na.rm = T), NA))

# plot only dispersers
g3 <- mov_track_days %>% 
  dplyr::filter(id %in% dispersers) %>% 
  ggplot() + 
  # geom_violin(aes(x = phase, y = daily_dist)) +
  # geom_boxplot(aes(x = phase, y = daily_dist)) +
  geom_density(aes(x = daily_dist, color = phase, fill = phase), alpha = 0.2, size = 1.3) +
  # geom_point(aes(x = phase, y = mean(daily_dist))) +
  # facet_wrap(~ id, scales = "free_y", ncol = 2) +
  theme_classic() +
  labs(x = "Step length (km/day)", y = "Density", color = "Movement phase") +
  guides(fill = "none")
g3

g4 <- mov_track_days %>% 
  dplyr::filter(id %in% dispersers) %>% 
  ggplot() + 
  geom_density(aes(x = mean_ta, color = phase), size = 1.3) +
  # facet_wrap(~id, scales = "free_y", ncol = 2) +
  theme_classic() + 
  labs(x = "Turning angle", y = "Density", color = "Movement phase")
g4

# plot all individuals
g5 <- mov_track_days %>% 
  ggplot() + 
  #geom_violin(aes(x = dispersal_behavior, y = step_length)) +
  geom_density(aes(x = daily_dist, color = phase), size = 1.3) +
  facet_wrap(~id, scales = "free_y", ncol = 3) +
  # coord_map() +
  theme_classic() +
  scale_color_manual(values = cols) +
  labs(x = "Movement rate (km/day)", y = "Density", color = "Movement phase")
g5
ggsave("cumulative_distance_per_day_each.png", plot = g5, path = "output/02_understand_dispersal/", 
       width = 15, height = 18, units = "cm",
       dpi = 300)

g6 <- mov_track_days %>% 
  ggplot() + 
  geom_density(aes(x = mean_ta, color = phase), size = 1.3) +
  facet_wrap(~id, scales = "free_y", ncol = 4) +
  theme_classic() + 
  labs(x = "Turning angle", y = "Density", color = "Movement phase")
g6

#---
# Analysis considering one mean position per day

# do the same with one point per day
# calculate sl and ta
mov_track_parms_1d <- mov_track %>% 
  tidyr::nest(data = c(x_, y_, t_, hDOP)) %>% 
  dplyr::mutate(
    data_rsp = map(data, function(x) {
      x %>% 
        amt::track_resample(rate = hours(24), tolerance = hours(3)) %>% 
        amt::filter_min_n_burst(min_n = 3) %>% 
        amt::steps_by_burst(keep_cols = "start") 
    })) %>% 
  tidyr::unnest(cols = data_rsp) %>% 
  dplyr::mutate(z = x1_ + 1i * y1_,
                day = as_date(t1_))

mov_track_parms_1d <- mov_track %>% 
  tidyr::nest(data = c(x_, y_, t_, hDOP)) %>% 
  dplyr::mutate(
    data_rsp = map(data, function(x) {
      x %>% 
        amt::track_resample(rate = hours(24), tolerance = hours(3)) %>% 
        amt::filter_min_n_burst(min_n = 3) %>% 
        amt::steps_by_burst(keep_cols = "start")
    })) %>% 
  tidyr::unnest(cols = data_rsp) %>% 
  dplyr::mutate(z = x1_ + 1i * y1_,
                day = as_date(t1_))

mov_track_parms_1d <- mov_track %>% 
  tidyr::nest(data = c(x_, y_, t_, hDOP)) %>% 
  dplyr::mutate(
    data_rsp = map(data, function(x) {
      x %>% 
        amt::track_resample(rate = hours(24), tolerance = hours(3)) %>% 
        amt::filter_min_n_burst(min_n = 3) %>% 
        amt::steps_by_burst(keep_cols = "start") %>% 
        amt::time_of_day(where = "start")
    })) %>% 
  tidyr::unnest(cols = data_rsp) %>% 
  dplyr::mutate(z = x1_ + 1i * y1_,
                day = as_date(t1_))

# cols <- c("springgreen4", "steelblue", "firebrick")
# cols <- c("#009E73", "#0072B2", "#D55E00")
cols <- c("#117733", "#88CCEE", "#CC6677")

# plot only dispersers
mov_track_parms_1d$translocated %>% unique()
mov_track_parms_1d$id <- ifelse(mov_track_parms_1d$translocated == "translocated",
                                paste(mov_track_parms_1d$id, "(translocated)"),
                                mov_track_parms_1d$id)
mov_track_parms_1d$id %>% unique()
dispersers[c(1,2,3,6,8)] <- paste(dispersers[c(1,2,3,6,8)], "(translocated)") 

# speed
g3 <- mov_track_parms_1d %>% 
  dplyr::filter(id %in% dispersers) %>% 
  ggplot() + 
  # geom_density(aes(x = sl_/1000, color = phase), size = 1.3) +
  # geom_violin(aes(x = phase, y = sl_/1000)) +
  geom_boxplot(aes(x = phase, color = phase, y = sl_/1000), show.legend = F) +
  facet_wrap(~ id, scales = "fixed", ncol = 2) +
  theme_classic() +
  scale_color_manual(values = cols) +
  # labs(x = "Step length (km/day)", y = "Density", color = "Movement phase")
  labs(y = "Movement rate (km/day)", y = "")
g3
ggsave("speed_dispersers_each.png", plot = g3, path = "output/02_understand_dispersal/", 
       width = 15, height = 18, units = "cm",
       dpi = 300)

# pooled translocated vs non-translocated dispersers
g3_trans_non_trans <- mov_track_parms_1d %>% 
  dplyr::filter(id %in% dispersers) %>% 
  ggplot() + 
  # geom_density(aes(x = sl_/1000, color = phase), size = 1.3) +
  # geom_violin(aes(x = phase, y = sl_/1000)) +
  geom_boxplot(aes(x = phase, col = phase, y = sl_/1000), show.legend = F) +
  theme_classic() +
  scale_color_manual(values = cols) +
  facet_wrap(~ translocated, scales = "fixed", ncol = 2) +
  # labs(x = "Step length (km/day)", y = "Density", color = "Movement phase")
  labs(y = "Movement rate (km/day)", y = "")
g3_trans_non_trans
ggsave("speed_dispersers_trans_vs_nontrans.png", plot = g3_trans_non_trans, path = "output/02_understand_dispersal/", 
       width = 18, height = 15, units = "cm",
       dpi = 300)

# pooled all dispersers
g3_pooled <- mov_track_parms_1d %>% 
  dplyr::filter(id %in% dispersers) %>% 
  ggplot() + 
  # geom_density(aes(x = sl_/1000, color = phase), size = 1.3) +
  # geom_violin(aes(x = phase, y = sl_/1000)) +
  geom_boxplot(aes(x = phase, col = phase, y = sl_/1000), show.legend = F) +
  theme_classic() +
  scale_color_manual(values = cols) +
  # labs(x = "Step length (km/day)", y = "Density", color = "Movement phase")
  labs(y = "Movement rate (km/day)", y = "")
g3_pooled
ggsave("speed_dispersers_polled.png", plot = g3_pooled, path = "output/02_understand_dispersal/", 
       width = 12, height = 12, units = "cm",
       dpi = 300)

# turning angles
mov_track_parms_1d_plot <- mov_track_parms_1d %>% 
  dplyr::filter(id %in% dispersers) %>% 
  dplyr::mutate(x = ta_/pi*180)
dens_df <- mov_track_parms_1d_plot %>%
  dplyr::filter(!(id == "Pepira" & phase == "dispersal")) %>% 
  dplyr::group_by(id, phase) %>%
  dplyr::summarise(density = list(density.circular(circular(x, units = "degrees"), bw = 20, 
                                            na.rm = TRUE, from = circular(-pi), to = circular(pi))), 
            .groups = "drop") %>%
  dplyr::mutate(density_df = lapply(density, function(d) data.frame(x = as.numeric(d$x), y = d$y))) %>%
  tidyr::unnest(density_df) %>% 
  dplyr::select(-density)
g4 <- mov_track_parms_1d_plot %>%
  ggplot() + 
  # geom_density(aes(x = ta_/pi*180, color = phase), size = 1.3) +
  geom_smooth(data = dens_df, aes(x = x, y = y, color = phase),
              formula = y ~ s(x, bs = "cc", k = 5), method = "gam", size = 1.3,
              se = FALSE) +
  facet_wrap(~id, scales = "free_y", ncol = 2) +
  theme_classic() + 
  scale_color_manual(values = cols) +
  labs(x = "Turning angle", y = "Density", color = "Movement phase")
g4
ggsave("turning_angle_dispersers_each.png", plot = g4, path = "output/02_understand_dispersal/", 
       width = 15, height = 18, units = "cm",
       dpi = 300)

# pooled translocated vs non-translocated
dens_df <- mov_track_parms_1d_plot %>%
  group_by(translocated, phase) %>%
  summarise(density = list(density.circular(circular(x, units = "degrees"), bw = 20, 
                                            na.rm = TRUE, from = circular(-pi), to = circular(pi))), 
            .groups = "drop") %>%
  mutate(density_df = lapply(density, function(d) data.frame(x = as.numeric(d$x), y = d$y))) %>%
  unnest(density_df) %>% 
  select(-density)
# dens_df %>% summary()
# dens_df$x
g4_trans_nontrans <- ggplot() +
  # geom_histogram(data = data = mov_track_parms_1d_plot, aes(x = x, y = ..density.., fill = phase), size = 1.3) +
  geom_smooth(data = dens_df, aes(x = x, y = y, color = phase),
              formula = y ~ s(x, bs = "cc", k = 10), method = "gam", size = 1.3,
              se = FALSE) +
  # geom_line(data = dens_df, aes(x = as.numeric(x), y = y, color = phase)) +
  # geom_density(aes(x = x, color = phase), size = 1.3) +
  facet_wrap(~translocated, scales = "free_y", ncol = 2) +
  theme_classic() + 
  scale_color_manual(values = cols) +
  labs(x = "Turning angle", y = "Density", color = "Movement phase")
g4_trans_nontrans
ggsave("turning_angle_dispersers_trans_vs_nontrans.png", plot = g4_trans_nontrans, path = "output/02_understand_dispersal/", 
       width = 15, height = 18, units = "cm",
       dpi = 300)

# pooled all
dens_df <- mov_track_parms_1d_plot %>%
  group_by(phase) %>%
  summarise(density = list(density.circular(circular(x, units = "degrees"), bw = 20, 
                                            na.rm = TRUE, from = circular(-pi), to = circular(pi))), 
            .groups = "drop") %>%
  mutate(density_df = lapply(density, function(d) data.frame(x = as.numeric(d$x), y = d$y))) %>%
  unnest(density_df) %>% 
  select(-density)
# dens_df %>% summary()
# dens_df$x
g4_pooled <- ggplot() +
  # geom_histogram(data = data = mov_track_parms_1d_plot, aes(x = x, y = ..density.., fill = phase), size = 1.3) +
  geom_smooth(data = dens_df, aes(x = x, y = y, color = phase),
              formula = y ~ s(x, bs = "cc", k = 10), method = "gam", size = 1.3,
              se = FALSE) +
  theme_classic() + 
  scale_color_manual(values = cols) +
  labs(x = "Turning angle", y = "Density", color = "Movement phase")
g4_pooled
ggsave("turning_angle_dispersers_polled.png", plot = g4_pooled, path = "output/02_understand_dispersal/", 
       width = 12, height = 12, units = "cm",
       dpi = 300)

# plot all individuals

# speed
g5 <- mov_track_parms_1d %>% 
  ggplot() + 
  # geom_density(aes(x = sl_/1000, color = phase), size = 1.3) +
  # geom_violin(aes(x = phase, y = sl_/1000)) +
  geom_boxplot(aes(x = phase, color = phase, y = sl_/1000), show.legend = F) +
  facet_wrap(~id, scales = "fixed", ncol = 3) +
  # coord_map() +
  theme_classic() +
  scale_color_manual(values = cols) +
  labs(y = "Movement rate (km/day)", x = "Movement phase")
g5
ggsave("speed_allind_each.png", plot = g5, path = "output/02_understand_dispersal/", 
       width = 18, height = 22, units = "cm",
       dpi = 300)

g5_trans_res_nontrans <- mov_track_parms_1d %>% 
  ggplot() + 
  # geom_density(aes(x = sl_/1000, color = phase), size = 1.3) +
  # geom_violin(aes(x = phase, y = sl_/1000)) +
  geom_boxplot(aes(x = phase, color = phase, y = sl_/1000), show.legend = F) +
  theme_classic() +
  facet_wrap(~dispersal_translocated, scales = "fixed", ncol = 3) +
  scale_color_manual(values = cols) +
  labs(y = "Movement rate (km/day)", x = "Movement phase")
g5_trans_res_nontrans
ggsave("speed_allind_trans_nontrans.png", plot = g5_trans_res_nontrans, path = "output/02_understand_dispersal/", 
       width = 12, height = 12, units = "cm",
       dpi = 300)

# all individuals pooled
g5_pooled <- mov_track_parms_1d %>% 
  ggplot() + 
  # geom_density(aes(x = sl_/1000, color = phase), size = 1.3) +
  # geom_violin(aes(x = phase, y = sl_/1000)) +
  geom_boxplot(aes(x = phase, color = phase, y = sl_/1000), show.legend = F) +
  theme_classic() +
  scale_color_manual(values = cols) +
  labs(y = "Movement rate (km/day)", x = "Movement phase")
g5_pooled
ggsave("speed_allind_pooled.png", plot = g5_pooled, path = "output/02_understand_dispersal/", 
       width = 12, height = 12, units = "cm",
       dpi = 300)

# turning angles
dens_df <- mov_track_parms_1d %>%
  dplyr::mutate(x = ta_/pi*180) %>% 
  dplyr::filter(!(id == "Pepira" & phase == "dispersal")) %>% 
  group_by(id, phase) %>%
  summarise(density = list(density.circular(circular(x, units = "degrees"), bw = 20, 
                                            na.rm = TRUE, from = circular(-pi), to = circular(pi))), 
            .groups = "drop") %>%
  mutate(density_df = lapply(density, function(d) data.frame(x = as.numeric(d$x), y = d$y))) %>%
  unnest(density_df) %>% 
  select(-density)
g6 <- ggplot() + 
  geom_smooth(data = dens_df, aes(x = x, y = y, color = phase),
              formula = y ~ s(x, bs = "cc", k = 5), method = "gam", size = 1.3,
              se = FALSE) +  
  facet_wrap(~id, scales = "free_y", ncol = 3) +
  theme_classic() + 
  scale_color_manual(values = cols) +
  labs(x = "Turning angle", y = "Density", color = "Movement phase")
g6
ggsave("turning_angle_allind_each.png", plot = g6, path = "output/02_understand_dispersal/", 
       width = 18, height = 22, units = "cm",
       dpi = 300)

dens_df <- mov_track_parms_1d %>%
  dplyr::mutate(x = ta_/pi*180) %>% 
  dplyr::filter(!(id == "Pepira" & phase == "dispersal")) %>% 
  group_by(dispersal_translocated, phase) %>%
  summarise(density = list(density.circular(circular(x, units = "degrees"), bw = 20, 
                                            na.rm = TRUE, from = circular(-pi), to = circular(pi))), 
            .groups = "drop") %>%
  mutate(density_df = lapply(density, function(d) data.frame(x = as.numeric(d$x), y = d$y))) %>%
  unnest(density_df) %>% 
  select(-density)
g6_trans_nontrans <- ggplot() + 
  # geom_density(aes(x = ta_/pi*180, color = phase), size = 1.3) +
  geom_smooth(data = dens_df, aes(x = x, y = y, color = phase),
              formula = y ~ s(x, bs = "cc", k = 5), method = "gam", size = 1.3,
              se = FALSE) +
  facet_wrap(~dispersal_translocated, scales = "fixed", ncol = 3) +
  theme_classic() + 
  scale_color_manual(values = cols) +
  labs(x = "Turning angle", y = "Density", color = "Movement phase")
g6_trans_nontrans
ggsave("turning_angle_trans_vs_nontrans.png", plot = g6_trans_nontrans, path = "output/02_understand_dispersal/", 
       width = 18, height = 22, units = "cm",
       dpi = 300)

dens_df <- mov_track_parms_1d %>%
  dplyr::mutate(x = ta_/pi*180) %>% 
  dplyr::filter(!(id == "Pepira" & phase == "dispersal")) %>% 
  group_by(phase) %>%
  summarise(density = list(density.circular(circular(x, units = "degrees"), bw = 20, 
                                            na.rm = TRUE, from = circular(-pi), to = circular(pi))), 
            .groups = "drop") %>%
  mutate(density_df = lapply(density, function(d) data.frame(x = as.numeric(d$x), y = d$y))) %>%
  unnest(density_df) %>% 
  select(-density)
g6_pooled <- ggplot() + 
  # geom_density(aes(x = ta_/pi*180, color = phase), size = 1.3) +
  geom_smooth(data = dens_df, aes(x = x, y = y, color = phase),
              formula = y ~ s(x, bs = "cc", k = 5), method = "gam", size = 1.3,
              se = FALSE) +
  theme_classic() + 
  scale_color_manual(values = cols) +
  labs(x = "Turning angle", y = "Density", color = "Movement phase")
g6_pooled
ggsave("turning_angle_allind_pooled.png", plot = g6_pooled, path = "output/02_understand_dispersal/", 
       width = 12, height = 12, units = "cm",
       dpi = 300)

# g6_pooled_pt <- mov_track_parms_1d %>% 
#   dplyr::mutate(fase = factor(phase, labels = c("pré-dispersão", "dispersão", "pós-dispersão"))) %>% 
#   ggplot() + 
#   geom_density(aes(x = ta_/pi*180, color = fase), size = 1.3) +
#   theme_classic() + 
#   scale_color_manual(values = cols) +
#   labs(x = "Ângulo de virada", y = "Densidade", color = "Fase do movimento")
# g6_pooled_pt
# ggsave("turning_angle_allind_pooled_pt.png", plot = g6_pooled_pt, path = "output/02_understand_dispersal/", 
#        width = 12, height = 12, units = "cm",
#        dpi = 300)

# combined figure for all inds combined
g_comb <- ggarrange(g5_pooled +
                      labs(x = ""),
            # theme(axis.title.x=element_blank(),
            #       axis.text.x=element_blank()), 
          g6_pooled,
          labels = c("A", "B"),
          common.legend = TRUE, legend = "bottom")
g_comb
ggsave("sl_ta_all_polled.png", plot = g_comb, path = "output/02_understand_dispersal/", 
       width = 16, height = 10, units = "cm",
       dpi = 300)

# combined figure for all trans vs non-trans
g_comb <- ggarrange(g5_trans_res_nontrans +
                      labs(x = ""),
                    # theme(axis.title.x=element_blank(),
                    #       axis.text.x=element_blank()), 
                    g6_trans_nontrans,
                    ncol = 1,
                    labels = c("A", "B"),
                    common.legend = TRUE, legend = "bottom")
g_comb
ggsave("sl_ta_all_polled_trans_vs_nontrans.png", plot = g_comb, path = "output/02_understand_dispersal/", 
       width = 20, height = 20, units = "cm",
       dpi = 300)

# models

# Movement rates

# adding 1 to sl to fit Gamma distribution
# working with km instead of m
mov_track_parms_1d <- mov_track_parms_1d %>% 
  dplyr::mutate(sl = (sl_ + 1)/1000)

# mov_track_parms <- mov_track_parms %>% 
#   dplyr::mutate(sl = (sl_ + 1)/1000)

# mov_track_days <- mov_track_days %>% 
#   dplyr::mutate(sl = daily_dist)

# only random effects
(m00 <- lme4::glmer(sl ~ (1|id), data = mov_track_parms_1d,
                   family = Gamma(link = "log")))
summary(m00)

# fixed effects only
(m0 <- glm(sl ~ phase - 1, data = mov_track_parms_1d, family = Gamma(link = "log")))
summary(m0)

# individuals as random intercepts
# relevel(phase, ref = "dispersal")
(m1 <- lme4::glmer(sl ~ phase - 1 + (1|id), data = mov_track_parms_1d,
                   family = Gamma(link = "log")))
summary(m1)

# fixed effects only - variation across dispersers, residents, and translocated
(m2 <- glm(sl ~ phase*dispersal_translocated - 1, data = mov_track_parms_1d, 
           family = Gamma(link = "log")))
summary(m2)

# individuals as random intercepts
# relevel(phase, ref = "dispersal")
(m3 <- lme4::glmer(sl ~ phase*dispersal_translocated - 1 + (1|id), data = mov_track_parms_1d,
                   family = Gamma(link = "log")))
summary(m3)

# random effects
bbmle::AICctab(m00, m0, m1, m2, m3, delta = T, base = T, weights = TRUE)

# average and sd
coef(m1)
broom.mixed::tidy(m1)
(fixed_params <- broom.mixed::tidy(m1, effects = "fixed")[,c("term", "estimate")])
(random_params <- broom.mixed::tidy(m1, effect = "ran_pars", ran_prefix = T))

# without accounting for random effects
pr <- ggpredict(m1, "phase")
# pr <- ggemmeans(m1, "phase")
pr
plot(pr)

# accounting for random effects
pr <- ggpredict(m1, terms =  c("phase"), type = "fixed", ci_level = 0.95, back.transform = T)
pr

g7 <- plot(pr, colors = "grey")# +
  # ylim(0, 8)
g7

g7 <- pr %>% 
  as.data.frame() %>% 
  ggplot() +
  geom_point(aes(x, predicted, color = x), show.legend = F, size = 3) +
  geom_point(data = mov_track_parms_1d, aes(jitter(as.numeric(phase), 1), sl_/1000), color = "grey",
             alpha = 0.2, inherit.aes = FALSE) +
  geom_point(aes(x, predicted, color = x), show.legend = F, size = 3) +
  geom_linerange(aes(x, ymin = conf.low, ymax = conf.high, color = x), 
                 show.legend = F, size = 1.2) +
  theme_classic() +
  scale_color_manual(values = cols) +
  ylim(0, quantile(mov_track_parms_1d$sl, 0.95)) +
  labs(y = "Movement rate (km/day)", x = "")
g7

# combined figure for all inds combined
g_comb2 <- ggarrange(g7 +
                      labs(x = ""),
                    # theme(axis.title.x=element_blank(),
                    #       axis.text.x=element_blank()), 
                    g6_pooled,
                    labels = c("A", "B"),
                    common.legend = TRUE, legend = "bottom")
g_comb2
ggsave("sl_ta_all_polled_estimated.png", plot = g_comb2, path = "output/02_understand_dispersal/", 
       width = 16, height = 10, units = "cm",
       dpi = 300)

#---------------
# just for comparison, if we want to compute the difference between translocated 
# and non-translocated animals

# accounting for random effects
pr_tr <- ggpredict(m3, terms =  c("dispersal_translocated", "phase"), type = "fixed", ci_level = 0.95, back.transform = T)
pr_tr

g7 <- plot(pr_tr, colors = "grey")# +
# ylim(0, 8)
g7

g7_trans <- pr_tr %>% 
  as.data.frame() %>% 
  dplyr::slice(-c(4,5)) %>% 
  ggplot() +
  geom_point(aes(x, predicted, color = x), show.legend = F, size = 3) +
  # geom_point(data = mov_track_parms_1d, aes(jitter(as.numeric(phase), 1), sl_/1000), color = "grey",
  #            alpha = 0.2, inherit.aes = FALSE) +
  geom_point(aes(x, predicted, color = x), show.legend = F, size = 3) +
  geom_linerange(aes(x, ymin = conf.low, ymax = conf.high, color = x), 
                 show.legend = F, size = 1.2) +
  theme_classic() +
  scale_color_manual(values = cols) +
  ylim(0, quantile(mov_track_parms_1d$sl, 0.95)) +
  labs(y = "Movement rate (km/day)", x = "") +
  facet_wrap(~group)
g7_trans

# # in pt
# g7_pt <- pr %>% 
#   dplyr::mutate(fase = factor(x, labels = c("pré-dispersão", "dispersão", "pós-dispersão"))) %>% 
#   as.data.frame() %>% 
#   ggplot() +
#   geom_point(aes(fase, predicted, color = fase), show.legend = F, size = 3) +
#   geom_point(data = mov_track_parms_1d, aes(jitter(as.numeric(phase), 1), sl_/1000), color = "grey",
#              alpha = 0.2, inherit.aes = FALSE) +
#   geom_point(aes(fase, predicted, color = fase), show.legend = F, size = 3) +
#   geom_linerange(aes(fase, ymin = conf.low, ymax = conf.high, color = fase), 
#                  show.legend = F, size = 1.2) +
#   theme_classic() +
#   scale_color_manual(values = cols) +
#   ylim(0, 6) +
#   labs(y = "Taxa de movimento (km/dia)", x = "")
# g7_pt

# g_comb2_pt <- ggarrange(g7_pt +
#                        labs(x = ""),
#                      # theme(axis.title.x=element_blank(),
#                      #       axis.text.x=element_blank()), 
#                      g6_pooled_pt,
#                      labels = c("A", "B"),
#                      common.legend = TRUE, legend = "bottom")
# g_comb2_pt
# ggsave("sl_ta_all_polled_estimated_pt.png", plot = g_comb2_pt, path = "output/02_understand_dispersal/", 
#        width = 16, height = 10, units = "cm",
#        dpi = 300)


# for the dispersers
(m0_disp <- lme4::glmer(sl ~ (1|id), 
                        data = filter(mov_track_parms_1d, id %in% dispersers),
                        family = "Gamma"))

(m1_disp <- lme4::glmer(sl ~ phase - 1 + (1|id), 
                        data = filter(mov_track_parms_1d, id %in% dispersers),
                        family = "Gamma"))

(m3_disp <- lme4::glmer(sl ~ phase*translocated - 1 + (1|id), 
                        data = filter(mov_track_parms_1d, id %in% dispersers),
                        family = "Gamma"))

bbmle::AICctab(m0_disp, m1_disp, m3_disp, delta = T, base = T)

ggpredict(m1_disp, c("phase"), type = "random", back.transform = T)
ggpredict(m3_disp, c("phase", "translocated"), type = "fixed", back.transform = T) %>% 
  plot()

pr_r <- ggpredict(m1_disp, c("phase", "id"), type = "random", back.transform = T)
pr_r

g8 <- plot(pr_r) +
  labs(y = "Movement rate (km/day)") +
  ylim(0,6)
g8

g7_r <- pr %>% 
  as.data.frame() %>% 
  ggplot() +
  geom_point(aes(x, predicted), size = 0) +
  geom_rect(mapping = aes(xmin = as.numeric(x) - .4, xmax = as.numeric(x) + .4, 
                          ymin = conf.low, ymax = conf.high), 
            fill = "grey90") +
  geom_segment(mapping = aes(x = as.numeric(x) - .4, xend = as.numeric(x) + .4,
                             y = predicted, yend = predicted), size = 1) +
  geom_pointrange(data = as.data.frame(pr_r), aes(jitter(as.numeric(x), 1), predicted,
                                                  ymin = conf.low, ymax = conf.high,
                                                  color = group)) +
  theme_classic() +
  #ylim(0, 6) +
  labs(y = "Movement rate (km/day)", x = "", color = "Individual")
g7_r
ggsave("sl_estimated_main_random_effects.png", plot = g7_r, path = "output/02_understand_dispersal/", 
       width = 15, height = 12, units = "cm",
       dpi = 300)

# simulate
# ggpredict(m1, c("phase"), type = "random", ci.lvl = 95)
ggpredict(m1, c("phase", "id"), type = "random", ci.lvl = 95)

#---
# for the circular parameters
library(circular)

vals <- list()
for(i in levels(mov_track_parms_1d$phase)) {
  dat <- mov_track_parms_1d %>% dplyr::filter(phase == i)
  
  mle_vals <- list()
  mle_vals$mu <- circular::mle.vonmises(dat$ta_)$mu
  mle_vals$mu.ci <- circular::mle.vonmises.bootstrap.ci(dat$ta_)$mu.ci
  
  mle_vals$kappa <- circular::mle.vonmises(dat$ta_)$kappa
  mle_vals$kappa.ci <- circular::mle.vonmises.bootstrap.ci(dat$ta_)$kappa.ci
  
  vals[[i]] <- mle_vals
}
as.data.frame(unlist(vals))

# png("output/turning_angles_fitted.png", width = 15, height = 15, units = "cm",
#     res = 300)
plot(0, 0, type = "n", xlim = c(-pi, pi), ylim = c(0.1, 0.23),
     xlab = "Turning anlge", ylab = "Density")
mapply(function(a, col) curve(dvonmises(x, mu = a$mu, kappa = a$kappa),
                         from = -pi, to = pi, col = col, add = T, lwd = 2),
       a = vals, col = cols)
legend("topright", legend = levels(mov_track_parms_1d$phase), col = cols, lwd = 2)
# dev.off()

(m0_a <- with(mov_track_parms_1d, 
            lm.circular(y = ta_, 
                        x = c(1,2,3)[phase], 
                        init = 0.5,
                        type = "c-l")))
summary(m0_a)
m0_a$mu
m0_a$coefficients
m0_a$mu + c(1,2,3)*m0_a$coefficients

curve(circular::dvonmises(x, mu = m0_a$mu + m0_a$coefficients, kappa = m0_a$kappa),
      from = -pi, to = pi)
curve(circular::dvonmises(x, mu = m0_a$mu + 2*m0_a$coefficients, kappa = m0_a$kappa),
      from = -pi, to = pi, add = T, col = 2)


# install.packages("circglmbayes", dep = T)
# library(circglmbayes)
# 
# llvonmises <- function(x, a, b, c, d, e, f, phase) {
#   m = c(a, b, c)[phase]
#   kapp = c(d, e, f)[phase]
#   
#   ll <- c()
#   for(i in 1:length(x)) {
#     ll1 <- circular::dvonmises(x[i], mu = m[i], kappa = kapp[i], log = T)
#     ll <- c(ll, ll1)
#   }
#   -sum(ll)
# }
# 
# mu_guess <- mle.vonmises(mov_track_parms_1d$ta_)$mu
# kappa_guess <- mle.vonmises(mov_track_parms_1d$ta_)$kappa

# bbmle::mle2(llvonmises, start = list(a = mu.guess, b = mu.guess, c = mu.guess, 
#                                      d = kappa.guess, e = kappa.guess, f = kappa.guess),
#             data = list(x = mov.track.parms.1d$ta_, phase = mov.track.parms.1d$phase),
#             method = "L-BFGS-B", lower = c(d = 0, e = 0, f = 0))