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

# Print options for this document
options(width = 165)
opts_knit$set(root_dir = rprojroot::find_rstudio_root_file()) # root project folder
opts_chunk$set(eval = F, error = F, message = F, warning = F, cache = F, echo = T, results = T)

# --------------- label=setup
# Set up 

# Clean everything before beginning
rm(list = ls())

# --------------- label=load_data_and_functions

# our data
load("output/02_understand_dispersal/individual_info_concise.rda")
individual_info_short |> 
  print(width = Inf)
ind_data <- individual_info_short

ind_data_sum <- ind_data |> 
  dplyr::filter(!is.na(euclidean_distance_km)) |> 
  dplyr::group_by(sex) |> 
  dplyr::summarise(individual = NA,
                   translocated = "yes",
                   mean_dispersal_age_months = mean(dispersal_age),
                   sd_dispersal_age_months = sd(dispersal_age),
                   min_dispersal_age_months = min(dispersal_age),
                   max_dispersal_age_months = max(dispersal_age),
                   dispersal_duration_days = mean(dispersal_duration),
                   sd_dispersal_duration_days = sd(dispersal_duration),
                   dispersal_mass_kg = mean(weight_kg),
                   dispersal_Euclidean_distance_km = mean(euclidean_distance_km),
                   sd_dispersal_Euc_distance_km = sd(euclidean_distance_km),
                   min_dispersal_Euc_distance_km = min(euclidean_distance_km),
                   max_dispersal_Euc_distance_km = max(euclidean_distance_km),
                   dispersal_total_distance_km = mean(total_distance_km, na.rm = TRUE),
                   min_dispersal_total_distance_km = min(total_distance_km, na.rm = TRUE),
                   max_dispersal_total_distance_km = max(total_distance_km, na.rm = TRUE),
                   ecoregion = "Atlantic Forest",
                   state = "SP",
                   country = "Brazil",
                   telemetry_method = "GPS",
                   n = n(),
                   data = "mean,range",
                   year = 2022, 
                   study = "Our study",
                   source = "Our study") |> 
  dplyr::mutate(sex = ifelse(sex == "Male", "M", "F"))
  
# literature
lit_disp <- readr::read_csv("text/literature_puma_data_dispersal.csv")
lit_disp |> print(width = Inf)

colnames(ind_data)
colnames(lit_disp)
# View(lit_disp)

# single
dat_all <- lit_disp |> 
  dplyr::filter(keep_analysis == 1) |> 
  dplyr::mutate(source = "literature",
                mean_dispersal_age_months = as.numeric(mean_dispersal_age_months),
                dispersal_duration_days = as.numeric(dispersal_duration_days)) |>
  dplyr::bind_rows(ind_data_sum) #|> 
  # dplyr::filter(dispersal_Euclidean_distance_km < 2000)

#---------
# General

# Number of studies
lit_disp |> 
  dplyr::count(study) |> 
  pull(study)
# 24 studies (28 total, 4 with two sampling sites)

# where
lit_disp |> 
  dplyr::group_by(study) |> 
  dplyr::summarise(country = unique(country)) |> 
  pull(country) |> 
  table()
# 24 studies/28 areas (4 Canada, 1 Chile, 23 areas US, 19 studies US)
23/24 # prop in North America

# monitoring method
(met <- lit_disp |> 
  dplyr::group_by(study) |> 
  dplyr::summarise(tel = unique(telemetry_method))) |> 
  print(n = 100)
table(met$tel)
table(met$tel)/(nrow(met)-4)

# method to identify dispersal
(met2 <- lit_disp |> 
    dplyr::group_by(study) |> 
    dplyr::summarise(met = unique(method_dispersal))) |> 
  print(n = 100)
table(met$tel)
table(met$tel)/(nrow(met)-4)

#---------
# Some plots

#-----
# Euclidean distance - hist

# mean
dat_all$data
ed_dat <- dat_all |> 
  dplyr::filter(!grepl("Hawley|Our study", study)) |> 
  dplyr::slice(grep("mean|raw", dat_all$data))

summary(ed_dat)
ed <- ind_data_sum$dispersal_Euclidean_distance_km
sum(ed_dat$dispersal_Euclidean_distance_km < ed)/(nrow(ed_dat))

ed_dat_m <- dat_all |> 
  dplyr::filter(!grepl("Hawley|Our study", study), sex == "M") |> 
  dplyr::slice(grep("mean|raw", dat_all$data))

summary(ed_dat_m)
ed <- ind_data_sum$dispersal_Euclidean_distance_km
sum(ed_dat_m$dispersal_Euclidean_distance_km < ed)/(nrow(ed_dat_m))

g_ed_mean <-  dat_all |> 
  dplyr::filter(!grepl("Hawley|Our study", study)) |> 
  dplyr::slice(grep("mean|raw", dat_all$data)) |> 
  ggplot(aes(x = dispersal_Euclidean_distance_km)) +
  geom_histogram() +
  geom_vline(aes(xintercept = dispersal_Euclidean_distance_km), data = ind_data_sum, 
             col = 2, size = 1.5) +
  # facet_wrap(~sex) +
  labs(x = "Mean Euclidean dispersal distance (km)",
       y = "Number of studies") +
  xlim(0, 1100) +
  ggeffects::theme_ggeffects()
g_ed_mean

# min
dat_all |> 
  dplyr::slice(grep("min|range|raw", dat_all$data)) |> 
  ggplot(aes(x = min_dispersal_Euc_distance_km)) +
  geom_histogram() +
  geom_vline(aes(xintercept = min_dispersal_Euc_distance_km), data = ind_data_sum, 
             col = 2, size = 1.5) +
  labs(x = "Minimum Euclidean dispersal distance (km)",
       y = "Number of studies") +
  ggeffects::theme_ggeffects()

# max
ed_dat_max <- dat_all |> 
  dplyr::filter(!grepl("Hawley|Our study", study)) |> 
  dplyr::slice(grep("max|range|raw", dat_all$data))

summary(ed_dat_max)
ed_max <- ind_data_sum$max_dispersal_Euc_distance_km
sum(ed_dat_max$max_dispersal_Euc_distance_km < ed_max, na.rm = TRUE)/(nrow(ed_dat_max))

ed_dat_max_m <- dat_all |> 
  dplyr::filter(!grepl("Hawley|Our study", study), sex == "M") |> 
  dplyr::slice(grep("max|range|raw", dat_all$data))

summary(ed_dat_max_m)
ed <- ind_data_sum$max_dispersal_Euc_distance_km
sum(ed_dat_max_m$max_dispersal_Euc_distance_km < ed, na.rm = TRUE)/(nrow(ed_dat_max_m))

g_ed_max <- dat_all |> 
  dplyr::filter(!grepl("Hawley|Our study", study)) |> 
  dplyr::slice(grep("max|range|raw", dat_all$data)) |> 
  ggplot(aes(x = max_dispersal_Euc_distance_km)) +
  geom_histogram() +
  geom_vline(aes(xintercept = max_dispersal_Euc_distance_km), data = ind_data_sum, 
             col = 2, size = 1.5) +
  labs(x = "Maximum Euclidean dispersal distance (km)",
       y = "Number of studies") +
  ylim(0, 7) +
  xlim(0, 1100) +
  ggeffects::theme_ggeffects()
g_ed_max

g_ed_hist <- ggpubr::ggarrange(g_ed_mean, g_ed_max + ylab(""), ncol = 2, labels = c("A", "B"))
g_ed_hist

ggsave("euc_dist_hist.png", plot = g_ed_hist, path = "output/05_literature/", 
       width = 20, height = 10, units = "cm", dpi = 300)

# range of values
range(ed_dat$dispersal_Euclidean_distance_km, na.rm = T)
range(ed_dat_max$max_dispersal_Euc_distance_km, na.rm = T)

# Only males
g_ed_mean_m <-  dat_all |> 
  dplyr::filter(!grepl("Hawley|Our study", study), sex == "M") |> 
  dplyr::slice(grep("mean|raw", dat_all$data)) |> 
  ggplot(aes(x = dispersal_Euclidean_distance_km)) +
  geom_histogram() +
  geom_vline(aes(xintercept = dispersal_Euclidean_distance_km), data = ind_data_sum, 
             col = 2, size = 1.5) +
  # facet_wrap(~sex) +
  labs(x = "Mean Euclidean dispersal distance (km)",
       y = "Number of studies") +
  xlim(0, 1100) +
  ggeffects::theme_ggeffects()
g_ed_mean_m

# max
g_ed_max_m <- dat_all |> 
  dplyr::filter(!grepl("Hawley|Our study", study), sex == "M") |> 
  dplyr::slice(grep("max|range|raw", dat_all$data)) |> 
  ggplot(aes(x = max_dispersal_Euc_distance_km)) +
  geom_histogram() +
  geom_vline(aes(xintercept = max_dispersal_Euc_distance_km), data = ind_data_sum, 
             col = 2, size = 1.5) +
  labs(x = "Maximum Euclidean dispersal distance (km)",
       y = "Number of studies") +
  ylim(0, 7) +
  xlim(0, 1100) +
  ggeffects::theme_ggeffects()
g_ed_max_m

g_ed_hist_m <- ggpubr::ggarrange(g_ed_mean_m, g_ed_max_m + ylab(""), ncol = 2, labels = c("A", "B"))
g_ed_hist_m

ggsave("euc_dist_hist_malesonly.png", plot = g_ed_hist, path = "output/05_literature/", 
       width = 20, height = 10, units = "cm", dpi = 300)

# boxplots
g_ed_sex_mean <- dat_all |> 
  dplyr::mutate(sex = ifelse(sex == "F", "Female", "Male")) |> 
  dplyr::filter(!grepl("Hawley|Our study", study)) |> 
  dplyr::slice(grep("mean|raw", dat_all$data)) |> 
  ggplot(aes(y = dispersal_Euclidean_distance_km)) +
  geom_boxplot(aes(x = sex)) +
  labs(x = "Sex", y = "Mean Euclidean dispersal distance (km)") +
  ggeffects::theme_ggeffects() +
  ylim(0, 1100)
g_ed_sex_mean

g_ed_sex_max <- dat_all |> 
  dplyr::mutate(sex = ifelse(sex == "F", "Female", "Male")) |>
  dplyr::filter(!grepl("Hawley|Our study", study), !is.na(sex)) |> 
  dplyr::slice(grep("max|range|raw", dat_all$data)) |> 
  ggplot(aes(y = max_dispersal_Euc_distance_km)) +
  geom_boxplot(aes(x = sex)) +
  labs(x = "Sex", y = "Maximum Euclidean dispersal distance (km)") +
  ggeffects::theme_ggeffects() +
  ylim(0, 1100)
g_ed_sex_max

g_ed_sex <- ggpubr::ggarrange(g_ed_sex_mean, g_ed_sex_max, ncol = 2, labels = c("A", "B"))
g_ed_sex

ggsave("euc_dist_sex_boxplot.png", plot = g_ed_sex, path = "output/05_literature/", 
       width = 20, height = 10, units = "cm", dpi = 300)


#-----
# Total distance - hist

# mean
dat_all$data
tot_dat <- dat_all |> 
  # dplyr::slice(grep("mean|raw|NA", dat_all$data)) |> 
  dplyr::filter(!is.na(dispersal_total_distance_km), !grepl("Hawley", study)) 

View(tot_dat)

# print
# tot_dat |> 
#   dplyr::select(study, country, year, sex, n, contains("Euclidean"), contains("total")) |> 
#   dplyr::arrange(sex, year) |> 
#   readr::write_csv("output/05_literature/table_total_disp_dist.csv")

# only four studies reported it
summary(tot_dat)
tot <- ind_data_sum$dispersal_total_distance_km
sum(tot_dat$dispersal_total_distance_km < tot, na.rm = TRUE)/(nrow(tot_dat))

g_tot_mean <-  tot_dat |> 
  dplyr::filter(!grepl("Our study", study)) |> 
  ggplot(aes(x = dispersal_total_distance_km)) +
  geom_histogram() +
  geom_vline(aes(xintercept = dispersal_total_distance_km), data = ind_data_sum, 
             col = 2, size = 1.5) +
  labs(x = "Mean total dispersal distance (km)",
       y = "Number of studies") +
  ggeffects::theme_ggeffects()
g_tot_mean

# max
tot_dat <- dat_all |> 
  # dplyr::slice(grep("mean|raw|NA", dat_all$data)) |> 
  dplyr::filter(!is.na(dispersal_total_distance_km), !grepl("Hawley", study)) 

summary(tot_dat_max)
tot_max <- ind_data_sum$max_dispersal_total_distance_km
sum(tot_dat_max$max_dispersal_total_distance_km < tot_max, na.rm = TRUE)/(nrow(tot_dat_max) - 1)

g_tot_max <- tot_dat |> 
  dplyr::filter(!grepl("Our study", study)) |> 
  ggplot(aes(x = max_dispersal_total_distance_km)) +
  geom_histogram() +
  geom_vline(aes(xintercept = max_dispersal_total_distance_km), data = ind_data_sum, 
             col = 2, size = 1.5) +
  labs(x = "Maximum total dispersal distance (km)",
       y = "Number of studies") +
  ggeffects::theme_ggeffects()
g_tot_max

# Very few studies reported dispersal distance
# Several used GPS data and reported dispersal events, but did not have a focus on estimating
# dispersal distances

#-----
# dispersal age - hist

# mean
dat_all$mean_dispersal_age_months
age_dat <- dat_all |> 
  dplyr::filter(source != "Our study") |> 
  dplyr::filter(!is.na(mean_dispersal_age_months))

summary(age_dat)
age_our <- ind_data_sum$mean_dispersal_age_months
sum(age_dat$mean_dispersal_age_months < age_our)/(nrow(age_dat))
range(age_dat$mean_dispersal_age_months)
mean(age_dat$mean_dispersal_age_months)

g_age_mean <-  dat_all |> 
  dplyr::filter(source != "Our study") |> 
  dplyr::filter(!is.na(mean_dispersal_age_months)) |>  
  ggplot(aes(x = mean_dispersal_age_months)) +
  geom_histogram() +
  geom_vline(aes(xintercept = mean_dispersal_age_months), data = ind_data_sum, 
             col = 2, size = 1.5) +
  labs(x = "Mean dispersal age (months)",
       y = "Number of studies") +
  ggeffects::theme_ggeffects()
g_age_mean

# min
g_age_min <-  dat_all |> 
  dplyr::filter(source != "Our study") |> 
  dplyr::filter(!is.na(min_dispersal_age_months)) |>  
  ggplot(aes(x = min_dispersal_age_months)) +
  geom_histogram() +
  geom_vline(aes(xintercept = min_dispersal_age_months), data = ind_data_sum, 
             col = 2, size = 1.5) +
  labs(x = "Minimum dispersal age (months)",
       y = "Number of studies") +
  ggeffects::theme_ggeffects()
g_age_min

# max
dat_all |> 
  dplyr::filter(source != "Our study") %>%
  hist(.$max_dispersal_age_months, na.rm = TRUE)

g_age_max <-  dat_all |> 
  dplyr::filter(source != "Our study") |> 
  dplyr::filter(!is.na(max_dispersal_age_months)) |>  
  ggplot(aes(x = max_dispersal_age_months)) +
  geom_histogram() +
  geom_vline(aes(xintercept = max_dispersal_age_months), data = ind_data_sum, 
             col = 2, size = 1.5) +
  labs(x = "Maximum dispersal age (months)",
       y = "Number of studies") +
  ggeffects::theme_ggeffects()
g_age_max

ggsave("disp_age_avg.png", plot = g_age_mean, path = "output/05_literature/", 
       width = 10, height = 10, units = "cm", dpi = 300)

#-------------------------------
# Organize dispersal table for paper

# install.packages("forestploter")
# library(forestploter)
# 
# library(grid)
# library(forestploter)
# 
# # Read provided sample example data
# dt <- read.csv(system.file("extdata", "example_data.csv", package = "forestploter"))
# 
# # Keep needed columns
# dt <- dt[,1:6]
# 
# # indent the subgroup if there is a number in the placebo column
# dt$Subgroup <- ifelse(is.na(dt$Placebo), 
#                       dt$Subgroup,
#                       paste0("   ", dt$Subgroup))
# 
# # NA to blank or NA will be transformed to carachter.
# dt$Treatment <- ifelse(is.na(dt$Treatment), "", dt$Treatment)
# dt$Placebo <- ifelse(is.na(dt$Placebo), "", dt$Placebo)
# dt$se <- (log(dt$hi) - log(dt$est))/1.96
# 
# # Add blank column for the forest plot to display CI.
# # Adjust the column width with space. 
# dt$` ` <- paste(rep(" ", 20), collapse = " ")
# 
# # Create confidence interval column to display
# dt$`HR (95% CI)` <- ifelse(is.na(dt$se), "",
#                            sprintf("%.2f (%.2f to %.2f)",
#                                    dt$est, dt$low, dt$hi))
# head(dt)
# 
# # plot
# p <- forest(dt[,c(1:3, 8:9)],
#             est = dt$est,
#             lower = dt$low, 
#             upper = dt$hi,
#             sizes = dt$se,
#             ci_column = 4,
#             ref_line = 1,
#             arrow_lab = c("Placebo Better", "Treatment Better"),
#             xlim = c(0, 4),
#             ticks_at = c(0.5, 1, 2, 3),
#             footnote = "This is the demo data. Please feel free to change\nanything you want.")
# 
# # Print plot
# plot(p)

#------------------
# Dispersal Euclidaean distance
colnames(dat_all)

dat_fp_ed <- dat_all |> 
  dplyr::arrange(sex) |> 
  dplyr::select(Study = study, sex, n, data, est = dispersal_Euclidean_distance_km, 
                lower = min_dispersal_Euc_distance_km, 
                upper = max_dispersal_Euc_distance_km,
                country, state, telemetry_method)
# dat_fp_ed$`Average distance [min - max]` <- ifelse(
#   is.na(dat_fp_ed$lower) & is.na(dat_fp_ed$upper), sprintf("%.0f", dat_fp_ed$est),
#   sprintf("%.0f [%.0f to %.0f]",
#           dat_fp_ed$est, dat_fp_ed$lower, dat_fp_ed$upper))
# dat_fp_ed$lower <- ifelse(is.na(dat_fp_ed$lower), dat_fp_ed$est, dat_fp_ed$lower)
# dat_fp_ed$upper <- ifelse(is.na(dat_fp_ed$upper), dat_fp_ed$est, dat_fp_ed$upper)
# 
# p_ed <- forest(dat_fp_ed[,c(2:3)],
#             est = dat_fp_ed$est,
#             lower = dat_fp_ed$lower, 
#             upper = dat_fp_ed$upper,
#             sizes = 1,
#             ci_column = 2, ,
#             ref_line = 500,
#             # arrow_lab = c("Placebo Better", "Treatment Better"),
#             xlim = c(0, 1000),
#             # ticks_at = c(0.5, 1, 2, 3),
#             # footnote = "This is the demo data. Please feel free to change\nanything you want.")
# )
# 
# # Print plot
# plot(p_ed)

lit_disp <- readr::read_csv("text/literature_puma_data_dispersal.csv")
lit_disp |> print(width = Inf)

colnames(ind_data)
colnames(lit_disp)
# View(lit_disp)

# single
dat_all <- lit_disp |> 
  dplyr::filter(keep_analysis == 1) |> 
  dplyr::mutate(source = "literature",
                mean_dispersal_age_months = as.numeric(mean_dispersal_age_months),
                dispersal_duration_days = as.numeric(dispersal_duration_days)) |>
  dplyr::bind_rows(ind_data_sum) #|> 

dat_fp_ed <- dat_all |> 
  dplyr::arrange(sex, year) |> 
  dplyr::select(study, sex, year, n, data, est = dispersal_Euclidean_distance_km, 
                lower = min_dispersal_Euc_distance_km, 
                upper = max_dispersal_Euc_distance_km,
                country, state, telemetry_method) |> 
  dplyr::mutate(lower = ifelse(is.na(lower), est, lower),
                upper = ifelse(is.na(upper), est, upper))

# View(dat_fp_ed)

dat_fp_ed <- dat_fp_ed |> 
  dplyr::filter(data != "raw") |> 
  dplyr::bind_rows(
    dat_fp_ed |> 
      dplyr::filter(data == "raw") |> 
      dplyr::group_by(study, sex) |> 
      dplyr::summarize(n = n(), 
                       year = unique(year),
                       lower = min(est),
                       upper = max(est),
                       est = mean(est))
  ) |> 
  dplyr::left_join(
    dat_fp_ed |> 
      dplyr::group_by(study) |> 
      dplyr::summarise(nF = ifelse(sex == "F", n, 0),
                       nM = ifelse(sex == "M", n, 0)),
    by = "study"
  )
 

# Hawley - outlier 2400 km
# Hemker - no sex
# Weaver - only maximum

g_fp_ed <- dat_fp_ed |> 
  dplyr::mutate(#study = paste0(study, " (", nF, "F, ", nM, "M)"),
                study = forcats::fct_reorder(as.factor(study), est, .desc = TRUE)) |> 
  dplyr::filter(!grepl("Hawley|Hemker|Weaver", study)) |> 
  dplyr::mutate(sex = ifelse(sex == "F", "Female", "Male")) |> 
  ggplot(aes(y = study)) +
  geom_point(aes(x = est), shape = 15, size = 3) +
  geom_linerange(aes(xmin = lower, xmax = upper)) +
  facet_wrap(~sex) +
  scale_x_log10(limits = c(1,4000)) +
  theme_minimal() +
  geom_text(aes(x = 3000, label = sprintf("n = %d", n)), 
            nudge_x = 0, nudge_y = -0, size = 3) +
  labs(x = "Euclidean dispersal distance (km)",
       y = "")
g_fp_ed

ggsave("euc_dist_forest_plot.png", plot = g_fp_ed, path = "output/05_literature/", 
       width = 20, height = 20, units = "cm", dpi = 300)
