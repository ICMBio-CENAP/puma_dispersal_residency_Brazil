#' ---
#' title: 'ARIMA fit for dispersing pumas - for several individuals'
#' author: Bernardo Niebuhr,
#' output: 
#'   html_document: default
#'   github_document: default
#' ---

# Load packages
if(!require(install.load)) install.packages("install.load"); library(install.load)

install.load::install_load("ezknitr", "knitr")
install.load::install_load("tidyverse", "lubridate", "purrr")
install.load::install_load("amt", "forecast", "disperser")

# Print options for this document
options(width = 165)
opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) # root project folder
opts_chunk$set(eval = F, error = F, message = F, warning = F, cache = F, echo = T, results = T)

# --------------- label=setup
# Set up 

# Clean everything before beginning
rm(list = ls())

# --------------- label=load_data_and_functions

# Load data
library(PardasIPC)
data("pumas")
pumas

# remove data for Kurupi in the process of dying
date_Kurupi_dying <- lubridate::ydm_h("20191211 01")
pumas <- pumas %>% 
  dplyr::filter(!(id == "Kurupi" & timestamp > date_Kurupi_dying))

# pumas %>% 
#   dplyr::group_by(id) %>% 
#   dplyr::summarize(last_date = max(timestamp))

# nest and resample data for 1d, 1h, and 6h fix rate
pumas_resamp <- pumas %>% 
  dplyr::select(-species, -country, -hDOP) %>% 
  tidyr::nest(data = X:timestamp) %>% 
  dplyr::mutate(
    data_1d = purrr::map(data, function(x) x %>% 
                           dplyr::arrange(timestamp) %>%
                           dplyr::mutate(Time = as.numeric(round(difftime(timestamp, timestamp[1], unit = "day")))) %>%
                           plyr::ddply("Time", summarize, X = mean(X), Y = mean(Y), Time = mean(Time), Date = mean(timestamp)) %>%
                           dplyr::as_tibble()),
    data_1h = purrr::map(data, function(x) x %>% 
                           dplyr::arrange(timestamp) %>% 
                           dplyr::mutate(Time = as.numeric(round(difftime(timestamp, timestamp[1], unit = "hour")))) %>% 
                           dplyr::select(X, Y, Time, Date = timestamp) %>% 
                           dplyr::as_tibble()),
    data_6h = purrr::map(data, function(x) x %>% 
                           amt::mk_track(X, Y, timestamp) %>% 
                           amt::track_resample(rate = hours(6), tolerance = hours(1)) %>%
                           dplyr::mutate(Time = row_number()) %>%
                           dplyr::select(X = x_, Y = y_, Time, Date = t_)))

# change to longer format to ease dealing with multiple fix rates per individual
pumas_resamp_long <- pumas_resamp %>% 
  dplyr::select(-data) %>% 
  tidyr::pivot_longer(contains("data"), names_to = "fix_rate", values_to = "dat")


# calculate dates and models for dispersal
pumas_resamp_long_depset <- pumas_resamp_long %>% 
  dplyr::mutate(
    # find departure and settling dates
    depart_settle = purrr::map(dat, function(x)
      with(x, find_date_departure_settlling(Time, X, Y, 
                                     event = c("departure", "settling"),
                                     condition.settling.to.departure = T))),
    # extract dates from object
    departure_date = purrr::map(depart_settle, function(x) x$departure$time),
    settling_date = purrr::map(depart_settle, function(x) x$settling$time),
    # compare models using the estimated departure and settling dates
    compared_models = purrr::pmap(
      list(id, dat, departure_date, settling_date),
      function(id, x, a, b) {
        print(id)
        with(x, compare_dispersal_models(Time, X, Y, cp = c(a, b)))}),
    # get the name of the maximum likelihood model
    max_lik_model = purrr::map(compared_models, function(x) names(x)[which.max(x)]),
    delta_ll_from_resident = purrr::map(compared_models, function(x) round(max(as.numeric(unlist(x)), na.rm = T) - x[[1]], digits = 2))) %>%
  tidyr::unnest(c(contains(".date"), max_lik_model, delta_ll_from_resident))


# pdfs with plots
str(pumas_resamp_long_depset, 1)
pumas_names <- unique(pumas_resamp_long_depset$id)

# for each individual
for(i in pumas_names) {
  
  # open output
  print(i)
  pdf(paste0("output/01_arima_fits/scan_track_", i, ".pdf"))
  
  ind <- pumas_resamp_long_depset %>% 
    dplyr::filter(id == i)
  
  # Plot
  with(ind$dat[[2]], scan_track(Time, X, Y))
  
  # Plot fits
  purrr::pmap(
    list(ind$dat, ind$departure_date, ind$settling_date, ind$max_lik_model),
    function(x, a, b, title)
      if(title == "settling") with(x, scan_track(Time, X, Y, settling.time = b, main = title)) else
        if(title == "departure") with(x, scan_track(Time, X, Y, departure.time = a, main = title)) else
          if(title == "depart-settle") with(x, scan_track(Time, X, Y, departure.time = a, settling.time = b, main = title)) else
            with(x, scan_track(Time, X, Y, main = title))
  )
  
  # close output
  dev.off()
}

# which model to consider for each individual, from checking each one visually and through the models
pumas_names
selected_model <- c("resident", "resident", "depart-settle", "settling", 
                    "depart-settle", "resident", "depart-settle", "depart-settle", 
                    "depart-settle", "resident", "depart-settle",
                    "resident", "depart-settle", "resident")
selected_scale <- c(NA, NA, "data_6h", "data_6h",
                    "data_6h", NA, "data_1d", "data_6h",
                    "data_1d", NA, "data_6h", 
                    NA, "data_6h", NA)

# just a summary
pumas_resamp_long_depset %>% 
  dplyr::select(id, dispersal_behavior, fix_rate, departure_date, settling_date, max_lik_model) %>% 
  print(n = 50)

# estimating dates (POSIXct) of departure and settling, when applicable
dates_depset <- pumas_resamp_long_depset %>% 
  dplyr::select(id, dispersal_behavior, fix_rate, dat, departure_date, settling_date, max_lik_model) %>% 
  dplyr::filter(paste0(id, fix_rate) %in% paste0(pumas_names, ifelse(is.na(selected_scale), "data_1d", selected_scale))) %>% 
  dplyr::mutate(selected_model = selected_model,
                departure_date_posix = purrr::pmap(list(dat, departure_date, selected_model),
                                                   function(x, a, b) 
                                                     return(if_else(b == "departure" | b == "depart-settle",
                                                                    x$Date[x$Time == floor(a)], lubridate::ymd_h(NA)))),#) %>% #,
                settling_date_posix = purrr::pmap(list(dat, settling_date, selected_model),
                                                  function(x, a, b)
                                                    return(if_else(b == "settling" | b == "depart-settle",
                                                                   x$Date[x$Time == floor(a)], lubridate::ymd_h(NA))))) %>%
  tidyr::unnest(contains("posix"))

# dates_depset$settling_date_posix <- lubridate::ymd_h(NA)
# for(i in 1:nrow(dates_depset)) {
#   if(dates_depset$selected_model[i] == "settling" | dates_depset$selected_model[i] == "depart-settle") {
#     x <- dates_depset$dat[[i]]
#     dates_depset$settling_date_posix[i] <- x$Date[x$Time == floor(dates_depset$settling_date[i])]
#   }
# }

# now we return to the original data and classify the behavior according to this
pumas_behav <- pumas %>% 
  dplyr::left_join(
    dates_depset %>%
      dplyr::select(id, selected_model, contains("posix")),
    by = "id") %>% 
  dplyr::mutate(departure_date_posix = if_else(is.na(departure_date_posix), ymd_h("20000101 01"), departure_date_posix),
                settling_date_posix = if_else(is.na(settling_date_posix), ymd_h("20500101 01"), settling_date_posix),
                phase = if_else(selected_model == "resident", "post-dispersal", "dispersal"),
                phase = if_else(timestamp < departure_date_posix, "pre-dispersal", phase),
                phase = if_else(timestamp > settling_date_posix, "post-dispersal", phase),
  )

# check    
ggplot(pumas_behav) +
  geom_point(aes(X, Y, color = phase)) +
  facet_wrap(~id, scales = "free")

# great!

# save results
save(pumas_resamp_long_depset, file = "output/mov_data_departure_settle_arima_output.rda")
save(dates_depset, file = "output/dates_departure_settling_arima.rda")
save(pumas_behav, file = "data/movement_data_pumas_dispersal_behavior.rda")
