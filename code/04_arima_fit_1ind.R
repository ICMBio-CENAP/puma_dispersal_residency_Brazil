#' ---
#' title: "Puma space use and dispersal in the tropics: fitting ARIMA functions to identify dispersal for single individuals"
#' author: Bernardo Niebuhr
#' date: ""
#' output: 
#'   html_document: default
#'   github_document: default
#' ---
#' 
#' Here we use some of the functionalities so far developed in the R package
#' `disperser` to estimate dates of departure and arrival for dispersal events
#' for pumas in Southeastern Brazil.

# --------------- label=setup, warning=FALSE, message=FALSE, echo=TRUE

# Load packages
if(!require(install.load)) install.packages("install.load"); library(install.load)

install.load::install_load("ezknitr", "knitr")
install.load::install_load("tidyverse", "lubridate")
install.load::install_load("amt", "forecast")

library(devtools)
devtools::install_github("EliGurarie/disperser")


# Print options for this document
options(width = 165)
opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) # root project folder
opts_chunk$set(error = F, message = F, warning = F, cache = F, echo = T, results = T)

# Set up 

# Clean everything before beginning
rm(list = ls())

#' We start by loading the data and selecting one individual, Piloto, to test and 
#' illustrate the procedure.

# --------------- label=load_data_and_functions

# Load data
library(PardasIPC)
data("pumas")

# ----- label=for_one_individual

selected_individual <- "Piloto"

# all this other individuals could be tested, since they seem to have dispersed
pumas %>% 
  dplyr::filter(dispersal_behavior == "dispersal") %>% 
  pull(id) %>% unique() %>% sort()

#' # Using `find_date_departure` and `find_date_settling` functions
#' 
#' First we resample the data to one point per day, which is enough for estimating 
#' dispersal events, since they take place at the time scales of dozens of days.
#' 

disperser_d <- pumas %>%
  dplyr::filter(id == selected_individual) %>%
  dplyr::arrange(timestamp) %>% 
  dplyr::mutate(Time = as.numeric(round(difftime(timestamp, timestamp[1], unit = "day")))) %>% 
  plyr::ddply("Time", summarize, X = mean(X), Y = mean(Y), Time = mean(Time), Date = mean(timestamp)) %>% 
  dplyr::as_tibble()

disperser_d

#' We then plot the trajectory to have an idea about how the movement for 
#' this individual happened.

# Plot
with(disperser_d, scan_track(Time, X, Y))

#' ## Departure
#' 
#' Now we start by estimating the departure time.

# departure day
(dep_time <- with(disperser_d, find_date_departure(Time, X, Y, method = 'ML')))
# Get the correspondent date-time
(dep_dttm <- disperser_d$Date[disperser_d$Time == floor(dep_time)])
# plot
with(disperser_d, scan_track(Time, X, Y, departure.time = dep_time))

#' ## Settling
#' 
#' Now we estimate the settling date, considering both all the data (from the first day) 
#' or only allowing the settling date to happen after the 
#' departure date estimated above.

# settling day, considering from the first day
(setl_time <- with(disperser_d, find_date_settling(Time, X, Y, method = 'ML')))
(setl_dttm <- disperser_d$Date[disperser_d$Time == floor(setl_time)])
with(disperser_d, scan_track(Time, X, Y, settling.time = setl_time))

# settling day, considering from the departure date estimated above
(setl_time <- with(disperser_d, find_date_settling(Time, X, Y, method = 'ML',
                                                   lower = dep_time)))
(setl_dttm <- disperser_d$Date[disperser_d$Time == floor(setl_time)])
with(disperser_d, scan_track(Time, X, Y, settling.time = setl_time))

#' WE can now plot both departure and settlement. 
with(disperser_d, scan_track(Time, X, Y, departure.time = dep_time, 
                             settling.time = setl_time))


#' # Re-doing the same thing using a short method: 
#' the `find_date_departure_settling` function
#' 

#' Now we test another function which already estimates both departure and settling.

# test finding both and selecting methods, 1d fix rate
dates_1d <- with(disperser_d, 
                 find_date_departure_settlling(Time, X, Y,
                                               event = c("departure", "settling"),
                                               condition.settling.to.departure = T))

dates_1d

with(disperser_d, 
     with(dates_1d, scan_track(Time, X, Y, departure.time = departure$time, 
                               settling.time = settling$time)))



#' # Now we do it for another individual
#' 
#' Now we do it for the individual, the puma  Mineiro. Again, for this individual the expected
#' days are 14-29 (departure) and 62 (settling).

selected_individual <- "Mineiro"

# Test for a dispersal individual, 1d fix rate
disperser_1d <- pumas %>%
  dplyr::filter(id == selected_individual) %>%
  dplyr::arrange(timestamp) %>% 
  dplyr::mutate(Time = as.numeric(round(difftime(timestamp, timestamp[1], unit = "day")))) %>% 
  plyr::ddply("Time", summarize, X = mean(X), Y = mean(Y), Time = mean(Time), Date = mean(timestamp)) %>% 
  dplyr::as_tibble()


# plot
with(disperser_1d, scan_track(Time, X, Y))

#' Testing both departure and settling:

# test finding both and selecting methods, 1d fix rate
dates_1d <- with(disperser_1d, 
                 find_date_departure_settlling(Time, X, Y,
                                               event = c("departure", "settling"),
                                               condition.settling.to.departure = F))

dates_1d

with(disperser_1d, 
     with(dates_1d, scan_track(Time, X, Y, departure.time = dates_1d$departure$time, 
                               settling.time = dates_1d$settling$time)))


#' This is strange - these are not good estimates... it gets even worse if we use
#' the argument `condition.settling.to.departure = TRUE`.
#' 
#' Let's try with a different input, with 1 position per 6h.

disperser_6h <- pumas %>%
  dplyr::filter(id == selected_individual) %>%
  amt::mk_track(X, Y, timestamp) %>% 
  amt::track_resample(rate = hours(6), tolerance = hours(1)) %>%
  dplyr::mutate(Time = row_number()) %>%
  dplyr::select(X = x_, Y = y_, Time, Date = t_)


dates_6h <- with(disperser_6h, find_date_departure_settlling(Time, X, Y, 
                                                      event = c("departure", "settling"),
                                                      condition.settling.to.departure = T))
dates_6h
# that seems good!
# remember that these estimates are in units of 6h now, not in days anymore.

with(disperser_6h, 
     with(dates_6h, scan_track(Time, X, Y, departure.time = departure$time, 
                               settling.time = settling$time)))

#' # Using the function `compare_dispersal_models`

#' Finally, the function below tries different models (residency, dispersal, 
#' departing, settling, or both departing and settling) and calculate their likelihood.
#' Let's try it. We use as input the candidate departure and arrival times 
#' estimated above.

with(disperser_6h, compare_dispersal_models(Time, X, Y, 
                                            cp = c(dates_6h$departure$time, 
                                                   dates_6h$settling$time)))

#' We see that, using these estimates for the departure and settling dates, there is 
#' higher support for the `depart-settle` model, which is reasonable, given the data.
