#' ---
#' title: 'Testing different methods to find dispersal'
#' author: Bernardo Niebuhr
#' ---
#' 

# Load packages
if(!require(install.load)) install.packages("install.load"); library(install.load)

install.load::install_load("ezknitr", "knitr")
install.load::install_load("tidyverse", "lubridate")
install.load::install_load("viridis")
install.load::install_load("sf")

# Print options for this document
options(width = 165)
opts_knit$set(root.dir = rprojroot::find_rstudio_root_file()) # root project folder
opts_chunk$set(error = F, message = F, warning = F, cache = F, echo = T, results = T)


#' 0) Read data

load("data/movement_data_pumas_tiete_legado.rda")

selected.individual <- "Piloto"
new.projection <- "+proj=aea +lat_1=-2 +lat_2=-22 +lat_0=-12 +lon_0=-54 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

# Test for a dispersal individual
disperser.d <- mov.data %>%
  dplyr::filter(name == selected.individual) %>%
  dplyr::arrange(timestamp) %>% 
  dplyr::mutate(Time = as.numeric(round(difftime(timestamp, timestamp[1], unit = "day")))) %>% 
  plyr::ddply("Time", summarize, X = mean(X), Y = mean(Y), Time = mean(Time), Date = mean(timestamp)) %>% 
  sf::st_as_sf(coords = c("X", "Y"), crs = 4326) %>% 
  sf::st_transform(crs = new.projection) %>% 
  dplyr::bind_cols(dplyr::as_tibble(sf::st_coordinates(.))) %>% 
  sf::st_drop_geometry()

#' # 1) FTP/residence time

library(adehabitatLT)
library(move)
library(recurse)

# ltraj object
disperser.ltraj <- with(disperser.d, as.ltraj(xy = cbind(X, Y), date = Date,
                                                     id = selected.individual,
                                                     proj4string = sp::CRS(new.projection)))
# disperser.ltraj         
plot(disperser.ltraj, xlab = "x", ylab = "y")

# move object
# transform the data into move object
move.data <- with(disperser.d, 
                  move(x = X, y = Y,
                  time = Date,
                  proj = sp::CRS(new.projection),
                  data = disperser.d, 
                  animal = selected.individual, 
                  sensor = 'GPS'))

# move.data
move::plot(move.data, pch = 19, xlab = 'Longitude', ylab = 'Latitude')

# Revisitation analysis - recurse

# visualization

# Calculate revisitation
mov.ind <- move.data # one individual
radius <- 20000 # in meters

revisitation <- getRecursions(mov.ind, radius = radius)

# plot

# time vs positions
lab_dates <- mov.ind$Date[seq(1, length(mov.ind), length.out = 10)]

g1 <- as(mov.ind, 'data.frame') %>% 
  ggplot(aes(x = X, y = Y, color = as.numeric(Date))) + 
  geom_point() +
  scale_color_viridis(breaks = as.numeric(lab_dates), 
                      labels = as_date(lab_dates)) +
  labs(x = '', y = '', color = "Time") +
  theme(legend.key.height = unit(2.5, "cm"))
g1

# plot number of visits
plot(revisitation, mov.ind, pch = 19, legendPos = c(685000, 7240000))
drawCircle(450000, -1200000, radius)
segments(450000, -1200000, 450000, -1200000+radius)

g2 <- as(mov.ind, 'data.frame') %>% 
  dplyr::mutate(residencetime = revisitation$residenceTime) %>% 
  ggplot(aes(x = X, y = Y, color = residencetime)) + 
  geom_point() +
  scale_colour_viridis() +
  labs(x = '', y = '', color = "Residence time (h)")
g2

# Revisitation analysis - adehabitatLT

mov.ind <- disperser.ltraj[c(1)] # one individual

# Plot distance vs time
plot(mov.ind)
plotltr(disperser.ltraj[c(1)])

# residence time analysis
radius <- 5000
plot(residenceTime(mov.ind, radius = radius, maxt = 1, units = "days")) # to visualize

radius <- 20000
plot(residenceTime(mov.ind, radius = radius, maxt = 1, units = "days")) # to visualize
res1 <- residenceTime(mov.ind, radius = radius, maxt = 1, addinfo = T, units = "days") # to analyze

# searching for clusters
# this functions performs segmentation of a time series
segmentation <- adehabitatLT::lavielle(res1,  which = paste0("RT.", radius), Kmax = 20, Lmin = 2)

# estimating the number of segments - k = number of groups, Jk is the number of positions per group
# chooseseg(segmentation) 

# now we find the limits of the segments, given a number of defined segments
pa <- findpath(segmentation, K = 5)
pa # see this object, different bursts

# select segment (given its number) and transform to dataframe
df <- ld(pa[c(2)]) # we are taking the second cluster, where the peak is

# take id, begin date and end date
id <- as.character(df$id[1])
date1  <- df$date[1]
date2  <- df$date[nrow(df)]
# calculate positions
x1 <- df$x[1]; y1 <- df$y[1]
x2 <- df$x[nrow(df)]; y2 <- df$y[nrow(df)]

(day1 <- difftime(date1, disperser.d$Date[1], "days"))
(day2 <- difftime(date2, disperser.d$Date[1], "days"))

# data.frame
df1 <- tibble(id, x1, y1, x2, y2, date1, date2, day1, day2, duration = date2 - date1)

# just to check
plot(mov.ind, xlab = "x", ylab = "y")
points(df1$x1, df1$y1, col = "purple", pch = 20, cex = 2)
points(df1$x2, df1$y2, col = "purple", pch = 20, cex = 2)

# initialize data frame for all individuals
(df.all <- df1)

#' # 2) smoove

library(smoove)

# visualize
with(disperser.d, scan_track(z = X + 1i*Y, time = Time))

# set window sweep
simSweep <- with(disperser.d, sweepRACVM(Z = X + 1i*Y, T = Time, 
                                      windowsize = 50, windowstep = 5, 
                                      model = "UCVM", progress=FALSE))

par(mfrow = c(1,1))
plotWindowSweep(simSweep)

# obtain candidate points
CP.all <- findCandidateChangePoints(windowsweep = simSweep, clusterwidth = 0)
CP.all %>% as.vector

CP.clustered <- findCandidateChangePoints(windowsweep = simSweep, clusterwidth = 20)
CP.clustered %>% as.vector

# select significant change points
getCPtable(CPs = CP.clustered, modelset = c("UCVM"), tidy = NULL)
getCPtable(CPs = CP.clustered, modelset = c("UCVM", "ACVM"), iterate = TRUE)
getCPtable(CPs = CP.clustered, modelset = c("UCVM", "ACVM"), criterion = "AIC")

# estimating the full model
simCP.table <- simSweep %>% findCandidateChangePoints(clusterwidth = 20) %>% 
  getCPtable(modelset = c("UCVM", "ACVM"))

# the results change in this case
# simCP.table <- simSweep %>% findCandidateChangePoints(clusterwidth = 20) %>% 
#   getCPtable(modelset = c("UCVM", "ACVM"), criterion = "AIC")

simPhaselist <- estimatePhases(simCP.table)

summarizePhases(simPhaselist)

# png("output/test_smoove.png", width = 15, height = 20, 
#     units = "cm", res = 300)

layout(c(1,1,1,2:6))
# extract locations and times
Z <- disperser.d$X + 1i*disperser.d$Y
T <- disperser.d$Time

# Note, these are also contained in the attributes of *any* of the 
# intermediate objects, e.g.
Z <- attributes(simPhaselist)$Z
T <- attributes(simPhaselist)$time

require(gplots) # for rich colors
cols <- rich.colors(length(simPhaselist))
T.cuts <- c(T[1], simCP.table$CP, T[length(T)])
Z.cols <- cols[cut(T, T.cuts, include.lowest = TRUE)]

phaseTable <- summarizePhases(simPhaselist)

plot(Z, asp=1, type="l", xpd=FALSE)
points(Z, col=Z.cols, pch=21, bg = alpha(Z.cols, 0.5), cex=0.8)
legend("topleft", legend = paste0(phaseTable$phase, ": ", phaseTable$model), 
       fill=cols, ncol=2, bty="n", title = "Phase: model")

par(mar=c(0,0,1,0), 
    xpd=NA)
plotPhaseParameter("tau", simPhaselist, ylab="", xaxt="n", xlab="", col=cols, log="y")
plotPhaseParameter("eta", simPhaselist,  ylab="", xaxt="n", xlab="", col=cols)
plotPhaseParameter("mu.x", simPhaselist,  ylab= "", xaxt="n", xlab="", col=cols)
plotPhaseParameter("mu.y", simPhaselist,  ylab= "", xaxt="n", xlab="", col=cols)
plotPhaseParameter("rms", simPhaselist,  ylab= "", xlab="time", col=cols)
# dev.off()

# phases 3, 4, 5 seem to be dispersal (with higher tau)

# take values
id <- selected.individual
dispersal.phase <- summarizePhases(simPhaselist)[3:5,]
day1 <- dispersal.phase$start[1]
day2 <- dispersal.phase$start[3]
# time
date1 <- disperser.d$Date[disperser.d$Time == floor(day1)]
date2 <- disperser.d$Date[disperser.d$Time == floor(day2)]
# coordinates
x1 <- disperser.d$X[disperser.d$Time == floor(day1)]
x2 <- disperser.d$X[disperser.d$Time == floor(day2)]
y1 <- disperser.d$Y[disperser.d$Time == floor(day1)]
y2 <- disperser.d$Y[disperser.d$Time == floor(day2)]

df1 <- tibble(id, x1, y1, x2, y2, date1, date2, day1, day2, duration = date2 - date1)

# check
plot(disperser.ltraj, xlab = "x", ylab = "y")
points(df1$x1, df1$y1, col = "purple", pch = 20, cex = 2)
points(df1$x2, df1$y2, col = "purple", pch = 20, cex = 2)

# initialize data frame for all individuals
df.all <- rbind(df.all, df1)

#' # 3) Net-squared displacement
disperser.d.nsd <- disperser.d %>% 
  # removing starting dispersal
  dplyr::slice(30:350) %>% 
  dplyr::mutate(nsd = sqrt((X - X[1])**2 + (Y - Y[1])**2)/1000)

with(disperser.d.nsd, plot(Time, nsd))

# dispersal functions
disp.fn <- function(t, dist, tim, dur) {
  dist / (1 + exp((tim - t)/dur))
}

# plot curve + parameters and meanings
dist = 3000
tim = 150
dur = 20
curve(disp.fn(x, dist = dist, tim = tim, dur = dur), from = 0, to = 365, xlab = 'Time', ylab = 'Distance',
      main = 'dist = 3000; tim = 150, dur = 20')
abline(h = dist, lty = 2); text(x = 20, y = dist-50, labels = expression(italic('dist'))); text(x = 20, y = dist + 50, labels = 'distance')
abline(h = dist/2, lty = 3)
abline(v = tim, lty = 3); text(x = tim + 10, y = 0, labels = expression(italic('tim')))
arrows(x0 = tim, y0 = 3*dist/4, x1 = tim + dur, y1 = 3*dist/4, code = 3, length = 0.1, angle = 90, lwd = 2); 
text(x = tim+20, y = 3*dist/4+100, labels = expression(italic('dur')))
# interpretation
abline(v = tim - 2*dur, lty = 4); text(x = tim - 2*dur +20, y = 100, labels = 'timing')
arrows(x0 = tim - 2*dur, y0 = 7*dist/8, x1 = tim + 2*dur, y1 = 7*dist/8, code = 3, length = 0.2, angle = 90, lwd = 2); 
text(x = tim, y = 7*dist/8+100, labels = 'duration')


# fit using nls
# functions/formulas
disp <- nsd ~ dist / (1 + exp((tim - Time)/dur))
disp.start.exp <- nsd ~ dist * exp(t/dur)

(m1 <- stats::nls(disp,
                  data = disperser.d.nsd,
                  start = list(dist = 80, tim = 50, dur = 10)))

curve(disp.fn(x, dist = coef(m1)[1], tim = coef(m1)[2], dur = coef(m1)[3]), 
      from = 0, to = max(disperser.d$Time), 
      add = T, col = 2, lwd = 2,
      xlab = 'Time', ylab = 'Distance (km)')
with(disperser.d.nsd, lines(Time, nsd, type = 'l'))
abline(v = coef(m1)[2] - 2*coef(m1)[3], col = 1)
abline(v = coef(m1)[2] + 2*coef(m1)[3], col = 1)

# take values
id <- selected.individual
day1 <- coef(m1)[2] - 2*coef(m1)[3]
day2 <-  coef(m1)[2] + 2*coef(m1)[3]
# time
date1 <- disperser.d$Date[disperser.d$Time == floor(day1)]
date2 <- disperser.d$Date[disperser.d$Time == floor(day2)]
# coordinates
x1 <- disperser.d$X[disperser.d$Time == floor(day1)]
x2 <- disperser.d$X[disperser.d$Time == floor(day2)]
y1 <- disperser.d$Y[disperser.d$Time == floor(day1)]
y2 <- disperser.d$Y[disperser.d$Time == floor(day2)]

df1 <- tibble(id, x1, y1, x2, y2, date1, date2, day1, day2, duration = date2 - date1)

# check
plot(disperser.ltraj, xlab = "x", ylab = "y")
points(df1$x1, df1$y1, col = "purple", pch = 20, cex = 2)
points(df1$x2, df1$y2, col = "purple", pch = 20, cex = 2)

# initialize data frame for all individuals
df.all <- rbind(df.all, df1)

#' # 4) Hiden Markov models

library(moveHMM)

disperser.d.hmm <- disperser.d %>% 
  dplyr::mutate(x = X/1000,
                y = Y/1000,
                ID = selected.individual)

hmm.data <- disperser.d.hmm %>% 
  as.data.frame %>% 
  prepData(type='UTM', coordNames = c('x', 'y'))
head(hmm.data)

# Ploting data and basic movement data statistics
plot(hmm.data)

#--- 
#  Fitting the model

# Let's fit a 2- and a 3-state model with no covariates
# Let's first fit the model for step lengths only

## 2-state model: initial parameters for gamma distribution
mu0 <- c(0.5, 1) # step mean (two parameters: one for each state)
sigma0 <- c(0.5, 1) # step SD
stepPar0 <- c(mu0, sigma0)

## call to fitting function
m1 <- fitHMM(data = hmm.data, nbStates = 2, stepPar0 = stepPar0, 
             angleDist = 'none', formula = ~1)
m1

## 3-state model: initial parameters for gamma distribution
mu1 <- c(0.5, 1, 2) # step mean (two parameters: one for each state)
sigma1 <- c(0.5, 1, 2) # step SD
stepPar1 <- c(mu1, sigma1)

m2 <- fitHMM(data = hmm.data, nbStates = 3, stepPar0 = stepPar1, 
             angleDist = 'none', formula = ~1)
m2

# Compare models
AIC(m1, m2)

# Plot the model
plot(m1)
plot(m2)

# but this would also make sense!

## 4-state model: initial parameters for gamma and von Mises distributions
mu1 <- c(0.5, 1, 2, 4) # step mean (two parameters: one for each state)
sigma1 <- c(0.5, 1, 2, 3) # step SD
stepPar1 <- c(mu1, sigma1)
angleMean1 <- c(0, pi/4, pi/2, 0) # angle mean
kappa1 <- c(0.1, 1, 5, 2) # angle concentration
anglePar1 <- c(angleMean1, kappa1)

m5 <- fitHMM(data = hmm.data, nbStates = 4, stepPar0 = stepPar1, 
             anglePar0 = anglePar1, formula = ~1)
AIC(m1, m2, m5)
plot(m5) # most parsimonious

# take values
id <- selected.individual

plotStates(m2)
viterbi(m2)
days <- which(viterbi(m2) == 3)
day1 <- days[1]
day2 <-  days[length(days)-3]
# time
date1 <- disperser.d$Date[disperser.d$Time == floor(day1)]
date2 <- disperser.d$Date[disperser.d$Time == floor(day2)]
# coordinates
x1 <- disperser.d$X[disperser.d$Time == floor(day1)]
x2 <- disperser.d$X[disperser.d$Time == floor(day2)]
y1 <- disperser.d$Y[disperser.d$Time == floor(day1)]
y2 <- disperser.d$Y[disperser.d$Time == floor(day2)]

df1 <- tibble(id, x1, y1, x2, y2, date1, date2, day1, day2, duration = date2 - date1)

# check
plot(disperser.ltraj, xlab = "x", ylab = "y")
points(df1$x1, df1$y1, col = "purple", pch = 20, cex = 2)
points(df1$x2, df1$y2, col = "purple", pch = 20, cex = 2)

# initialize data frame for all individuals
df.all <- rbind(df.all, df1)

#' # 5) Our approach

source("code/source_arima_function_dispersal_residency.R")
library(amt)

with(disperser.d.nsd, ScanTrack(Time, X, Y, unit = "days"))

# test finding both and selecting methods, 1d fix rate
dates.1d <- with(disperser.d.nsd, findDepartSettlingDate(Time, X, Y, Date,
                                                  event = c("departure", "settling"),
                                                  condition.settling.to.departure = F))
dates.1d
# bad dates

# Test for a dispersal individual
disperser.6h <- mov.data %>%
  dplyr::filter(name == selected.individual) %>%
  amt::mk_track(X, Y, timestamp) %>%
  amt::track_resample(rate = hours(6), tolerance = hours(1)) %>%
  dplyr::mutate(Time = row_number()) %>%
  dplyr::select(X = x_, Y = y_, Time, Date = t_) %>%
  sf::st_as_sf(coords = c("X", "Y"), crs = 4326) %>%
  sf::st_transform(crs = new.projection) %>%
  dplyr::bind_cols(dplyr::as_tibble(sf::st_coordinates(.))) %>%
  sf::st_drop_geometry()

# test finding both and selecting methods, 6h fix rate
dates.6h <- with(disperser.6h, findDepartSettlingDate(Time, X, Y, Date,
                                                      event = c("departure", "settling"),
                                                      condition.settling.to.departure = F))
dates.6h
# that seems good!

with(disperser.6h,
     with(dates.6h, ScanTrack(Time, X, Y, dep.time = departure$time,
                              setl.time = settling$time)))
# that's good!
with(disperser.6h, compare.dispersal.models(Time, X, Y,
                                            cp = c(dates.6h$departure$time,
                                                   dates.6h$settling$time)))


# take values
id <- selected.individual
# time
date1 <- dates.6h$departure$date
date2 <- dates.6h$settling$date
day1 <- difftime(date1, disperser.6h$Date[1], "days")
day2 <- difftime(date2, disperser.6h$Date[1], "days")

# coordinates
x1 <- disperser.d$X[disperser.d$Time == floor(day1)]
x2 <- disperser.d$X[disperser.d$Time == floor(day2)]
y1 <- disperser.d$Y[disperser.d$Time == floor(day1)]
y2 <- disperser.d$Y[disperser.d$Time == floor(day2)]

df1 <- tibble(id, x1, y1, x2, y2, date1, date2, day1, day2, duration = date2 - date1)

# check
par(mfrow = c(1,1))
plot(disperser.ltraj, xlab = "x", ylab = "y")
points(df1$x1, df1$y1, col = "purple", pch = 20, cex = 2)
points(df1$x2, df1$y2, col = "purple", pch = 20, cex = 2)

# initialize data frame for all individuals
df.all <- rbind(df.all, df1)

#' # Compare all approaches
df.all %>%
  knitr::kable()

df.all %>%
  dplyr::select(-c(2:5)) %>% 
  knitr::kable()



library(marcher)
marcher.estimate <- with(disperser.d, estimate_shift(T = Time, X = X, Y = Y))

plot(marcher.estimate)
summary(marcher.estimate)
