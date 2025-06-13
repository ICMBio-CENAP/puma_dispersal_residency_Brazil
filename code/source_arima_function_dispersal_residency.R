# problems: we set number of parameters k manualy
library(forecast)

#' Fit residency (ranging?) and dispersal models to data
#' 
#' Fits ranging and dispersal models to animal movement data using arima models.
#' 
#' @param cp vector of change points between behaviors (in day counts from the first day of monitoring)
#' @param Time day of monitoring (could be hour or just timestamp?)
#' @param X x or longitude coordinate
#' @param Y y or latitude coordinate
#' @param event type of behaviro or transitions between behaviors to be fit.
#' May be any of "settling" (from dispersal to residency), "departure" (from residency
#' to dispersal), "depart-settle" (residency-dispersal-residency), "resident", or "dispersal"
#' @param method.arima method to fit arima models; any of "CSS-ML", "ML", or "CSS. For more
#' details, see `?arima`.
#' @param ... other arguments for the `arima` function.
#' 
#' @return maximum likelihood of the model
ll.finddate <- function(Time, X, Y, cp, 
                        event = c("settling", "departure", "depart-settle", "resident", "dispersing")[1], 
                        method.arima = c("CSS-ML", "ML", "CSS")[1],
                        ...){
  
  if(event == "settling"){ d1 = 1; d2 = 0 } else
    if(event == "departure") { d1 = 0; d2 = 1 } else 
      if(event == "depart-settle") { d1 = 0; d2 = 1 } else
        if(event == "resident") { d1 = 0 } else
          if(event == "dispersing") { d1 = 1 } else
            stop("Event must be one of `resident`, `dispersing`, `settling`, `departure`, or `depart-settle`.")
  # if(event == "settling"){ order1 <- c(1,1,0); order2 <- c(1,0,0); k = 5 } else
  #   if(event == "departure") { order1 <- c(1,0,0); order2 <- c(1,1,0); k = 5} else 
  #     if(event == "depart-settle") { order1 <- c(1,0,0); order2 <- c(1,1,0); order3 <- c(1,0,0); k = 8} else
  #       if(event == "resident") { order1 <- c(1,0,0); k = 2 } else
  #         if(event == "dispersing") { order1 <- c(1,1,0); k = 2} else
  #           stop("Event must be one of `resident`, `dispersing`, `settling`, `departure`, or `depart-settle`.")
  
  
  
  if(event == "settling" | event == "departure") {
    fits <- list(
      X.fit1 = auto.arima(X[Time < cp[1]], d = d1, D = 0, max.p = 2, max.q = 2), 
      X.fit2 = auto.arima(X[Time >= cp[1]], d = d2, D = 0, max.p = 2, max.q = 2),
      Y.fit1 = auto.arima(Y[Time < cp[1]], d = d1, D = 0, max.p = 2, max.q = 2), 
      Y.fit2 = auto.arima(Y[Time >= cp[1]], d = d2, D = 0, max.p = 2, max.q = 2))
  } else {
    if(event == "depart-settle") {
      fits <- list(
        X.fit1 = auto.arima(X[Time < cp[1]], d = d1, D = 0, max.p = 2, max.q = 2),
        X.fit2 = auto.arima(X[Time >= cp[1] & Time < cp[2]], d = d2, D = 0, max.p = 2, max.q = 2),
        X.fit3 = auto.arima(X[Time >= cp[2]], d = d1, D = 0, max.p = 2, max.q = 2),
        Y.fit1 = auto.arima(Y[Time < cp[1]], d = d1, D = 0, max.p = 2, max.q = 2), 
        Y.fit1 = auto.arima(Y[Time >= cp[1] & Time < cp[2]], d = d2, D = 0, max.p = 2, max.q = 2),
        Y.fit3 = auto.arima(Y[Time >= cp[2]], d = d1, D = 0, max.p = 2, max.q = 2))
    } else {
      if(event == 'resident' | event == 'dispersing') {
        fits <- list(
          X.fit1 = auto.arima(X, d = d1, D = 0, max.p = 2, max.q = 2),
          Y.fit1 = auto.arima(Y, d = d1, D = 0, max.p = 2, max.q = 2))
      }
    }
  }
  
  # if(event == "settling" | event == "departure") {
  #   fits <- list(
  #     X.fit1 = arima(X[Time < cp[1]], order = order1, method = method.arima, ...),
  #     X.fit2 = arima(X[Time >= cp[1]], order = order2, method = method.arima,  ...),
  #     Y.fit1 = arima(Y[Time < cp[1]], order = order1, method = method.arima,  ...),
  #     Y.fit2 = arima(Y[Time >= cp[1]], order = order2, method = method.arima,  ...))
  # } else {
  #   if(event == "depart-settle") {
  #     fits <- list(
  #       X.fit1 = arima(X[Time < cp[1]], order = order1, method = method.arima,  ...),
  #       X.fit2 = arima(X[Time >= cp[1] & Time < cp[2]], order = order2, method = method.arima,  ...),
  #       X.fit3 = arima(X[Time >= cp[2]], order = order3, method = method.arima,  ...),
  #       Y.fit1 = arima(Y[Time < cp[1]], order = order1, method = method.arima,  ...),
  #       Y.fit2 = arima(Y[Time >= cp[1] & Time < cp[2]], order = order2, method = method.arima,  ...),
  #       Y.fit3 = arima(Y[Time >= cp[2]], order = order3, method = method.arima,  ...))
  #   } else {
  #     if(event == 'resident' | event == 'dispersing') {
  #       fits <- list(
  #         X.fit1 = arima(X, order = order1, method = method.arima,  ...),
  #         Y.fit1 = arima(Y, order = order1, method = method.arima,  ...))
  #     }
  #   }
  # }
  
  # How should we define/calculate the number of parameters?
  # Now we are defining k mannualy, but it is not the best solution...
  # For the time being, I comment this part and we return only the loglikelihood
  # ll <- ll = sum(sapply(fits, function(fit) fit$loglik))
  # AIC <- -2*ll + 2*k
  # 
  # return(list(ll = ll, k = k, AIC = AIC))
  
  return(ll = sum(sapply(fits, function(fit) fit$loglik)))
}

#' Fit residency (ranging?) and dispersal models to data - v2
#' 
#' Fits ranging and dispersal models to animal movement data using arima models.
#' 
#' This is exactly the same function as the one above, but parameterized differently.
#' Instead of providing cp as a vector of values - c(depart, settle) - we use two
#' values as different function arguments, cp and cp2. 
#' This is just to test because I found it easier to use this way with some of 
#' the optimization methods.
#' 
# ll.finddate2 <- function(Time, X, Y,
#                         cp = NULL, cp2 = NULL,
#                         event = c("settling","departure","depart-settle","resident","dispersing")[1],
#                         method.arima = c("CSS-ML", "ML", "CSS")[1],
#                         ...){
#   if(event == "settling"){ order1 <- c(1,1,0); order2 <- c(1,0,0); k = 5 } else
#     if(event == "departure") { order1 <- c(1,0,0); order2 <- c(1,1,0); k = 5} else
#       if(event == "depart-settle") { order1 <- c(1,0,0); order2 <- c(1,1,0); order3 <- c(1,0,0); k = 8} else
#         if(event == "resident") { order1 <- c(1,0,0); k = 2 } else
#           if(event == "dispersing") { order1 <- c(1,1,0); k = 2} else
#             stop("Event must be one of `resident`, `dispersing`, `settling`, `departure`, or `depart-settle`.")
# 
#   if(event == "settling" | event == "departure") {
#     fits <- list(
#       X.fit1 = arima(X[Time < cp], order = order1, method = method.arima,  ...),
#       X.fit2 = arima(X[Time >= cp], order = order2, method = method.arima,  ...),
#       Y.fit1 = arima(Y[Time < cp], order = order1, method = method.arima,  ...),
#       Y.fit2 = arima(Y[Time >= cp], order = order2, method = method.arima,  ...))
#   } else {
#     if(event == "depart-settle") {
#       fits <- list(
#         X.fit1 = arima(X[Time < cp], order = order1, method = method.arima,  ...),
#         X.fit2 = arima(X[Time >= cp & Time < cp2], order = order2, method = method.arima, ...),
#         X.fit3 = arima(X[Time >= cp2], order = order3, method = method.arima,  ...),
#         Y.fit1 = arima(Y[Time < cp], order = order1, method = method.arima,  ...),
#         Y.fit2 = arima(Y[Time >= cp & Time < cp2], order = order2, method = method.arima,  ...),
#         Y.fit3 = arima(Y[Time >= cp2], order = order3, method = method.arima,  ...))
#     } else {
#       if(event == 'resident' | event == 'dispersing') {
#         fits <- list(
#           X.fit1 = arima(X, order = order1, method = method.arima,  ...),
#           Y.fit1 = arima(Y, order = order1, method = method.arima,  ...))
#       }
#     }
#   }
# 
#   # ll <- ll = sum(sapply(fits, function(fit) fit$loglik))
#   # AIC <- -2*ll + 2*k
#   #
#   # return(list(ll = ll, k = k, AIC = AIC))
# 
#   return(ll = sum(sapply(fits, function(fit) fit$loglik)))
# }

# with(disperser.d, optimize(ll.finddate, Time=days, X = X, Y = Y, maximum = TRUE,
#       lower = min(days), upper = max(days), event = "settling", method.arima = "ML"))


#' Find best candidate dates for departure and settling
#' 
#' These functions estimate the maximum likelihood departure and settling dates
#' for a dispersal process, based on the ll.finddate function.
#' 
#' @param Time day of monitoring (could be hour of monitoring?)
#' @param X x or longitude coordinate
#' @param Y y or latitude coordinate
#' @param ... additional arguments for the `optimize` function.
#' 
#' @return the date that corresponds to the maximum likelihood model for departure 
#' or settling
findDepartureDate <- function(Time, X, Y, lower = min(Time), upper = max(Time), ...)
  optimize(ll.finddate, Time = Time, X = X, Y = Y, maximum = TRUE, 
           lower = lower, upper = upper, event = "departure", ...)$maximum

findSettlingDate <- function(Time, X, Y, lower = min(Time), upper = max(Time), ...) 
  optimize(ll.finddate, Time = Time, X = X, Y = Y, maximum = TRUE, 
           lower = lower, upper = upper, event = "settling", ...)$maximum

#' ## Obs
#' 
#' Later on, this function could return also the date-time value for this day
#' (e.g. day 132, 2020-04-02)


#' Finding both departure and settling times
#' 
#' This function does the same as the ones above, but it may find either the
#' departure and the settling dates (or both), and it can use either one (or all)
#' methods to fit the arima models (CSS-ML, ML, CSS). This is done so that the estimation
#' method that do not produce errors or warnings may be select, for a more reliable 
#' estimation.
#' 
#' @param Time day of monitoring (could be hour or other periods of time)
#' @param X x or longitude coordinate
#' @param Y y or latitude coordinate
#' @param Date date and time (POSIXct) of monitoring
#' @param event type of behavioral transitions to be fit.
#' May be any either "settling" (from dispersal to residency), "departure" (from residency
#' to dispersal), or both
#' @param method.arima a vector of methods to fit arima models; any of "CSS-ML", "ML", or "CSS. For more
#' details, see `?arima`. By default all of them are fit.
#' @param condition.settling.to.departure logical. If TRUE, the departure date is set as the lower limit
#' for estimating the settling date. This is used only if both "departure" and "settling" dates 
#' are estimated.
#' @param ... other arguments for the `arima` function.
#' 
#' @return a list with two elements, one for "departure" and other for "settling", based on the estimation
#' of a arima fitting method that does not raise errors or warnings during model fit. Each of these elements
#' is itself a list with three values:  
#' -  time: time of the departure or settling, in units of hours, days, or whichever is the unit for the 
#' argument `Time`;
#' - loglik: the loglikelihood of the fit from ll.finddate;
#' - date: the date, in POSIX.ct format, corresponding to the time of departure or settling
findDepartSettlingDate <- function(Time, X, Y, Date = NULL, 
                                   lower = min(Time), upper = max(Time),
                                   event = c("departure", "settling")[1], 
                                   methods.arima = c("CSS-ML", "ML", "CSS"),
                                   condition.settling.to.departure = FALSE,
                                   ...) {
  # check if the input arguments are right
  
  # initialize output list
  depart.settle <- list(departure = NULL, settling = NULL)
  
  # loop for one or both "departure" and "settling" events
  for(event.type in event) {
    
    # look for departure date using the select methods
    # I tried putting this thing to handle errors and warnings below into a function but it did not work
    fitted.list <- list()
    
    # find the date using optimize using each arima fitting method
    fitted.list <-
      #with(disperser.d, 
      lapply(methods.arima,
             function(method, ...) {
               warn <- err <- NULL
               
               # check whether the initial date should be conditioned to the departure date, for settling
               if(event.type == "settling" & length(event) == 2 & condition.settling.to.departure == TRUE) {
                 lower.date <- depart.settle$departure$time
               } else {
                 lower.date <- lower
               }
               
               res <- withCallingHandlers(
                 tryCatch(optimize(ll.finddate, Time = Time, X = X, Y = Y, maximum = TRUE,
                                   lower = lower.date, upper = upper, event = event.type,
                                   method.arima = method, ...), 
                          error=function(e) {
                            err <<- conditionMessage(e)
                            NULL
                          }), 
                 warning=function(w) {
                   warn <<- append(warn, conditionMessage(w))
                   invokeRestart("muffleWarning")
                 })
               list(res = res, warn = warn, err = err)
             })
      #)
    
    # get errors
    error.list <- lapply(fitted.list, function(x) x[["err"]])
    no.errors <- sapply(error.list, is.null)
    
    # get warnings
    warning.list <- lapply(fitted.list, function(x) x[["warn"]])
    no.warnings <- sapply(warning.list, is.null)
    
    # get the possibility with no errors or warnings, if possible
    # if not possible, take the first one without errors
    index <- which(no.errors & no.warnings)
    if(length(index) > 0) {
      index <- index[1]
    } else {
      index <- which(no.errors)
      if(length(index) > 0) {
        index <- index[1]
      } else {
        stop(paste("All three methods ('CSS-ML', 'ML', and 'CSS') gave errors for the event", event.type))
      }
    }
    
    # take results
    depart.settle[[event.type]] <- fitted.list[[index]]$res
    # rename the output to ease interpretation
    names(depart.settle[[event.type]]) <- c("time", "loglik")
    
    # get exact date, if this is present in the arguments
    if(!is.null(Date)) depart.settle[[event.type]]$date <- Date[Time == floor(depart.settle[[event.type]]$time)] else
      depart.settle[[event.type]]$date <- NULL
  }
  
  return(depart.settle) 
}
  
  


#' Tentative function to find both settling and departure dates using optim
#' It did not work, I am not sure how to use `optim` yet
findSettlingDepartDate.optim <- function(Time, X, Y, 
                                   event = c("settling", "departure")[1], 
                                   ...) 
  optim(par = 1, ll.finddate, Time = Time, X = X, Y = Y,
        lower = min(Time), upper = max(Time), event = event, ...)

# Not used anymore; I incorporated the possibility to plot
# fitted values also in the same function, below
#
# ScanTrack <- function(time, x, y = NULL, ...)
# {
#   if(is.null(y)) if(is.complex(x)){y <- Im(x); x <- Re(x)} else if(ncol(x) == 2){y <- x[,2]; x <- x[,1]}
#   
#   layout(rbind(c(1,2), c(1,3)))
#   plot(x,y,asp=1, type="o", pch=19, col=rgb(0,0,0,.5), cex=0.5, ...)
#   plot(time,x, type="o", pch=19, col=rgb(0,0,0,.5), xaxt="n", xlab="", cex=0.5, ...)
#   plot(time,y, type="o", pch=19, col=rgb(0,0,0,.5), cex=0.5, ...)
# }

#' Scan a track and plot it
#' 
#' This function plots x,y coorinates of a track, as well as x and y as function
#' of time. It also distinguishes locations between dispersal (red) and ranging 
#' (black) phases and identifies departure and settlement dates, 
#' if these arguments are used.
#' 
#' @param time day/time of monitoring (could be hour of monitoring or any counter of periods, given a fix rate?)
#' @param x vector of x or longitude coordinates; can also be a `data.frame` 
#' with (x,y) coordinates in the columns or an imaginary number in the 
#' format z = x + 1i*y.
#' @param y vector of y or latitude coordinates, if x represent only the x coordinates.
#' @param dep.time departure day/time, in case this was estimated; the track will be
#' divided into before (ranging) and after departure (dispersal) phases.
#' @param setl.time settling day/time, in case this was estimated; the track will be
#' divided into before (dispersal) and after settling (ranging) phases. In case
#' `dep.time` is also given, the dispersal phase will be plotted between the 
#' departure and settling days.
#' @param ... additional arguments for the `plot` function.
#' 
#' @return maximum likelihood of the model 
ScanTrack <- function(time, x, y = NULL, dep.time = NULL, setl.time = NULL, 
                      unit = c("", "hours", "days")[1], ...)
{
  if(is.null(y)) 
    if(is.complex(x)) { y <- Im(x); x <- Re(x)} else 
      if(ncol(x) == 2){y <- x[,2]; x <- x[,1]} else
        stop("y is NULL, please provide coordinates to the x argument either as 
             a vector of imaginary numbers or as a data.frame with (x, y) in its columns")
  
  layout(rbind(c(1,2), c(1,3)))
  plot(x, y, asp = 1, type = "o", pch=19, col=rgb(0,0,0,.5), cex=0.5, ...)
  if(!is.null(dep.time)) {
    if(is.null(setl.time)) { 
      points(x, y, type = "o", asp = 1, pch = 19, col = (time > dep.time) + 1, cex=0.5, ...) 
    } else {
      points(x, y, type = "o", asp = 1, pch = 19, col = (time < setl.time & time > dep.time) + 1, 
             cex=0.5, ...) 
    }
  } else {
    if(!is.null(setl.time)) points(x, y, type = "o", asp = 1, pch = 19, 
                                  col = (time < setl.time) + 1, cex=0.5, ...)
  }
  # points(x, y, asp = 1, pch = 19, col = (time > setl.time) + (time < dep.time) + 1, cex=0.5, ...)
  # plot(x,y,asp=1, type="o", pch=19, col=(time > setl.time) + (time < dep.time) + 1, cex=0.5, ...)
  # lines(x,y,asp=1, col = 1, col=(time > setl.time) + (time < dep.time) + 1)
  par(mar = c(c(5, 4, 1, 1) + 0.1))
  plot(time, x, type = "o", pch = 19, col = rgb(0,0,0,.5), xaxt = "n", xlab = "", cex = 0.5, ...)
  if(!is.null(dep.time)) abline(v = dep.time, col = 2)
  if(!is.null(setl.time)) abline(v = setl.time, col = 2)
  lab <- ifelse(unit == "", "Time", paste0("Time (", unit, ")"))
  plot(time, y, type = "o", xlab = lab, pch = 19, col = rgb(0,0,0,.5), cex = 0.5, ...)
  if(!is.null(dep.time)) abline(v = dep.time, col = 2)
  if(!is.null(setl.time)) abline(v = setl.time, col = 2)
}

#' Compare models of residency and dispersal
#' 
#' This function fits residency and dispersal models using `arima` fits, through the function
#' ll.finddate, calculate their maximum likelihood, and compare the linkihood of the different
#' models, given the movement data and values for departure and settling. The models fit are:
#' (i) residency, (ii) dispersal, (iii) departure (a transition from residency to dispersal),
#' (iv) settling (a transition from dispersal to residency), and (v) depart-settle (three phases,
#' residency-dispersal-residency, separated by departure and settling events).
#' 
#' @param Time time or day of monitoring
#' @param X vector of x or longitude coordinates
#' @param Y vector of y or latitude coordinates
#' @param cp vector of change points between behaviors, in the format c(dep, setl) 
#' (in days or time units counted from the first day of monitoring). 
#' When event is "resident" or "dispersal", these values are ignored, if given. 
#' For "departure" and "settling", only the first or second values are considered, respectively.
#' For event == "depart-settle", both values are considered.
#' @param event type of behaviour or transitions between behaviors to be fit.
#' May be any of "settling" (from dispersal to residency), "departure" (from residency
#' to dispersal), "depart-settle" (residency-dispersal-residency), "resident", or "dispersal".
#' By default all of them are fit.
#' @param method.arima a vector of methods to fit arima models; any of "CSS-ML", "ML", or "CSS. For more
#' details, see `?arima`. By default all of them are fit.
#' @param ... other arguments for the `arima` function.
#' 
#' @return a list with the maximum likelihood of each of the models (argument "event") tested.
compare.dispersal.models <- function(Time, X, Y, cp, 
                                     event = c("settling", "departure", "depart-settle", "resident", "dispersing"), 
                                     method.arima = c("CSS-ML", "ML", "CSS"),
                                     ...) {
  
  # possible events
  possible.events <- c("resident", "dispersing", "departure", "settling", "depart-settle")
  
  # check if the events are valid
  if(any(!sapply(event, function(x) x %in% possible.events))) 
    stop(paste0("Wrong event. The 'event's should be any of those options: '", 
               paste(possible.events, collapse = "','"), "'."))
  
  # list of loglikelihoods for each model
  ll.list <- vector(mode = "list", length = length(possible.events))
  names(ll.list) <- possible.events
  
  # loop for one or both "departure" and "settling" events
  for(event.type in event) {
    
    if(event.type == "settling") change.point <- cp[2] else
      if(event.type == "departure") change.point <- cp[1] else
        change.point <- cp
    
    # look for departure date using the select methods
    # I tried putting this thing to handle errors and warnings below into a function but it did not work
    fitted.list <- list()
    
    # find the date using optimize using each arima fitting method
    fitted.list <-
      #with(disperser.6h,
      lapply(method.arima,
             function(method, ...) {
               warn <- err <- NULL
               
               res <- withCallingHandlers(
                 tryCatch(ll.finddate(Time, X, Y, change.point, 
                                      event = event.type, method.arima = method), 
                          error=function(e) {
                            err <<- conditionMessage(e)
                            NULL
                          }), 
                 warning=function(w) {
                   warn <<- append(warn, conditionMessage(w))
                   invokeRestart("muffleWarning")
                 })
               list(res = res, warn = warn, err = err)
             })
    #)
    
    # get errors
    error.list <- lapply(fitted.list, function(x) x[["err"]])
    no.errors <- sapply(error.list, is.null)
    
    # get warnings
    warning.list <- lapply(fitted.list, function(x) x[["warn"]])
    no.warnings <- sapply(warning.list, is.null)
    
    # get results
    res <- sapply(fitted.list, function(x) x[["res"]])
    
    # get the possibility with no errors or warnings, if possible
    # if not possible, take the first one without errors
    index <- which(no.errors & no.warnings)
    if(length(index) > 0) {
      index <- index[which.max(res[index])]
    } else {
      index <- which(no.errors)
      if(length(index) > 0) {
        index <- index[which.max(res[index])]
      } else {
        warning(paste("All three methods ('CSS-ML', 'ML', and 'CSS') gave errors for the event", event.type))
      }
    }
    
    # take results
    if(length(index) > 0) {
      ll.list[[event.type]] <- fitted.list[[index]]$res
    } else {
      ll.list[[event.type]] <- "ERROR"
    }
  
  }
    
  return(ll.list)
}


# from an old script, I preferred to keep it here
#
# x <- arima.sim(list(ar = .9), 100)
# plot(x, type = "o")
# 
# x.int <- cumsum(x)
# plot(x.int, type = "o")
# 
# arima(x + 1e3, order = c(1,0,0))
# arima(x.int, order = c(1,1,0))
# arima(diff(x.int), order = c(1,0,0))
# 
# require(smoove)
# 
# (cvm.fit <- with(disperser.d,
#                 estimateUCVM(Z = X + 1i*Y, T = Time,
#                              time.units = 'day',
#                              method = 'zLike')))
# cvm.fit
# 
# with(disperser.d[100:349,],
#      findSingleBreakPoint(Z = X + 1i*Y, T = Time, method = "sweep"))

# vignette("smoove")



