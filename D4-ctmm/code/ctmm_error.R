
###########################################################################
#' Day 3. Continuous-time movement models
#'
#' Error modelling
#' Example workflow using the 'ctmm' R package
#' See: https://www.biorxiv.org/content/10.1101/2020.06.12.130195v2.full
###########################################################################

library(ctmm)
library(here)

# Workflow: --------------------------------------------------------------

#' STEP A: Do you need to model error?
#' How do the scales of error compare to the scales of movement?

#' STEP B: Do you have calibrated data or calibration data?
#' Calibration data can be collected or opportunistic.
#' Without calibration data, you should supply a prior.

#' IF YOUR DATA NEEDS CALIBRATION:
#' STEP C: What columns do you have in your data?
  #' DOP values? (HDOP, VDOP, PDOP, GDOP, ...)
  #' location classes?
  #' number of satellites?
  #' time-to-fix timeout?
#' ctmm will try to pick out the best data

#' STEP D: Error model selection

#' STEP E: Run error-informed movement analyses

## 1. Read data: ----------------------------------------------------------

data(turtle) # Wood turtle (Glyptemys insculpta) GPS tracking data
# tracked in Virginia, USA

names(turtle) # first two are calibration data - not turtles

# Calibration data:
head(turtle[[1]])
# HDOP: horizontal dilution or precision -- proportional to RMS error
# Location class: 3D
plot(turtle[1:2], col = rainbow(2))

# Animal location data:
plot(turtle$F231, main = "F231 data, no error model")
uere(turtle)


# 2. Fit error parameters to calibration data: ---------------------------

help("uere.fit") # used to estimate the RMS uere from calibration data

head(turtle[1:2])
uere <- uere.fit(turtle[1:2]) 
# do not run uere.fit on tracking data!

# Estimated error model parameters:
summary(uere)
# tells us how much the recorded position is likely to
# deviate from the true location

# Apply error model to data:
uere(turtle) <- uere
plot(turtle$F231, main = "F231 data, error model")


# 3. What if we are not sure about the error data? ------------------------

## 3.1. Are the HDOP and location class values informative? ---------------

# Load un-calibrated data:
data(turtle)
uereList <- list() # make a list to store error models

# First attempt: let's use everything!
uereList$all <- uere.fit(turtle[1:2]) 

# Second attempt: let's drop the location class information
turtle_test <- turtle[1:2] # copy calibration data
# delete location class column
turtle_test[[1]]$class <- NULL
turtle_test[[2]]$class <- NULL
uere(turtle_test) <- NULL
uereList$HDOP <- uere.fit(turtle_test)

# Third attempt: let's further drop the HDOP values
# delete HDOP column
turtle_test[[1]]$HDOP <- NULL
turtle_test[[2]]$HDOP <- NULL
uereList$none <- uere.fit(turtle_test)

summary(uereList$all) # error-model fit 1 (with all available columns)
summary(uereList$HDOP) # error-model fit 2 (HDOP only)
summary(uereList$none) # error-model fit 2 (none)

# Compare models:
summary(uereList)
# AICc: super-fancy AIC values
# reduced Z-squared statistic (goodness of fit)
# compare to reduced chi-squared statistic (~1 is good)


## 3.2. Are these GPS tags identical? -------------------------------------
# i.e., do the two calibration tags warrant separate error models?

# Calculate tag UEREs:
gpsList <- lapply(turtle[1:2], uere.fit)

# compare calibration parameters
summary(uereList$all) # joint model
summary(gpsList[[1]])
summary(gpsList[[2]])

summary(list(joint = uereList$all, individual = gpsList))

# 4. Error calibration: ---------------------------------------------------

# Calibrate turtle data with best error model:
uere(turtle) <- uereList$all
head(turtle[[3]]) # error columns now in data
# VAR.xy: estimated variance of positional error

# Calculate residuals of calibration data w.r.t best error model
res_all <- lapply(turtle[1:2], residuals)

# Calculate residuals of calibration data w.r.t. worst error model
uere(turtle_test) <- uereList$none
res_none <- lapply(turtle_test, residuals)

# Plot residuals
plot(res_all, ext = extent(RES_none))
plot(res_none, ext = extent(RES_none))

# 5. Error-informed movement analyses: ------------------------------------

names(turtle)
data_turtle <- turtle$F231
plot(data_turtle)

help("outlie") # check for outliers

outliers <- outlie(data_turtle)
plot(outliers)
head(outliers)

# Good location estimates:
VALID <- outliers$speed < 0.05
# biological threshold for this species (wood turtle)
data_turtle <- data_turtle[VALID, ]

# Re-check:
plot(data_turtle)
outliers <- outlie(data_turtle)
plot(outliers)

# Create guesstimate interactively:
ctmm.guess(data_turtle)
# * check the error box

# Create guesstimate non-interactively:
guess <- ctmm.guess(data_turtle,
                    CTMM = ctmm(error = TRUE), interactive = FALSE)

# Fit movement models:
fitList <- ctmm.select(data_turtle, guess, verbose = TRUE, trace = 2)
# Note: argument 'verbose = TRUE' returns all candidate models
# save(fitList, file = here::here("data","fits_turtle.rda"))
load(here::here("D4-ctmm", "data", "fits_turtle.rda"))

# Look at all models:
summary(fitList)

# Look only at best-fit model:
summary(fitList[[1]]) # always the first element on the list

# Compare to prior model:
summary(uere(data_turtle))

# Compare to movement model without error model:
fit <- fitList[[1]]
fit$error <- FALSE
fit_no_error <- ctmm.fit(data_turtle, fit, trace = 2)

summary(fit)$CI
summary(fit_no_error)$CI

hr <- akde(data_turtle, fit)
hr_no_error <- akde(data_turtle, fit_no_error)

# Plotting AKDE with and without error model:

ext <- extent(list(hr, hr_no_error))

par(mfrow = c(1, 2))
plot(data_turtle, error = FALSE, UD = hr_no_error, ext = ext)
title(expression("Without error"))
plot(data_turtle, error = FALSE, UD = hr, ext = ext)
title(expression("With error"))
par(mfrow = c(1, 1))

summary(hr_no_error)
summary(hr)


# Extra. What if you do not have calibration data? ------------------------
# you will need to suply a prior!

# Load un-calibrated data once again:
data(turtle)

# You will need to match the class structure (2D, 3D here):
summary(uere(turtle))

# Supply point estimates:
# 20-meter 2D error at HDOP = 1
# 10-meter 3D error at HDOP = 1
uere(turtle) <- c(20, 10)

# Extract calibration object:
prior <- uere(turtle)
summary(prior) # the default uncertainty is currently zero
prior$DOF

# Set DOF for wide credible intervals:
prior$DOF[] <- 2 # very small effective sample size!
summary(prior)

# Assign prior to data:
uere(turtle) <- prior

# Fit movement model:
guess <- ctmm.guess(turtle[[3]],
                    CTMM = ctmm(error = TRUE), interactive = FALSE)
fit_prior <- ctmm.select(turtle[[3]], guess, trace = 2)
# save(fit_prior, file = here::here("data", "fit_turtle-prior.rda"))
load(here::here("D4-ctmm", "data", "fit_turtle-prior.rda"))

summary(fit_prior)
summary(prior) # compare update to prior

