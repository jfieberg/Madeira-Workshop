
###########################################################################
#' Day 3. Continuous-time movement models
#'
#' Speed, distance, diffusion
#' Example workflow using the 'ctmm' R package
#' See: Noonan et al. 2019. Movement Ecology, 7(1), 1–15.
###########################################################################

library(ctmm)
library(here)

# Data preparation: -------------------------------------------------------

data(buffalo)

# North-up projection:
projection(buffalo) <- median(buffalo)

# Select only the first buffalo:
data_buffalo <- buffalo[[1]]

# Units operator:
?`%#%`

1 %#% 'day' # day in seconds
1 %#% 'year' # year in seconds

# Consider only the first week of data:
data_buffalo <- data_buffalo[
  data_buffalo$t <= data_buffalo$t[1] + 1 %#% "week", ]
plot(data_buffalo, col = color(data_buffalo, by = "time"), error = FALSE)

# Select best fit:
guess <- ctmm.guess(data_buffalo, interactive = FALSE)
fit <- ctmm.select(data_buffalo, guess, trace = 2)
# or load model fits from ctmm.select():
load(here::here("D4-ctmm", "data", "fit_buffalo_cilla.rda"))
fit <- fitList[[1]]

# Speed estimate here is RMS Gaussian:
summary(fit)

# Gaussian (regular speed - not RMS):
speed(fit)

# Non-parametric speed estimation:
ctsd <- speed(data_buffalo, fit)
ctsd

# Impact of coarsening the data: ------------------------------------------

SUB <- data_buffalo
fit.SUB <- fit

# Removing every other time:
SUB <- SUB[as.logical(1:nrow(SUB) %% 2), ]
plot(SUB, col = color(SUB, by = "time"), error = FALSE)
fit.SUB <- ctmm.select(SUB, fit.SUB, trace = 2)
# RMS Gaussian:
summary(fit)
summary(fit.SUB)
# Gaussian (regular speed - not RMS):
speed(fit)
speed(fit.SUB)
# Non-parametric speed estimation:
ctsd
speed(SUB, fit.SUB)
# repeat until data become too coarse

# keep in mind the stationary assumption of the model
# see the appendix of Noonan et al.

# Population-level inferences: --------------------------------------------

help("meta")

# Load in the fitted movement models:
load(here::here("D4-ctmm", "data", "fits_buffalo.rda"))

# Estimate mean spead for each animal:
ctsdList <- list()
for (i in seq_along(length(buffalo))) {
  ctsdList[[i]] <- speed(buffalo[[i]], fitList[[i]])
}
names(ctsdList) <- names(buffalo)
# save(ctsdList, file = here::here("data", "speeds_buffalo.rda"))
load(here::here("D4-ctmm", "data", "speeds_buffalo.rda"))

meta(ctsdList, sort = TRUE)

# Instantaneous speeds: ---------------------------------------------------

inst_speeds <- speeds(buffalo[[1]], fitList[[1]])
head(inst_speeds) # standard units (m/s)
