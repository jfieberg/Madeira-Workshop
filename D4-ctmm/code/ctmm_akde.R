
###########################################################################
#' Day 3. Continuous-time movement models
#'
#' Home range estimation
#' Example workflow using the 'ctmm' R package
#' See: Silva et al. 2021 Methods in Ecology and Evolution 13(3), 534-544
###########################################################################

library(ctmm)
library(here)

# Load pre-run objects:
load(here::here("D4-ctmm", "data", "hr.rda"))

# Data preparation: -------------------------------------------------------

# Loading tracking datasets:
data(buffalo)
data(gazelle)

data_buffalo <- buffalo$Pepper # or buffalo[[4]]
head(data_buffalo)

data_gazelle <- gazelle[[11]]
head(data_gazelle)

# Plotting locations:
plot(data_buffalo, col = "red", lwd = 3)
plot(data_gazelle, col = "blue", lwd = 3)

nrow(data_buffalo)
nrow(data_gazelle)

# Range residency assumption: ---------------------------------------------

level <- 0.95 # we want to display 95% confidence intervals
xlim <- c(0, 1 %#% "day") # to create a window of one day

# Checking for the range residency assumption:
svf1 <- variogram(data_buffalo)
par(mfrow = c(1,2))
plot(svf1, fraction = 0.5, level = level)
abline(v = 1, col = "red", lty = 2) # adding a line at ~one month
plot(svf1, xlim = xlim, level = level)
par(mfrow = c(1,1))

# Model selection: --------------------------------------------------------
# Selecting the best-fit movement model through model selection:

# Calculate an automated model guesstimate:
guess1 <- ctmm.guess(data_buffalo, interactive = FALSE)
summary(guess1)

# Automated model selection, starting from guesstimate:
start_time <- Sys.time()
fit1 <- ctmm.select(data_buffalo, guess1,
                    method = "pHREML", verbose = TRUE)
## reminder: it will default to pHREML if no method is specified.
Sys.time() - start_time # Time difference of 2.435679 mins
summary(fit1)

plot(svf1, CTMM = fit1[[1]],
     units = TRUE, fraction = 0.5, level = c(0.95, 0.50), 
     col = "black", col.CTMM = "red")

summary(fit1[[1]])

# Home range estimator: ---------------------------------------------------
# Feeding a movement model into the home range estimator

# Run an area-corrected AKDE (default):
AKDE1 <- akde(data_buffalo, fit1[[1]], debias = TRUE)

# (For comparison purposes) assuming IID:
fit1_iid <- ctmm.fit(data_buffalo)
KDE1 <- akde(data_buffalo, fit1_iid)
plot(svf1, CTMM = fit1_iid,
     units = TRUE, fraction = 0.5, level = c(0.95, 0.50), 
     col = "black", col.CTMM = "red")

summary(AKDE1, level.UD = 0.95)$CI # 95% AKDE area
summary(KDE1, level.UD = 0.95)$CI # 95% KDE area

( 1 - summary(KDE1)$CI[1,2] / summary(AKDE1)$CI[1,2] ) * 100
# ~39% underestimation by KDE vs. AKDE

# Creating an extent that includes both UDs at the 95% CI level:
ext <- extent(list(AKDE1, KDE1))

# Plotting KDE and AKDE side-by-side:
par(mfrow = c(1, 2))
plot(data_buffalo, UD = KDE1, ext = ext)
title(expression("KDEc"))
plot(data_buffalo, UD = AKDE1, ext = ext)
title(expression("AKDEc"))
par(mfrow = c(1, 1))

# Mitigation measures: ----------------------------------------------------
# Evaluating additional biases, applying mitigation measures

## Irregular representation in time: --------------------------------------

plot(data_buffalo, lwd = 3)

# Sample sizes:
summary(AKDE1)$DOF["area"] # effective sample size of Pepper
nrow(data_buffalo) # absolute sample size of Pepper

# Plot all sampling intervals:
dt.plot(data_buffalo) # Pepper (buffalo[[4]])
abline(h = 2 %#% "hours", col = "red")

# Minimum sampling interval:
"minutes" %#% min(diff(data_buffalo$t))

# Calculate wAKDE:

start_time <- Sys.time()
AKDE1_weighted <- akde(data_buffalo, fit1[[1]], weights = TRUE)
Sys.time() - start_time # Time difference of 1.298815 mins

summary(AKDE1_weighted)$CI # 95% home range area (weighted)

ext <- extent(list(AKDE1_ML, AKDE1, AKDE1_weighted), level = 0.95)

# Plotting pHREML (with and without weights) side-by-side:
par(mfrow = c(1,2))
plot(data_buffalo, UD = AKDE1, ext = ext)
title(expression("pHREML AKDE"["C"]))
plot(data_buffalo, UD = AKDE1_weighted, ext = ext)
title(expression("pHREML wAKDE"["C"]))
par(mfrow = c(1,1))

## Low sample sizes: ------------------------------------------------------

plot(data_gazelle, lwd = 3)
svf2 <- variogram(data_gazelle)
plot(svf2)

guess2 <- ctmm.guess(data_gazelle, interactive = FALSE)
fit2 <- ctmm.select(data_gazelle, guess2, method = "pHREML")
summary(fit2)

AKDE2 <- akde(data_gazelle, fit2)
summary(AKDE2)$CI

# Sample sizes:
summary(AKDE2)$DOF["area"] # effective sample size
nrow(data_gazelle) # absolute sample size

1/summary(fit2)$DOF["area"]^2 # expected order of pHREML bias

start_time <- Sys.time()
fit2_boot <- ctmm.boot(data_gazelle, fit2, error = 0.01, trace = 2)
# save(fit2_boot, file = here::here("data", "hr_bootstrap.rda"))
#! Note: this function takes a long time to run!
( total_time <- Sys.time() - start_time )
# Time difference of 43.93944 mins

load(here::here("D4-ctmm", "data", "hr_bootstrap.rda"))
summary(fit2_boot)

1/summary(fit2_boot)$DOF["area"]^3 # expected order of bias

AKDE2_boot <- akde(data_gazelle, fit2_boot, weights = TRUE)
summary(AKDE2)$CI
summary(AKDE2_boot)$CI

( 1 - summary(AKDE2_boot)$CI[1,2] / summary(AKDE2)$CI[1,2] ) * 100
# AKDE (without bootstrapping) underestimates by ~6%

ext <- extent(list(AKDE2, AKDE2_boot), level = 0.95)

# Plotting pHREML and bootstrapped-pHREML side-by-side:
par(mfrow = c(1, 2))
plot(data_gazelle, UD = AKDE2, ext = ext)
title(expression("pHREML AKDE"["C"]))
plot(data_gazelle, UD = AKDE2_boot, ext = ext)
title(expression("Bootstrapped pHREML wAKDE"["C"]))
par(mfrow = c(1, 1))

# save(fit1,
#      fit2,
#      fit2_boot,
#      KDE1,
#      AKDE1,
#      AKDE2,
#      AKDE1_weighted,
#      AKDE2_boot,
#      file = here::here(data", "hr.rda"))
# save(fit2_boot, file = here::here("data", "bootstrap.rda"))

# Population-level inferences: --------------------------------------------

# Load pre-run objects:
load(here::here("D4-ctmm", "data", "meta.RData"))

# Fit movement models for all individuals:
start_time <- Sys.time()
fitList <- list()
for(i in seq_along(buffalo)) {
  guess <- ctmm.guess(buffalo[[i]], interactive = FALSE)
  fitList[[i]] <- ctmm.select(buffalo[[i]], guess, trace = 2)
}
Sys.time() - start_time

# Calculate AKDEs on a consistent grid:
hrList <- akde(buffalo, fitList, trace = 2)

# Plot AKDEs:
pal <- color(hrList, by = "individual")
plot(hrList,
     col.DF = pal,
     col.level = pal,
     col.grid = NA,
     level = NA)

# What is the mean home range area of an average individual:

meta(hrList,
     col = c(pal, "black"), 
     verbose = TRUE,
     sort = TRUE) 

# Compare mea home range between groups:

meta(list(south = hrList[1:3],
          north = hrList[4:6]),
     plot = TRUE, 
     verbose = TRUE) 

# What is the population range?

MEAN <- mean(hrList) # distribution of the sample
# (assumes no "untracked" individuals)
plot(buffalo, MEAN)

PKDE <- pkde(buffalo, hrList) # distribution of the population
# (extrapolating to a larger population; including "untracked")
plot(buffalo, PKDE)

ext <- extent(list(MEAN, PKDE))
pal <- c("red", "black", "blue")

par(mfrow = c(1, 2))
plot(buffalo, MEAN, col = pal, ext = ext)
title("mean()")
plot(buffalo, PKDE, col = pal, ext = ext)
title("pkde()")
par(mfrow = c(1, 1))

summary(MEAN)$CI
summary(PKDE)$CI

# save(fitList,
#      hrList,
#      MEAN,
#      PKDE,
#      file = here::here("data", "meta.rda"))

# Extra: ------------------------------------------------------------------

# Fit multiple movement models with parallelization:

library(doParallel)
library(foreach)

?detectCores()
cl <- parallel::makeCluster(detectCores() - 1)
doParallel::registerDoParallel(cl)

# Create a "fitting" function:

fitting_model <- function(i) {
  guess <- ctmm.guess(dat[[i]], interactive = FALSE)
  ctmm.select(dat[[i]], guess, verbose = TRUE, trace = 2)
}

# Use the fitting function in a foreach loop and
# parallel backend by using %dopar%:

start_time <- Sys.time()
fitList <- foreach(i = 1:length(dat), 
                   .packages = "ctmm") %dopar% {
                     fitting_model(i)
                   }
Sys.time() - start_time

parallel::stopCluster(cl)
