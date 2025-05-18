
###########################################################################
#' Day 3. Continuous-time movement models
#'
#' Resource selection functions
#' Example workflow using the 'ctmm' R package
#' See: Alston et al. 2023 Methods in Ecology and Evolution 4:2 643-654
###########################################################################

library(ctmm)
library(here)

## 1. Read data: ----------------------------------------------------------

# Load lowland tapir movement dat and tree cover raster:
load(here::here("D3-ctmm", "data", "tapir.rda"))
length(tapir) # 29 individuals

# Select first individual:
i <- 1
dat <- tapir[[i]]

# Plot location dat and tree cover:
plot(dat, error = 2, R = treecover,
     main = "Lowland tapir under tree cover")

# Environmental covariates must be raster objects and in a named list:
R <- list(tree = treecover)
# see raster::as.factor() for categorical variables

# 2. Fit an autocorrelated movement model: --------------------------------

# Estimate starting model (isotropic, includes error):
guess <- ctmm.guess(dat,
                    CTMM = ctmm(error = TRUE,
                                isotropic = TRUE),
                    interactive = FALSE)

# Select best-fitting autocorrelated movement model:
fit <- ctmm.select(dat, guess, trace = 2)
summary(fit)

# Save fitted model:
# save(fit, file = here::here("data", "fit_tapir-iso.rda"))
load(here::here("D3-ctmm", "data", "fit_tapir-iso.rda"))

# 3. Run Autocorrelated Kernel Density Estimate (AKDE): -------------------

AKDE <- akde(dat, fit, weights = TRUE)
plot(dat, error = 2, UD = AKDE, R = treecover,
     col.grid = NA, main = "AKDE")

# Compare with IID model (no autocorrelation) [comparison only]
iid <- ctmm.fit(dat, CTMM = ctmm(isotropic = TRUE))
KDE <- akde(dat, iid)

# 4. RSF (IID vs autocorrelated): -----------------------------------------

help("rsf.fit")

# Check assigned weights without autocorrelation:
plot(dat$timestamp,
     mean(KDE$DOF.area) * KDE$weights,
     xlab = "Time", ylab = "Weight")

# iRSF without autocorrelation:
RSF.IID <- rsf.fit(dat, KDE, R = R)
# (iterates until the 1% error threshold)

# Check assigned weights *with* autocorrelation:
plot(dat$timestamp,
     mean(AKDE$DOF.area) * AKDE$weights,
     xlab = "Time", ylab = "Weight")

# iRSF with autocorrelation:
RSF <- rsf.fit(dat, AKDE, R = R)
# save(RSF, file = here::here("data", "rsf_tapir.rda"))
load(here::here("D3-ctmm", "data", "rsf_tapir.rda"))
summary(RSF)

# with faster integration method (if you don’t need time-dependent model)
RSF <- rsf.fit(dat, AKDE, R = R, integrator = "Riemann")
summary(RSF)

# 5. RSF model selection: -------------------------------------------------

RSFs <- rsf.select(
  dat, AKDE, R = R,
  formula = ~ I(sqrt(tree)) + tree + I(tree^2),
  integrator = "Riemann",
  verbose = TRUE,
  trace = TRUE
)
summary(RSFs)

# View the top selected model:
RSF <- RSFs[[1]]
summary(RSF)

# with multiple individuals:
help("mean.ctmm")
# mean() will give you population-level information

# Interpreting effect sizes:

treecover # scaled between 0–1
# relative selection between no tree cover (0) and full tree cover (1)
exp( summary(RSF)$CI[1,] * (sqrt(1) - sqrt(0)) )

# 6. Generate iRSF-informed distributions: --------------------------------

# iRSF distribution that was fit:
help("agde")

AGDE <- agde(dat, RSF, R = R)
plot(dat, AGDE, main = "iRSF")

# 7. Generate suitability maps: -------------------------------------------

help("suitability")

SUIT <- suitability(dat, CTMM = RSF, R = R, grid = AKDE)
names(SUIT) # brick with 3 layers (lower, point estimate, upper)
plot(dat, error = 2, R = SUIT[["est"]],
     col.grid = NA, main = "Suitability")

# 8. Generate RSF-informed AKDE: ------------------------------------------

# Re-estimate AKDE while incorporating RSF
RAKDE <- akde(dat, RSF, R = R, weights = TRUE)
plot(dat, error = 2, UD = RAKDE, col.grid = NA, main = "iRSF-AKDE")
