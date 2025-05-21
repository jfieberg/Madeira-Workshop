
###########################################################################
#' Day 3. Continuous-time movement models
#'
#' Occurrence distribution - Kriging
#' Example workflow using the 'ctmm' R package
#' See: Fleming et al. 2016. Ecology, 97(3), 576-582.
###########################################################################

library(ctmm)
library(here)

load(here::here("D4-ctmm", "data", "occurrence.rda"))

# Data preparation: -------------------------------------------------------

# Loading tracking datasets:
data(buffalo)
data_buffalo <- buffalo[[1]]
plot(data_buffalo, lwd = 3)

# Fit movement model:
guess <- ctmm.guess(data_buffalo, interactive = FALSE)
fitList <- ctmm.select(data_buffalo, guess, verbose = TRUE)

# Q: What is the occurrence distribution?
# A: Given a random time *in the sampling period*, where was the animal

# Q: What is the range distribution?
# A: Where the animal will be at some time in the future/past
#    *under the same behaviors*
# A: Long-term space use *for continuing behaviors*

# Include Brownian motion models
fitList[["BM"]] <- ctmm.fit(data_buffalo,
                            ctmm(tau = Inf, isotropic = TRUE))
# this one is not as commonly used, but let's throw it in
fitList[["BM anisotropic"]] <- ctmm.fit(data_buffalo,
                                        ctmm(tau = Inf))

# You can't compare stationary (IID, OU, OUF)
# and conditionally stationary (BM,IOU) models with likelihood
summary(fitList)
# but you can compare within:
summary(fitList[c("BM", "BM anisotropic")])

svf <- variogram(data_buffalo, CI = "Gauss")
zoom(svf, fitList[[1]]) # the selected model looks okay

# The Brownian motion model looks...
zoom(svf, fitList$BM)

# Create occurrence distribution (OD) using best-fit model:
cilla_od <- occurrence(data_buffalo, fitList[[1]])
plot(cilla_od)

# Conventional (non-dynamic) Brownian bridge
cilla_bbmm <- occurrence(data_buffalo, fitList$BM)

# Compare the two:
ext <- extent(list(cilla_od, cilla_bbmm))
par(mfrow = c(1, 2))
plot(cilla_od, ext = ext)
title("Kriging")
plot(cilla_bbmm, ext = ext)
title("BBMM")
par(mfrow = c(1, 1))

# Range versus occurrence: ------------------------------------------------

# AKDE:
cilla_hr <- akde(data_buffalo, fitList[[1]])

# Comparison plot of home range and occurrence region:
ext <- extent(list(cilla_hr, cilla_od))
par(mfrow = c(1, 2))
plot(data_buffalo, cilla_hr, ext = ext, col.grid = NA)
title("Range")
plot(cilla_od, ext = ext, col.DF = "blue", col.level = NA)
title("Occurrence")
par(mfrow = c(1, 1))

SUB <- SUB2 <- data_buffalo

# Manipulation 1 - coarsen the data:
SUB <- SUB[as.logical(1:nrow(SUB) %% 2),]
cilla_hr.SUB <- akde(SUB, fitList[[1]])
cilla_od.SUB <- occurrence(SUB, fitList[[1]])
par(mfrow = c(1, 2))
# plot(cilla_hr, ext = ext, col.grid = NA, col.level = NA)
# title("Range UD")
# plot(cilla_od, ext = ext, col.level = NA)
# title("Occurrence UD")
plot(cilla_hr.SUB, ext = ext, col.grid = NA, col.level = NA)
title("Range UD (Subset)")
plot(cilla_od.SUB, ext = ext, col.level = NA)
title("Occurrence UD (Subset)")
#' Rerun this block to coarsen the data further.

# Manipulation 2 - truncate the data:
SUB2 <- SUB2[1:(nrow(SUB2)/2),]
cilla_hr.SUB2 <- akde(SUB2, fitList[[1]])
cilla_od.SUB2 <- occurrence(SUB2, fitList[[1]])
par(mfrow = c(1, 2))
# plot(cilla_hr, ext = ext, col.grid = NA, col.level = NA)
# title("Range UD")
# plot(cilla_od, ext = ext, col.level = NA)
# title("Occurrence UD")
plot(cilla_hr.SUB2, ext = ext, col.grid = NA, col.level = NA)
title("Range UD (Subset 2)")
plot(cilla_od.SUB2, ext = ext, col.level = NA)
title("Occurrence UD (Subset 2)")
par(mfrow = c(1, 1))
# Rerun this block to truncate data further.

# save(fitList,
#      cilla_od,
#      cilla_bbmm,
#      cilla_hr,
#      file = here::here("D4-ctmm", "data", "occurrence.rda"))

# range area = predicted space use, given the same behaviors
# occurrence area = uncertainty (sampling dependent and limited to the sampling period)
# neither estimate the amount of space used during the sampling period!
