
###########################################################################
#' Day 3. Continuous-time movement models
#'
#' Introduction to continuous-time movement models
#' Example workflow using the 'ctmm' R package
#' 
###########################################################################

library(ctmm)
library(ggplot2)

# Help files & vignettes:
help(package = "ctmm")
browseVignettes(package = "ctmm")
help("ctmm-FAQ", package = "ctmm")
browseURL("https://groups.google.com/g/ctmm-user")

# Development branch of ctmm (more recent than CRAN):
# remotes::install_github("ctmm-initiative/ctmm")

# 'ctmm' point-and-click app - if you know anyone that doesn't use R
# remotes::install_github("ctmm-initiative/ctmmweb")
# ctmmweb::app()

# 'ctmm' MoveApps
browseURL("https://www.moveapps.org/apps/browser?q=ctmm")

# 'movedesign' app - for study design of animal tracking projects
# remotes::install_github("ecoisilva/movedesign")
# movedesign::run_app()

# Import and visualize: ---------------------------------------------------

# You can either get data through MoveBank, or
# import data with as.telemetry()
help("as.telemetry")

# Loading data from Movebank CSV (which can be compressed)
# dat <- as.telemetry('file.zip')
# you can also import from a move object, data.frame, etc.

# Load buffalo data:
data(buffalo)
?buffalo

# List of buffalo telemetry objects:
class(buffalo)

# Number of individuals:
length(buffalo)

class(buffalo[[1]])
head(buffalo[[1]])

# Names of individuals:
names(buffalo)

# Summary of buffalo dataset:
summary(buffalo)

help("plot.telemetry")

# Plot all individuals:
plot(buffalo) # but all the same color
plot(buffalo, col = color(buffalo, by = "individual"))

# Projections: ------------------------------------------------------------

projection(buffalo)
# you want a projection that is locally flat to minimize distortion

# By default, as.telemetry() will choose a two-point equidistant projection,
# which is safer for migratory species, but does not preserve North = up.
# The algorithm can be found in:
ctmm:::median_longlat
# and automates the estimation of k=2 geometric median (robust) clusters

# Show north on plot:
compass()

# Center the projection on the geometric median of the data:
projection(buffalo) <- median(buffalo)
projection(buffalo)

# now North=up, which is fine for this dataset
plot(buffalo, col = pal, main = "Azimuthal-equidistant projection")
compass()

# Variogram: --------------------------------------------------------------

names(buffalo)
data_buffalo <- buffalo$Cilla # selecting buffalo named Cilla

# Plot telemetry object:
plot(data_buffalo, main = "Cilla")

# Color by time:
pal <- color(data_buffalo, by = "time")
plot(data_buffalo, col = pal)
# easier to see migrations/dispersals

ggplot(data_buffalo) +
  geom_point(aes(x, y, color = timestamp)) +
  scale_color_viridis("Timestamp", trans = "time")

# Calculate a variogram object from the telemetry object:
help("variogram")
svf <- variogram(data_buffalo)
plot(svf, main = "Variogram")
# on average how far apart (in distance^2) given a time lag
# between any two points

# More accurate CIs (too slow for larger datasets):
svf <- variogram(data_buffalo, CI = "Gauss")
# frequently you want to zoom in to the beginning of the variogram
# so you can plot with zoom slider:
zoom(svf)

# Model selection: --------------------------------------------------------

# Guesstimate model parameters:
help("ctmm.guess")
# variogram will be calculated automatically (with default arguments)

# this is interactive mode
ctmm.guess(data_buffalo, variogram = svf)

# this is noninteractive mode
guess <- ctmm.guess(data_buffalo, variogram = svf, interactive = FALSE)
summary(guess)

# Automated model selection:
help("ctmm.select")

fitList <- ctmm.select(data_buffalo, guess, trace = 2, verbose = TRUE)
# current candidate models: OUF, OUf, OUΩ, IOU, BM, IID
# save(fitList, file = here::here("data", "fit_buffalo_cilla.rda"))
load(here::here("D3-ctmm", "data", "fit_buffalo_cilla.rda"))

# Summarize results:
summary(fitList)

# IID was not attempted because the
# nested-model hierarchy is OUF -> OU -> IID
# (we can include the IID model manually)
fitList[["IID anisotropic"]] <- ctmm.fit(data_buffalo)
fitList[["IID"]] <- ctmm.fit(data_buffalo, ctmm(isotropic = TRUE))

# Summarize results (now including IID model):
summary(fitList)

# Looking at IID anisotropic model only:
summary(fitList$`IID anisotropic`)

# Compare mean and covariance to data:
plot(data_buffalo,
     fitList$`IID anisotropic`,
     main = "IID Gaussian Distribution")

# Compare empirical variogram to that of model:
zoom(svf, fitList$`IID anisotropic`, main = "IID Variogram")

# Calculate residuals:
res_iid <- residuals(data_buffalo, fitList$`IID anisotropic`)
plot(res_iid, main = "IID Residuals")

# Calculate correlogram of residuals:
ACF_iid <- correlogram(res_iid, fast = FALSE)
zoom(ACF_iid, main = "ACF of 'IID' Residuals")

summary(fitList)
summary(fitList[[1]]) # the first element is the best-fit model
# area here is Gaussian area
# speed here is Gaussian RMS speed

summary(data_buffalo)
(4.967566 %#% 'months') / (7.505372 %#% 'days')
help("%#%")

summary(fitList[[1]])
sigfig(summary(fitList[[1]])$CI[2,] * 3)

plot(data_buffalo, fitList[[1]],
     main = "Anisotropic Gaussian") # anisotropic
plot(data_buffalo, fitList[[2]],
     main = "Isotropic Gaussian") # isotropic

zoom(svf, fitList[[1]], main = "OUF Variogram")
# not perfect, but much better

# Plot residuals (for best fit model):
res <- residuals(data_buffalo, fitList[[1]])
plot(res, main = "OUF Residuals")

#  Calculate correlogram of residuals (for best fit model):
ACF <- correlogram(res, fast = FALSE)
zoom(ACF, main = "ACF of 'OUF' Residuals")

# Non-resident models: ----------------------------------------------------

guess <- ctmm.guess(data_buffalo, ctmm(range = FALSE), interactive = FALSE)
fitList_nonrange <- ctmm.select(data_buffalo, guess,
                                verbose = TRUE, trace = 2)

summary(fitList_nonrange)
zoom(svf, fitList_nonrange[["IOU anisotropic"]])
zoom(svf, fitList_nonrange[["BM anisotropic"]])

# Likelihoods/AICs cannot be compared
summary(c(fitList, fitList_nonrange))
