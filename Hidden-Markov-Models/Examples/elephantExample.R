library(momentuHMM)
library(readr)

URL <- "https://datarepository.movebank.org/server/api/core/bitstreams/77e38af2-b5fc-4c10-9e6d-f36b2adf9ea9/content"
rawData <- readr::read_csv(url(URL))

# select and rename relevant columns
rawData <- rawData[,c(11,3,4,5,6)]
colnames(rawData) <- c("ID","time","lon","lat","temp")

# only keep first track
rawData <- subset(rawData,ID==unique(ID)[1])

# convert times from factors to POSIX
rawData$time <- as.POSIXct(rawData$time,tz="GMT")

# project to UTM coordinates using package sp
library(sp)
llcoord <- SpatialPoints(rawData[,3:4], 
                         proj4string=CRS("+proj=longlat +datum=WGS84"))
utmcoord <- spTransform(llcoord,CRS("+proj=utm +zone=30 ellps=WGS84"))

# add UTM locations to data frame
rawData$x <- attr(utmcoord,"coords")[,1]
rawData$y <- attr(utmcoord,"coords")[,2]

crwOut<-crawlWrap(rawData,timeStep="hour",theta=c(6.855, -0.007),fixPar=c(NA,NA))
plot(crwOut,ask=FALSE)

# create momentuHMMData object from crwData object
elephantData <- momentuHMM::prepData(data=crwOut, covNames="temp")

# add cosinor covariate based on hour of day
elephantData$hour <- as.integer(strftime(elephantData$time, format = "%H", tz="GMT"))

# acf plot of step lengths
acf(elephantData$step[!is.na(elephantData$step)],lag.max=300)

# label states
stateNames <- c("encamped","exploratory")
# distributions for observation processes
dist = list(step = "gamma", angle = "wrpcauchy")

# initial parameters
Par0_m1 <- list(step=c(100,500,100,200),angle=c(0.3,0.7))

# fit model
m1 <- momentuHMM::fitHMM(data = elephantData, nbStates = 2, dist = dist, Par0 = Par0_m1, 
             estAngleMean = list(angle=FALSE), stateNames = stateNames)


# formula for transition probabilities
formula <- ~ temp * cosinor(hour, period = 24)

# initial parameters (obtained from nested model m1)
Par0_m2 <- getPar0(model=m1, formula=formula)

# fit model
m2 <- momentuHMM::fitHMM(data = elephantData, nbStates = 2, dist = dist, Par0 = Par0_m2$Par, 
             beta0=Par0_m2$beta, stateNames = stateNames, formula=formula)


# formulas for parameters of state-dependent observation distributions
DM <- list(step = list(mean = ~ temp * cosinor(hour, period = 24),
                       sd = ~ temp * cosinor(hour, period = 24)),
           angle = list(concentration = ~ temp))

# initial parameters (obtained from nested model m2)
Par0_m3 <- getPar0(model=m2, formula=formula, DM=DM)

# fit model
m3 <- momentuHMM::fitHMM(data = elephantData, nbStates = 2, dist = dist, Par0 = Par0_m3$Par, 
             beta0 = Par0_m3$beta, DM = DM, stateNames = stateNames,
             formula = formula)

# calculate and compare AIC
AIC(m1,m2,m3)

# decode most likely state sequence
states <- momentuHMM::viterbi(m3)
# derive percentage of time spent in each state
timeInStates(m3)

# plot results for model m3
plot(m3, plotCI = TRUE, covs = data.frame(hour=12))

# plot satellite image for model m3
# commented out because Google now requires an API key for Google Maps
# plotSat(m3,zoom=8,col=c("firebrick3","seagreen4"),projargs = proj4string(utmcoord),ask=FALSE)

# compute pseudo-residuals for the steps and the angles
pr <- momentuHMM::pseudoRes(m3)

# plot the ACF of step pseudo-residuals
acf(pr$stepRes[!is.na(pr$stepRes)],lag.max = 300)

save.image("Results/elephantExample.RData")
