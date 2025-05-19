library(moveHMM)

# read data from Morales et al 2004
trackData <- read.table("https://esapubs.org/archive/ecol/E085/072/elk_data.txt",sep="\t",header=TRUE)
# remove blank rows
trackData <- trackData[1:735,]

# simplify covariate names
covNames <- colnames(trackData)[7:ncol(trackData)]
covNames <- gsub("meters","",gsub(".","",covNames,fixed=TRUE))
colnames(trackData)[7:ncol(trackData)] <- covNames

#convert locations and covariates from meters to km
trackData$Easting <- trackData$Easting/1000
trackData$Northing <- trackData$Northing/1000
trackData[,covNames]<-trackData[,covNames]/1000 

# prepare data for analysis
colnames(trackData)[1] <- "ID"
elkData <- moveHMM::prepData(trackData,type="UTM",coordNames=c("Easting","Northing"))
plot(elkData, ask=FALSE)

# standardize covariate values
for(i in covNames){
  elkData[[i]] <-scale(elkData[[i]])
}

# specify initial values for optimization
shape <- c(1,2) # step shape (two parameters: one for each state)
scale <- c(0.5,5) # step scale
zeromass0 <- c(0.01,0.001) # step zero-mass (needed to account for zero step lengths)
stepPar0 <- c(shape,scale,zeromass0)
angleMean0 <- c(pi,0) # angle mean
rho <- c(0.2,0.75) # angle concentration
anglePar0 <- c(angleMean0,rho)

doubleSwitch <- moveHMM::fitHMM(data=elkData,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0,stepDist="weibull",angleDist="wrpcauchy")

# t.p.m. formula
formula <- ~ dist_water+dist_swamp+dist_otw+dist_openfor+dist_ntw+dist_mixfor+dist_devel+dist_ddf+dist_conifer+dist_alvar

# fit "double switch with covariates" model
doubleSwitchCovs <- moveHMM::fitHMM(data=elkData,nbStates=2,stepPar0=stepPar0,anglePar0=anglePar0,formula=formula,stepDist="weibull",angleDist="wrpcauchy")

AIC(doubleSwitch,doubleSwitchCovs)

plot(doubleSwitchCovs,ask=FALSE)
moveHMM::plotPR(doubleSwitchCovs)

save.image("Results/elkExample.RData")
