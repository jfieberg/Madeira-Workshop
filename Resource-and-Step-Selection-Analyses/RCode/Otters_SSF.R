#' ##  Habitat selection of otters
#'  Coded example for an SSF analysis  
#'  
#'  Authors: S. Muff, J. Signer, J. Fieberg
#'
#'  This code replicates the analysis presented in Section 4.2 of Muff, Signer, Fieberg (2022). Accounting for individual‐specific variation in habitat‐selection studies: Efficient estimation of mixed‐effects models using Bayesian or frequentist computation, Journal of Animal Ecology, 89, p.80-92, https://doi.org/10.1111/1365-2656.13087
#'
#'  The code is adapted from our data and code repository: https://conservancy.umn.edu/handle/11299/204737
#'
#'  The otter data are taken from Weinberger et al (2016). Flexible habitat selection paves the way for a recovery of otter populations in the European Alps, Biological Conservation 199, p. 88-95, https://doi.org/10.1016/j.biocon.2016.04.017
#'
#' Load libraries and read in data
#+ warning=FALSE, message=FALSE
library(survival)
library(TwoStepCLogit)
library(glmmTMB)
library(tictoc)
options(width=160)# setting output width in html report
dat <-  read.csv(here("Resource-and-Step-Selection-Analyses/Data", "d_otter.csv"))

head(dat)
str(dat)

########################################################################################
### Part 1: Data wrangling 
########################################################################################

#' NAT1, REST1 and STAU1 are the three factor levels of the factor variable "habitat type", encoded as dummy variables, where
#'
#' - NAT1: natural habitat (reference category)
#' - REST1: residual water
#' - STAU1: a reservoir

#' Further, the two continuous variables in the model are:
 
#' - Sohlbrei: the river width
#' - Breaks_Dis: step length
 
#' Finally, `Loc` is the binary response variable that indicates if a habitat point was used (1) or available (0).

#' Some data manipulation:

#' Add numerical variable for animals:
dat$ANIMAL_ID <- as.numeric(as.factor(dat$NA_ANIMAL))

#' Stratum ID is given as "NA_ID" in the data; 
#' It is easier to have sequential enumeration, so let's generate a new 
#' stratum-ID variable str_ID:
d.map <- data.frame(NA_ID=unique(dat$NA_ID),str_ID=1:length(unique(dat$NA_ID)))
dat$str_ID <- d.map[match(dat$NA_ID,d.map$NA_ID),"str_ID"]
dat <- dat[order(dat$str_ID),]

#' Scale and center the two continuous variables river width (Sohlenbrei) and step length (Breaks_Dis)
dat$Sohlenbrei <- scale(dat$Sohlenbrei)
dat$Breaks_Dis <- scale(dat$Breaks_Dis)

#' Getting to know the data better:
str(dat)
#' Note that str_ID is the stratum ID. In each stratum there is exactly one realized (Loc=1) and 9 available locations (Loc=0):
dat[1:30,]



########################################################################################
### Random effects models
########################################################################################

###  Random effects SSFs using glmmTMB

#' Set up the model without fitting it. Independent random effects are added as  (0 + STAU1 | ANIMAL_ID), for animal-specific slopes of "STAU1", for example:

TMBStruc = glmmTMB(Loc ~ -1 + STAU1 + REST1 + Sohlenbrei +  
                     Breaks_Dis +  (1|str_ID) + 
                     (0 + STAU1 | ANIMAL_ID) + 
                     (0 + REST1 | ANIMAL_ID)   + 
                     (0 + Sohlenbrei | ANIMAL_ID), 
                   family=poisson, data=dat, doFit=FALSE) 

#' Variance parameters in the model (there are 4, the first we want to fix to a large value)
TMBStruc$parameters$theta
nvarparm<-length(TMBStruc$parameters$theta)

#' Since (1|str_ID) is the first parameter in the formula, it is the one we have to fix.
#' Therefore, set the value of the (log of the) standard deviation of the first random effect (here (1|str_ID)) to the fixed value:
TMBStruc$parameters$theta[1] = log(1e3) 

#' Tell glmmTMB not to change the first standard deviation (NA), all other values are freely estimated (and are different from each other)
TMBStruc$mapArg = list(theta=factor(c(NA,1:(nvarparm-1))))

#' Fit the model and look at the summary:
glmm.TMB.random <- glmmTMB:::fitTMB(TMBStruc)
summary(glmm.TMB.random)

#' 95\% CIs for fixed and random effects (standard deviations) are obtained via the confint() function:
confint(glmm.TMB.random)



#'There is a one-step way to carry out the regression with newer versions of glmmTMB:

glmm.TMB.random = glmmTMB(Loc ~ -1 + STAU1 + REST1 + Sohlenbrei +  
                            Breaks_Dis +  (1|str_ID) + 
                            (0 + STAU1 | ANIMAL_ID) + 
                            (0 + REST1 | ANIMAL_ID)   + 
                            (0 + Sohlenbrei | ANIMAL_ID),
                          family=poisson, data=dat,
                          map=list(theta=factor(c(NA,1:3))),
                          start=list(theta=c(log(1e3),0,0,0))
)
 
######################################################
# A note for those using an older version of glmmTMB
######################################################


#' It it has been a problem in a previous version of glmmTMB that the confint() function only showed a table of the length equal to the number of parameters estimated. As the variance for str_ID was not estimated but fixed to 10^6, but was still listed, the last variance component (here the one for (0 + Sohlenbrei | ANIMAL_ID)) was not shown. If you face this problem, you can solve the issue by moving the component (1|str_ID) to the last position in the formula of the above code, and then use 
#'
#' `TMBStruc$parameters$theta[4] = log(1e3)` 
#' 
#' `TMBStruc$mapArg = list(theta=factor(c(1:3, NA)))`
