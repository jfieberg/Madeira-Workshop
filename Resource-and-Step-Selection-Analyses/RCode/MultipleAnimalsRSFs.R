#' ---
#' title: "Multiple RSF's with Fisher data"
#' author: "John Fieberg & Johannes Signer"
#' date: ""
#' output:
#'  html_document:
#'    toc: yes
#'    toc_float: true
#' ---

#' **Purpose**:  demonstrate methods for analyzing data from multiple animals, allowing for multiple random slopes. 
#' 
#' **NOTE** the mixed models here take considerable time to run (you may want to use caching to speed up if you plan to run multiple times). Change the following line, below, to enable caching:
#' opts_chunk$set(fig.width=12,fig.height=4.5, error=TRUE,cache = FALSE).
#' 
#' #### Preamble
#' 
#' Load libraries
#+ warning=FALSE, message=FALSE
library(glmmTMB)
library(tidyverse) 
library(survival)
library(TwoStepCLogit)
library(amt)
library(here)
library(broom) 
library(tictoc)  
options(width=165)
knitr::opts_chunk$set(fig.width=12,fig.height=4.5, error=TRUE)
set.seed(09081940)

#' Record time for running all code
tic("total")


#' Load the fisher data  
#+ eval=TRUE
dat <- read_csv(here("Resource-and-Step-Selection-Analyses/Data", "fisher_data.csv")) %>% 
  filter(id %in% c(1465, 1466, 1072, 1078, 1016, 1469)) 
 
#' We will make use of nested data frames in R, which will allow us to apply functions
#' separately to data from each individual using the map functions in the Purr pacakge.
#' If you are unfamiliar with `purrr` syntax, may want to view one or more of the tutorials, 
#' below, or make use of the [purrr cheat sheet](https://github.com/rstudio/cheatsheets/blob/master/purrr.pdf).
#'
#' - http://www.rebeccabarter.com/blog/2019-08-19_purrr/
#' - https://www.r-bloggers.com/2020/05/one-stop-tutorial-on-purrr-package-in-r/
#' - https://jennybc.github.io/purrr-tutorial/index.html
dat_all <- dat %>% nest(data=c(x,y,t)) 
dat_all$sex <- c("f", "f", "f", "m", "m", "m")

#' Include sex of each animal and create tracks with an appropriate coordinate reference system
#' using the amt package
dat_all <- dat_all %>% 
  mutate(trk = map(data, function(d) {
    make_track(d, x, y, t, crs = 4326) %>% 
      transform_coords(crs_to = 5070)
  }))
dat_all

#' Summarize sampling rates, 
dat_all %>% mutate(sr = map(trk, summarize_sampling_rate)) %>% 
  dplyr::select(id, sr) %>% unnest(cols=c(sr))

#' 10 minutes seems to appropriate for all animals.
#' Resample the track to 10 minutes with a tolerance of 2 minutes.

dat1 <- dat_all %>% mutate(dat_clean = map(trk, ~ {
  .x %>% track_resample(rate = minutes(10), tolerance = seconds(120))
}))

#' Read in data layers
elevation <- rast(here("Resource-and-Step-Selection-Analyses/Data/elevation/", "ASTER ASTGTM2 Elevation-20100101120000000-0-0.tif"))
landuse <- rast(here("Resource-and-Step-Selection-Analyses/Data/landuse", "landuse_study_area.tif"))
popden <- rast(here("Resource-and-Step-Selection-Analyses/Data/pop_den", "pop_den.tif"))

# Center and scale elevation and popden
elevation[] <- (elevation[] - mean(elevation[], na.rm = TRUE))/sd(elevation[], na.rm = TRUE)
popden[] <- (popden[] - mean(popden[], na.rm = TRUE))/sd(popden[], na.rm = TRUE)

# Reproject rasters to EPSG:5070
landuse <- project(landuse, "EPSG:5070")
elevation <- project(elevation, "EPSG:5070")
popden <- project(popden, "EPSG:5070")

# Create a binary layer where landuse classes 41 to 43 represent forest
forest <- landuse %in% 41:43
names(elevation)<-"elevation"
names(forest)<-"forest"
names(popden)<-"popden"

#' # Resource Selection Functions (RSF)
#' 
#' ## Data development for RSF
#' 
#' Now start with an RSF by creating random points per animal within the animal's MCP
dat_rsf <- dat1 %>% mutate(rp = map(dat_clean, ~ .x %>% random_points() %>%
                                      extract_covariates(forest) %>%
                                      extract_covariates(elevation) %>% 
                                      extract_covariates(popden))) %>%
                                      unnest(rp) %>%
                                      dplyr::select(id, sex, case_, elevation, popden, forest)

                           

#' Change id column, to 1:6
dat_rsf$id <- as.numeric(factor(dat_rsf$id))

#' Make response numeric (required for INLA)
dat_rsf$y <- as.numeric(dat_rsf$case_)
dat_rsf$forest <- as.numeric(dat_rsf$forest)

#' We use a weighted likelihood for to fit the RSF. To this end, we need
#'  to create a variable for the weights, where used points (`case_ = TRUE`)
#'   keep weight 1, and available points (`case_ = FALSE`) obtain a large weight $W$ (here $W=1000$):
#+ echo=TRUE, message=FALSE
dat_rsf$w <- 1000^(1 - dat_rsf$case_)


#' ### Individual fits
#'
#' Lets start by fitting "full" models to individuals and then look
#'  to see if coeficients are correlated (high values of one coefficient 
#'  correspond to high/low values of another coefficient)
rsffits <- dat_rsf %>% nest(data=!c(id,sex)) %>% 
  mutate(mod = map(data, function(x){glm(case_ ~  forest + elevation + popden,
                                         data = x, weight=w,family = binomial)})) 


#' Now, summarize fits to individual animals
rsffits2 <- rsffits %>%
  dplyr::mutate(tidy = purrr::map(mod, broom::tidy),
                n = purrr::map(data, nrow) %>% simplify()) %>%
  unnest(tidy) %>% dplyr::select(-c(data, mod))
rsffits2

#' We can do statistics on statistics now
se<-function(x){sd(x)/sqrt(length(x))}
statsonstats<-rsffits2 %>% group_by(sex, term) %>% 
  filter(term %in% c("elevation", "forest", "popden"))%>%
  summarize(mean=mean(estimate), se=se(estimate), n = n()) %>%
  mutate(up=mean+qt(0.025, df = n-1)*se, low=mean+qt(0.975, df = n-1)*se)

ggplot(statsonstats, aes(x=sex, y=mean))+
  geom_point(size=2) + 
  geom_errorbar(aes(ymin=low, ymax=up),  width=0.2, size=1)+
  facet_wrap(~term, scales="free")

#' We can look at whether the coefficients tend to covary (e.g., do individuals
#' that select for high population density also select/avoid forest?).
#' With only 6 individuals, we have little information here but this plot
#' could be useful with more individuals to see if we should allow random 
#' effects to covary. 
rsf_coefs_wide <- rsffits2 %>% 
  dplyr::select(-(std.error:p.value)) %>%
  pivot_wider(., names_from=term, values_from=estimate)

#+ fig.height=8, fig.width=8
rsf_coefs_wide %>% dplyr::select(elevation, popden, forest) %>%    
   GGally::ggpairs(.) 

#' ### glmmTMB independent random slopes model
#' 
#' Process: 
#' 
#' 1. Set up model, but do not fit
#' 2. Set random intercept variance to large fixed value, set other variance components to 0
#' 3. Fit the model
#'

#' Here, we are assuming that the coefficient for elevation, popden, and forest
#' all vary independent of one another.  We could allow them to covary if we 
#' switched the first bullet to the second:
#' 
#' -  (0 + elevation|id) + (0 + popden |id) +(0 +forest | id)
#' -  (0 + elevation + popden + forest | id)
#' 
#' But, this would require 3 variance parameters and 3 covariance parameters.
#' That would be too complex a model given we only have 6 individuals!
tic("mixed mod")
fisher.tmp <- glmmTMB(case_ ~ elevation + forest + popden -1 + (1|id) + 
                        (0 + elevation|id) + (0 + popden |id) +(0 +forest | id),
                      family=binomial(), data = dat_rsf,
                      doFit=FALSE, weights = w)

#' Set variance of random intercept to 10^6
fisher.tmp$parameters$theta[1] <- log(1e3)
nvarparm<-length(fisher.tmp$parameters$theta)
fisher.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
fisher.rsf.ind <- glmmTMB:::fitTMB(fisher.tmp)
summary(fisher.rsf.ind)
toc()

#' ### Comparisons
#' 
#' Compare fixed effects (above) to mean coefficient from individual fits. These are quite similar.
# MEAN COEFFICIENT FROM INDIVIDUAL FITS 
rsffits2 %>%   
  filter(term!="(Intercept)") %>% 
  group_by(term) %>% 
  summarize(mean=mean(estimate), se.mean=sd(estimate)/sqrt(n()))

# POPULATION MEAN COEFFICIENT FROM MIXED MODEL
summary(fisher.rsf.ind)$coef$cond

#' We can also compare the variance of individual estimates to the estimated
#' variance components from mixed model. Note, the variance of the individual
#' coefficients will be larger (too big) since it will capture both
#' true variability and sampling variability.

# VARIANCE OF COEFFICIENTS FROM INDIVIDUAL FITS
rsffits2 %>% filter(term!="(Intercept)") %>% 
  group_by(term) %>% 
  summarize(var_elev=var(estimate), sd_elev=sd(estimate))

# VARIANCE COMPONENT FROM MIXED MODEL
summary(fisher.rsf.ind)$varcor

#' Look at individual coefficients. Create a data set with the individual
#' coefficients and asympototic normal-based CI for parameters from the
#' individual fits
rsf_coefs.se<- rsffits2  %>% filter(term!="(Intercept)") %>%
  mutate(conf.high = estimate + 1.96*std.error, 
         conf.low = estimate - 1.96*std.error,
         method="indiv_fits") %>%
  rename(Estimate=estimate) %>%  
  dplyr::select(id, term, Estimate, method, conf.high, conf.low) 
class(rsf_coefs.se)<-"data.frame" 


# COEFFICIENTS FROM MIXED MODEL
mixed_coefs <- coef(fisher.rsf.ind)$cond$id[,-1]
mixed_coefs$id <- rownames(mixed_coefs)
mixed_coefs <- mixed_coefs %>% 
  pivot_longer(c('elevation', 'forest', 'popden'), 
               names_to="term", 
               values_to="Estimate") %>% 
  mutate(method = "Mixed")
mixed_coefs$conf.high <- mixed_coefs$conf.low <- NA

#' ### Plot results

#' Individual estimates
allests<-rbind(mixed_coefs, rsf_coefs.se)

#' Dataset with mean coefficients
betas1<-rsffits2 %>% 
  filter(term!="(Intercept)") %>% 
  group_by(term) %>% 
  summarize(mean=mean(estimate)) %>%
  mutate(method = "indiv_fits")

# POPULATION MEAN COEFFICIENT FROM MIXED MODEL
betas2<-data.frame(mean=summary(fisher.rsf.ind)$coef$cond[,1])
betas2$method<-"Mixed"
betas2$term<-c("elevation", "forest", "popden")


#' Plot everything with nice colors
#' Note that because of the large sample size, individual coefficients show very little shrinkage towards the overall mean.
#+ fig.width=14, fig.height=6
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(data=allests, 
       aes(x=id, y=Estimate, col=method))+
  geom_point(size=3.5)+
  xlab("")+ylab(expression(hat(beta)))+facet_wrap(~term, scales="free") +
  geom_hline(aes(yintercept=mean, col=method), betas1)+
  geom_hline(aes(yintercept=mean, col=method), betas2)+
  scale_colour_manual(values=cbp1)+
  theme(text = element_text(size=20))+
  geom_errorbar( aes(x=id, ymin=conf.low, ymax=conf.high),
               width=0.2, size=1) 

    
 
#' Total Elapsed time      
toc()


#' ## Document Footer	
   
#' Session Information:	
#' 	
sessionInfo()	  
   