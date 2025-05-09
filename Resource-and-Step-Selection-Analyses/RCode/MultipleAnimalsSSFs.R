#' ---
#' title: "Multiple SSF's with Fisher data"
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
 

#' ## Data prep for fitting step-selection functions
dat1 <- dat_all %>% mutate(dat_clean = map(trk, ~ {
  .x %>% track_resample(rate = minutes(10), tolerance = seconds(120))
}))

dat_ssf <- dat1 %>% 
  mutate(stps = map(dat_clean, ~ .x %>% 
                      steps_by_burst() %>% 
                      random_steps() %>%
                      extract_covariates(forest) %>%
                      extract_covariates(elevation) %>% 
                      extract_covariates(popden))) %>%
  unnest(stps) %>%
  mutate(
    y = as.numeric(case_),
    id = as.numeric(factor(id)), 
    step_id = paste(id, step_id_, sep = ":")) %>%
    dplyr::select(id, sex, y, step_id, elevation, popden, forest)
  dat_ssf
  dat_ssf$forest <- as.numeric(dat_ssf$forest)

#' ### glmmTMB: random slopes model
#' 
#' The process is the same as for the mixed RSF:
#' 
#' 1. Set up model, but do not fit
#' 2. Set random intercept variance to large fixed value, set other variance components to 0
#' 3. Fit the model
#'
#' However, there are a few differences:
#' 
#' - we will use a Poisson likelihood rather than logistic 
#' - we won't need weights
#' - we will include fixed intercepts for each step_id
#' 
tic("mixed ssf")
fisher.tmp <- glmmTMB(y ~ elevation + popden + forest -1 + (1|step_id) + (0 + elevation|id) +
                       (0 + popden|id) + (0+ forest |id), family=poisson(), data = dat_ssf,
                      doFit=FALSE)

#' Set variance of random intercept to 10^6
fisher.tmp$parameters$theta[1] <- log(1e3)
nvarparm<-length(fisher.tmp$parameters$theta)
fisher.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
fisher.ssf <- glmmTMB:::fitTMB(fisher.tmp)
summary(fisher.ssf)
toc()

#' ### Using Ts.Estim
library(TwoStepCLogit)

#' This will use a formal two-step approach to estimating the variances of
#' the random coefficients, accounting for sampling variability.  By defualt,
#' the method assumes the coefficients vary independent of each other, but this
#' can be changed by adding the argument D = "UN"
tic("two-step")
twostep<-Ts.estim(formula = y ~  elevation+ popden+forest+
           strata(step_id) +cluster(id), 
           data = dat_ssf, 
           random = ~ elevation + popden + forest, all.m.1=F) 
twostep
toc()

#' ### Individual fits 
ssffits <- dat_ssf %>% group_by(id) %>% nest %>% 
  mutate(mod = map(data, function(x) (fit_issf(y ~ elevation + popden + forest + strata(step_id), data = x))))


#' Now, summarize fits to individual animals
ssffits2 <- ssffits %>%
  mutate(tidy = map(mod, ~ broom::tidy(.x$model)),
         n = map_int(data, nrow)) %>%
  unnest(tidy) %>% 
  dplyr::select(-c(data, mod)) %>%
  filter(term!="(Intercept)")  

 
#' ### Comparisons
#' 
#' Compare fixed effects from mixed model and TS.Estim (above) to mean coefficient from individual fits
# MEAN COEFFICIENT FROM INDIVIDUAL FITS 
se<-function(x){sd(x)/sqrt(length(x))}
ssffits2 %>% filter(term!="(Intercept)") %>% 
  group_by(term) %>%   
  summarize(mean=mean(estimate), se=se(estimate))

# POPULATION MEAN COEFFICIENT FROM MIXED MODEL
summary(fisher.ssf)$coef$cond

# POPULATION MEAN COEFFICIENT FROM TWO-STEP (TS.Estim)
twostep$beta

#'Compare variance of individual esimtates to variance component from mixed model.
# VARIANCE OF COEFFICIENTs FROM INDIVIDUAL FITS
ssffits2 %>% filter(term!="(Intercept)") %>% 
  group_by(term) %>% 
  summarize(var_elev=var(estimate), sd_elev=sd(estimate))

# VARIANCE OF COEFFICIENTs FROM MIXED MODEL
summary(fisher.ssf)$varcor #var

# VARIANCE OF COEFFICIENTs FROM TWO-STEP (TS.Estim)
sqrt(twostep$D) # sd

#' Look at individual coefficients.  Append asympototic normal-based CI 
#' for parameters in models fit to individuals
ssf_coefs.se<-ssffits2 %>%  
   mutate(conf.high = estimate + 1.96*std.error) %>%
   mutate(conf.low = estimate - 1.96*std.error) %>%
   mutate(method = "indiv_fits") %>%
  dplyr::select(id, term, estimate, method, conf.high, conf.low)
           
# COEFFICIENTs FROM MIXED MODEL
mixed_coefs <- coef(fisher.ssf)$cond$id[,-1]
mixed_coefs$id <- rownames(mixed_coefs)
mixed_coefs <- mixed_coefs %>% 
  pivot_longer(c('elevation', 'forest', 'popden'), 
               names_to="term", values_to="estimate") %>%
  mutate(method = "Mixed") %>%
  mutate(conf.high = NA) %>%
  mutate(conf.low = NA)

# COEFFICIENTS FROM TWO-STEP (TS.Estim)
twostepcoefs<-matrix(twostep$beta,
                     nrow(twostep$r.effect), 
                     ncol(twostep$r.effect), byrow=TRUE)+twostep$r.effect
twostepcoefs<-as.data.frame(twostepcoefs)
twostepcoefs$method<-"two_step"
twostepcoefs$id<-rownames(twostep$r.effect)
twostepcoefs<-twostepcoefs %>% 
  pivot_longer(c('elevation', 'popden','forest'), names_to="term", values_to="estimate")%>%
  dplyr::select(c(id, term, estimate, method))
twostepcoefs$conf.high <- twostepcoefs$conf.low <- NA

#' ### Plot results
#' 
#' Combine individual estimates
allests<-rbind(mixed_coefs, ssf_coefs.se, twostepcoefs)

#' Dataset with mean coefficients
betas1<-ssffits2 %>% 
  filter(term!="(Intercept)") %>% 
  group_by(term) %>%   
  summarize(mean=mean(estimate)) %>%
  mutate(method="indiv_fits")

# POPULATION MEAN COEFFICIENT FROM MIXED MODEL
betas2<-data.frame(mean=summary(fisher.ssf)$coef$cond[,1])
betas2$method<-"Mixed"
betas2$term<-c("elevation", "popden", "forest")


# POPULATION MEAN COEFFICIENT FROM TWO-STEP (TS.Estim)
betas3<-data.frame(mean=twostep$beta, method=rep("two_step",3), 
                   term=c("elevation", "popden", "forest" ))




#' Note here that a few coefficients for elevation and popD show considerable shrinkage. These are coefficients with large SEs as seen by inspecting the individual fits.
#+ fig.width=14, fig.height=6
ggplot(data=allests, 
       aes(x=id, y=estimate, col=method))+
  geom_point(size=3.5, position=position_dodge(width=0.3))+
  xlab("")+ylab(expression(hat(beta)))+facet_wrap(~term, scales="free")


#' With SEs and dropping the two-step method
#+ fig.width=14, fig.height=6
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
"#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(data=allests[allests$method!="two_step",], 
       aes(x=id, y=estimate, col=method))+
  geom_point(size=3.5)+
  xlab("")+ylab(expression(hat(beta)))+facet_wrap(~term, scales="free")+
  geom_hline(aes(yintercept=mean, col=method), betas1)+
  geom_hline(aes(yintercept=mean, col=method), betas2)+
  scale_colour_manual(values=cbp1)+
  theme(text = element_text(size=20))+
  geom_errorbar( aes(x=id, ymin=conf.low, ymax=conf.high),
                 width=0.2, size=1, data=ssf_coefs.se)

#' Total Elapsed time      
toc()


#' ## Document Footer	
   
#' Session Information:	
#' 	
sessionInfo()	  
   