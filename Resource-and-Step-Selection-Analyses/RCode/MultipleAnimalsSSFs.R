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

#' 
#' # Setup
#' 
#' Load libraries
#+ warning=FALSE, message=FALSE
library(glmmTMB)
library(tidyverse) 
library(survival)
library(TwoStepCLogit)
library(amt)
library(terra)
library(here)
library(broom) 
library(tictoc)  
options(width=165)
knitr::opts_chunk$set(fig.width=12,fig.height=4.5, error=TRUE)
set.seed(09081940)

#' Record time for running all code
tic("total")

#' ## Environmental data
#' 
#' Read in and process the environmental data layers (elevation, population density. forest)
# Read in data
elevation <- rast(here("Resource-and-Step-Selection-Analyses/Data/elevation/", "ASTER ASTGTM2 Elevation-20100101120000000-0-0.tif"))
landuse <- rast(here("Resource-and-Step-Selection-Analyses/Data/landuse", "landuse_study_area.tif"))
popden <- rast(here("Resource-and-Step-Selection-Analyses/Data/pop_den", "pop_den.tif"))

# Center and scale elevation and popden
elevation[] <- (elevation[] - mean(elevation[], na.rm = TRUE))/sd(elevation[], na.rm = TRUE)
popden[] <- (popden[] - mean(popden[], na.rm = TRUE))/sd(popden[], na.rm = TRUE)

# Reproject rasters to EPSG:5070
landuse <- project(landuse, "EPSG:5070")
elevation <- project(elevation, landuse)
popden <- project(popden, landuse)

# Create a binary layer where landuse classes 41 to 43 represent forest
forest <- landuse %in% 41:43
names(elevation)<-"elevation"
names(forest)<-"forest"
names(popden)<-"popden"

#' Plot the layers
plot(landuse)
plot(elevation)
plot(popden)

#' ## Fisher data
#' 
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

#' Add sex as a column of the nested data frame
dat_all$sex <- c("f", "f", "f", "m", "m", "m")

#' Create tracks with an appropriate coordinate reference system
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
#'

#' # Prepare data for fitting step-selection functions
#' 
#' First, we will create a semi-regular trajectory by subsampling the track 
#' to include locations that are separated by 10 +/- 2 minutes.
dat1 <- dat_all %>% mutate(dat_clean = map(trk, ~ {
  .x %>% track_resample(rate = minutes(10), tolerance = seconds(120))
}))

#' Rather than fit tentative step length and turn angle distributions to 
#' each individual separately, we will pool the data and fit a single gamma
#' distribution to the **pooled** step length data and a single von-Mises
#' distribution to the pooled turn angle data.  
#' 

#' We will then generate random steps using these distributions and include
#' step length, log(step length), and cos(turn angle) in the model as 
#' both fixed and random effects.  That will allow each individual to have its
#' own movement kernel with the fixed effects used to update the movement kernel
#' for a "typical" individual (one with all random effects set to 0).  The random
#' effects will  be used to update the movement kernel for all observed
#' individuals in the data set. These updates can be accomplished using the 
#' [mixedSSA package](https://github.com/smthfrmn/mixedSSA). Here, we will largely
#' focus on the habitat selection parameters and the mean correction factors
#' for the movement parameters (needed to estimate the movement kernel for a 
#' "typical" individual).
#' 
#' 
#' First, create a steps version of the data set using the steps_by_burst 
#' function applied to the data from each individual. We will unnest the
#' data so that we can fit a single set of distributions to the pooled step
#' length and turn angle data. 
dat_ssf <- dat1 %>% 
  mutate(stps = map(dat_clean, ~ .x %>% 
                      steps_by_burst())) %>%
  select(id, sex, stps) %>%
  unnest(stps)

#' Unfortunately, some of the classes were lost in the process, so we have to 
#' add some of the classes back on to the object.
class(dat_ssf) <- c("bursted_steps_xyt", "steps_xyt", "steps_xy", class(dat_ssf)) 

#' Now, we can pool the steps across individuals to fit a single tentative 
#' distribution to the step lengths and turn angles, and use this distribution
#' to generate random steps.  Also, we will extract the covariates at the 
#' end of the observed and random steps.
dat_ssf2 <- dat_ssf %>% random_steps() %>%
  extract_covariates(elevation) %>%
  extract_covariates(popden) %>%
  extract_covariates(forest) %>%
  mutate(
    y = as.numeric(case_),
    id = as.numeric(factor(id)), 
    cos_ta_= cos(ta_),
    log_sl_ = log(sl_),
    step_id = paste(id, step_id_, sep = ":")) %>%
    dplyr::select(id, sex, y, step_id, elevation, popden, forest,
                  cos_ta_, log_sl_, sl_, ta_)
 
#' create forest indicator variable (1 if forest, 0 otherwise)
dat_ssf2$forest <- as.numeric(dat_ssf2$forest)

#' # Models
#' 
#' ## glmmTMB: random slopes model
#' 
#' The process is the same as for the mixed RSF:
#' 
#' 1. Set up model, but do not fit it
#' 2. Set random intercept variance to large fixed value, set other variance components to 0
#' 3. Fit the model
#'
#' However, there are a few differences:
#' 
#' - we will use a Poisson likelihood rather than logistic 
#' - we won't need weights
#' - we will include fixed intercepts for each step_id
#' 
#' In our model below, we will assume all coefficients vary independently 
#' since we have so few individuals. With more individuals, we could estimate
#' covariances between the random effects.  The "(0 + x | id)" is used to 
#' let R know that we want to model random coefficients for x and that these
#' coefficients should vary indepdently of the intercept or other random effects
#' included in the model.
tic("mixed ssf")
fisher.tmp <- glmmTMB(y ~ elevation + popden + forest  + 
                        (1|step_id) +  
                        sl_ + log_sl_ + cos_ta_ +
                        (0+ sl_|id) + (0+log_sl_|id) + (0+cos_ta_|id)+
                        (0 + elevation|id) + (0 + popden|id) + (0+ forest |id)
                      , family=poisson(), data = dat_ssf2,
                      doFit=FALSE)

#' Set variance of random intercept to 10^6 and fit the model
fisher.tmp$parameters$theta[1] <- log(1e3)
nvarparm<-length(fisher.tmp$parameters$theta)
fisher.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
fisher.ssf <- glmmTMB:::fitTMB(fisher.tmp)
summary(fisher.ssf)
toc()

#' ## Using Ts.Estim (a formal two-step approach)
library(TwoStepCLogit)

#' This will use a formal two-step approach to estimating the variances of
#' the random coefficients, accounting for sampling variability.  By defualt,
#' the method assumes the coefficients vary independent of each other, but this
#' can be changed by adding the argument D = "UN"
tic("two-step")
twostep<-Ts.estim(formula = y ~  elevation+ popden+forest+ 
           sl_ + log_sl_ + cos_ta_ +
           strata(step_id) +cluster(id), 
           data = dat_ssf2, 
           random = ~ elevation + popden + forest + sl_ + log_sl_ + cos_ta_, all.m.1=F) 
twostep
toc()

#' ## Individual fits 
#' 
#' Next, we will fit a model separately to each individual and then capture and 
#' process the coefficients in a second stage.
ssffits <- dat_ssf2 %>%  nest(data = -c(id, sex)) %>% 
  mutate(mod = map(data, function(x) (fit_issf(y ~ elevation + popden + forest +
                                        sl_ + log_sl_ + cos_ta_ +strata(step_id), 
                                        data = x))))


#' Now, summarize fits to individual animals
ssffits2 <- ssffits %>%
  mutate(tidy = map(mod, ~ broom::tidy(.x$model)),
         n = map_int(data, nrow)) %>%
  unnest(tidy) %>% 
  dplyr::select(-c(data, mod)) %>%
  filter(term!="(Intercept)")  

# Mean coefficient and SE
se<-function(x){sd(x)/sqrt(length(x))}
ssffits2 %>% filter(term!="(Intercept)") %>% 
  group_by(term) %>%   
  summarize(mean=mean(estimate), se=se(estimate))

#' # Graphical Comparison of Results
#' 
#' ## Mean coefficients
#'  
# Dataset with mean coefficients from individual fits
betas1<-ssffits2 %>% 
  filter(term!="(Intercept)") %>% 
  group_by(term) %>%   
  summarize(mean=mean(estimate), se=se(estimate)) %>%
  mutate(method="indiv_fits") %>%
  select(mean, se, method, term)

# Mean coefficients from mixed model (dropping the intercept)
betas2<-data.frame(mean=summary(fisher.ssf)$coef$cond[-1,1],
                   se = summary(fisher.ssf)$coef$cond[-1,2])
betas2$method<-"Mixed"
betas2$term<-c("elevation", "popden", "forest", "sl_", "log_sl_", "cos_ta_")


# Mean coefficients from TWO-STEP (TS.Estim)
betas3<-data.frame(mean=twostep$beta, method=rep("two_step",3), 
                   se = twostep$se,
                   term=c("elevation", "popden", "forest", "sl_", "log_sl_", "cos_ta_"))

all_means <- rbind(betas1, betas2, betas3)

#'  Plot the results.  Estimates of mean coefficients are similar using all 
#'  3 methods, though the point estimates differ slightly for elevation and popden
#'  due to a few individuals with poorly estimated coefficients that are also 
#'  far from the mean.
#+ fig.width=14, fig.height=6
ggplot(data = all_means, aes(x = method, y = mean)) +
  geom_point() +
  geom_errorbar( aes(x = method, ymin = mean - 1.96*se, ymax = mean + 1.96*se),
                 width=0.2, linewidth=1) + 
  xlab("") + ylab(expression(hat(beta))) +
    facet_wrap(~term, scales="free")

#' ## Variance parameters
#' 
#'Compare variance of individual estimates to variance component from mixed model.
# Variance of the coefficients from the individual fits (this will be biased high)
varpar1 <- ssffits2 %>% filter(term!="(Intercept)") %>% 
  group_by(term) %>% 
  summarize(sd=sd(estimate)) %>%
  mutate(method = "Individual fits") %>%
  select(sd, method, term)

# Variance estimated by the mixed model
varpar2 <- data.frame(sd = sqrt(as.numeric(VarCorr(fisher.ssf)$cond[-1]))) %>% 
  mutate(method = "Mixed",
         term = c("sl_", "log_sl_", "cos_ta_","elevation", "popden", "forest")) 


# Variance using the formal two step appraoch
varpar3 <- data.frame(sd = diag(sqrt(twostep$D))) %>%
  mutate(method = "two_step",
         term=c("elevation", "popden", "forest", "sl_", "log_sl_", "cos_ta_")) 

# Combine into a single data frame for plotting
allvars <- rbind(varpar1, varpar2, varpar3)

#' Plot the estimated variance parameters. Here we see that var(individual fits) > var(two step) > var(mixed model).
#+ fig.width=14, fig.height=6
ggplot(data = allvars, aes(x = method, y = sd)) +
  geom_point() +
  xlab("") + ylab("Variance") +
  facet_wrap(~term, scales="free")
#' 

#' ## Individual Coefficients
#' 
#' Look at individual coefficients.
#  Individual fits (we will also append asympototic normal-based CI 
# for parameters in models fit to individuals)
ssf_coefs.se<-ssffits2 %>%  
   mutate(conf.high = estimate + 1.96*std.error) %>%
   mutate(conf.low = estimate - 1.96*std.error) %>%
   mutate(method = "indiv_fits") %>%
  dplyr::select(id, term, estimate, method, conf.high, conf.low)
           
# Predicted coefficients from the mixed model
mixed_coefs <- coef(fisher.ssf)$cond$id[,-1] 
mixed_coefs$id <- rownames(mixed_coefs)
mixed_coefs <- mixed_coefs %>% 
  pivot_longer(c('elevation', 'forest', 'popden', "sl_", "log_sl_", "cos_ta_"), 
               names_to="term", values_to="estimate") %>%
  mutate(method = "Mixed") %>%
  mutate(conf.high = NA) %>%
  mutate(conf.low = NA)

# Coefficients from the two-step approach (TS.Estim)
twostepcoefs<-matrix(twostep$beta,
                     nrow(twostep$r.effect), 
                     ncol(twostep$r.effect), byrow=TRUE)+twostep$r.effect
twostepcoefs<-as.data.frame(twostepcoefs)
twostepcoefs$method<-"two_step"
twostepcoefs$id<-rownames(twostep$r.effect)
twostepcoefs<-twostepcoefs %>% 
  pivot_longer(c('elevation', 'popden','forest', 'sl_', 'log_sl_', 'cos_ta_'), 
               names_to="term", values_to="estimate")%>%
  dplyr::select(c(id, term, estimate, method))
twostepcoefs$conf.high <- twostepcoefs$conf.low <- NA

#' Combine individual estimates
allests<-rbind(mixed_coefs, ssf_coefs.se, twostepcoefs)

#' ### Plot results
 
#' Note here that a few coefficients for elevation and popdenshow considerable shrinkage. These are coefficients with large SEs as seen by inspecting the individual fits.
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
   