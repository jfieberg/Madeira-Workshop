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
library(raster)
library(survival)
library(TwoStepCLogit)
library(amt)
library(here)
library(broom) 
library(tictoc) 
library(glmmTMB) 
options(width=165)
knitr::opts_chunk$set(fig.width=12,fig.height=4.5, error=TRUE)


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
elevation <- raster(here("Resource-and-Step-Selection-Analyses/Data/elevation/", "ASTER ASTGTM2 Elevation-20100101120000000-0-0.tif"))
landuse <- raster(here("Resource-and-Step-Selection-Analyses/Data/landuse", "landuse_study_area.tif"))
popden <- raster(here("Resource-and-Step-Selection-Analyses/Data/pop_den", "pop_den.tif"))
elevation[] <- (elevation[] - mean(elevation[], na.rm = TRUE))/sd(elevation[], na.rm = TRUE)
popden[] <- (popden[] - mean(popden[], na.rm = TRUE))/sd(popden[], na.rm = TRUE)
 
landuse <- raster::projectRaster(landuse, crs = "+init=epsg:5070")
elevation <- raster::projectRaster(elevation, crs = "+init=epsg:5070")
popden <- raster::projectRaster(popden, crs = "+init=epsg:5070")

forest <- raster::match(landuse, 41:43, nomatch = 0) > 0
raster::plot(forest)
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
se<-function(x){sqrt(var(x))/length(x)}
statsonstats<-rsffits2 %>% group_by(sex, term) %>% 
  filter(term %in% c("elevation", "forest", "popden"))%>%
  summarize(mean=mean(estimate), se=se(estimate)) %>%
  mutate(up=mean+1.96*se, low=mean-1.96*se)

ggplot(statsonstats, aes(x=sex, y=mean))+
  geom_point(size=2) + 
  geom_errorbar(aes(ymin=low, ymax=up),  width=0.2, size=1)+
  facet_wrap(~term, scales="free")

#' Keep just the estimates for later comparisons
rsf_coefs<-rsffits2 %>% 
  dplyr::select(-(std.error:p.value))

rsf_coefs_wide <- rsf_coefs %>% 
  pivot_wider(., names_from=term, values_from=estimate)

#' Explore coefficients using a pair plot
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p
}

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

#' THis is the model we'd like to fit (where coefficients may be correlated), 
#' but you will wait a long time. And, it requires too many variance parameters 
#' for the small number of individuals

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
rsf_coefs %>% filter(term!="(Intercept)") %>% 
  group_by(term) %>% 
  summarize(mean=mean(estimate), se.mean=sd(estimate)/sqrt(n()))

# POPULATION MEAN COEFFICIENT FROM MIXED MODEL
summary(fisher.rsf.ind)$coef$cond

#'Compare variance of individual esimtates to variance component from mixed model. Note, the variance of the individual coefficients will be larger (too big) since it includes true variability and sampling variability.

# VARIANCE OF COEFFICIENTS FROM INDIVIDUAL FITS
rsf_coefs %>% filter(term!="(Intercept)") %>% 
  group_by(term) %>% 
  summarize(var_elev=var(estimate), sd_elev=sd(estimate))

# VARIANCE COMPONENT FROM MIXED MODEL
summary(fisher.rsf.ind)$varcor

#' Look at individual coefficients
# COEFFICIENTS FROM INDIVIDUAL FITS
rsf_coefs2 <- rsf_coefs %>% filter(term!="(Intercept)") %>% arrange(term, id)
rsf_coefs2$method<-"indiv_fits"   
rsf_coefs2<-rsf_coefs2 %>%
  rename(Estimate=estimate)%>%
  dplyr::select(c(id, term, Estimate, method))

# COEFFICIENTS FROM MIXED MODEL
mixed_coefs <- coef(fisher.rsf.ind)$cond$id[,-1]
mixed_coefs$id <- rownames(mixed_coefs)
mixed_coefs <- mixed_coefs %>% 
  pivot_longer(c('elevation', 'forest', 'popden'), names_to="term", values_to="Estimate")
mixed_coefs$method<-"Mixed"

#' ### Plot results

#' Individual estimates
allests<-rbind(mixed_coefs, rsf_coefs2)

#' Dataset with mean coefficients
betas1<-rsf_coefs %>% 
  filter(term!="(Intercept)") %>% 
  group_by(term) %>% 
  summarize(mean=mean(estimate))

betas1$method="indiv_fits"

# POPULATION MEAN COEFFICIENT FROM MIXED MODEL
betas2<-data.frame(mean=summary(fisher.rsf.ind)$coef$cond[,1])
betas2$method<-"Mixed"
betas2$term<-c("elevation", "forest", "popden")

#' Now, get asympototic normal-based CI for individual fits
rsf_coefs.se<-rsffits2  %>% filter(term!="(Intercept)") %>%
  mutate(conf.high = estimate + 1.96*std.error, 
         conf.low = estimate - 1.96*std.error,
         method="indiv_fits") %>%
  rename(Estimate=estimate) %>%  
  dplyr::select(id, term, Estimate, method, conf.high, conf.low) 
class(rsf_coefs.se)<-"data.frame" 
 

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
               width=0.2, size=1, data=rsf_coefs.se)

    

#' # SSF model

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
    step_id = paste0(id, step_id_, sep = "-")) %>%
    dplyr::select(id, sex, y, step_id, elevation, popden, forest)
  dat_ssf


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
tic("two-step")
twostep<-Ts.estim(formula = y ~  elevation+ popden+forest+
           strata(step_id) +cluster(id),data = dat_ssf, random = ~ elevation+popden+forest, all.m.1=F) 
twostep
toc()

#' ### Individual fits 
ssffits <- dat_ssf %>% group_by(id) %>% nest %>% 
  mutate(mod = map(data, function(x) (fit_issf(y ~ elevation + popden + forest + strata(step_id), data = x))))


#' Now, summarize fits to individual animals
ssffits <- ssffits %>%
  mutate(tidy = map(mod, ~ broom::tidy(.x$model)),
         n = map_int(data, nrow))

ssffits %>% unnest(tidy) %>% filter(term!="(Intercept)") %>% print(., n=1000)

#' Keep just coefficients for comparisons
ssf_coefs<-ssffits %>% unnest(tidy) %>% 
  dplyr::select(-(std.error:p.value)) 

#' ### Comparisons
#' 
#' Compare fixed effects from mixed model and TS.Estim (above) to mean coefficient from individual fits
# MEAN COEFFICIENT FROM INDIVIDUAL FITS 
se<-function(x){sd(x)/sqrt(length(x))}
ssf_coefs %>% filter(term!="(Intercept)") %>% 
  group_by(term) %>%   
  summarize(mean=mean(estimate), se=se(estimate))

# POPULATION MEAN COEFFICIENT FROM MIXED MODEL
summary(fisher.ssf)$coef$cond

# POPULATION MEAN COEFFICIENT FROM TWO-STEP (TS.Estim)
twostep$beta

#'Compare variance of individual esimtates to variance component from mixed model.
# VARIANCE OF COEFFICIENTs FROM INDIVIDUAL FITS
ssf_coefs %>% filter(term!="(Intercept)") %>% 
  group_by(term) %>% 
  summarize(var_elev=var(estimate), sd_elev=sd(estimate))

# VARIANCE OF COEFFICIENTs FROM MIXED MODEL
summary(fisher.ssf)$varcor #var

# VARIANCE OF COEFFICIENTs FROM TWO-STEP (TS.Estim)
sqrt(twostep$D) # sd

#' Look at individual coefficients
# COEFFICIENTs FROM INDIVIDUAL FITS
ssf_coefs2 <- ssf_coefs %>% filter(term!="(Intercept)") %>% arrange(term, id)
ssf_coefs2$method<-"indiv_fits"                       
names(ssf_coefs2)[5]<-"Estimate"
ssf_coefs2<-ssf_coefs2 %>%
  dplyr::select(c(id, term, Estimate, method))
          
# COEFFICIENTs FROM MIXED MODEL
mixed_coefs <- coef(fisher.ssf)$cond$id[,-1]
mixed_coefs$id <- rownames(mixed_coefs)
mixed_coefs <- mixed_coefs %>% 
  pivot_longer(c('elevation', 'forest', 'popden'), names_to="term", values_to="Estimate")
mixed_coefs$method<-"Mixed"


# COEFFICIENTS FROM TWO-STEP (TS.Estim)
twostepcoefs<-matrix(twostep$beta, nrow(twostep$r.effect), ncol(twostep$r.effect), byrow=TRUE)+twostep$r.effect
twostepcoefs<-as.data.frame(twostepcoefs)
twostepcoefs$method<-"two_step"
twostepcoefs$id<-rownames(twostep$r.effect)
twostepcoefs<-twostepcoefs %>% 
  pivot_longer(c('elevation', 'popden','forest'), names_to="term", values_to="Estimate")%>%
  dplyr::select(c(id, term, Estimate, method))

#' ### Plot results
#' 
#' Combine individual estimates
allests<-rbind(mixed_coefs, ssf_coefs2,twostepcoefs)

#' Dataset with mean coefficients
betas1<-ssf_coefs %>% filter(term!="(Intercept)") %>% group_by(term) %>%   summarize(mean=mean(estimate))
betas1$method="indiv_fits"

# POPULATION MEAN COEFFICIENT FROM MIXED MODEL
betas2<-data.frame(mean=summary(fisher.ssf)$coef$cond[,1])
betas2$method<-"Mixed"
betas2$term<-c("elevation", "popden", "forest")


# POPULATION MEAN COEFFICIENT FROM TWO-STEP (TS.Estim)
betas3<-data.frame(mean=twostep$beta, method=rep("two_step",3), term=c("elevation", "popden", "forest" ))


#' Now, get asympototic normal-based CI for individual fits
ssf_coefs.se<-ssffits %>%   unnest(tidy) 
ssf_coefs.se<-ssf_coefs.se %>% dplyr::select(-data) %>% dplyr::select(-mod)
ssf_coefs.se$conf.high<-ssf_coefs.se$estimate + 1.96*ssf_coefs.se$std.error
ssf_coefs.se$conf.low<-ssf_coefs.se$estimate - 1.96*ssf_coefs.se$std.error
class(ssf_coefs.se)<-"data.frame"
names(ssf_coefs.se)[4]<-"Estimate"
ssf_coefs.se<-ssf_coefs.se %>% filter(term!="(Intercept)")
ssf_coefs.se$method="indiv_fits"

#' Note here that a few coefficients for elevation and popD show considerable shrinkage. These are coefficients with large SEs as seen by inspecting the individual fits.
#+ fig.width=14, fig.height=6
ggplot(data=allests, 
       aes(x=id, y=Estimate, col=method))+
  geom_point(size=3.5, position=position_dodge(width=0.3))+
  xlab("")+ylab(expression(hat(beta)))+facet_wrap(~term, scales="free")


#' With SEs and dropping the two-step method
#+ fig.width=14, fig.height=6
ggplot(data=allests[allests$method!="two_step",], 
       aes(x=id, y=Estimate, col=method))+
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
   