library(spOccupancy)
library(tidyverse)

# We fit a model similar to the imperfect detection one implemented in Tobler et
# el. 2019, as implemented by the sfMsPGOcc function in the spOccupancy package
# This is an occupancy model with a detection component, so we need repeat
# occasions at each site. We exclude the sites with less than 4 weeks of
# sampling, and divide the detections into monthly counts.

#### CT data ####
allsites_pts <- terra::vect("Data/combined_survey_pts.geojson")
allsites_cov_db <- as.data.frame(allsites_pts)
# campts_db <- filter(as.data.frame(allsites_pts), has_cam==1)

# camvect <- vect("Data/ct_covars.geojson")
# campts_db <- as.data.frame(camvect) %>% mutate(site = gsub("#", "", deployment_id), days = as.numeric(days))
data_raw <- read_csv("Data/ct_records.csv",
                     col_types = cols_only(deployment_id = 'c', common_name = 'c', timestamp = col_datetime("%Y-%m-%d %H:%M") ))
# include only cameras with more than 12 weeks (~3 months)
cams_include <- summarise(data_raw, .by = deployment_id, t=max(timestamp)-min(timestamp)) %>% 
  filter(t>12*7) %>% 
  mutate(site = sub("#","",deployment_id)) %>% 
  select(site) 
DAT <- read_csv("Data/records_wild.csv") %>% inner_join(cams_include)

# for the occupancy model the order of species in the array matters.reorder
# species to put the most common first, and then others that will behave
# differently. This should help identifiability and convergence times.
sp_order <- c("central_american_agouti", "oncilla", "baird's_tapir","greater_grison", "white_lipped_peccary", "ocelot","white_nosed_coati", 
              "guan_black", "tinamou_great", "margay", "puma", "spotted_paca", "tayra","northern_tamandua", "coyote", "crab_eating_raccoon",
              "jaguar", "tinamou_little", "collared_peccary", "guan_crested", "striped_hog_nosed_skunk", "central_american_red_brocket", "common_opossum",
              "jaguarundi", "cottontail_dices", "opossum_four_eyed", "nine_banded_armadillo", "great_curassow")

# model data
y_df <- filter(DAT, !grepl("uid", common_name)) %>% 
  left_join(select(allsites_cov_db,site)) %>% 
  mutate(smplwk = as.numeric(1+floor(difftime(timestamp,min(timestamp),units = "weeks"))), .by = site) %>% 
  mutate(smplmon = 1+smplwk%/%4) %>% 
  filter(max(smplwk)>4, .by = site) %>% 
  count(smplmon, site, common_name) %>% 
  complete(nesting(smplmon, site), common_name, fill = list(n=0)) %>% # make undetected species explicit (only existing sampling time) fill with 0s
  complete(smplmon, site, common_name) %>% # fill missing times/species with NA
  arrange(smplmon,  site, common_name) # make all have the same order (for array creation)

nsp <- length(unique(y_df$common_name))
nst <- length(unique(y_df$site))
nocc <- max(y_df$smplmon)
nsp*nst*nocc
mody <- array(data = as.numeric(y_df$n>0), dim = c(nsp, nst, nocc), 
              dimnames = list(species = unique(y_df$common_name), sites = unique(y_df$site), occ = 1:nocc))
# reorder y according to order established previously
mody <- mody[sp_order,,]
covs_df <- distinct(y_df,site) %>% left_join(campts_db) %>% 
  select(site, easting, northing, ghm1, lc2,alt1, lfdistance_2) %>% 
  rename(x = easting, y = northing, forest = lc2, alt = alt1, ghm = ghm1, Fdist = lfdistance_2) %>% 
  mutate(forest = factor(forest, labels = c("Dense", "Open"))) 
moddata <- list(y = mody, 
                occ.covs = covs_df[,c("forest","alt", "ghm", "Fdist")] ,
                # det.covs = covs_df[,"site"],
                coords = as.matrix(covs_df[,c("x", "y")])
)

# Priors
sitedists <- dist(moddata$coords)# distance between sites
mindist <- min(sitedists)
maxdist <- max(sitedists)
prior.list <- list(phi.unif = list(3/maxdist, 3/mindist))

# initial values
inits.list <- list(beta.comm = 0, alpha.comm = 0, # Community-level occurrence (beta.comm) and detection (alpha.comm) regression coefficients
                   beta = 0, alpha = 0, # species-level occ and det (could be length nsp)
                   tau.sq.beta = 1, tau.sq.alpha = 1, # comm=level variance for occ and det
                   phi = 3/1000, # spatial decay parameter. 3/dist at which corr drops to 0.05
                   sigma.sq.psi = 1, sigma.sq.p = 1, # random effect variances for occ and det formulas
                   z = apply(mody, c(1,2), max, na.rm = T)
)


# Run model
# Null model, intercept only
sfmod_null <- sfMsPGOcc(occ.formula = ~1, 
                        det.formula = ~1,
                        data = moddata, 
                        inits = inits.list,
                        priors = prior.list,
                        tuning = list(phi = 1),
                        n.factors = 3, # number of latent factors
                        verbose = T, n.report = 50,
                        n.batch = 240, batch.length = 25, # total n of samples is n.batch*batch.length
                        n.chains = 3, n.burn = 4000, n.thin = 2,n.omp.threads = 3,
                        k.fold = 4, k.fold.threads = 4
)
summary(sfmod_null)
waicOcc(sfmod_null) 
# Null model wAIC: 8052.4. The convergence for the spatial covariance terms is
# quite poor (rhats are 1.06, 3.05, 2.51). Maybe try longer chains (6000->30000)
# with more burn (4000->10000) thinning (2->50), or more factors (3->4)
sfmod_null_updated <- update(sfmod_null, 
                             n.factors = 4, n.batch = 1200, n.burn = 10000, n.thin = 50)
summary(sfmod_null_updated) 
sfmod_null_updated$k.fold.deviance
# The R-hat values for the occupancy and detection params are now better,
# although there are still issues with the spatial covariance estimates (Rhat = 1.3,1.0,1.1, 3.1).
sfmod_null$k.fold.deviance # there isn't a huge difference in the k-fold deviance, some higher, some lower
waicOcc(sfmod_null_updated) # 8026.6, much lower than the non-updated model
summary(ppcOcc(sfmod_null_updated, "freeman-tukey", group = 1)) # the model seems to have low Bayesian p-values for common spp (e.g. p = 0 for agoutis)

# Create alternative models with covars
sfmodfull <- update(sfmod_null, occ.formula = ~forest+scale(alt)+I(scale(alt)^2)+scale(ghm),
                    n.factors = 5, n.batch = 600, n.burn = 5000, n.thin = 25)
# the R hat is still very high for the spatial covariance latent factors, with
# very limited eff. sample size (ESS)
sfmodalt <- update(sfmodfull, occ.formula = ~scale(alt))
sfmodalt2 <- update(sfmodfull, occ.formula = ~scale(alt)+I(scale(alt)^2))
sfmodalt2ghm <- update(sfmodfull, occ.formula = ~scale(alt)+I(scale(alt)^2)+scale(ghm))
# summary: better Rhat values than the full model.
sfmodaltghm <- update(sfmodfull, occ.formula = ~scale(alt)+scale(ghm))
sfmodaltlc <- update(sfmodfull, occ.formula = ~forest+scale(alt))
sfmodalt2lc <- update(sfmodfull, occ.formula = ~forest+scale(alt)+I(scale(alt)^2))

# Compare models based on wAIC
waicOcc(sfmodfull) # waic = 8062 -> 7995.03
waicOcc(sfmodalt)  # wAIC = 8079
waicOcc(sfmodalt2) # wAIC = 8058 -> 7995.71
waicOcc(sfmodalt2ghm) # wAIC = 8054 -> 7997.57
waicOcc(sfmodaltghm) # wAIC = 8079
waicOcc(sfmodaltlc) # wAIC = 8091
waicOcc(sfmodalt2lc) # wAIC = 8064 -> 8002.06
# lowest wAIC is given by full model. The spatial covariance is not properly estimated though.


# Model summary
summary(sfmodalt2ghm)
# Phi estimates have rhat>>1.1 (2,4). Intercept and alt2 also have Rhat>1.1. 
plot(sfmodalt2ghm, 'beta', dens = F)

# Models without spatial covariates.
lfmodfull <- lfMsPGOcc(occ.formula = ~forest+scale(alt)+I(scale(alt)^2)+scale(ghm), 
                   det.formula = ~1,
                   data = moddata, 
                   inits = inits.list, 
                   n.factors = 4, # number of latent factors
                   n.samples = 15000, 
                   n.chains = 3, n.burn = 5000, n.thin = 25
) #working
waicOcc(update(lfmodfull, occ.formula = ~1)) # 8120
waicOcc(update(lfmodfull, occ.formula = ~forest+scale(alt))) # 8075
waicOcc(lfmodfull) # 8026 vs. 7995 for the spatial full model.
waicOcc(lfmodalt2) # 8072
waicOcc(lfmod) # 8080 for the original model.

# model with no species correlations
basemodalt2 <- msPGOcc(occ.formula = ~scale(alt)+I(scale(alt)^2)+scale(ghm), 
                   det.formula = ~1,
                   data = moddata, 
                   inits = inits.list, 
                   n.samples = 5000, 
                   n.chains = 3, n.burn = 3000, n.thin = 2)
# compare the two
waicOcc(basemod) # 8202 for the original model
waicOcc(basemodalt2) # 8220, among the worst models

# the model with the species correlations has lower wAIC, and therefore a better fit
# diagnostics
ppc.out <- ppcOcc(sfmodfull, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out)
# The model's Bayesian p value is 0.237, so not terrible, but rather low (0.5
# is good, <0.1 or >0.9 is poor). The fit is terrible for agoutis (p=0), and p<0.1 for collared
# peccaries, common opossums, rabbits, curassows, pumas, pacas, and coatis. 

### Plots ####

##### Coef plot ####
betas <- sfmodfull$beta.samples
# plot occupancy vs predictor values for every sp
data.frame(coef = colnames(betas),
           mean = colMeans(betas),
           lq = apply(betas,2,quantile, prob = 0.025),
           uq = apply(betas,2,quantile, prob = 0.975)
) %>% separate_wider_delim(coef, "-", names = c("parameter", "species")) %>% 
  # filter(parameter!="(Intercept)") %>% 
  mutate(parameter = gsub("scale(|)", "", parameter)) %>% 
  ggplot(aes(y = fct_reorder(species, mean, first), x = mean))+
  geom_vline(xintercept = 0, linetype=2)+
  geom_linerange(aes(xmin = lq, xmax = uq))+
  geom_point()+
  facet_wrap(~parameter)

##### Predicted occupancy ####
#We can use the model predictions to see how the expected richness compares
#against observed richness 
modpreds <- predict(sfmodalt2ghm, data.frame(int = 1, alt = scale(moddata$occ.covs$alt), alt2 = scale(moddata$occ.covs$alt)^2, ghm = scale(moddata$occ.covs$ghm)),moddata$coords, type = 'occupancy')
predrich <- apply(modpreds$z.0.samples,c(1,3),sum)|>apply(2,quantile,p = c(0.25,0.5,0.75))
rownames(predrich) <- c('lq','med','uq')
predrich <- t(predrich)

data.frame(site=rownames(comm),nobs=specnumber(comm)) %>% right_join(covs_df) %>% cbind(predrich) %>% 
  ggplot(aes(x = nobs))+
  geom_abline(slope = 1)+
  geom_pointrange(aes(y = med, ymin = lq,ymax = uq))
# generally model predicts greater than observed, with the greatest difference
# at the lower end. Is this evidence of species missing, or just poor model fit?

# Plot predicted occupancy vs alt.
cent_alt <- attr(scale(moddata$occ.covs$alt), "scaled:center")
scal_alt <- attr(scale(moddata$occ.covs$alt), "scaled:scale")

# newalts <- seq(min(moddata$occ.covs$alt), max(moddata$occ.covs$alt), length.out=100)
newalts <- c(100, 1200, 2000)
# newforest <- c(0, 1)
altpredvars <- scale(newalts, center = cent_alt,scale = scal_alt)
# forpredvars <- scale(newforest, center = cent_forest,scale = scal_forest)
# prediction covars. vary only the altitude, leave other covariates at 0
X.0 <- data.frame(intercept = 1, forest = 1, alt = altpredvars, alt2 = altpredvars^2, ghm = rep(c(-2,2),each = 3))

# get mean coordinates 
coords.0 <- matrix(colMeans(moddata$coords),ncol=2,nrow=nrow(X.0), byrow = T) 
modpreds <- predict(sfmodfull, X.0 = X.0, coords.0 = coords.0)
modpreds2occ <- modpreds2$psi.0.samples
dim(modpreds2occ)

# Plot of occupancy vs alt (L,M,H) for example spp
preddb <- data.frame(expand.grid(1:3000, sp = sfmodfull$sp.names,alt = altpredvars, ghm = c(-2,2)),psi=as.matrix(modpreds2occ)) 

# data.frame(expand.grid(sp = sfmodalt2$sp.names, alt=factor(newalts))) %>% 
# c("white_nosed_coati", "collared_peccary", "tinamou_great", "ocelot", "jaguar", "baird's_tapir", "spotted_paca")
filter(preddb, sp %in% c("tinamou_great", "jaguar", "ocelot",
                   "baird's_tapir","white_lipped_peccary","white_nosed_coati",
                   "jaguarundi","margay","opossum_four_eyed")) %>% # keep only a subset of species to represent the different patterns
  mutate(sp = factor(sp, 
                     levels = c("baird's_tapir","white_lipped_peccary","jaguar", 
                                "ocelot","margay","jaguarundi",
                                "tinamou_great", "white_nosed_coati","opossum_four_eyed"),
                     labels = c("Tapir", "White-lipped peccary", "Jaguar",  
                                "Ocelot", "Margay", "Jaguarundi", 
                                "Great tinamou","Coati", "Four-eyed opossum"))) %>% 
  ggplot(aes(factor(alt),psi, fill = factor(ghm)))+
  # geom_violin(scale = "width", fill = "seagreen")+
  # geom_violin(scale = 'width', show.legend = F)+
  stat_summary(aes(color = factor(ghm)), fun = "median", fun.min = \(x) quantile(x, p = 0.25), fun.max = \(x) quantile(x, p = 0.75), position = position_dodge(0.9), show.legend = F)+
  facet_wrap(~sp)+
  theme_classic(base_size = 12)+
  labs(x = "Altitude (masl)", y = "Predicted occupancy")+
  scale_x_discrete(labels = newalts)+
  scale_fill_manual(values = alpha(c("darkgreen", "goldenrod"), alpha = 0.5))+
  scale_color_manual(values = c("darkgreen", "darkgoldenrod"))+
  theme(strip.background = element_blank(), panel.border = element_rect(fill = NA))


#### Pitfall data ####
DATb <- read.csv("Data/dung_beetles_prc.csv") %>% 
  select(-Date) %>% 
  # join traps from same station
  summarise(.by = Plot, across(where(is.numeric), sum)) %>%
  inner_join(select(allsites_cov_db, site), by = join_by("Plot"=="site"))
# create sp matrix
beet_sp_mat <- DATb %>% 
  column_to_rownames("Plot") %>% 
  mutate(across(1:ncol(.), \(x) as.numeric(x>0))) %>% # make binary
  select(where(\(x) sum(x)>3)) %>% # only species present at multiple sites
  filter(rowSums(.)>0) %>% # filter empty sites
  as.matrix() %>% 
  t() # transpose

# create env matrix. This approach does not allow for factors like land cover.
beet_covs <- tibble(site = colnames(beet_sp_mat)) %>% left_join(allsites_cov_db) %>% 
  select(site, easting, northing, ghm1, lc2, alt1, lfdistance_2) %>% 
  mutate(lc2 = lc2==10) %>% 
  rename(x = easting, y = northing, alt = alt1, ghm = ghm1, forest = lc2, Fdist = lfdistance_2) %>% 
  column_to_rownames("site") %>% 
  as.matrix()

beet.data <- list(y = beet_sp_mat, 
                  covs = beet_covs[,c("forest", "alt", "ghm", "Fdist")], 
                  coords = beet_covs[,c("x", "y")])
### run JSDM model
beetjsdmnull <- lfJSDM(formula = ~1, 
                       data = beet.data, 
                       n.samples = 2000, n.thin = 2, n.chains = 3, n.factors = 4, n.burn = 1000)
beetjsdmfull <- lfJSDM(formula = ~forest+scale(alt)+I(scale(alt)^2)+scale(ghm),
                   data = beet.data, 
                   n.samples = 2000, n.thin = 2, n.chains = 3, n.factors = 4, n.burn = 1000)
sapply(list(beetjsdmfull, beetjsdmnull),waicOcc)
summary(beetjsdmfull,"community")
plot(beetjsdmfull, "beta.comm")
# The full model has lower wAIC than the null model. The model seems to perform
# better than the one for vertebrates. The overall community occupancy is
# positively linked with altitude and gHM, and negatively influenced by alt^2
# and distance to forest. The effect of forest is the least clear, likely due to
# collinearity. ESSs are low, run for longer
beetjsdmalt <- update(beetjsdmnull, formula = ~scale(alt))
beetjsdmalt2 <- update(beetjsdmnull, formula = ~scale(alt)+I(scale(alt)^2))
beetjsdmalt2ghm <- update(beetjsdmnull, formula = ~scale(alt)+I(scale(alt)^2)+scale(ghm))
beetjsdmalt2lc <- update(beetjsdmnull, formula = ~forest+scale(alt)+I(scale(alt)^2))
beetjsdmaltlc <- update(beetjsdmnull, formula = ~forest+scale(alt))
beetjsdmaltghm <- update(beetjsdmnull, formula = ~scale(alt)+scale(ghm))

sapply(list(beetjsdmalt, beetjsdmalt2, beetjsdmfull, beetjsdmalt2lc, 
            beetjsdmaltlc, beetjsdmaltghm, beetjsdmalt2ghm),waicOcc)
# The full has the lowest wAIC, but it is comparable to the mode model with alt, alt^2, and ghm, and the one with alt and lc. 
summary(beetjsdmalt2ghm, "community")
summary(beetjsdmaltlc, "community")

### Models with spatial covars
sfbeetJSDMnull <- sfJSDM(~1, 
                         data = beet.data, 
                         inits = list(phi = 3/10000),
                         priors = prior.list,
                         n.factors = 2, # number of latent factors
                         n.batch = 200, batch.length = 25, # total n of samples is n.batch*batch.length
                         n.chains = 3, n.burn = 3000, n.thin = 4)

sfbeetJSDMfull <- update(sfbeetJSDMnull, formula = ~forest+scale(alt)+I(scale(alt)^2)+scale(ghm))
sapply(list(sfbeetJSDMfull, sfbeetJSDMnull), waicOcc)
# Here too the full model performs better

sfbeetJSDMalt2 <- update(sfbeetJSDMfull, formula = ~scale(alt)+I(scale(alt)^2))
sfbeetJSDMalt <- update(sfbeetJSDMfull, formula = ~scale(alt))
sfbeetJSDMaltghm <- update(sfbeetJSDMfull, formula = ~scale(alt)+scale(ghm))
sfbeetJSDMalt2ghm <- update(sfbeetJSDMfull, formula = ~scale(alt)+I(scale(alt)^2)+scale(ghm))
sfbeetJSDMaltlc <- update(sfbeetJSDMfull, formula = ~forest+scale(alt))
sfbeetJSDMalt2lc <- update(sfbeetJSDMfull, formula = ~forest+scale(alt)+I(scale(alt)^2))
sfbeetJSDMaltlcghm <- update(sfbeetJSDMfull, formula = ~forest+scale(alt)+scale(ghm))

# AIC ranking
sapply(list(sfbeetJSDMfull,   sfbeetJSDMaltghm, sfbeetJSDMaltlcghm, sfbeetJSDMalt2ghm,
            sfbeetJSDMalt2lc), waicOcc) 
# the model with the lowest AIC is the full model (2402). All
# these models, however, have higher wAIC than the model with no spatial
# covariates (wAIC = 2376)
rm(sfbeetJSDMalt, sfbeetJSDMalt2, sfbeetJSDMaltlc)

##### occupancy model with bait as detection covar #####
beet_sp_mat <- read.csv("Data/dung_beetles_prc.csv") %>% 
  filter(Bait !="") %>% # remove stations with no bait info
  select(-c(Date,Plot)) %>% 
  inner_join(select(allsites_cov_db, site)) %>% 
  mutate(across(where(is.numeric), \(x) as.numeric(x>0))) %>% # make binary
  pivot_longer(where(is.numeric), names_to = "sp", values_to = "pres") %>% 
  filter(sum(pres)>0, .by = site) %>% # filter empty sites
  # filter(sum(pres)>5, .by = sp) %>% # only species present at multiple (>5) sites
  # complete(nesting(site,Bait), sp, fill = list(pres = 0)) %>% 
  complete(site,Bait,sp) %>% 
  arrange(Bait,site,sp)
nsp <- length(unique(beet_sp_mat$sp))
nst <- length(unique(beet_sp_mat$site))
nocc <- length(unique(beet_sp_mat$Bait))
nsp*nst*nocc
beetmody <- array(data = beet_sp_mat$pres, dim = c(nsp, nst, nocc),
              dimnames = list(species = unique(beet_sp_mat$sp), 
                              sites = unique(beet_sp_mat$site),
                              bait = unique(beet_sp_mat$Bait)))
beetmody <- beetmody[,apply(beetmody,2,sum,na.rm=T)>5,]
# create env matrix
beet_covs <- tibble(site = dimnames(beetmody)[[2]]) %>% left_join(allsites_cov_db) %>% 
  select(site, easting, northing, ghm1, lc2,alt1,  lfdistance_2) %>% 
  rename(x = easting, y = northing, forest = lc2, alt = alt1, ghm = ghm1, pa_dist = distance, lgfor_dist = lfdistance) %>% 
  mutate(forest = factor(forest, labels = c("Dense","Open"))) %>% 
  column_to_rownames("site") %>% 
  as.matrix()
baitcov <- read.csv("Data/dung_beetles_prc.csv") %>% 
  filter(Bait !="") %>% # remove stations with no bait info
  select(-c(Date,Plot)) %>% 
  inner_join(select(campts_db, site)) %>% 
  mutate(across(where(is.numeric), \(x) as.numeric(x>0))) %>% # make binary
  pivot_longer(where(is.numeric), names_to = "sp", values_to = "pres") %>% 
  filter(sum(pres)>0, .by = site) %>% # filter empty sites
  filter(sum(pres)>5, .by = sp) %>% # only species present at multiple (>5) sites
  distinct(site,Bait) %>% 
  pivot_wider(values_from = Bait, names_from = Bait) %>% 
  arrange(site)

beet.data.occ <- list(y = beetmody, 
                  occ.covs = beet_covs[,c("landcov", "alt", "ghm", "pa_dist", "lgfor_dist")], 
                  det.covs = list(bait = baitcov),
                  coords = beet_covs[,c("x", "y")])

lfmodbeetnull <- lfMsPGOcc(occ.formula = ~1,
                       det.formula = ~1,
                       data = beet.data.occ, 
                       n.samples = 1000, n.thin = 2, n.chains = 3, n.factors = 4)

lfmodbeetfull <- update(lfmodbeetnull, occ.formula = ~scale(alt)+scale(landcov)+scale(ghm))


#### Plots ####
##### Parameter plots ####
beetbetas <- beetjsdmfull$beta.samples
data.frame(coef = colnames(beetbetas),
           mean = colMeans(beetbetas),
           lq = apply(beetbetas,2,quantile, prob = 0.025),
           uq = apply(beetbetas,2,quantile, prob = 0.975)
) %>% separate_wider_delim(coef, "-", names = c("parameter", "species")) %>% 
  # filter(parameter!="(Intercept)") %>% 
  mutate(parameter = gsub("scale(|)", "", parameter)) %>% 
  ggplot(aes(y = fct_reorder(species,mean, first), x = mean))+
  geom_vline(xintercept = 0, linetype=2)+
  geom_linerange(aes(xmin = lq, xmax = uq))+
  geom_point()+
  facet_wrap(~parameter)

#####--- Predicted occupancy ####
cent_alt <- attr(scale(beet.data$covs[,'alt']), "scaled:center")
scal_alt <- attr(scale(beet.data$covs[,'alt']), "scaled:scale")

newalts <- c(100, 500, 1500)
altpredvars <- as.numeric(scale(newalts, center = cent_alt,scale = scal_alt))
newghm <- c(-1,1)

# prediction covars. 
X.0 <- data.frame(intercept = 1, forest = 1, alt = altpredvars, alt2 = altpredvars^2, ghm = rep(newghm,each = length(newalts)))
X.0 <- as.matrix(X.0)
# get mean coordinates 
coords.0 <- matrix(colMeans(beet.data$coords), ncol=2, nrow=nrow(X.0), byrow = T) 
modpreds <- predict(beetjsdmfull, X.0 = X.0, coords.0 = coords.0) # 
modpredsocc <- modpreds$psi.0.samples
dim(modpredsocc)
# dim(out$psi.0.samples)

# Plot of occupancy vs alt (L,M,H) for example spp
preddb <- data.frame(expand.grid(1:1500, sp = beetjsdmfull$sp.names,
                                 alt = altpredvars, ghm = c(-1,1)),
                     psi=as.matrix(modpredsocc))

# data.frame(expand.grid(sp = sfmodalt2$sp.names, alt=factor(newalts))) %>% 
filter(preddb, sp %in% c("Onthophagus_acuminatus", "Coprophanaeus_telamon", "Coprophanaeus_pecki",
                         "Canthon_angustatus","Onthophagus_praecellens","Onthophagus_haematopus")) %>% # keep only a subset of species to represent the different patterns
  ggplot(aes(factor(alt),psi, fill = factor(ghm)))+
  geom_violin(scale = 'width', show.legend = F)+
  stat_summary(aes(color = factor(ghm)), fun = "median", fun.min = \(x) quantile(x, p = 0.25), fun.max = \(x) quantile(x, p = 0.75), position = position_dodge(0.9), show.legend = F)+
  facet_wrap(~sp)+
  theme_classic(base_size = 12)+
  labs(x = "Altitude (masl)", y = "Predicted occupancy")+
  scale_x_discrete(labels = newalts)+
  scale_fill_manual(values = alpha(c("darkgreen", "goldenrod"), alpha = 0.5))+
  scale_color_manual(values = c("darkgreen", "darkgoldenrod"))+
  theme(strip.background = element_blank(), panel.border = element_rect(fill = NA))
