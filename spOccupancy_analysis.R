library(spOccupancy)
library(tidyverse)

# We fit a model similar to the imperfect detection one implemented in Tobler et
# el. 2019, as implemented by the sfMsPGOcc function in the spOccupancy package
# This is an occupancy model with a detection component, so we need repeat
# occasions at each site. We exclude the sites with less than 4 weeks of
# sampling, and divide the detections into monthly counts.

#### CT data ####

camvect <- terra::vect("Data/ct_covars.geojson")
campts_db <- as.data.frame(camvect) %>% mutate(site = gsub("#", "", deployment_id), days = as.numeric(days))

data_raw <- read_csv("Data/ct_records.csv",
                     col_types = cols_only(deployment_id = 'c', common_name = 'c', timestamp = col_datetime("%Y-%m-%d %H:%M") ))
# include only cameras with more than 12 weeks (~3 months)
cams_include <- summarise(data_raw, .by = deployment_id, t=max(timestamp)-min(timestamp)) %>% 
  filter(t>12*7) %>% 
  mutate(site = sub("#","",deployment_id)) %>% 
  select(site) 
DAT <- read_csv("Data/records_wild.csv") %>% inner_join(cams_include)

# model data
y_df <- filter(DAT, !grepl("uid", common_name)) %>% 
  left_join(select(campts_db,site)) %>% 
  mutate(smplwk = as.numeric(1+floor(difftime(timestamp,min(timestamp),units = "weeks"))), .by = site) %>% 
  mutate(smplmon = 1+smplwk%/%4) %>% 
  filter(max(smplwk)>4, .by = site) %>% 
  count(smplmon, site, common_name) %>% 
  complete(nesting(smplmon, site), common_name, fill = list(n=0)) %>% 
  complete(smplmon, site, common_name) %>% 
  arrange(smplmon,  site, common_name)
nsp <- length(unique(y_df$common_name))
nst <- length(unique(y_df$site))
nocc <- max(y_df$smplmon)
nsp*nst*nocc
mody <- array(data = as.numeric(y_df$n>0), dim = c(nsp, nst, nocc), 
              dimnames = list(species = unique(y_df$common_name), sites = unique(y_df$site), occ = 1:nocc))
covs_df <- distinct(y_df,site) %>% left_join(campts_db) %>% 
  select(site, easting, northing, ghm1, lc1,alt1, evi_mean1, lfdistance, distance) %>% 
  rename(x = easting, y = northing, forest = lc1, alt = alt1, ghm = ghm1, pa_dist = distance, lgfor_dist = lfdistance) %>% 
  mutate(forest = ifelse(forest==12,1,0)) 
moddata <- list(y = mody, 
                occ.covs = covs_df[,c("forest","alt", "ghm", "pa_dist", "lgfor_dist")] ,
                # det.covs = covs_df[,"site"],
                coords = as.matrix(covs_df[,c("x", "y")])
)



# Priors
prior.list <- list(phi.unif = list(3/100))

# initial values
inits.list <- list(beta.comm = 0, alpha.comm = 0, # Community-level occurrence (beta.comm) and detection (alpha.comm) regression coefficients
                   beta = 0, alpha = 0, # species-level occ and det (could be length nsp)
                   tau.sq.beta = 1, tau.sq.alpha = 1, # comm=level variance for occ and det
                   phi = 3/100, # spatial decay parameter. 3/phi is dist at which corr drops to 0.05
                   sigma.sq.psi = 1, sigma.sq.p = 1, # random effect variances for occ and det formulas
                   z = apply(mody, c(1,2), max, na.rm = T)
)


# Run model
sfmod <- sfMsPGOcc(occ.formula = ~scale(forest)+scale(alt)+scale(ghm), 
                   det.formula = ~1,
                   data = moddata, 
                   # inits = inits.list, 
                   # priors = prior.list, 
                   n.factors = 2, # number of latent factors
                   n.batch = 200, batch.length = 25, # total n of samples is n.batch*batch.length
                   n.chains = 3, n.burn = 3000, n.thin = 2
)


lfmod <- lfMsPGOcc(occ.formula = ~scale(forest)+scale(alt)+scale(ghm), 
                   det.formula = ~1,
                   data = moddata, 
                   inits = inits.list, 
                   # priors = prior.list, 
                   n.factors = 2, # number of latent factors
                   n.samples = 5000, 
                   n.chains = 3, n.burn = 3000, n.thin = 2
) #working

# model with no species correlations
basemod <- msPGOcc(occ.formula = ~scale(forest)+scale(alt)+scale(ghm), 
                   det.formula = ~1,
                   data = moddata, 
                   inits = inits.list, 
                   # priors = prior.list, 
                   n.samples = 5000, 
                   n.chains = 3, n.burn = 3000, n.thin = 2)
# compare the two
waicOcc(lfmod)
waicOcc(basemod)
# the model with the species correlations has lower wAIC, and therefore a better fit
# diagnostics
ppc.out <- ppcOcc(lfmod, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out)
# The model's Bayesian p value is 0.244, so not terrible, but rather low (0.5
# is good, <0.1 or >0.9 is poor). The fit is quite poor for agoutis, collared
# peccaries, rabbits, curassows, pumas, pacas, and tayras. We can try to fit a
# different model with other occ or det covariates.
lfmod2 <- lfMsPGOcc(occ.formula = ~scale(forest)+scale(alt)+scale(pa_dist), 
                    det.formula = ~1,
                    data = moddata, 
                    inits = inits.list, 
                    # priors = prior.list, 
                    n.factors = 2, # number of latent factors
                    n.samples = 5000, 
                    n.chains = 3, n.burn = 3000, n.thin = 2)

summary(lfmod2)
lfmod3 <- lfMsPGOcc(occ.formula = ~scale(forest)+scale(alt)+scale(ghm)+scale(pa_dist==0), 
                    det.formula = ~1,
                    data = moddata, 
                    inits = inits.list, 
                    # priors = prior.list, 
                    n.factors = 2, # number of latent factors
                    n.samples = 5000, 
                    n.chains = 3, n.burn = 3000, n.thin = 2)
waicOcc(lfmod2)
waicOcc(lfmod3)
waicOcc(lfmod)
# the first model with ghm is still better, although not substantially better
# than the model with four covariates
summary(lfmod3)
plot(lfmod3, 'beta', density = F) # plots look good
ppc.out3 <- ppcOcc(lfmod3, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out3)

#### Plots ####

# extract coefficients
betas <- sfmod$beta.samples
# plot occupancy vs predictor values for every sp
data.frame(coef = colnames(betas),
           mean = colMeans(betas),
           lq = apply(betas,2,quantile, prob = 0.025),
           uq = apply(betas,2,quantile, prob = 0.975)
) %>% separate_wider_delim(coef, "-", names = c("parameter", "species")) %>% 
  filter(parameter!="(Intercept)") %>% 
  mutate(parameter = gsub("scale(|)", "", parameter)) %>% 
  ggplot(aes(y = species, x = mean))+
  geom_vline(xintercept = 0, linetype=2)+
  geom_linerange(aes(xmin = lq, xmax = uq))+
  geom_point()+
  facet_wrap(~parameter)

# Plot predicted occupancy vs alt.
cent_alt <- attr(scale(moddata$occ.covs$alt), "scaled:center")
scal_alt <- attr(scale(moddata$occ.covs$alt), "scaled:scale")
cent_forest <- attr(scale(moddata$occ.covs$forest), "scaled:center")
scal_forest <- attr(scale(moddata$occ.covs$forest), "scaled:scale")

# newalts <- seq(min(moddata$occ.covs$alt), max(moddata$occ.covs$alt), length.out=100)
newalts <- c(100, 1200, 2000)
newforest <- c(0, 1)
altpredvars <- scale(newalts, center = cent_alt,scale = scal_alt)
forpredvars <- scale(newforest, center = cent_forest,scale = scal_forest)
# prediction covars. vary only the altitude, leave other covariates at 0
X.0 <- expand.grid(intercept = 1, forest = forpredvars, alt = altpredvars, ghm = 0)
# get mean coordinates 
coords.0 <- matrix(apply(sfmod$coords, 2, mean),ncol=2,nrow=nrow(X.0), byrow = T) 
modpreds <- predict(sfmod, X.0 = X.0, coords.0 = coords.0)

dim(modpreds$psi.0.samples)
predmean <- as.numeric(apply(modpreds$psi.0.samples, c(2,3), median))
predql <- as.numeric(apply(modpreds$psi.0.samples, c(2,3), quantile, probs = 0.25))
predqu <- as.numeric(apply(modpreds$psi.0.samples, c(2,3), quantile, probs = 0.75))

# Plot of occupancy vs alt (L,M,H) and forest cover (1,0) for example spp
data.frame(expand.grid(sp = sfmod$sp.names,forest = factor(newforest), alt=factor(newalts)),psi=predmean,psiu=predqu,psil=predql) %>% 
# data.frame(expand.grid(sp = sfmod$sp.names, alt = newalts, scalealt = altpredvars,), 
#            psi = predmean, psiu = predqu, psil = predql) %>% 
  filter(sp %in% c("spotted_paca", "central_american_agouti", "cottontail_dices", 
                   "tinamou_great", "jaguar", "common_opossum", "ocelot",
                   "white_nosed_coati", "opossum_four_eyed")) %>% # keep only a subset of species to represent the different patterns
  ggplot(aes(alt, psi,color = forest,shape=forest))+
  geom_pointrange(aes(ymin = psil,ymax=psiu), position = position_dodge(width = 0.5), show.legend = F)+
  facet_wrap(~sp)+
  theme_classic(base_size = 16)+
  scale_color_manual(values = c("goldenrod","darkgreen"))+
  scale_shape_manual(values = c(17, 19))+
  labs(x = "Altitude (masl)", y = "Predicted occupancy")+
  theme(strip.background = element_blank(), panel.border = element_rect(fill = NA))


#### Pitfall data ####
DATb <- read.csv("Data/dung_beetles_prc.csv") %>% 
  select(-Date) %>% 
  # join traps from same station
  summarise(.by = Plot, across(where(is.numeric), sum)) %>%
  inner_join(select(campts_db, site), by = join_by("Plot"=="site"))
# create sp matrix
beet_sp_mat <- DATb %>% 
  column_to_rownames("Plot") %>% 
  mutate(across(1:ncol(.), \(x) as.numeric(x>0))) %>% # make binary
  select(where(\(x) sum(x)>5)) %>% # only species present at multiple (>5) sites
  filter(rowSums(.)>0) %>% # filter empty sites
  as.matrix() %>% 
  t() # transpose
# create env matrix
beet_covs <- tibble(site = rownames(beet_sp_mat)) %>% left_join(campts_db) %>% 
  select(site, easting, northing, ghm1, lc1,alt1, evi_mean1, lfdistance, distance) %>% 
  rename(x = easting, y = northing, forest = lc1, alt = alt1, ghm = ghm1, pa_dist = distance, lgfor_dist = lfdistance) %>% 
  mutate(landcov = ifelse(forest==12,1,0)) %>% 
  column_to_rownames("site") %>% 
  as.matrix()

beet.data <- list(y = beet_sp_mat, 
                  covs = beet_covs[,c("forest", "alt", "ghm", "pa_dist", "lgfor_dist")], 
                  coords = beet_covs[,c("x", "y")])
# run model
beetjsdm <- lfJSDM(formula = ~scale(forest)+scale(alt)+scale(ghm),
                   data = beet.data, 
                   n.samples = 1000, n.thin = 2, n.chains = 3, n.factors = 4)


# occupancy model with bait as detection covar
beet_sp_mat <- read.csv("Data/dung_beetles_prc.csv") %>% 
  filter(Bait !="") %>% # remove stations with no bait info
  select(-c(Date,Plot)) %>% 
  inner_join(select(campts_db, site)) %>% 
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
beet_covs <- tibble(site = dimnames(beetmody)[[2]]) %>% left_join(campts_db) %>% 
  select(site, easting, northing, ghm1, lc1,alt1, evi_mean1, lfdistance, distance) %>% 
  rename(x = easting, y = northing, forest = lc1, alt = alt1, ghm = ghm1, pa_dist = distance, lgfor_dist = lfdistance) %>% 
  mutate(landcov = ifelse(forest==12,1,0)) %>% 
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

lfmodbeet <- lfMsPGOcc(occ.formula = ~scale(alt)+scale(landcov)+scale(ghm),
                       det.formula = ~1,
                   data = beet.data.occ, 
                   n.samples = 1000, n.thin = 2, n.chains = 3, n.factors = 4)

# Plot predicted occupancy vs alt.
cent_alt <- attr(scale(beet.data.occ$occ.covs[,'alt']), "scaled:center")
scal_alt <- attr(scale(beet.data.occ$occ.covs[,'alt']), "scaled:scale")
cent_forest <- attr(scale(beet.data.occ$occ.covs[,'landcov']), "scaled:center")
scal_forest <- attr(scale(beet.data.occ$occ.covs[,'landcov']), "scaled:scale")

altpredvars <- scale(c(100,1000, 1500), center = cent_alt,scale = scal_alt)
forpredvars <- scale(newforest, center = cent_forest,scale = scal_forest)
# prediction covars. vary only the altitude, leave other covariates at 0
X.0 <- expand.grid(intercept = 1, alt = altpredvars, landcov = forpredvars, ghm = 0)
# get mean coordinates 
coords.0 <- matrix(apply(lfmodbeet$coords, 2, mean),ncol=2,nrow=nrow(X.0), byrow = T) 
modpreds <- predict(lfmodbeet, X.0 = X.0, coords.0 = coords.0)

dim(modpreds$psi.0.samples)
predmean <- as.numeric(apply(modpreds$psi.0.samples, c(2,3), median))
predql <- as.numeric(apply(modpreds$psi.0.samples, c(2,3), quantile, probs = 0.25))
predqu <- as.numeric(apply(modpreds$psi.0.samples, c(2,3), quantile, probs = 0.75))

# Plot of occupancy vs alt (L,M,H) and forest cover (1,0) for example spp
data.frame(expand.grid(sp = sfmod$sp.names,forest = factor(newforest), alt=factor(newalts)),psi=predmean,psiu=predqu,psil=predql) %>% 
  # data.frame(expand.grid(sp = sfmod$sp.names, alt = newalts, scalealt = altpredvars,), 
  #            psi = predmean, psiu = predqu, psil = predql) %>% 
  filter(sp %in% c("spotted_paca", "central_american_agouti", "cottontail_dices", 
                   "tinamou_great", "jaguar", "common_opossum", "ocelot",
                   "white_nosed_coati", "opossum_four_eyed")) %>% # keep only a subset of species to represent the different patterns
  ggplot(aes(alt, psi,color = forest,shape=forest))+
  geom_pointrange(aes(ymin = psil,ymax=psiu), position = position_dodge(width = 0.5), show.legend = F)+
  facet_wrap(~sp)+
  theme_classic(base_size = 16)+
  scale_color_manual(values = c("goldenrod","darkgreen"))+
  scale_shape_manual(values = c(17, 19))+
  labs(x = "Altitude (masl)", y = "Predicted occupancy")+
  theme(strip.background = element_blank(), panel.border = element_rect(fill = NA))

#### Plots ####
beetbetas <- beetjsdm$beta.samples
data.frame(coef = colnames(beetbetas),
           mean = colMeans(beetbetas),
           lq = apply(beetbetas,2,quantile, prob = 0.025),
           uq = apply(beetbetas,2,quantile, prob = 0.975)
) %>% separate_wider_delim(coef, "-", names = c("parameter", "species")) %>% 
  filter(parameter!="(Intercept)") %>% 
  mutate(parameter = gsub("scale(|)", "", parameter)) %>% 
  ggplot(aes(y = species, x = mean))+
  geom_vline(xintercept = 0, linetype=2)+
  geom_linerange(aes(xmin = lq, xmax = uq))+
  geom_point()+
  facet_wrap(~parameter)

beetbetasocc <- lfmodbeet$beta.samples
data.frame(coef = colnames(beetbetasocc),
           mean = colMeans(beetbetasocc),
           lq = apply(beetbetasocc,2,quantile, prob = 0.025),
           uq = apply(beetbetasocc,2,quantile, prob = 0.975)
) %>% separate_wider_delim(coef, "-", names = c("parameter", "species")) %>% 
  filter(parameter!="(Intercept)") %>% 
  mutate(parameter = gsub("scale(|)", "", parameter)) %>% 
  ggplot(aes(y = species, x = mean))+
  geom_vline(xintercept = 0, linetype=2)+
  geom_linerange(aes(xmin = lq, xmax = uq))+
  geom_point()+
  facet_wrap(~parameter)
