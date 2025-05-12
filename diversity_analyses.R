# Community diversity analyses. 

# Author: Juan S. Vargas
# Date: 25/4/25

# Data were obtained from camera-traps and pitfall traps placed along an
# altitudinal gradient in southern Costa Rica

# Load libraries
library(vegan)
library(tidyverse)
library(terra)
library(mgcv)
library(ape)
library(AICcPermanova)
library(gratia)

#### Import data ####
# Import camera-trap data as matrix for vegan
comm <- read.csv("Data/vert_comm_mat.csv", row.names = 1)

# Import pitfall trap data
commb <- read.csv("Data/beet_comm_mat.csv", row.names = 1)

# Import station data #
allsites_pts <- vect("Data/combined_survey_pts.geojson")
allsites_cov_db <- as.data.frame(allsites_pts)

#### Community composition ####
# Calculate distances
vert_dissim <- vegdist(comm)
beet_dissim <- vegdist(commb)
vert_covs_raw <- tibble(site = rownames(comm)) %>% left_join(allsites_cov_db) %>% 
  select(site, alt1, lc2, temp_mean1, ghm1, lfdistance_2) 
vert_covs <- mutate(vert_covs_raw, lc = factor(lc2, labels = c("Dense","Open")),
         alt = scale(alt1),
         temp = scale(temp_mean1),
         ghm = scale(ghm1),
         fdist = scale(log1p(lfdistance_2))) %>% 
  select(lc,alt,temp,ghm,fdist)
  
beet_covs_raw <- tibble(site = rownames(commb)) %>% left_join(allsites_cov_db) %>% 
  select(site, alt1, lc2, temp_mean1, ghm1, lfdistance_2) 
beet_covs <- mutate(beet_covs_raw, lc = factor(lc2, labels = c("Dense","Open")),
         alt = scale(alt1),
         temp = scale(temp_mean1),
         ghm = scale(ghm1),
         fdist = scale(log1p(lfdistance_2))) %>% 
  select(lc,alt,temp,ghm,fdist)

##### PERMANOVA #####

# Analysis to see how the different variables affect the dissimilarity across
# sites
permanova_comp <- fit_models(make_models(vars = c("lc", "alt", "ghm", "fdist", "temp")), 
                             veg_data = comm, env_data = vert_covs)
select_models(permanova_comp)
# We compared different covariate combinations, and the one with the lowest AICc
# included altitude, distance to forest, and temperature. However, there are
# seven other models that have similar AICc. Among the candidates, the third
# best one also has the lowest VIF. This one includes altitude and distance to
# forest. The VIF for the best model is not that high (3.2) so we could ignore
# it.
permanova_comp_beet <- fit_models(make_models(vars = c("lc", "alt", "ghm", "fdist", "temp")), 
                             veg_data = commb, env_data = beet_covs)
select_models(permanova_comp_beet)
# In the case of beetles, all 15 top models rank similarly, although the full one is the worst.
vert_permanova <- adonis2(comm~alt+fdist+temp, data = vert_covs, method = 'bray')
beet_permanova <- adonis2(commb~alt+ghm+fdist, data = beet_covs, method = 'bray')
vert_permanova
beet_permanova

# These results suggest that differences in species composition for vertebrates
# are influenced mostly by differences in altitude and distance to forest.
# Temperature and ghm were included in the models for verts and beetles,
# respectively. Results are similar to those for beetles, except
# for temperature, which did not have a strong effect on beetle community
# composition. Altitude and land cover have the strongest effects for both

##### PERMANOVA viz #####
# Run PCoA for visualization of PERMANOVA
vert_pcoa <- pcoa(vert_dissim)
beet_pcoa <- pcoa(beet_dissim)

# Fit environ variables
vert_envfit_pcoa <- envfit(vert_pcoa$vectors, vert_covs)
beet_envfit_pcoa <- envfit(beet_pcoa$vectors, beet_covs)

# test if there is more variability in open vs dense forest sites for vertebrates
anova(betadisper(vert_dissim, vert_covs$lc))
plot(betadisper(vert_dissim, vert_covs$lc), hull = F, ellipse = T, col = c("darkgreen", "goldenrod"), segments = F)
# There seems to be greater distance to the centroid in open forest sites
# (betadisper: F_1,115 = 5.13, p = 0.025). It is important to keep this in mind
# as it may be driving apparent differences due to other covariates

anova(betadisper(beet_dissim, beet_covs$lc))
plot(betadisper(beet_dissim, beet_covs$lc), 
     hull = F, ellipse=T, col = c("darkgreen", "goldenrod"), segments = F)
# For beetles there is no difference in variability. F_1,70 = 0.20, p = 0.65.

# Plot PCoA
plot(vert_pcoa$vectors[,1:2], las=1, asp = 1, xlab = "PCoA 1", ylab = "PCoA 2", type = 'n')
ordiellipse(vert_pcoa$vectors, vert_covs$lc, draw = "polygon", col = c("darkgreen", "goldenrod"))
points(vert_pcoa$vectors[vert_covs$lc=="Dense",1:2], pch = 16, col = "darkgreen")
points(vert_pcoa$vectors[vert_covs$lc=="Open",1:2], pch = 17, col = "goldenrod")
ordisurf(vert_pcoa$vectors, vert_covs_raw$lfdistance_2, add = T, col = "gray50")
plot(vert_envfit_pcoa, labels = list(factors = c("Dense", "Open")), bg='white', col = "gray10")

# PCoA for beetles
plot(beet_pcoa$vectors[,1:2], type = 'n', las=1, asp = 1, xlab = "PCoA 1", ylab = "PCoA 2")
ordiellipse(beet_pcoa$vectors, beet_covs$lc, draw = "polygon", col = c("darkgreen", "goldenrod"))
points(beet_pcoa$vectors[beet_covs$lc=="Dense",1:2], pch = 16, col = "darkgreen")
points(beet_pcoa$vectors[beet_covs$lc=="Open",1:2], pch = 17, col = "goldenrod")
ordisurf(beet_pcoa$vectors, beet_covs_raw$alt1, add = T, col = "gray50")
plot(beet_envfit_pcoa, labels = list(factors = c("Dense", "Open")), bg='white', col = "gray10")

#### Community metrics ####
par(mfrow = c(2,2))
# plot richness (observed) against covariates individually
plot(specnumber(comm)~vert_covs$alt)
plot(specnumber(comm)~vert_covs$temp)
plot(specnumber(comm)~vert_covs$ghm)
plot(specnumber(comm)~vert_covs$fdist)

plot(specnumber(commb)~beet_covs$alt)
plot(specnumber(commb)~beet_covs$temp)
plot(specnumber(commb)~beet_covs$ghm)
plot(specnumber(commb)~beet_covs$fdist)

# plot diversity (Shannon) against covariates individually
plot(diversity(comm)~vert_covs$alt)
plot(diversity(comm)~vert_covs$temp)
plot(diversity(comm)~vert_covs$ghm)
plot(diversity(comm)~vert_covs$fdist)
# beetles
plot(diversity(commb)~beet_covs$alt)
plot(diversity(commb)~beet_covs$temp)
plot(diversity(commb)~beet_covs$ghm)
plot(diversity(commb)~beet_covs$fdist)

##### Vert richness #####
vert_gam_db <- data.frame(s = diversity(comm), n = specnumber(comm), vert_covs)
beet_gam_db <- data.frame(s = diversity(commb), n = specnumber(commb), beet_covs)
plot(vert_gam_db)
plot(beet_gam_db)

# Sites with higer richness also have higher diversity for both verts and
# beetles
# Vertebrate richness linear model
vert_rich_lm <- gam(n~lc+ghm+temp+alt+fdist, family = poisson, data=vert_gam_db)
# Full GAM, with smooths for all covars
vert_rich_gam <- gam(n~lc+s(ghm)+s(temp)+s(alt)+s(fdist), 
                     family = poisson, data = vert_gam_db, method = 'REML')
summary(vert_rich_gam)
# The full model reduces the ghm to a linear effect.
concurvity(vert_rich_gam, full = F)
# There is high concurvity between temp and alt. Similarly, ghm is confounded by temp, alt and fdist. 
# Let's remove temperature, and make ghm explicitly linear.
vert_rich_gam2 <- update(vert_rich_gam, .~. - s(ghm)+ghm-s(temp))
summary(vert_rich_gam2)
concurvity(vert_rich_gam2, full = F)
# This model has significant smooths for altitude and fdist. The linear terms
# are not significant and can be removed. There is also still some confusion
# between altitude and fdist (~0.5). I will leave only altitude
vert_rich_gam3 <- update(vert_rich_gam2, .~.-lc-ghm-s(fdist))
summary(vert_rich_gam3)
anova(vert_rich_gam, vert_rich_gam2, vert_rich_gam3, vert_rich_lm, test = "Chisq")
# All models are different. The one with the lowest residual deviance is the
# full one. 
AIC(vert_rich_gam, vert_rich_gam2, vert_rich_gam3, vert_rich_lm)
# The lowest AIC is given by the full model, followed by model 2, which deals
# with concurvity.
vert_rich_gam4 <- update(vert_rich_gam2, .~.-lc-ghm)
draw(vert_rich_gam4, residuals = T)
appraise(vert_rich_gam4)
anova(vert_rich_gam, vert_rich_gam2, vert_rich_gam3, vert_rich_gam4, vert_rich_lm, test='Chisq')
AIC(vert_rich_gam, vert_rich_gam2, vert_rich_gam3, vert_rich_gam4, vert_rich_lm)

# I will report the results of model 4, which has smooth terms for distance to
# forest and altitude, and no linear terms. This model has comparable AIC to the
# full model, but does not have problems of concurvity.
summary(vert_rich_gam4)


##### Beet richness #####

# linear model
beet_rich_lm <- gam(n~lc+ghm+temp+alt+fdist, data = beet_gam_db)
summary(beet_rich_lm) # no sig effects
# beetle richness GAM
beet_rich_gam <- gam(n~lc+s(ghm)+s(temp)+s(alt)+s(fdist), 
                     family = poisson, data = beet_gam_db)
summary(beet_rich_gam)
concurvity(beet_rich_gam, full = F)
# no significant smooths, nor linear effects. Furthermore, there is noticeable
# concurvity, highest in the full model is between ghm and altitude. I'll remove
# ghm first.
beet_rich_gam2 <- update(beet_rich_gam, .~.-s(ghm))
summary(beet_rich_gam2)
concurvity(beet_rich_gam2, full = F)
# This model has a significant smooth for altitude, and temp is reduced to a
# linear effect. There is still concurvity of temp and fdist, these depend on
# altitude. Let's first remove the smooth of temp and make it explicitly linear.
beet_rich_gam3 <- update(beet_rich_gam2, .~.-s(temp)+temp)
summary(beet_rich_gam3)
concurvity(beet_rich_gam3)
# There is still  concurvity issues, I'll remove both separately
beet_rich_gam4 <- update(beet_rich_gam3, .~.-s(fdist))
beet_rich_gam5 <- update(beet_rich_gam3, .~.-s(alt))
summary(beet_rich_gam4)
summary(beet_rich_gam5)
# The model with only distance to forest is much poorer.
# compare vs linear model
anova(beet_rich_lm,beet_rich_gam, beet_rich_gam2, beet_rich_gam3, beet_rich_gam4, test = "Chisq") # the nonlinear model fits the data much better
anova(beet_rich_gam, beet_rich_gam3, test = "Chisq")
AIC(beet_rich_lm,beet_rich_gam, beet_rich_gam2, beet_rich_gam3, beet_rich_gam4)
# The full model is not significantly different from model 3. GAMS 1, 2, and 3
# all perform similarly.
summary(beet_rich_gam3)
# Removing lc reduces AIC, albeit only slightly.
summary(update(beet_rich_gam3, .~.-lc))
# I will report the values of model 3, which includes smooth terms for altitude
# and distance to forest, and linear terms for temperature and land cover.
appraise(beet_rich_gam3)
draw(beet_rich_gam3, residuals = T)

##### Vert diversity #####
vert_div_lm <- gam(s~lc+ghm+temp+fdist+alt, data = vert_gam_db, method = "REML")
vert_div_gam <- gam(s~lc+s(ghm)+s(temp)+s(fdist)+s(alt), data = vert_gam_db, method = "REML")
summary(vert_div_lm)
summary(vert_div_gam)
# All the smooth effects are reduced to linear terms in the GAM. Let's check
# collinearity in the lm.
performance::check_collinearity(vert_div_lm)
# temperature has moderate correlation, let's remove it.
vert_div_lm2 <- update(vert_div_lm, .~.-temp)
summary(vert_div_lm2)
AIC(vert_div_lm, vert_div_lm2)
anova(vert_div_lm, vert_div_lm2, test = "F")
# This model has significant linear effects of fdist and alt. You could remove the terms for lc and ghm.
AIC(vert_div_lm, vert_div_lm2, update(vert_div_lm2, .~.-lc), update(vert_div_lm2, .~.-lc-ghm))
anova(vert_div_lm, vert_div_lm2, update(vert_div_lm2, .~.-lc), update(vert_div_lm2, .~.-lc-ghm), test = "F")
vert_div_lm3 <- update(vert_div_lm2, .~.-lc-ghm)
summary(vert_div_lm3)
# The better performing model includes only linear effects, for fdist, and alt.
# The model is quite poor, with R2 = 0.144.
# This model has significant negative effect of log distance to forest, and a
# positive effect of altitude. Sites at higher altitude would thus have higher
# vertebrate diversity, and sites further away from forest would have lower
# diversity.


##### Beet diversity #####
beet_div_lm <- gam(s~lc+ghm+alt+temp+fdist, data=beet_gam_db, method = "REML")
summary(beet_div_lm) # no effects of anything
beet_div_gam <- gam(s~lc+s(ghm)+s(alt)+s(temp)+s(fdist), data=beet_gam_db, method = "REML") # in the full model the ghm and temp are both reduced to linear effects. There is also no sig effect of ghm or lc.
summary(beet_div_gam)
# The ghm and temp are reduced to linear terms, neither of which is significant. 
beet_div_gam2 <- update(beet_div_gam, .~.-s(ghm)-s(temp))
summary(beet_div_gam2)
concurvity(beet_div_gam2, full = F)
# removing these terms gives significant smooths of altitude and fdist. There is
# some concurvity among them, let's remove fdist.
beet_div_gam3 <- update(beet_div_gam2, .~.-s(fdist))
summary(beet_div_gam3)

AIC(beet_div_gam, beet_div_gam2, beet_div_gam3, beet_div_lm)
anova(beet_div_lm, beet_div_gam2, test="F")
# The lowest AIC is given by model 2, followed by the full model. THese are not
# significantly different in terms of how well they fit the data.
AIC(beet_div_gam2, update(beet_div_gam2, .~.-s(fdist)))
#' The model includes smooth effects of altitude and fdist. Overall it's hard to tell any real effect of these covariates on richness or diversity
draw(beet_div_gam2, residuals = T)

##### MFC ####
acous_gam_db <- left_join(summarised_MFC, allsites_cov_db) %>% 
  mutate(site = factor(site), lc = factor(lc2, labels = c("Dense","Open")), alt = scale(alt1), fdist = scale(log1p(lfdistance_2)), ghm = scale(ghm1), temp = scale(temp_mean1))
acous_gam <- gam(median_MFC~lc+s(alt)+s(ghm)+s(temp)+s(fdist)+s(site, bs = "re"),data = acous_gam_db, method = 'REML') # this code does not have explicit import of this dataset yet
appraise(acous_gam)
draw(acous_gam, residuals = T, parametric = T)
concurvity(acous_gam, full = F)
summary(acous_gam)
# remove fdist, which shows up as a linear effect
acous_gam2 <- update(acous_gam, .~.-s(fdist)+fdist)
concurvity(acous_gam2, full = F)
summary(acous_gam2)
draw(acous_gam2, residuals = T)
# remove ghm
acous_gam3 <- update(acous_gam, .~.-s(fdist)+fdist-s(ghm)+ghm)
concurvity(acous_gam3, full = F)
summary(acous_gam3)
draw(acous_gam3, residuals = T)
# remove alt
acous_gam4 <- update(acous_gam3, .~.-s(alt)+alt)
summary(acous_gam4)

AIC(acous_gam, acous_gam2, acous_gam3, acous_gam4)
# All models are nearly equivalent in terms of AIC.
performance::check_collinearity(acous_gam4)
# all linear terms are highly correlated, let's remove lc first
performance::check_collinearity(update(acous_gam4, .~.-lc-s(temp)-ghm-alt-fdist))
# All covariates have high VIF
acous_nullgamm <- gam(median_MFC~s(site, bs = "re"),data = acous_gam_db, method = 'REML')
summary(acous_nullgamm)
AIC(acous_nullgamm, acous_gam4, acous_gam3, acous_gam2, acous_gam)
#The null model has the lowest AIC, though still comparable. There is therefore
#no point in the smooths or linear predictors, all the variability can be
#explained by the site. 

### Similarity between taxa ####
# Mantel correlation between dissimilarity matrices of verts and beetles
# Find sites in common, subset data frames
commonsites <- intersect(rownames(comm), rownames(commb))
commb_commonsites <- commb[commonsites,]
comm_commonsites <- comm[commonsites,]
# calculate distances and mantel correlation
subsetdistvert <- dist(comm_commonsites)
subsetdistbeet <- dist(commb_commonsites)
mantel(subsetdistbeet, subsetdistvert)
# There is no statistically significant correlation. r = 0.0133, p = 0.283.

# Effect of distance
commonsitedist <- tibble(site = commonsites) %>% left_join(allsites_cov_db) %>% 
  select(easting, northing) %>% dist()
mantel(subsetdistbeet, commonsitedist)
mantel(subsetdistvert, commonsitedist)
# In the case of beetles, there is a significant, through not so strong  between
# the distance across sites and their dissimilarity. r = 0.169, p = 0.001.
plot(commonsitedist[upper.tri(commonsitedist)],subsetdistbeet[upper.tri(subsetdistbeet)])

####Things to check ####
