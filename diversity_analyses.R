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
# Standardize matrices to account for different
#effort across sites. We divide by effort
comm_std <- comm/effort
# Calculate distances
vert_dissim <- vegdist(comm_std)
beet_dissim <- vegdist(commb)
vert_covs_raw <- data.frame(site = rownames(comm)) %>% left_join(allsites_cov_db) %>% 
  select(site, alt1, lc2, ghm1, lfdistance_2) 
vert_covs <- mutate(vert_covs_raw, lc = factor(lc2, labels = c("Dense","Open")),
         alt = scale(alt1),
         # temp = scale(temp_mean1),
         ghm = scale(ghm1),
         fdist = scale(log1p(lfdistance_2))) %>% 
  select(lc,alt,ghm,fdist)
  
beet_covs_raw <- data.frame(site = rownames(commb)) %>% left_join(allsites_cov_db) %>% 
  select(site, alt1, lc2, ghm1, lfdistance_2) 
beet_covs <- mutate(beet_covs_raw, lc = factor(lc2, labels = c("Dense","Open")),
         alt = scale(alt1),
         # temp = scale(temp_mean1),
         ghm = scale(ghm1),
         fdist = scale(log1p(lfdistance_2))) %>% 
  select(lc,alt,ghm,fdist)

##### PERMANOVA #####

# Analysis to see how the different variables affect the dissimilarity across
# sites
permanova_comp <- fit_models(make_models(vars = c("lc", "alt", "ghm", "fdist")), 
                             veg_data = comm_std, env_data = vert_covs)
select_models(permanova_comp)
# We compared different covariate combinations, and the one with the lowest AICc
# included land cover, altitude, and ghm. However, there are
# five other models that have similar AICc. Among the candidates, the second
# best one also has the lowest VIF. This one includes altitude and land cover.
# The VIF for the best model is not that high (1.67) so we could ignore
# it.
permanova_comp_beet <- fit_models(make_models(vars = c("lc", "alt", "ghm", "fdist")), 
                             veg_data = commb, env_data = beet_covs)
select_models(permanova_comp_beet)
# In the case of beetles, all 6 top models rank similarly. The best one includes
# altitude, ghm, and distance to forest. THe VIF for this one is 2.01.
vert_permanova <- adonis2(comm_std~lc+alt+ghm, data = vert_covs, method = 'bray')
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
# (betadisper: F_1,115 = 4.30, p = 0.040). It is important to keep this in mind
# as it may be driving apparent differences due to other covariates

anova(betadisper(beet_dissim, beet_covs$lc))
plot(betadisper(beet_dissim, beet_covs$lc), 
     hull = F, ellipse=T, col = c("darkgreen", "goldenrod"), segments = F)
# For beetles there is no difference in variability. F_1,70 = 0.79, p = 0.38.

# Plot PCoA
plot(vert_pcoa$vectors[,1:2], las=1, asp = 1, xlab = "PCoA 1", ylab = "PCoA 2", type = 'n')
ordiellipse(vert_pcoa$vectors, vert_covs$lc, draw = "polygon", col = c("darkgreen", "goldenrod"))
points(vert_pcoa$vectors[vert_covs$lc=="Dense",1:2], pch = 16, col = "darkgreen")
points(vert_pcoa$vectors[vert_covs$lc=="Open",1:2], pch = 17, col = "goldenrod")
ordisurf(vert_pcoa$vectors, vert_covs_raw$alt1, add = T, col = "gray50")
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
vert_rich_lm <- gam(n~lc+ghm+alt+fdist, family = poisson, data=vert_gam_db)
# Full GAM, with smooths for all covars
vert_rich_gam <- gam(n~lc+s(ghm)+s(alt)+s(fdist), 
                     family = poisson, data = vert_gam_db, method = 'REML')
summary(vert_rich_gam)
# The full model reduces the ghm and fdist to linear effects.
check_concurvity(vert_rich_gam)
# There is moderate concurvity for fdist, followed by alt.
# Let's make ghm and fdist explictily linear.
vert_rich_gam2 <- update(vert_rich_gam, .~. - s(ghm)+ghm-s(fdist)+fdist)
summary(vert_rich_gam2)
check_concurvity(vert_rich_gam2)
check_collinearity(vert_rich_gam2)
# This model has significant smooth for altitude. The linear terms are not
# significant. lc and fdist could be correlated, I'll remove lc first since it
# has the highest VIF (4.57)
vert_rich_gam3 <- update(vert_rich_gam2, .~.-lc)
summary(vert_rich_gam3) # sig smooth of alt (X2 = 18.04, p = 0.0038), sig neg effect of fdist (beta = -0.18, z = -3.937, p<0.001)
vert_rich_gam4 <- update(vert_rich_gam3, .~.-ghm) # ghm was not significant
vert_rich_gam5 <- update(vert_rich_gam2, .~.-fdist-ghm)

AIC(vert_rich_gam, vert_rich_gam2, vert_rich_gam3, vert_rich_lm, vert_rich_gam4, vert_rich_gam5)
# The lowest AIC is given by model 5 (smooth alt, lc), but it is comparable to
# model 4 (smooth alt, fdist)
draw(vert_rich_gam5, residuals = T)
appraise(vert_rich_gam5)
anova(vert_rich_gam4, vert_rich_gam5, test = 'Chisq')
# I will report the results of model 5, which has smooth terms for altitude and
# a fixed effect of land cover
summary(vert_rich_gam5) # smooth is sig (X2 = 15.17, p = 0.0116), sig lc (beta = -0.408, z = -4.456, p<0.001)

##### Beet richness #####

# linear model
beet_rich_lm <- gam(n~lc+ghm+alt+fdist, data = beet_gam_db, method = 'REML')
summary(beet_rich_lm) # no sig effects
check_collinearity(beet_rich_lm) # land cover and distance to forest are correlated.
# This model has moderate correlation between lc and fdist, let's remove lc.
beet_rich_lm2 <- update(beet_rich_lm, .~.-lc)

# beetle richness GAM
beet_rich_gam <- gam(n~lc+s(ghm)+s(alt)+s(fdist), 
                     family = poisson, data = beet_gam_db, method = 'REML')
summary(beet_rich_gam)
check_concurvity(beet_rich_gam)
# no significant smooths, nor linear effects. Fdist is reduced to a linear
# effect. Furthermore, there is noticeable concurvity, highest in the full model
# is for distance to forest.
beet_rich_gam2 <- update(beet_rich_gam, .~.-s(fdist)+fdist)
summary(beet_rich_gam2)
check_concurvity(beet_rich_gam2)
check_collinearity(beet_rich_gam2)
#' This model has no significant smooths, but lc does have a low p = 0.07
#' Concurvity is now high among the parametric terms. LC has the highest VIF

beet_rich_gam3 <- update(beet_rich_gam2, .~.-lc)
summary(beet_rich_gam3) # significant altitude smooth (X2 = 17.13, p = 0.00356), lc no longer sig.
check_concurvity(beet_rich_gam3)
# Concurvity is now low
beet_rich_gam4 <- update(beet_rich_gam3, .~.-s(ghm)+ghm)
summary(beet_rich_gam4)
check_collinearity(beet_rich_gam4)
beet_rich_gam5 <- update(beet_rich_gam4, .~.-ghm)
summary(beet_rich_gam5) # sig smooth alt, fdist not sig.
# Rank models
AIC(beet_rich_lm,beet_rich_lm2,beet_rich_gam, beet_rich_gam2, beet_rich_gam3, beet_rich_gam4, beet_rich_gam5)
# Most models are comparable, the lowest is model 5 (s(alt)+fdist). 

anova(beet_rich_gam5,beet_rich_lm2, test = "Chisq") # res deviance is much lower in the GAM
anova(beet_rich_gam5, beet_rich_lm, test = "Chisq")
anova(beet_rich_gam5, beet_rich_gam, test = "Chisq")
anova(beet_rich_gam5, beet_rich_gam2, test = "Chisq")
# Model 5 is not significantly different from the full model or gam 2, though
# they are different from the linear models
summary(beet_rich_gam5)
# I will report the values of model 5, which includes smooth terms for altitude
# and a linear term of fdist.
appraise(beet_rich_gam5)
draw(beet_rich_gam5, residuals = T)

##### Vert diversity #####
# linear model
vert_div_lm <- gam(s~lc+ghm+fdist+alt, data = vert_gam_db, method = "REML")
summary(vert_div_lm)
check_collinearity(vert_div_lm) # VIF generally low, highest is 4.85 for lc

# GAM
vert_div_gam <- gam(s~lc+s(ghm)+s(fdist)+s(alt), data = vert_gam_db, method = "REML")
summary(vert_div_gam)
# All the smooth effects are reduced to linear terms in the GAM.
# Checking collinearity all VIF are low, the higest is for lc followed by fdist.
vert_div_lm2 <- update(vert_div_lm, .~.-lc)
summary(vert_div_lm2) # sig fdist and alt.
vert_div_lm3 <- update(vert_div_lm2, .~.-ghm)
summary(vert_div_lm3)# sig fdist and alt. Coefs are nearly identical.

AIC(vert_div_lm, vert_div_lm2, vert_div_lm3) # lm3 has the lowest AIC
anova(vert_div_lm3, vert_div_lm2, test = "F") # lm2 and lm3 are not sig. diff.
# This model has significant linear effects of fdist and alt. You could remove the terms for lc and ghm.
summary(vert_div_lm3)
# The better performing model includes only linear effects, for fdist, and alt.
# The model is quite poor, with R2 = 0.144. This model has significant negative
# effect of distance to forest (beta = -0.17, t = -4.318, p<0.001), and a
# positive effect of altitude (beta = 0.082, t = 2.055, p = 0.0422). Sites at
# higher altitude would thus have higher vertebrate diversity, and sites further
# away from forest would have lower diversity.

##### Beet diversity #####
beet_div_lm <- gam(s~lc+ghm+alt+fdist, data=beet_gam_db, method = "REML")
summary(beet_div_lm) # sig ghm and alt

beet_div_gam <- gam(s~lc+s(ghm)+s(alt)+s(fdist), data=beet_gam_db, method = "REML") # in the full model the ghm and temp are both reduced to linear effects. There is also no sig effect of ghm or lc.
summary(beet_div_gam)
# All terms are reduced to linear effects.
check_collinearity(beet_div_lm) # moderate correlation of lc and fdist
beet_div_lm2 <- update(beet_div_lm, .~.-lc)
summary(beet_div_lm2)
check_collinearity(beet_div_lm2)
# Correlations are low. 
AIC(beet_div_lm, beet_div_lm2)
anova(beet_div_lm, beet_div_lm2, test='F')
# model 2 has lower AIC, but not significantly different

##### MFC ####
acous_gam_db <- left_join(summarised_MFC, allsites_cov_db) %>% 
  mutate(site = factor(site), lc = factor(lc2, labels = c("Dense","Open")), alt = scale(alt1), fdist = scale(log1p(lfdistance_2)), ghm = scale(ghm1), temp = scale(temp_mean1))
acous_glmm <- gam(median_MFC~lc+alt+ghm+fdist+s(site, bs = 're'), data = acous_gam_db, method = "REML")
summary(acous_glmm) # sig. random effect of site.
# GAM
acous_gam <- gam(median_MFC~lc+s(alt)+s(ghm)+s(fdist)+s(site, bs = "re"),data = acous_gam_db, method = 'REML') # this code does not have explicit import of this dataset yet
summary(acous_gam) # ghm and fdist are linear
acous_gam <- gam(median_MFC~lc+s(alt)+ghm+fdist+s(site, bs = "re"),data = acous_gam_db, method = 'REML')
summary(acous_gam)
check_concurvity(acous_gam) # high concurvity in parametric terms
check_collinearity(acous_gam) # extremely high VIF, starting with fdist
acous_gam2 <- update(acous_gam, .~.-fdist)
summary(acous_gam2)
check_concurvity(acous_gam2) # still high concurvity in parametric
check_collinearity(acous_gam2)
acous_gam3 <- update(acous_gam2, .~.-ghm)
summary(acous_gam3)
check_concurvity(acous_gam3)
check_collinearity(acous_gam3)
summary(acous_gam3)
# remove alt
acous_gam4 <- update(acous_gam3, .~.-s(alt)+alt)
summary(acous_gam4)

AIC(acous_gam, acous_gam2, acous_gam3, acous_gam4)
# All models are nearly equivalent in terms of AIC.
check_collinearity(acous_gam4) # still high VIF

acous_nullgamm <- gam(median_MFC~s(site, bs = "re"),data = acous_gam_db, method = 'REML')
summary(acous_nullgamm)
AIC(acous_nullgamm, acous_gam4, acous_gam3, acous_gam2, acous_gam)
# GAM 3 has the lowest AIC, but all are comparable, and it is almost identical to models 1 and 2. 
draw(acous_gam3, residuals = T)

### Visualize diversity ####
# Vert richness
draw(vert_rich_gam4, residuals = T)
# Beetle richness
draw(beet_rich_gam3, residuals = T, parametric = T)
# Vert diversity

### Similarity between taxa ####
# Mantel correlation between dissimilarity matrices of verts and beetles
# Find sites in common, subset data frames
commonsites <- intersect(rownames(comm), rownames(commb))
commb_commonsites <- commb[commonsites,]
comm_commonsites <- comm_std[commonsites,]
# calculate distances and mantel correlation
subsetdistvert <- vegdist(comm_commonsites)
subsetdistbeet <- vegdist(commb_commonsites)
mantel(subsetdistbeet, subsetdistvert)
# There is a statistically significant correlation. r = 0.175, p = 0.001.
plot(subsetdistbeet, subsetdistvert)

# Effect of distance
commonsitedist <- tibble(site = commonsites) %>% left_join(allsites_cov_db) %>% 
  select(easting, northing) %>% dist()
mantel(subsetdistbeet, commonsitedist)
mantel(subsetdistvert, commonsitedist)
# Geographical distance is correlated with community dissimilarity for both
# vertebrates (r = 0.1977, p = 0.001) and beetles (r = 0.321, p = 0.001). There
# is a lot of noise in this relationship though.

####Things to check ####
