# Community diversity analyses. 

# Author: Juan S. Vargas
# Date: 25/4/25

# Data were obtained from camera-traps and pitfall traps placed along an
# altitudinal gradient in southern Costa Rica

# Load libraries
library(vegan)
library(tidyverse)
library(mgcv)
library(ape)
library(AICcPermanova)
library(gratia)
library(performance)
library(ggpubr)
library(datawizard)

#### Import data ####
# Import camera-trap data as matrix for vegan
comm <- read.csv("Data/vert_comm_mat.csv", row.names = 1)

# Import pitfall trap data
commb <- read.csv("Data/beet_comm_mat.csv", row.names = 1)  

# Import station data #
allsites_pts <- terra::vect("Data/combined_survey_pts.geojson")
allsites_cov_db <- as.data.frame(allsites_pts)

#### Community composition #### 
# Standardize matrices to account for different
#effort across sites. We divide by effort
comm_std <- comm/effort
# Calculate disimilarities
vert_dissim <- vegdist(comm_std)
beet_dissim <- vegdist(commb)

vert_covs_raw <- data.frame(site = rownames(comm)) %>% left_join(allsites_cov_db) %>% 
  select(site, alt1, lc2, ghm1, lfdistance_2) 
vert_covs <- mutate(vert_covs_raw, 
                    lc = factor(lc2, labels = c("Dense","Open")),
                    alt = standardise(alt1),
                    alt2 = standardise(alt1)^2,
                    # temp = scale(temp_mean1),
                    ghm = standardise(ghm1),
                    fdist = standardise(log1p(lfdistance_2))) %>% 
  select(lc, alt, alt2, ghm, fdist)
  
beet_covs_raw <- data.frame(site = rownames(commb)) %>% left_join(allsites_cov_db) %>% 
  select(site, alt1, lc2, ghm1, lfdistance_2) 
beet_covs <- mutate(beet_covs_raw, lc = factor(lc2, labels = c("Dense","Open")),
         alt = standardise(alt1),
         alt2 = standardise(alt1)^2,
         # temp = scale(temp_mean1),
         ghm = standardise(ghm1),
         fdist = standardise(log1p(lfdistance_2))) %>% 
  select(lc,alt,alt2, ghm,fdist)

##### PERMANOVA #####

# Analysis to see how the different variables affect the dissimilarity across
# sites
permanova_comp <- fit_models(make_models(vars = c("lc", "alt", "alt2","ghm", "fdist")), 
                             veg_data = comm_std, env_data = vert_covs)
select_models(permanova_comp)
# We compared different covariate combinations, and the one with the lowest AICc
# included land cover, altitude, and alt^2. However, there are 3 other models
# that have similar AICc. The VIF for the best model is not that high (2.35) so
# we could ignore it.
permanova_comp_beet <- fit_models(make_models(vars = c("lc", "alt","alt2",  "ghm", "fdist")), 
                             veg_data = commb, env_data = beet_covs)
select_models(permanova_comp_beet)
# In the case of beetles, all 6 top models rank similarly. The best one includes
# altitude, ghm, and lc. THe VIF for this one is 1.98. It performs similar to 
# the same model with lc instead of ghm, and the other variables combined with alt.
vert_permanova <- adonis2(comm_std~lc+alt+alt2, data = vert_covs, method = 'bray')
beet_permanova <- adonis2(commb~lc+alt+ghm, data = beet_covs, method = 'bray')
vert_permanova
beet_permanova

  # These results suggest that differences in species composition for vertebrates
# are influenced mostly by differences in altitude and land cover, with some
# smaller effect of disturbance.


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
# (betadisper: F_1,115 = 3.77, p = 0.055). It is important to keep this in mind
# as it may be driving apparent differences due to other covariates

anova(betadisper(beet_dissim, beet_covs$lc))
plot(betadisper(beet_dissim, beet_covs$lc), 
     hull = F, ellipse=T, col = c("darkgreen", "goldenrod"), segments = F)
# For beetles there is no difference in variability. F_1,70 = 0.581, p = 0.449.


# Plot PCoA
par(mfrow = c(1,2))
plot(vert_pcoa$vectors[,1:2], las=1, asp = 1, xlab = "PCoA 1", ylab = "PCoA 2", type = 'n', main = "Vertebrates")
ordiellipse(vert_pcoa$vectors, vert_covs$lc, draw = "polygon", col = c("darkgreen", "goldenrod"))
points(vert_pcoa$vectors[vert_covs$lc=="Dense",1:2], pch = 16, col = "darkgreen")
points(vert_pcoa$vectors[vert_covs$lc=="Open",1:2], pch = 17, col = "goldenrod")
ordisurf(vert_pcoa$vectors, vert_covs_raw$alt1, add = T, col = "gray50")
plot(vert_envfit_pcoa, labels = list(factors = c("Dense", "Open"), vectors = c("Alt", "Alt^2", "gHM", "Forest dist")), bg='white', col = "gray10")

# PCoA for beetles
plot(beet_pcoa$vectors[,1:2], type = 'n', las=1, asp = 1, xlab = "PCoA 1", ylab = "PCoA 2", main = "Beetles")
ordiellipse(beet_pcoa$vectors, beet_covs$lc, draw = "polygon", col = c("darkgreen", "goldenrod"))
points(beet_pcoa$vectors[beet_covs$lc=="Dense",1:2], pch = 16, col = "darkgreen")
points(beet_pcoa$vectors[beet_covs$lc=="Open",1:2], pch = 17, col = "goldenrod")
ordisurf(beet_pcoa$vectors, beet_covs_raw$alt1, add = T, col = "gray50")
plot(beet_envfit_pcoa, labels = list(factors = c("Dense", "Open"), vectors = c("Alt", "Alt^2", "gHM", "Forest dist")), bg='white', col = "gray10")
par(mfrow = c(1,1))

## NMDS approach ##
vert_NMDS <- metaMDS(comm_std, trymax = 100)
beet_NMDS <- metaMDS(commb)

# Fit environ variables
vert_envfit_nmds <- envfit(vert_NMDS, vert_covs)
beet_envfit_nmds <- envfit(beet_NMDS, beet_covs)



#### Community metrics ####

# plot richness (observed) against covariates individually
plot(specnumber(comm)~vert_covs$alt)
plot(specnumber(comm)~vert_covs$ghm)
plot(specnumber(comm)~vert_covs$fdist)

plot(specnumber(commb)~beet_covs$alt)
plot(specnumber(commb)~beet_covs$ghm)
plot(specnumber(commb)~beet_covs$fdist)

# plot diversity (Shannon) against covariates individually
plot(diversity(comm)~vert_covs$alt)
plot(diversity(comm)~vert_covs$ghm)
plot(diversity(comm)~vert_covs$fdist)
# beetles
plot(diversity(commb)~beet_covs$alt)
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
# The full model reduces ghm to a linear effects
vert_rich_gam <- update(vert_rich_gam, .~.-s(ghm)+ghm)
check_concurvity(vert_rich_gam)
# There is high concurvity for fdist, moderate for the parametric terms.
vert_rich_gam2 <- update(vert_rich_gam, .~. -s(fdist)+fdist)
summary(vert_rich_gam2)
check_concurvity(vert_rich_gam2)
check_collinearity(vert_rich_gam2)
# This model has significant smooth for altitude. The linear terms are not
# significant. lc and fdist could be correlated, I'll remove fdist first since it
# has the highest VIF (6.75)
vert_rich_gam3 <- update(vert_rich_gam2, .~.-fdist)
summary(vert_rich_gam3) # sig smooth of alt (X2 = 16.98, p = 0.0060), sig neg effect of open lc (beta = -0.41, z = -3.929 p<0.001)
vert_rich_gam4 <- update(vert_rich_gam3, .~.-ghm) # ghm was not significant
vert_rich_gam5 <- update(vert_rich_gam3, .~.-lc-ghm) # alternative including fdist instead of lc

AIC(vert_rich_gam, vert_rich_gam2, vert_rich_gam3, vert_rich_lm, vert_rich_gam4, vert_rich_gam5)
# The lowest AIC is given by model 4 (smooth alt, lc), but it is comparable to
# model 3 (smooth alt, lc, )
draw(vert_rich_gam4, residuals = T)
appraise(vert_rich_gam4)
anova(vert_rich_gam4, vert_rich_gam3, test = 'Chisq')
# I will report the results of model 4, which has smooth terms for altitude and
# a fixed effect of land cover
summary(vert_rich_gam4) # smooth is sig (X2 = 17.0, p = 0.0060), sig lc (beta = -0.410, z = -4.247, p<0.001)

##### Beet richness #####

# linear model
beet_rich_lm <- gam(n~lc+ghm+alt+fdist, data = beet_gam_db, method = 'REML')
summary(beet_rich_lm) # no sig effects
check_collinearity(beet_rich_lm) # lfdist has the highest VIF (10), and lc has moderate (VIF = 9.47)
# This model has moderate correlation between lc and fdist, let's remove lc.
beet_rich_lm2 <- update(beet_rich_lm, .~.-fdist)
check_collinearity(beet_rich_lm2) # removing fdist reduces all VIF 

# beetle richness GAM
beet_rich_gam <- gam(n~lc+s(ghm)+s(alt)+s(fdist), 
                     family = poisson, data = beet_gam_db, method = 'REML')
summary(beet_rich_gam)
# no significant smooths, nor linear effects. ghm and Fdist reduced to linear
# effect. 
beet_rich_gam <- update(beet_rich_gam, .~.-s(ghm)+ghm-s(fdist)+fdist)
check_concurvity(beet_rich_gam)
# There is moderate concurvity in the parametric terms. Highest collinearity is for fdist
check_collinearity(beet_rich_gam)
beet_rich_gam2 <- update(beet_rich_gam, .~.-fdist)
summary(beet_rich_gam2) # sig effect of lc (beta = -0.29, z = -2.209, p = 0.0272)
check_concurvity(beet_rich_gam2) # low concurvity
check_collinearity(beet_rich_gam2) # low correlations in parametric terms

beet_rich_gam3 <- update(beet_rich_gam2, .~.-ghm) # remove non sig ghm term
summary(beet_rich_gam3) # significant altitude smooth (X2 = 18.22, p = 0.00263), sig lc (-0.284, z = -2.125, p = 0.0336)

# Rank models
AIC(beet_rich_lm,beet_rich_lm2,beet_rich_gam, beet_rich_gam2, beet_rich_gam3)
# Lowest AIC given by model 3 (s(alt) + lc). Comparable to model 2 though (+ ghm)
anova(beet_rich_gam3, beet_rich_gam2, test = "Chisq")

# I will report the values of model 3, which includes smooth terms for altitude
# and a linear term of lc
appraise(beet_rich_gam3)
draw(beet_rich_gam3, residuals = T)

##### Vert diversity #####
# linear model
vert_div_lm <- gam(s~lc+ghm+fdist+alt, data = vert_gam_db, method = "REML")
summary(vert_div_lm)
check_collinearity(vert_div_lm) # moderate VIF for fdist (8.48) and lc (7.95). I will remove fdist
vert_div_lm2 <- update(vert_div_lm, .~.-fdist)
# GAM
vert_div_gam <- gam(s~lc+s(ghm)+s(fdist)+s(alt), data = vert_gam_db, method = "REML")
summary(vert_div_gam)
# All the smooth effects are reduced to linear terms in the GAM. So I'll stick with lm
summary(vert_div_lm2) # sig lc (-0.414, t = -3.637, p = 0.0004) and alt (0.080, t = 1.989, p = 0.049).
vert_div_lm3 <- update(vert_div_lm2, .~.-ghm) # remove nonsig ghm
summary(vert_div_lm3)# sig fdist and alt. Coefs are nearly identical. lc (-0.408, t = -4.431, p<0.001) and alt (0.080, t = 1.996, p = 0.0483).

AIC(vert_div_lm, vert_div_lm2, vert_div_lm3) # lm3 has the lowest AIC
anova(vert_div_lm3, vert_div_lm2, test = "F") # lm2 and lm3 are not sig. diff.
# This model has significant linear effects of lc and alt. The model is quite
# poor, with R2 = 0.15. Sites at higher altitude would thus have higher
# vertebrate diversity, and sites in open forest sites would have lower
# diversity.

##### Beet diversity #####
beet_div_lm <- gam(s~lc+ghm+alt+fdist, data=beet_gam_db, method = "REML")
summary(beet_div_lm) # sig ghm and alt

beet_div_gam <- gam(s~lc+s(ghm)+s(alt)+s(fdist), data=beet_gam_db, method = "REML") # in the full model the ghm and temp are both reduced to linear effects. There is also no sig effect of ghm or lc.
summary(beet_div_gam)
# All terms are reduced to linear effects.
check_collinearity(beet_div_lm) # High corr of fdist
beet_div_lm2 <- update(beet_div_lm, .~.-fdist)
summary(beet_div_lm2)
check_collinearity(beet_div_lm2) # Correlations are low. 
# remove nonsig lc term
beet_div_lm3 <- update(beet_div_lm2, .~.-lc)
summary(beet_div_lm3) # ghm no longer sig
beet_div_lm4 <- update(beet_div_lm3, .~.-ghm)
summary(beet_div_lm4)
beet_div_lm5 <- update(beet_div_lm2, .~.-ghm)
summary(beet_div_lm5)
AIC(beet_div_lm, beet_div_lm2, beet_div_lm3, beet_div_lm4, beet_div_lm5) 
anova(beet_div_lm3, beet_div_lm2, test='F')
# all models have comparable AIC, the top two (2 and 3) are nearly identical. I
# will stick with model 2 (alt, ghm, lc).

##### MFC ####
acous_gam_db <- left_join(summarised_MFC, allsites_cov_db) %>% 
  mutate(site = factor(site), lc = factor(lc2, labels = c("Dense","Open")), alt = scale(alt1), fdist = scale(log1p(lfdistance_2)), ghm = scale(ghm1), temp = scale(temp_mean1))
acous_glmm <- gam(median_MFC~lc+alt+ghm+fdist+s(site, bs = 're'), data = acous_gam_db, method = "REML")
summary(acous_glmm) # sig. random effect of site. No sig linear terms.
check_collinearity(acous_glmm) # very high corr, highest is for lc
# GAM
acous_gam <- gam(median_MFC~lc+s(alt)+s(ghm)+s(fdist)+s(site, bs = "re"),data = acous_gam_db, method = 'REML') # this code does not have explicit import of this dataset yet
summary(acous_gam) # ghm and fdist are linear
acous_gam <- gam(median_MFC~lc+s(alt)+ghm+fdist+s(site, bs = "re"),data = acous_gam_db, method = 'REML')
summary(acous_gam)
check_concurvity(acous_gam) # high concurvity in parametric terms
check_collinearity(acous_gam) # extremely high VIF, starting with lc
acous_gam2 <- update(acous_gam, .~.-lc)
summary(acous_gam2)
check_collinearity(acous_gam2)
acous_gam3 <- update(acous_gam2, .~.-ghm)
summary(acous_gam3) # alt reduced to nearly linear
acous_gam4 <- update(acous_gam3, .~.-s(alt)+alt)
summary(acous_gam4)
check_collinearity(acous_gam4) # High corr between fdist and alt.
acous_gam5 <- update(acous_gam4, .~.-fdist)
summary(acous_gam5)
AIC(acous_gam, acous_gam2, acous_gam3, acous_gam4, acous_gam5)
# model 5 (alt+ranef) has the lowest AIC.
acous_nullgamm <- gam(median_MFC~s(site, bs = "re"),data = acous_gam_db, method = 'REML')
acous_nullgamm2 <- gam(median_MFC~1,data = acous_gam_db, method = 'REML')
summary(acous_nullgamm)
summary(acous_nullgamm2)
AIC(acous_nullgamm, acous_nullgamm2, acous_gam5)
# The null model with only random effects has comparable AIC as model 5


##### Visualize diversity ####
# Vert richness
pred_alts_raw <- seq(min(vert_covs_raw$alt1), max(vert_covs_raw$alt1),length.out = 100)
newdata_raw <- expand.grid(intercept = 1, lc = levels(vert_covs$lc), alt = pred_alts_raw)
pred_alts <- seq(min(vert_covs$alt), max(vert_covs$alt),length.out = 100)
newdata <- expand.grid(intercept = 1, lc = levels(vert_covs$lc), alt = pred_alts)
pred_vert_rich <- predict.gam(vert_rich_gam4, newdata = newdata, type = 'response', se.fit = T)
vert_rich_pred_db <- data.frame(newdata_raw,n = pred_vert_rich$fit, se = pred_vert_rich$se.fit)
pvertrich <- data.frame(n = specnumber(comm),vert_covs_raw) %>% 
  mutate(lc = factor(lc2, labels = c("Dense", "Open"))) %>% 
  ggplot(aes(alt1,n))+
  geom_ribbon(aes(x = alt, ymin = n-se, ymax = n+se, fill = lc), data = vert_rich_pred_db)+
  geom_line(aes(alt,n, color = lc),vert_rich_pred_db)+
  geom_point(aes(color = lc, shape = lc))+
  labs(x = "Altitude (masl)", y = "Species richness", color = "Land cover", fill = "Land cover", shape = "Land cover")+
  scale_fill_manual(values = alpha(c("darkgreen", "goldenrod"), alpha = 0.3))+
  scale_color_manual(values = c("darkgreen", "goldenrod"))+
  theme_pubr(base_size = 14)

ggplot(vert_gam_db, aes(alt,n))+geom_point()+geom_line(aes(y = fitted(vert_rich_gam4),col=lc))
# Beetle richness
pred_alts_raw <- seq(min(beet_covs_raw$alt1), max(beet_covs_raw$alt1),length.out = 100)
newdata_raw <- expand.grid(intercept = 1, lc = levels(beet_covs$lc), alt = pred_alts_raw)
pred_alts <- seq(min(beet_covs$alt), max(beet_covs$alt),length.out = 100)
newdata <- expand.grid(intercept = 1, lc = levels(beet_covs$lc), alt = pred_alts)
pred_beet_rich <- predict.gam(beet_rich_gam3, newdata = newdata, type = 'response', se.fit = T)
beet_rich_pred_db <- data.frame(newdata_raw,n = pred_beet_rich$fit, se = pred_beet_rich$se.fit)
pbeetrich <- data.frame(n = specnumber(commb),beet_covs_raw) %>% 
  mutate(lc = factor(lc2, labels = c("Dense", "Open"))) %>% 
  ggplot(aes(alt1,n))+
  geom_ribbon(aes(x = alt, ymin = n-se, ymax = n+se, fill = lc), data = beet_rich_pred_db)+
  geom_line(aes(alt,n, color = lc),beet_rich_pred_db)+
  geom_point(aes(color = lc, shape = lc))+
  labs(x = "Altitude (masl)", y = "Species richness", color = "Land cover", fill = "Land cover", shape = "Land cover")+
  scale_fill_manual(values = alpha(c("darkgreen", "goldenrod"), alpha = 0.3))+
  scale_color_manual(values = c("darkgreen", "goldenrod"))+
  theme_pubr(base_size = 14)
# Vert diversity
pred_alts_raw <- seq(min(vert_covs_raw$alt1), max(vert_covs_raw$alt1),length.out = 100)
newdata_raw <- expand.grid(intercept = 1, lc = levels(vert_covs$lc), alt = pred_alts_raw)
pred_alts <- seq(min(vert_covs$alt), max(vert_covs$alt),length.out = 100)
newdata <- expand.grid(intercept = 1, lc = levels(vert_covs$lc), alt = pred_alts)
pred_vert_div <- predict(vert_div_lm3, newdata = newdata, type = 'response', se.fit = T)
vert_div_pred_db <- data.frame(newdata_raw,n = pred_vert_div$fit, se = pred_vert_div$se.fit)
pvertdiv <- data.frame(n = diversity(comm),vert_covs_raw) %>% 
  mutate(lc = factor(lc2, labels = c("Dense", "Open"))) %>% 
  ggplot(aes(alt1,n))+
  geom_ribbon(aes(x = alt, ymin = n-se, ymax = n+se, fill = lc), data = vert_div_pred_db)+
  geom_line(aes(alt,n, color = lc),vert_div_pred_db)+
  geom_point(aes(color = lc, shape = lc))+
  labs(x = "Altitude (masl)", y = "Diversity", color = "Land cover", fill = "Land cover", shape = "Land cover")+
  scale_fill_manual(values = alpha(c("darkgreen", "goldenrod"), alpha = 0.3))+
  scale_color_manual(values = c("darkgreen", "goldenrod"))+
  theme_pubr(base_size = 14)

# Beet diversity
pred_alts_raw <- seq(min(beet_covs_raw$alt1), max(beet_covs_raw$alt1),length.out = 100)
newdata_raw <- expand.grid(intercept = 1, lc = levels(beet_covs$lc), alt = pred_alts_raw)
pred_alts <- seq(min(beet_covs$alt), max(beet_covs$alt),length.out = 100)
newdata <- expand.grid(intercept = 1, lc = levels(beet_covs$lc), ghm = 0, alt = pred_alts)
pred_beet_div <- predict(beet_div_lm2, newdata = newdata, type = 'response', se.fit = T)
beet_div_pred_db <- data.frame(newdata_raw,n = pred_beet_div$fit, se = pred_beet_div$se.fit)
pbeetdiv <- data.frame(n = diversity(commb),beet_covs_raw) %>% 
  mutate(lc = factor(lc2, labels = c("Dense", "Open"))) %>% 
  ggplot(aes(alt1,n))+
  geom_ribbon(aes(x = alt, ymin = n-se, ymax = n+se, fill = lc), data = beet_div_pred_db)+
  geom_line(aes(alt,n, color = lc),beet_div_pred_db)+
  geom_point(aes(color = lc, shape = lc))+
  labs(x = "Altitude (masl)", y = "Diversity", color = "Land cover", fill = "Land cover", shape = "Land cover")+
  scale_fill_manual(values = alpha(c("darkgreen", "goldenrod"), alpha = 0.3))+
  scale_color_manual(values = c("darkgreen", "goldenrod"))+
  theme_pubr(base_size = 14)

ggarrange(pvertrich,pbeetrich,pvertdiv,pbeetdiv, labels = 'auto', common.legend = T, align = 'hv')
### Similarity between taxa ####
# Mantel correlation between dissimilarity matrices of verts and beetles
# Find sites in common, subset data frames
commonsites <- intersect(rownames(comm), rownames(commb))
commb_commonsites <- commb[commonsites,]
comm_commonsites <- comm_std[commonsites,]
# calculate betadiversity as dissimilarity and mantel correlation
subsetdistvert <- betadiver(comm_commonsites, method = "w")
subsetdistbeet <- betadiver(commb_commonsites, method = "w")

mantel(subsetdistbeet, subsetdistvert)
# There is a statistically significant correlation. r = 0.175, p = 0.001. (with Bray distance)
# r = 0.06569, p = 0.111 (for betadiver metric - Sorensen, presence/absence)
plot(subsetdistbeet, subsetdistvert)

# Include Effect of distance, calculate partial statistic
commonsitedist <- tibble(site = commonsites) %>% left_join(allsites_cov_db) %>% 
  select(easting, northing) %>% dist()
mantel.partial(subsetdistbeet, subsetdistvert, commonsitedist)
mantel(subsetdistvert, commonsitedist)
# Geographical distance is correlated with community dissimilarity for both
# vertebrates (r = 0.1977, p = 0.001) and beetles (r = 0.321, p = 0.001). There
# is a lot of noise in this relationship though.

# Compare and plot rich, div, for both taxa, and compare with MFC
# aggregate the index calculated into weekly values, and only taking sunrise and sunset values
comb_metrics_db <- acoustic_jointdb %>% filter(hour(dttm)>5) %>% # keep only daytime recordings
  summarise(mean_MFC=mean(MFC), .by = site) %>% 
  select(site,mean_MFC) %>%
  full_join(data.frame(site = rownames(comm), n_vert = specnumber(comm), s_vert = diversity(comm))) %>% 
  full_join(data.frame(site = rownames(commb), n_beet = specnumber(commb), s_beet = diversity(commb))) %>% 
  left_join(allsites_cov_db)

p1 <- ggplot(comb_metrics_db, aes(n_vert,n_beet))+geom_point()+theme_pubr(base_size = 14)+
  labs(x = "Vertebrate richness", y = "Beetle richness")
  # stat_smooth(method = 'lm', color = "gray50")
p2 <- ggplot(comb_metrics_db, aes(n_vert,mean_MFC))+geom_point()+theme_pubr(base_size = 14)+
  labs(x = "Vertebrate richness", y = "MFC")
  # stat_smooth(method = 'lm')
p3 <- ggplot(comb_metrics_db, aes(mean_MFC, n_beet))+geom_point()+theme_pubr(base_size = 14)+
  labs(x = "Beetle richness", y = "MFC")
  # stat_smooth(method = 'lm')
p4 <- ggplot(comb_metrics_db, aes(s_vert,s_beet))+geom_point()+theme_pubr(base_size = 14)+
  labs(x = "Vertebrate diversity", y = "Beetle diversity")
# stat_smooth(method = 'lm')

ggarrange(p1,p2,p4,p3, labels='auto')
####Things to check ####
# Maybe analyze associations between large herbivores and specific dung beetle species/richness?