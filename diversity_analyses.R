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

# Import camera-trap data
comm <- read.csv("Data/vert_comm_mat.csv", row.names = 1)

# Import pitfall trap data
commb <- read.csv("Data/beet_comm_mat.csv", row.names = 1)

# Import station data #
allsites_pts <- vect("Data/combined_survey_pts.geojson")
allsites_cov_db <- as.data.frame(allsites_pts)

#### Community composition exploration ####
# Calculate distances
vert_dissim <- vegdist(comm)
beet_dissim <- vegdist(commb)
vert_covs <- tibble(site = rownames(comm)) %>% left_join(allsites_cov_db) %>% 
  select(site, alt1, lc2, temp_mean1, ghm1, lfdistance_2) %>% 
  mutate(lc = factor(lc2, labels = c("Dense","Open")),
         alt = scale(alt1),
         temp = scale(temp_mean1),
         ghm = scale(ghm1),
         fdist = scale(log1p(lfdistance_2))) %>% 
  select(lc,alt,temp,ghm,fdist)
  
beet_covs <- tibble(site = rownames(commb)) %>% left_join(allsites_cov_db) %>% 
  select(site, alt1, lc2, temp_mean1, ghm1, lfdistance_2) %>% 
  mutate(lc = factor(lc2, labels = c("Dense","Open")),
         alt = scale(alt1),
         temp = scale(temp_mean1),
         ghm = scale(ghm1),
         fdist = scale(log1p(lfdistance_2))) %>% 
  select(lc,alt,temp,ghm,fdist)
# Create and plot NMDS
vert_NMDS <- metaMDS(vert_dissim)
beet_NMDS <- metaMDS(beet_dissim)

# PCoA instead
vert_pcoa <- pcoa(vert_dissim)
beet_pcoa <- pcoa(beet_dissim)

# Fit environ variables
vert_envfit <- envfit(vert_NMDS, vert_covs)
beet_envfit <- envfit(beet_NMDS, beet_covs)

vert_envfit_pcoa <- envfit(vert_pcoa$vectors, vert_covs)
beet_envfit_pcoa <- envfit(beet_pcoa$vectors, beet_covs)

# GAM plots (ordisurf), showing contour plots of how the covariates are
# distributed among sites
par(mfrow = c(2,2))
ordisurf(vert_NMDS, vert_covs$alt, select=T)
ordisurf(vert_NMDS, vert_covs$ghm, select=T)
ordisurf(vert_NMDS, vert_covs$temp, select=T)
ordisurf(vert_NMDS, vert_covs$fdist, select=T)
par(mfrow=c(1,1))

# PERMANOVA to see how the different variables affect the dissimilarity across
# sites
vert_permanova <- adonis2(vert_dissim~lc+alt+ghm+temp+fdist, data = vert_covs)
beet_permanova <- adonis2(beet_dissim~lc+alt+ghm+temp+fdist, data = beet_covs)
vert_permanova
beet_permanova

# These results suggest that differences in species composition for vertebrates
# are influenced by land cover, altitude, ghm, temperature and distance to large
# forest patches. Results are similar to those for beetles, except for
# temperature, which did not have a strong effect on beetle community
# composition. Altitude and land cover have the strongest effects for both

# Visualize using NMDS
plot(vert_NMDS, type='n')
ordiellipse(vert_NMDS, vert_covs$lc, draw = "polygon", col = c("darkgreen", "goldenrod"))
points(vert_NMDS, display = "sites",select = vert_covs$lc=="Dense",col = "darkgreen", pch = 16)
points(vert_NMDS, display = "sites",select = vert_covs$lc=="Open",col = "goldenrod", pch = 16)
plot(vert_envfit, labels = list(factors = c("Dense", "Open")), bg='white', col = "gray10")
# text(vert_NMDS, display = "species", cex=0.7)


# test if there is more variability in open vs dense forest sites for vertebrates
anova(betadisper(vert_dissim, vert_covs$lc))
plot(betadisper(vert_dissim, vert_covs$lc), hull = F, ellipse = T, col = c("darkgreen", "goldenrod"), segments = F)
# There seems to be greater distance to the centroid in open forest sites
# (betadisper: F_1,115 = 5.13, p = 0.025). It is important to keep this in mind
# as it may be driving apparent differences due to other covariates

plot(beet_NMDS, type='n')
ordiellipse(beet_NMDS, beet_covs$lc, draw = "polygon", col = c("darkgreen", "goldenrod"))
points(beet_NMDS, display = "sites",select = beet_covs$lc=="Dense",col = "darkgreen", pch = 16)
points(beet_NMDS, display = "sites",select = beet_covs$lc=="Open",col = "goldenrod", pch = 16)
plot(beet_envfit, labels = list(factors = c("Dense", "Open")), bg='white', col = "gray10")

anova(betadisper(beet_dissim, beet_covs$lc))
plot(betadisper(beet_dissim, beet_covs$lc), 
     hull = F, ellipse=T, col = c("darkgreen", "goldenrod"), segments = F)
# For beetles there is no difference in variability. F_1,70 = 0.20, p = 0.65.

#### Richness and diversity ####
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

plot(diversity(commb)~beet_covs$alt)
plot(diversity(commb)~beet_covs$temp)
plot(diversity(commb)~beet_covs$ghm)
plot(diversity(commb)~beet_covs$fdist)

#### Diversity models ####
vert_gam_db <- data.frame(s = diversity(comm), n = specnumber(comm), vert_covs)
beet_gam_db <- data.frame(s = diversity(commb), n = specnumber(commb), beet_covs)
plot(vert_gam_db)
plot(beet_gam_db)

# Sites with higer richness also have higher diversity for both verts and
# beetles
vert_rich_lm <- gam(n~lc+ghm+temp+alt+fdist, family = poisson, data=vert_gam_db)
vert_rich_gam <- gam(n~lc+ghm+s(temp)+s(alt)+s(fdist), 
                     family = poisson, data = vert_gam_db)
summary(vert_rich_gam)
# reduced model with only 
vert_rich_gam2 <- update(vert_rich_gam, .~.-lc-ghm-s(temp))
# The model that has the highest explained variance (39%, R^2 = 0.28) includes
# smooth terms for temperature, distance to forest, and altitude. However only
# temperature and distance to forest are statistically meaningful. temp: F_7.09
# = 1.604, p = 0.042, fdist: F_3.47 = 1.19, p = 0.0133; alt: F_4.79 = 1.02, p =
# 0.081. The linear effects (disturbance, )
summary(vert_rich_gam)
gam.check(vert_rich_gam)
# This model seems to be overestimating richness at the low predicted values.
plot(vert_rich_gam)

# compare models
anova(vert_rich_gam, vert_rich_lm, vert_rich_gam2, test = "Chisq")
#Comparing the models shows the nonlinear model with all covariates is
#significantly better than the linear model and the reduced gam

## Diversity model
vert_div_lm <- gam(s~lc+ghm+temp+fdist+alt, data = vert_gam_db)
vert_div_gam <- gam(s~lc+ghm+temp+s(fdist)+alt, data = vert_gam_db)
summary(vert_div_gam)
plot(vert_div_gam)
# The GAM for diversity is even poorer. The better performing model includes
# only linear effects, for fdist, and alt.
anova(vert_div_lm, vert_div_gam, test = "F")
# In fact the linear model provides a better fit. Furthermore, the model with
# only fdist and alt performs similarly to the one with all covariates
vert_div_lm_red <- gam(s~fdist+alt, data = vert_gam_db)
summary(vert_div_lm_red)
# This model has significant negative effect of log distance to forest, and a
# positive effect of altitude. Sites at higher altitude would thus have higher
# vertebrate diversity, and sites further away from forest would have lower
# diversity.

## Beetles models
# Richness linear model
beet_rich_lm <- gam(n~lc+ghm+temp+alt+fdist, data = beet_gam_db)
summary(beet_rich_lm) # no sig effects
# richness GAM
beet_rich_gam <- gam(n~lc+s(ghm)+s(temp)+s(alt)+s(fdist), 
                     family = poisson, data = beet_gam_db)
summary(beet_rich_gam)
beet_rich_gam2 <- update(beet_rich_gam, n~s(alt)+s(fdist))

# compare vs linear model
anova(beet_rich_lm,beet_rich_gam, test = "Chisq") # the nonlinear model fits the data much better
# I have tested most combinations, and the best fit is provided by the model
# that has smooth effects of altitude and distance to forest, and a linear effect of temperature
summary(beet_rich_gam2) # sig effect for alt: X2 = 13.10, p = 0.02. R2 = 0.14
plot(beet_rich_gam2)
gam.check(beet_rich_gam2)
# the predicted values are more spread out than would be expected

beet_div_lm <- gam(s~lc+ghm+alt+temp+fdist, data=beet_gam_db)
summary(beet_div_lm) # no effects of anything
beet_div_gam <- gam(s~lc+s(ghm)+s(alt)+s(temp)+s(fdist), data=beet_gam_db) # in the full model the ghm and temp are both reduced to linear effects. There is also no sig effect of ghm or lc.
beet_div_gam <- gam(s~alt, data=beet_gam_db)
summary(beet_div_gam)# no covars
plot(beet_div_gam)

#' The model includes only a smooth effect of altitude, but this is not
#significantly different from including just a linear positive effect of altitude.
# Overall it's hard to tell any real effect of these covariates on richness or diversity

######Things to check ######
