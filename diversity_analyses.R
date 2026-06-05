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
rownames(comm) <- gsub("#", "", rownames(comm))
# Import pitfall trap data
commb <- read.csv("Data/beet_comm_mat.csv", row.names = 1)
rownames(commb) <- gsub("#", "", rownames(commb))
# Import station data #
allsites_pts <- terra::vect("Data/combined_survey_pts.geojson")
allsites_cov_db <- as.data.frame(allsites_pts)

### Import temp data
temp_data_day <- read.csv("Data/LST_2020_allsites_day.csv",
                          col.names = c("date", "LST", "site")) %>%
  mutate(date = mdy(date), month = month(date))
# summarize by station
daytemp_sum <- summarise(temp_data_day, max_T = max(LST), min_T = min(LST),
                         mean_T = mean(LST), sd_T = sd(LST), .by = site)

#### Community composition #### 
# Standardize matrices to account for different
#effort across sites. We divide by effort
effort <- tibble(site = rownames(comm)) %>% 
  left_join(read_csv("Data/ct_deployments.csv"),
            by = join_by(site == deployment_id)) %>% 
  pull(days) %>% as.numeric()

comm_std <- comm / effort

# Subset stations with available temp data
comm_std <- comm_std[rownames(comm_std) %in% daytemp_sum$site,]
commb <- commb[rownames(commb) %in% daytemp_sum$site, ]

##### Calculate disimilarities ####
# distance with the Jaccard turnover index
vert_dissim_jac <- vegdist(comm_std, "jaccard", binary = TRUE)
beet_dissim_jac <- vegdist(commb, "jaccard", binary = TRUE)


# distances with the Simpson nestedness index
vert_dissim_sim <- designdist(comm_std, method = "J / pmin(A, B)",
                            terms = "binary", name = "Simpson")
beet_dissim_sim <- designdist(commb, method = "J / pmin(A, B)",
                            terms = "binary", name = "Simpson")

# covariates
vert_covs_raw <- tibble(site = rownames(comm_std)) %>%
  left_join(allsites_cov_db) %>%
  inner_join(daytemp_sum) %>%
  select(site, easting, northing, 
         alt1, ghm1, lfdistance_2, 
         mean_T, treecov_gfw_100) 
# standardised
vert_covs <- mutate(vert_covs_raw,
                    alt = standardise(alt1),
                    alt2 = standardise(alt1)^2,
                    ghm = ghm1,
                    fdist = log1p(lfdistance_2),
                    tmean = standardise(mean_T),
                    treecov = standardise(treecov_gfw_100)) %>%
  select(alt, alt2, ghm, fdist, tmean,
         treecov, easting, northing)
# could also standardize using decostand
beet_covs_raw <- tibble(site = rownames(commb)) %>%
  left_join(allsites_cov_db) %>%
  inner_join(daytemp_sum) %>%
  select(site, easting, northing,
         alt1, ghm1, lfdistance_2,
         mean_T, treecov_gfw_100)
beet_covs <- mutate(beet_covs_raw, 
                    alt = standardise(alt1),
                    alt2 = standardise(alt1)^2,
                    ghm = ghm1,
                    fdist = log1p(lfdistance_2),
                    tmean = standardise(mean_T),
                    treecov = treecov_gfw_100) %>%
  select(easting, northing, alt, alt2, ghm, fdist,
         tmean, treecov)

##### PERMANOVA #####

# Analysis to see how the different variables affect the dissimilarity across
# sites
permanova_comp <- fit_models(make_models(vars = c("alt", "alt2","ghm", "fdist", "tmean", "tsd", "treecov"),
                                         k = 5), 
                             com_data = comm_std, env_data = vert_covs_raw)
select_models(permanova_comp)
# We compared different covariate combinations, and the one with the lowest AICc
# included land cover, altitude, and alt^2. However, there are 12 other models
# that have similar AICc. The VIF for the best model is not that high (2.35) so
# we could ignore it. The second best model includes mean temp instead of
# altitude.

# doing the same but removing altitude
all_perm_mods <- make_models(vars = c("ghm", "fdist", "tmean", "treecov"),
                             ncores = 4,
                             k = 4)
permanova_comp2 <- fit_models(all_perm_mods, 
                             com_data = comm_std,
                             env_data = vert_covs)
select_models(permanova_comp2)
# we obtain AICc values ~ -139, with the best model including ghm, distance to
# forest, and mean temperature, with a VIF of 1.82

permanova_comp_beet <- fit_models(make_models(vars = c("ghm", "fdist", "tmean", "treecov")), 
                             com_data = commb, env_data = beet_covs)
select_models(permanova_comp_beet)
# In the case of beetles, there are 18 models that rank similarly. The top one
# included forest distance, and mean temp, and has a VIF of 1.21


### Same but using Jaccard binary index
# doing the same but removing some covariates (lc and altitude)
permanova_comp_jac <- fit_models(make_models(vars = c("ghm", "fdist",
                                                      "tmean", "treecov"),
                                          k = 4), 
                              com_data = vert_dissim_jac, env_data = vert_covs)

permanova_comp_beet_jac <- fit_models(make_models(vars = c("ghm", "fdist",
                                                           "tmean", "treecov"),
                                                  k = 4), 
                                  com_data = beet_dissim_jac, env_data = beet_covs)

### Same but using Simpson binary index
permanova_comp_sim <- fit_models(make_models(vars = c("ghm", "fdist", "tmean",
                                                      "treecov"),
                                             k = 4), 
                                 com_data = vert_dissim_sim, env_data = vert_covs)

permanova_comp_beet_sim <- fit_models(make_models(vars = c("ghm", "fdist", "tmean",
                                                           "treecov"),
                                                  k = 4), 
                                      com_data = beet_dissim_sim, env_data = beet_covs)
select_models(permanova_comp_jac)
select_models(permanova_comp_beet_jac)
select_models(permanova_comp_sim)
select_models(permanova_comp_beet_sim)


## Run PERMANOVA
vert_permanova <- adonis2(comm_std ~ fdist + tmean, data = vert_covs, method = 'bray', by = 'margin')
beet_permanova <- adonis2(commb ~ fdist + tmean, data = beet_covs, method = 'bray', by = 'margin')
vert_permanova
beet_permanova

# using jaccard index of turnover
vert_permanova_jac <- adonis2(vert_dissim_jac ~ fdist + tmean, data = vert_covs, by = 'margin')
beet_permanova_jac <- adonis2(beet_dissim_jac ~ fdist + tmean, data = beet_covs, by = 'margin')
vert_permanova_jac
beet_permanova_jac

# These results suggest that differences in species turnover for vertebrates are
# influenced mostly by differences in temperature and distance to forest, with
# some smaller effect of disturbance. These variables however explain only a
# small proportion of the variance, especially for beetles. The covariates do
# not seem to affect nestedness at all.


#####------ PCoA #####
# Run PCoA for visualization of PERMANOVA
vert_pcoa <- pcoa(vert_dissim)
beet_pcoa <- pcoa(beet_dissim)

vert_pcoa_jac <- pcoa(vert_dissim_jac)
beet_pcoa_jac <- pcoa(beet_dissim_jac)

vert_pcoa_sim <- pcoa(vert_dissim_sim)
beet_pcoa_sim <- pcoa(beet_dissim_sim)

# Fit environ variables
vert_envfit_pcoa <- envfit(vert_pcoa$vectors ~ tmean + fdist, data = vert_covs)
beet_envfit_pcoa <- envfit(beet_pcoa$vectors ~ tmean + fdist, data = beet_covs)

vert_envfit_pcoa_jac <- envfit(vert_pcoa_jac$vectors ~ tmean + fdist, data = vert_covs)
beet_envfit_pcoa_jac <- envfit(beet_pcoa_jac$vectors ~ tmean + fdist, data = beet_covs)

# scores
scores(vert_envfit_pcoa_jac, display = "vectors")
# Plot PCoA
taxa_cols <- hcl.colors(2, "Geyser", rev = TRUE)
# vertebrate PcOA
plot(vert_pcoa$vectors[, 1:2], las = 1, asp = 1,
     xlab = "Axis 1", ylab = "Axis 2", type = "n",
     main = "Vertebrates")
points(vert_pcoa$vectors, pch = 21, bg = 'turquoise')
ordisurf(vert_pcoa$vectors, vert_covs_raw$mean_T-273.15,
         add = T, col = "gray50")
plot(vert_envfit_pcoa,
     labels = list(vectors = c("Temp", "F dist", "gHM")),
     bg = rgb(1, 1, 1, 0.5), col = "gray30", p.max = 0.05)
# plot(vert_envfit_pcoa, labels = list(factors = c("Dense", "Open"),
# vectors = c("Alt", "Alt^2", "gHM", "Forest dist")),
# bg='white', col = "gray10")

# PCoA for beetles
plot(beet_pcoa$vectors[, 1:2],
     type = "n", las = 1, asp = 1,
     xlab = "Axis 1", ylab = "Axis 2", main = "Beetles",
     pty = "s")
# ordiellipse(beet_pcoa$vectors, beet_covs$lc, draw = "polygon", col = c("darkgreen", "goldenrod"))
# points(beet_pcoa$vectors[beet_covs$lc=="Dense",1:2], pch = 16, col = "darkgreen")
# points(beet_pcoa$vectors[beet_covs$lc=="Open",1:2], pch = 17, col = "goldenrod")
points(beet_pcoa$vectors, pch = 21, bg = "magenta")
ordisurf(beet_pcoa$vectors, beet_covs_raw$mean_T-273.15,
         add = T, col = "gray50")
plot(beet_envfit_pcoa,
     labels = list(vectors = c("Temp", "F dist")),
     bg = rgb(1, 1, 1, 0.5), col = "gray30")

# Plot with result of Jaccard dissimilarity
png("figures/PCoA.png", width = 8, height = 4, units = "in", res = 300)
par(mfrow = c(1,2), mar = c(3, 3, 1, 1))
###
plot(vert_pcoa_jac$vectors[, 1:2], las = 1, asp = 1,
     xlab = "Axis 1", ylab = "Axis 2", type = "n")
points(vert_pcoa_jac$vectors, pch = 17, col = taxa_cols[2])
ordisurf(vert_pcoa_jac$vectors, vert_covs_raw$lfdistance_2,
         add = T, col = "gray50")
plot(vert_envfit_pcoa_jac,
     labels = list(vectors = c("Temp", "Dist")),
     bg = rgb(1, 1, 1, 0.5), col = "gray30", p.max = 0.05, arrow.mul = 0.4)
mtext("a", side = 2, line = 2, padj = -12,  las = 2, font = 2)

# PCoA for beetles
plot(beet_pcoa_jac$vectors[, 1:2],
     type = "n", las = 1, asp = 1,
     xlab = "Axis 1", ylab = "Axis 2",
     pty = "s")
points(beet_pcoa_jac$vectors, pch = 16, col = taxa_cols[1])
ordisurf(beet_pcoa_jac$vectors, beet_covs_raw$mean_T,
         add = T, col = "gray50")
plot(beet_envfit_pcoa_jac,
     labels = list(vectors = c("Temp", "Dist")),
     bg = rgb(1, 1, 1, 0.5), col = "gray30", arrow.mul = 0.6)
mtext("b", side = 2, line = 2, padj = -12, las = 2, font = 2)
dev.off()
#### Diversity and richness ####

##### ggplots #####
# richness vs alt.


# richness vs T
p1 <- data.frame(n = specnumber(comm), vert_covs_raw) %>%
  ggplot(aes(mean_T - 273.15, n)) +
  geom_point(shape = 16, size = 2, color = rgb(0, 0, 0, 0.5))+
  labs(x = expression(paste("Mean temperature (",degree*C, ")")),
       y = "Species richness") +
  theme_pubr()

p5 <- data.frame(n = specnumber(commb),beet_covs_raw) %>%
  ggplot(aes(mean_T - 273.15, n)) +
  geom_point(shape = 16, size = 2, color = rgb(0, 0, 0, 0.5)) +
  labs(x = expression(paste("Mean temperature (",degree*C, ")")),
       y = "Species richness") +
  theme_pubr()

# diversity vs T
p2 <- data.frame(n = diversity(comm), vert_covs_raw) %>%
  ggplot(aes(mean_T - 273.15, n)) +
  geom_point(shape = 16, size = 2, color = rgb(0, 0, 0, 0.5)) +
  labs(x = expression(paste("Mean temperature (",degree*C, ")")),
  y = "Diversity") +
  theme_pubr()


p6 <- data.frame(n = diversity(commb), beet_covs_raw) %>%
  ggplot(aes(mean_T - 273.15, n))+
  geom_point(shape = 16, size = 2, color = rgb(0, 0, 0, 0.5)) +
  labs(x = expression(paste("Mean temperature (", degree*C, ")")),
       y = "Diversity") +
  theme_pubr()

# new plots
pvertrich2 <- data.frame(n = specnumber(comm),vert_covs_raw) %>%
  ggplot(aes(alt1, n)) +
  geom_point(aes(color = lfdistance_2 + 1), show.legend = FALSE) +
  labs(x = "Altitude (masl)", y = "Species richness",
       color = "Distance to forest") +
  theme_pubr() +
  theme(palette.color.continuous = hcl.colors(3, "Greens")) +
  coord_cartesian(xlim = c(0, 2550), ylim = c(0, 16))
pvertrich2

pbeetrich2 <- data.frame(n = specnumber(commb), beet_covs_raw) %>%
  ggplot(aes(alt1, n)) +
  geom_point(aes(color = lfdistance_2 + 1)) +
  labs(x = "Altitude (masl)", y = "Species richness",
       color = "Distance to forest") +
  theme_pubr() +
  theme(palette.color.continuous = hcl.colors(3, "Greens")) +
  coord_cartesian(xlim = c(0, 2550), ylim = c(0, 16))
pbeetrich2

# Join
ggarrange(pvertrich2, pbeetrich2, 
          labels = "auto",
          legend.grob = get_legend(pbeetrich2, position = "top"), legend = "top")
ggsave("figures/richness_vs_alt_fdist.png", dpi = 300, width = 7, height = 4)



#### Richness models ####
##### Vertebrates #####
vert_gam_db <- data.frame(n = specnumber(comm), vert_covs_raw)

# GAM of richness vs alt
gam_rich_vert <- gam(n ~ s(alt1, bs = "cr", k = 6),
                     family = poisson, data = vert_gam_db, method = "REML")
gam_rich_vert_null <- update(gam_rich_vert, . ~ 1) # intercept-only
gam_rich_vert_null2 <- update(gam_rich_vert, . ~ alt1) # linear fixed effect of altitude
appraise(gam_rich_vert)

AIC(gam_rich_vert, gam_rich_vert_null, gam_rich_vert_null2)
# The lowest AIC is given by the model with the smooth term
anova(gam_rich_vert, gam_rich_vert_null)
anova(gam_rich_vert, gam_rich_vert_null2)
# Both null models perform worse and are significantly different.

##### Beetles #####
beet_gam_db <- data.frame(n = specnumber(commb), beet_covs_raw)

# GAM of richness vs alt
gam_rich_beet <- gam(n ~ s(alt1, bs = "cr", k = 6),
                     family = poisson, data = beet_gam_db, method = "REML")
gam_rich_beet_null <- update(gam_rich_beet, . ~ 1) # intercept-only
gam_rich_beet_null2 <- update(gam_rich_beet, . ~ alt1) # linear fixed effect of altitude

AIC(gam_rich_beet, gam_rich_beet_null, gam_rich_beet_null2)
# In the case of beetles the intercept-only model had the lowest AIC, but they
# were all comparable (dAIC < 2)
anova(gam_rich_beet_null, gam_rich_beet)
anova(gam_rich_beet_null, gam_rich_beet_null2)

##### Prediction plots ####

# Vert GAM predicted richness vs alt
pred_alts_vert <- seq(min(vert_covs_raw$alt1),
                     max(vert_covs_raw$alt1), length.out = 100
)
newdata_vert <- expand.grid(intercept = 1,
                           alt1 = pred_alts_vert
)
pred_vert_rich <- predict.gam(gam_rich_vert,
                              newdata = newdata_vert, 
                              type = 'response',
                              se.fit = TRUE
)
vert_rich_pred_db <- data.frame(newdata_raw,
                                n = pred_vert_rich$fit,
                                se = pred_vert_rich$se.fit
)
# beetles
pred_alts_beet <- seq(min(beet_covs_raw$alt1),
                     max(beet_covs_raw$alt1), length.out = 100
)
newdata_beet<- expand.grid(intercept = 1, 
                           alt1 = pred_alts_beet
)

pred_beet_rich <- predict.gam(gam_rich_beet_null,
                              newdata = newdata_beet, 
                              type = 'response', se.fit = T
)
beet_rich_pred_db <- data.frame(newdata_beet,
                                n = pred_beet_rich$fit, se = pred_beet_rich$se.fit
)
## single plot with vert and beet richness vs alt
rbind(data.frame(n = specnumber(comm), vert_covs_raw, sp = "vert"),
      data.frame(n = specnumber(commb), beet_covs_raw, sp = "beet")) %>% 
  ggplot() +
  geom_point(aes(alt1, n, shape = sp, color = sp), show.legend = FALSE) +
  labs(x = "Altitude (masl)", y = "Species richness") +
  theme_pubr(base_size = 16) +
  theme(palette.color.discrete = taxa_cols)

# same but with predictions from GAM
rbind(data.frame(n = specnumber(comm), vert_covs_raw, sp = "vert"),
      data.frame(n = specnumber(commb), beet_covs_raw, sp = "beet")) %>% 
  ggplot() +
  geom_ribbon(aes(x = alt1, ymin = n - se, ymax = n + se),
              data = vert_rich_pred_db, fill = alpha(taxa_cols[2], 0.5)) +
  geom_line(aes(alt1, n),
            data = vert_rich_pred_db, color = taxa_cols[2], linetype = 2) +
  geom_ribbon(aes(x = alt1, ymin = n - se, ymax = n + se),
              data = beet_rich_pred_db, fill = alpha(taxa_cols[1], 0.5)) +
  geom_line(aes(alt1, n),
            data = beet_rich_pred_db, color = taxa_cols[1], linetype = 2) +
  geom_point(aes(alt1, n, shape = sp, color = sp), show.legend = FALSE) +
  labs(x = "Altitude (masl)", y = "Species richness") +
  theme_pubr(base_size = 16) +
  theme(palette.color.discrete = taxa_cols)

# export
ggsave("figures/richness_vs_alt.png", dpi = 300, width = 5, height = 4)


#### Vert diversity #####
# linear model
vert_div_lm <- gam(s ~ lc + ghm + fdist + alt, data = vert_gam_db, method = "REML")
summary(vert_div_lm)
check_collinearity(vert_div_lm) # moderate VIF for fdist (8.48) and lc (7.95). I will remove fdist
vert_div_lm2 <- update(vert_div_lm, . ~ . - fdist)
# GAM
vert_div_gam <- gam(s ~ lc + s(ghm) + s(fdist) + s(alt),
                    data = vert_gam_db, method = "REML")
summary(vert_div_gam)
# All the smooth effects are reduced to linear terms in the GAM. So I'll stick with lm
summary(vert_div_lm2) # sig lc (-0.414, t = -3.637, p = 0.0004) and alt (0.080, t = 1.989, p = 0.049).
vert_div_lm3 <- update(vert_div_lm2, .~. - ghm) # remove nonsig ghm
summary(vert_div_lm3)# sig fdist and alt. Coefs are nearly identical. lc (-0.408, t = -4.431, p<0.001) and alt (0.080, t = 1.996, p = 0.0483).

AIC(vert_div_lm, vert_div_lm2, vert_div_lm3) # lm3 has the lowest AIC
anova(vert_div_lm3, vert_div_lm2, test = "F") # lm2 and lm3 are not sig. diff.
# This model has significant linear effects of lc and alt. The model is quite
# poor, with R2 = 0.15. Sites at higher altitude would thus have higher
# vertebrate diversity, and sites in open forest sites would have lower
# diversity.

##### Beet diversity #####
beet_div_lm <- gam(s ~ lc + ghm + alt + fdist, data = beet_gam_db, method = "REML")
summary(beet_div_lm) # sig ghm and alt

beet_div_gam <- gam(s ~ lc + s(ghm) + s(alt) + s(fdist),
  data=beet_gam_db, method = "REML"
)
# in the full model the ghm and temp are both reduced to linear effects.
# There is also no sig effect of ghm or lc.
summary(beet_div_gam)
# All terms are reduced to linear effects.
check_collinearity(beet_div_lm) # High corr of fdist
beet_div_lm2 <- update(beet_div_lm, .~. - fdist)
summary(beet_div_lm2)
check_collinearity(beet_div_lm2) # Correlations are low.
# remove nonsig lc term
beet_div_lm3 <- update(beet_div_lm2, .~. - lc)
summary(beet_div_lm3) # ghm no longer sig
beet_div_lm4 <- update(beet_div_lm3, .~. - ghm)
summary(beet_div_lm4)
beet_div_lm5 <- update(beet_div_lm2, .~. - ghm)
summary(beet_div_lm5)
AIC(beet_div_lm, beet_div_lm2, beet_div_lm3, beet_div_lm4, beet_div_lm5)
anova(beet_div_lm3, beet_div_lm2, test = "F")
# all models have comparable AIC, the top two (2 and 3) are nearly identical. I
# will stick with model 2 (alt, ghm, lc).




##### Visualize diversity ####
# Vert richness
pred_alts_raw <- seq(min(vert_covs_raw$alt1),
  max(vert_covs_raw$alt1), length.out = 100
)
newdata_raw <- expand.grid(intercept = 1,
  lc = levels(vert_covs$lc), alt = pred_alts_raw
)
pred_alts <- seq(min(vert_covs$alt), max(vert_covs$alt), length.out = 100)
newdata <- expand.grid(intercept = 1, lc = levels(vert_covs$lc),
  alt = pred_alts
)
pred_vert_rich <- predict.gam(vert_rich_gam4,
  newdata = newdata, type = 'response', se.fit = TRUE
)
vert_rich_pred_db <- data.frame(newdata_raw,
  n = pred_vert_rich$fit, se = pred_vert_rich$se.fit
)
pvertrich <- data.frame(n = specnumber(comm),vert_covs_raw) %>%
  mutate(lc = factor(lc2, labels = c("Dense", "Open"))) %>%
  ggplot(aes(alt1, n)) +
  geom_ribbon(aes(x = alt, ymin = n-se, ymax = n+se, fill = lc),
    data = vert_rich_pred_db
  ) +
  geom_line(aes(alt, n, color = lc), vert_rich_pred_db) +
  geom_point(aes(color = lc, shape = lc)) +
  labs(x = "Altitude (masl)", y = "Species richness",
       color = "Land cover", fill = "Land cover", shape = "Land cover") +
  scale_fill_manual(values = alpha(c("darkgreen", "goldenrod"), alpha = 0.3)) +
  scale_color_manual(values = c("darkgreen", "goldenrod")) +
  theme_pubr(base_size = 14)

ggplot(vert_gam_db, aes(alt,n)) + geom_point()+
  geom_line(aes(y = fitted(vert_rich_gam4), col=lc))
# Beetle richness
pred_alts_raw <- seq(min(beet_covs_raw$alt1),
  max(beet_covs_raw$alt1), length.out = 100
)
newdata_raw <- expand.grid(intercept = 1, lc = levels(beet_covs$lc),
  alt = pred_alts_raw
)
pred_alts <- seq(min(beet_covs$alt), max(beet_covs$alt), length.out = 100)
newdata <- expand.grid(intercept = 1, lc = levels(beet_covs$lc),
  alt = pred_alts
)
pred_beet_rich <- predict.gam(beet_rich_gam3,
  newdata = newdata, type = 'response', se.fit = T
)
beet_rich_pred_db <- data.frame(newdata_raw,
  n = pred_beet_rich$fit, se = pred_beet_rich$se.fit
)
pbeetrich <- data.frame(n = specnumber(commb), beet_covs_raw) %>%
  mutate(lc = factor(lc2, labels = c("Dense", "Open"))) %>%
  ggplot(aes(alt1, n)) +
  geom_ribbon(aes(x = alt, ymin = n - se, ymax = n + se, fill = lc),
              data = beet_rich_pred_db) +
  geom_line(aes(alt, n, color = lc), beet_rich_pred_db) +
  geom_point(aes(color = lc, shape = lc)) +
  labs(x = "Altitude (masl)", y = "Species richness",
       color = "Land cover", fill = "Land cover", shape = "Land cover") +
  scale_fill_manual(values = alpha(c("darkgreen", "goldenrod"), alpha = 0.3)) +
  scale_color_manual(values = c("darkgreen", "goldenrod")) +
  theme_pubr(base_size = 14)
# Vert diversity
pred_alts_raw <- seq(min(vert_covs_raw$alt1), max(vert_covs_raw$alt1),
                     length.out = 100)
newdata_raw <- expand.grid(intercept = 1, lc = levels(vert_covs$lc),
                           alt = pred_alts_raw)
pred_alts <- seq(min(vert_covs$alt), max(vert_covs$alt), length.out = 100)
newdata <- expand.grid(intercept = 1, lc = levels(vert_covs$lc),
                       alt = pred_alts)
pred_vert_div <- predict(vert_div_lm3, newdata = newdata,
                         type = 'response', se.fit = T)
vert_div_pred_db <- data.frame(newdata_raw,
  n = pred_vert_div$fit, se = pred_vert_div$se.fit
)
pvertdiv <- data.frame(n = diversity(comm),vert_covs_raw) %>%
  mutate(lc = factor(lc2, labels = c("Dense", "Open"))) %>%
  ggplot(aes(alt1, n)) +
  geom_ribbon(aes(x = alt, ymin = n-se, ymax = n+se, fill = lc),
              data = vert_div_pred_db)+
  geom_line(aes(alt,n, color = lc), vert_div_pred_db) +
  geom_point(aes(color = lc, shape = lc)) +
  labs(x = "Altitude (masl)", y = "Diversity", color = "Land cover",
       fill = "Land cover", shape = "Land cover")+
  scale_fill_manual(values = alpha(c("darkgreen", "goldenrod"), alpha = 0.3)) +
  scale_color_manual(values = c("darkgreen", "goldenrod")) +
  theme_pubr(base_size = 14)

# Beet diversity
pred_alts_raw <- seq(min(beet_covs_raw$alt1), max(beet_covs_raw$alt1),
                     length.out = 100)
newdata_raw <- expand.grid(intercept = 1, lc = levels(beet_covs$lc),
                           alt = pred_alts_raw)
pred_alts <- seq(min(beet_covs$alt), max(beet_covs$alt), length.out = 100)
newdata <- expand.grid(intercept = 1, lc = levels(beet_covs$lc),
                       ghm = 0, alt = pred_alts)
pred_beet_div <- predict(beet_div_lm2, newdata = newdata,
                         type = 'response', se.fit = T)
beet_div_pred_db <- data.frame(newdata_raw,n = pred_beet_div$fit,
                               se = pred_beet_div$se.fit)
pbeetdiv <- data.frame(n = diversity(commb), beet_covs_raw) %>%
  mutate(lc = factor(lc2, labels = c("Dense", "Open"))) %>%
  ggplot(aes(alt1, n)) +
  geom_ribbon(aes(x = alt, ymin = n - se, ymax = n + se, fill = lc),
              data = beet_div_pred_db) +
  geom_line(aes(alt, n, color = lc), beet_div_pred_db) +
  geom_point(aes(color = lc, shape = lc)) +
  labs(x = "Altitude (masl)", y = "Diversity", color = "Land cover",
       fill = "Land cover", shape = "Land cover") +
  scale_fill_manual(values = alpha(c("darkgreen", "goldenrod"), alpha = 0.3)) +
  scale_color_manual(values = c("darkgreen", "goldenrod")) +
  theme_pubr(base_size = 14)

ggarrange(pvertrich, pbeetrich, pvertdiv, pbeetdiv, labels = 'auto',
          common.legend = T, align = 'hv')


####Things to check ####
# Maybe analyze associations between large herbivores and specific dung beetle species/richness?