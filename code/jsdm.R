library(terra)
library(sjSDM)
library(tidyverse)

camvect <- vect("Data/ct_covars.geojson")
campts_db <- as.data.frame(camvect) %>% mutate(site = gsub("#", "", deployment_id), days = as.numeric(days))

data_raw <- read_csv("Data/ct_records.csv",
                     col_types = cols_only(deployment_id = 'c', common_name = 'c', timestamp = col_datetime("%Y-%m-%d %H:%M") ))
# include only cameras with more than 12 weeks (~3 months)
cams_include <- summarise(data_raw, .by = deployment_id, t=max(timestamp)-min(timestamp)) %>% 
  filter(t>12*7) %>% 
  mutate(site = sub("#","",deployment_id)) %>% 
  select(site) 
DAT <- read_csv("Data/records_wild.csv") %>% inner_join(cams_include)
DAT_wcovs <- left_join(DAT, campts_db)
# create sp matrix
occ_data_bin <- DAT %>% count(site, common_name) %>% 
  filter(!grepl("uid",common_name)) %>% # remove 'uid' species (felids, opossums, peccaries, raccoons)
  # mutate(n = (n>0)) %>% 
  pivot_wider(names_from = common_name, values_from = n, values_fill = 0) %>% 
  column_to_rownames("site") %>% 
  as.matrix()
# create env matrix
env_mat <- DAT %>% distinct(site) %>% left_join(campts_db) %>% 
  select(site, easting, northing, ghm1, lc1,alt1, evi_mean1, lfdistance, distance) %>% 
  rename(x = easting, y = northing, landcov = lc1, alt = alt1, ghm = ghm1, pa_dist = distance, lgfor_dist = lfdistance) %>% 
  mutate(landcov = ifelse(landcov==12,1,0)) %>% 
  column_to_rownames("site") %>% 
  as.matrix()


occ_models <- list()

# fit a basic model with no biotic interactions or spatial autocorrelation.
occ_models[["base"]] <- sjSDM(Y = occ_data_bin, env = linear(data = scale(env_mat), formula = ~landcov+alt+ghm), # scales env. covariates: land cover (forest/grass), altitude, ghm, dist to pa
                               se = TRUE, family=binomial("probit"), device = "cpu")
SPeigen <- generateSpatialEV(env_mat[,c(1,2)])
occ_models[["spatEV"]] <- sjSDM(Y = occ_data_bin, env = linear(data = scale(env_mat), formula = ~landcov+alt+ghm), # scales env. covariates: land cover (forest/grass), altitude, ghm, dist to pa
                                spatial = linear(SPeigen, ~0+.),
                                se = TRUE, family=binomial("probit"), device = "cpu")
occ_models[["spatlin"]] <- sjSDM(Y = occ_data_bin, env = linear(data = scale(env_mat), formula = ~landcov+alt+ghm), # scales env. covariates: land cover (forest/grass), altitude, ghm, dist to pa
                                spatial = linear(env_mat, ~0+poly(x,y,degree=2)),
                                se = TRUE, family=binomial("probit"), device = "cpu")

anovas <- lapply(occ_models, anova)
anovas
# The ANOVA table for the models show that, for the base model, the abiotic
# component has higher deviance than the biotic component, explaining 19% of the
# explained variance, while biotic associations explain 12%. The full variance
# explained is 36%.
lapply(occ_models, logLik)
# The highest log-likelihood is obtained with the base model, which has no
# spatial component.

ggobj <- plot(occ_models$base)
ggobj+theme_classic(base_size = 16)
# This model clearly shows the influence of altitude, species like tinamous,
# tamanduas, armadillos, curassows, raccoons, and agoutis were associated with
# lowlands, while oncillas, rabbits, jaguars, guans, and brocket were associated with
# highlands.


#### Beetles ####
# probably need to remove rare species that appeared in a single site.
DATb <- read.csv("Data/dung_beetles_prc.csv") %>% 
  # join traps from same station
  summarise(.by = Plot, across(where(is.numeric), sum)) %>% 
  inner_join(select(campts_db, site), by = join_by("Plot"=="site"))
# create sp matrix
occ_data_binb <- DATb %>% 
  column_to_rownames("Plot") %>% 
  select(where(\(x) sum(x)>0)) %>% 
  filter(rowSums(.)>0) %>% 
  as.matrix()
# create env matrix
env_matb <- tibble(site = rownames(occ_data_binb)) %>% left_join(campts_db) %>% 
  select(site, easting, northing, ghm1, lc1,alt1, evi_mean1, lfdistance, distance) %>% 
  rename(x = easting, y = northing, landcov = lc1, alt = alt1, ghm = ghm1, pa_dist = distance, lgfor_dist = lfdistance) %>% 
  mutate(landcov = ifelse(landcov==12,1,0)) %>% 
  column_to_rownames("site") %>% 
  as.matrix()

occ_modelsb <- list()

# fit a basic model with no biotic interactions or spatial autocorrelation.
occ_modelsb[["base"]] <- sjSDM(Y = occ_data_binb, env = linear(data = scale(env_matb), formula = ~landcov+alt+ghm), # scales env. covariates: land cover (forest/grass), altitude, ghm, dist to pa
                              se = TRUE, family=binomial("probit"), device = "cpu")
SPeigenb <- generateSpatialEV(env_matb[,c(1,2)])
occ_modelsb[["spatEV"]] <- sjSDM(Y = occ_data_binb, env = linear(data = scale(env_matb), formula = ~landcov+alt+ghm), # scales env. covariates: land cover (forest/grass), altitude, ghm, dist to pa
                                spatial = linear(SPeigenb, ~0+.),
                                se = TRUE, family=binomial("probit"), device = "cpu")
occ_modelsb[["spatlin"]] <- sjSDM(Y = occ_data_binb, env = linear(data = scale(env_matb), formula = ~landcov+alt+ghm), # scales env. covariates: land cover (forest/grass), altitude, ghm, dist to pa
                                 spatial = linear(env_matb, ~0+poly(x,y,degree=2)),
                                 se = TRUE, family=binomial("probit"), device = "cpu")

anovasb <- lapply(occ_modelsb, anova)
anovasb
# plot correlations
sp_hclust_order <- hclust(dist(cov2cor(getCov(occ_modelsb$base))))$order
corrmat <- matrix(nrow = length(sp_hclust_order), ncol = length(sp_hclust_order))
for (i in seq_along(sp_hclust_order)) {
  for (j in seq_along(sp_hclust_order)) {
    corrmat[i,j] <- cov2cor(getCov(occ_modelsb$base))[sp_hclust_order[i],sp_hclust_order[j]]
  }
}
corr_ggdf <- expand.grid(occ_modelsb$base$species[sp_hclust_order], occ_modelsb$base$species[sp_hclust_order])
for (i in seq_along(corrmat)) {
  corr_ggdf[i,3] <- corrmat[i]
}

names(corr_ggdf) <- c("sp1", "sp2", "corr")
ggplot(corr_ggdf, aes(sp1,sp2))+
  geom_tile(aes(fill = corr))+
  scale_fill_gradientn(limits=c(-1,1),
                       colours = col2(200))+
  labs(x='',y='', fill='')+
  theme(axis.text.y=element_text(hjust = 0.5), 
        axis.text.x=element_text(angle = 45, hjust = 0),
        legend.key.height = unit(3, 'lines'))+
  scale_x_discrete(position = 'top', limits = rev)

#### Plots ####
# plot
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))

bestmod <- occ_models$base
bestmod_summary <- summary(bestmod)

bestmodcoefs <- as.data.frame(t(bestmod_summary$coefs[-1,])) %>% rownames_to_column("sp") %>% 
  pivot_longer(2:4, names_to = "coef", values_to = "value")

sigLabel <- as.data.frame(t((bestmod_summary$P<0.05)*bestmod_summary$coefs)[,-1])%>%   rownames_to_column("sp") %>% 
  mutate(across(2:4, \(x) case_when(x>0~"+",
                                    x<0~"-"))) %>% 
  pivot_longer(cols = 2:4, names_to = "coef", values_to = "text")


# arrange species by their dissimilarity for correlation matrix
sp_hclust_order <- hclust(dist(cov2cor(getCov(bestmod))))$order

ggplot(bestmodcoefs, aes(x=coef, y=sp))+
  geom_tile(aes(fill = value), show.legend = F)+
  scale_fill_gradientn(limits=c(-2,2),
                       colours = col2(200))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 0), 
        # legend.position = 'none',
        axis.title = element_blank())+
  scale_x_discrete(position = 'top', labels = c("Altitude", "Human pressure", "Land cover"))+
  geom_text(data=left_join(bestmodcoefs,sigLabel), aes(label=text), vjust=0.5) -> jsdm_coef_plot
jsdm_coef_plot


# # import common names database
# common <- read.csv("Data/species_commonName.csv")
# common_hclust <- data.frame(common, order = sp_hclust_order)

corrmat <- matrix(nrow = length(sp_hclust_order), ncol = length(sp_hclust_order))
for (i in seq_along(sp_hclust_order)) {
  for (j in seq_along(sp_hclust_order)) {
    corrmat[i,j] <- cov2cor(getCov(bestmod))[sp_hclust_order[i],sp_hclust_order[j]]
  }
}

# colnames(corrmat) <- rownames(corrmat) <- 
  
#   # base R plot
#   par(mar = c(2,15,10,5), las=2)
# plot(seq_along(sp_hclust_order), seq_along(sp_hclust_order), 
#      type = "n", 
#      xaxt = "n", yaxt = "n",bty="n",
#      xlab = "", ylab = "",asp=1)
# image(1:31,1:31,corrmat[ncol(corrmat):1,], add = T, col = hcl.colors(100,"Blue-Red 3", rev=T), zlim=c(-1,1))
# axis(2, at = 1:31, labels = bestmod$species[sp_hclust_order])
# axis(3, at = 31:1, labels = bestmod$species[sp_hclust_order])

## ggplot
corr_ggdf <- expand.grid(bestmod$species[sp_hclust_order], bestmod$species[sp_hclust_order])
for (i in seq_along(corrmat)) {
  corr_ggdf[i,3] <- corrmat[i]
}

names(corr_ggdf) <- c("sp1", "sp2", "corr")
ggplot(corr_ggdf, aes(sp1,sp2))+
  geom_tile(aes(fill = corr))+
  scale_fill_gradientn(limits=c(-1,1),
                       colours = col2(200))+
  labs(x='',y='', fill='')+
  theme(axis.text.y=element_text(hjust = 0.5), 
        axis.text.x=element_text(angle = 45, hjust = 0),
        legend.key.height = unit(3, 'lines'))+
  scale_x_discrete(position = 'top', limits = rev) -> jsdm_corr_plot
jsdm_corr_plot+theme(aspect.ratio = 1)


cowplot::plot_grid(jsdm_coef_plot, jsdm_corr_plot, align = 'h', rel_widths = c(1,3))
# need to check species name order

#### spAbundance approach ####
library(spOccupancy)
# We fit a model similar to the imperfect detection one implemented in Tobler et
# el. 2019, as implemented by the sfMsPGOcc function in the spOccupancy package
# This is an occupancy model with a detection component, so we need repeat
# occasions at each site. We exclude the sites with less than 4 weeks of
# sampling, and divide the detections into monthly counts.

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

# model formulas
occform <- ~1+scale(forest)+scale(alt)+scale(ghm)# formula for the abundance portion of the model
detform <- ~1 

# Run model
sfmod <- sfMsPGOcc(occ.formula = ~scale(forest)+scale(alt)+scale(ghm), 
                   det.formula = ~1,
                  data = moddata, 
                  inits = inits.list, 
                  # priors = prior.list, 
                  n.factors = 2, # number of latent factors
                  n.batch = 200, batch.length = 25, # total n of samples is n.batch*batch.length
                  n.chains = 3, n.burn = 3000, n.thin = 2
                  )
# this model hasn't run. Below is a non-spatial version
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
# The model's Bayesian p value is 0.2364, so not terrible, but rather low (0.5
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
ppc.df <- data.frame(fit = ppc.out3$fit.y, 
                     fit.rep = ppc.out3$fit.y.rep,
                     color = 'lightskyblue')
ppc.df$color[ppc.df$fit.rep>ppc.df$fit] <- 'lightsalmon'
plot(ppc.df$fit, ppc.df$fit.rep, bg=ppc.df$color, pch = 21, ylab = "Fit", xlab= ' True')
lines(ppc.df$fit, ppc.df$fit, col = "black")

# Plot
betas <- lfmod$beta.samples
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

