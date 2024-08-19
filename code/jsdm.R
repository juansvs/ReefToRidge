library(terra)
library(sjSDM)
library(tidyverse)

camvect <- vect("Data/ct_covars.geojson")
campts_db <- as.data.frame(camvect) %>% mutate(site = gsub("#", "", deployment_id), days = as.numeric(days))

data_raw <- read_csv("Data/ct_records.csv",
                     col_types = cols_only(deployment_id = 'c', common_name = 'c', timestamp = col_datetime("%Y-%m-%d %H:%M") ))
# include only cameras with more than we weeks (~3 months)
cams_include <- summarise(data_raw, .by = deployment_id, t=max(timestamp)-min(timestamp)) %>% 
  filter(t>28*3) %>% 
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
  as.matrix()
# create env matrix
env_matb <- DATb %>% left_join(campts_db, by = join_by("Plot"=="site")) %>% 
  select(Plot, easting, northing, ghm1, lc1,alt1, evi_mean1, lfdistance, distance) %>% 
  rename(x = easting, y = northing, landcov = lc1, alt = alt1, ghm = ghm1, pa_dist = distance, lgfor_dist = lfdistance) %>% 
  mutate(landcov = ifelse(landcov==12,1,0)) %>% 
  column_to_rownames("Plot") %>% 
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
library(spAbundance)
# We fit a model similar to the imperfect detection one implemented in TObler et
# el. 2019, as implemented by the sfMsDS function in the spAbundance package
abundform <- ~landcov+alt+ghm+(1|species)# formula for the abundance portion of the model
detform <- ~(1|site)
