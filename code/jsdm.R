library(terra)
library(sjSDM)
library(tidyverse)

camvect <- vect("Data/ct_covars.geojson")
campts_db <- as.data.frame(camvect) %>% mutate(site = gsub("#", "", deployment_id), days = as.numeric(days))

DAT <- read_csv("Data/records_wild.csv")
DAT_wcovs <- left_join(DAT, campts_db)
# create sp matrix
occ_data_bin <- DAT %>% count(site, common_name) %>% 
  mutate(n = (n>0)) %>% 
  pivot_wider(names_from = common_name, values_from = n, values_fill = 0) %>% 
  column_to_rownames("site") %>% 
  as.matrix()
# create env matrix
env_mat <- DAT %>% distinct(site) %>% left_join(campts_db) %>% 
  select(site, easting, northing, ghm1, lc1,alt1, evi_mean1, lfdistance, distance) %>% 
  rename(x = easting, y = northing, landcov = lc1, alt = alt1, ghm = ghm1, pa_dist = distance, lgfor_dist = lfdistance) %>% 
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
# The ANOVA table for the models show that, for the base model, the abiotic
# component has higher deviance than the biotic component, explaining 18% of the
# explained variance, while biotic associations explain 12%. The full variance
# explained is the sum of the two, 33%.
plot(anovas[[1]])
# The highest log-likelihood is obtained with the base model, which has no
# spatial component.
plot(occ_models$base)
# This model clearly shows the influence of altitude, species like tinamous,
# tamanduas, armadillos, curassows, raccoons, and agoutis were associated with
# lowlands, while oncillas, jaguars, guans, and brocket were associated with
# highlands.



## Species correlations
spcor <- getCor(occ_models[[1]])
image(spcor, asp=1, col = hcl.colors(12, 'Blue-Red2', rev=T), zlim = c(-1,1), 
      bty='n',xaxt='n',yaxt='n') # remove box and axes
mtext(occ_models[[1]]$species, 1,at = (0:30)/30, adj = 0, srt = )

