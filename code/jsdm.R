library(terra)
library(sjSDM)
camvect <- vect("Data/GIS/CT_wcovar.geojson")

DAT <- read_csv("../Data/records_wild.csv")
occ_data_bin <- DAT %>% 
  rename(species = common_name, site = deployment_id) %>%
  distinct(site, species) %>%
  mutate(n=1) %>% 
  pivot_wider(names_from = species, values_from = n, values_fill = 0) %>% 
  column_to_rownames("site") %>% 
  as.matrix()


env_mat <- tibble(site = rownames(occ_data_bin)) %>%
  left_join(as.data.frame(camvect)) %>% 
  column_to_rownames("site") %>% 
  select(easting, northing,LC_1,alt_1,ghm_1,distance) %>% 
  rename(x = easting, y = northing, landcov = LC_1, alt = alt_1, ghm = ghm_1) %>% 
  mutate(landcov = case_match(landcov, 12~1,21~0)) %>% # 1 is forest, 0 is grassland
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

sapply(occ_models, logLik)
# The highest log-likelihood is obtained with the base model, which has no spatial component.
plot(occ_models$base)
# This model clearly shows the influence of altitude, species like tinamous,
# tamanduas, armadillos, curassows, raccoons, and agoutis were associated with
# lowlands, while oncillas, jaguars, guans, and brocket were associated with
# highlands.