library(vegan)
library(tidyverse)

camvect <- terra::vect("Data/ct_covars.geojson")
campts_db <- as.data.frame(camvect) %>% mutate(site = gsub("#", "", deployment_id), days = as.numeric(days))

# create species matrix
DAT <- read_csv("Data/records_wild.csv") %>% rename(species = common_name) %>% 
  filter(!grepl("uid", species))
sp_mat <- DAT %>% count(site, species) %>% 
  pivot_wider(names_from = species, values_from = n, values_fill = 0) %>% 
  column_to_rownames("site") 
# create env matrix
env_mat <- tibble(site = rownames(sp_mat)) %>% 
  left_join(select(campts_db,site, lfdistance, ghm1, alt1, evi_mean1, temp_mean1, area)) %>% 
  rename(fdist = lfdistance, ghm = ghm1, alt = alt1, evi = evi_mean1, temp = temp_mean1) %>% 
  column_to_rownames("site")

# Ordination
mds <- metaMDS(sp_mat, trymax = 40)

# fit environmental variables
ordi_covar <- envfit(mds, env_mat)
ordi_covar

# see if there are significant differences across areas using ANOSIM
vertanosim <- anosim(sp_mat, env_mat$area)
summary(vertanosim)
# plot the two together, with ellipses for sites near and far (>1000m) from
# large forest patches
ordiplot(mds, type="n")
points(mds,"sites", pch=as.integer(as.factor(env_mat$area)))
# ordiellipse(mds, groups = env_mat$fdist<1000, col = ellcols, draw = "polygon")
ordiellipse(mds, groups = env_mat$area, draw = "polygon", col = hcl.colors(3))
text(mds, "species", col="darkblue", cex=0.9) 
plot(ordi_covar, col = "black", cex = 1.2)
# ellcols <- c("darkgreen", "gold")

### Beetles ---------
DATb <- read.csv("Data/dung_beetles_prc.csv") %>% 
  # join traps from same station
  summarise(.by = Plot, across(where(is.numeric), sum)) %>% 
  inner_join(select(campts_db, site), by = join_by("Plot"=="site"))
# create sp matrix
sp_mat_beet <- DATb %>% 
  column_to_rownames("Plot") %>% 
  select(where(\(x) sum(x)>5)) %>% # only species present at multiple (>5) sites
  filter(rowSums(.)>0) %>% # filter empty sites
  as.matrix()

# create env matrix
env_mat_beet <- tibble(site = rownames(occ_data_binb)) %>% left_join(campts_db) %>% 
  select(site, lfdistance, ghm1, alt1, evi_mean1, temp_mean1, area) %>% 
  rename(alt = alt1, ghm = ghm1, fdist = lfdistance) %>% 
  # mutate(forest = ifelse(landcov==12,"Y","N")) %>% 
  column_to_rownames("site")

# Ordination
mds_beet <- metaMDS(sp_mat_beet)

# ANoSIM
beetanosim <- anosim(sp_mat_beet, env_mat_beet$area)
beetanosim
# fit environmental variables
ordi_covar_beet <- envfit(mds_beet, env_mat_beet)
# plot the two together, with ellipses for sites near and far (>1000m) from
# large forest patches
ordiplot(mds_beet, type="n")
points(mds_beet,"sites", pch=as.integer(as.factor(env_mat_beet$area)))
ordiellipse(mds_beet, groups = env_mat_beet$area, draw = "polygon", col = hcl.colors(3))
text(mds_beet, "species", col="blue", cex=0.8) 
plot(ordi_covar_beet, col = "black")
