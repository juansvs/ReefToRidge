#' Juan S. Vargas
#' February 2024
#'
#' Run this code to generate a cleaner database from the raw camera-trap data.
#' The output is a new csv file with only the station, species, and time. In it
#' we only keep indep events (dt>5 mins) of wild terrestrial vertebrates.
#' 

# load libraries
library(tidyverse)

#### Read in data ####
data_raw <- read_csv("Data/ct_records.csv",
                     col_types = cols_only(deployment_id = 'c', common_name = 'c', timestamp = col_datetime("%Y-%m-%d %H:%M") ))
stn_data_raw <- read_csv("Data/ct_deployments.csv", 
                         col_select = c(deployment_id, placename, longitude, latitude,start_date, end_date))
data_raw
stn_data_raw

#### check times and dates ####
# check times
ggplot(data_raw, aes(hour(timestamp)))+geom_histogram(binwidth = 1)
# dates
ggplot(data_raw, aes(timestamp, deployment_id))+stat_summary(fun.min = min,fun.max = max,geom = "linerange")
# there is a wrong date, that has year 2035, let's find and correct.
data_raw %>% slice_max(timestamp,n=10)
# The problem is three records from station 129, possibly due to a malfunction.
# The records are of nothing so can be ignored
filter(data_raw, common_name!="Nothing") %>% 
  ggplot(aes(timestamp, deployment_id))+stat_summary(fun.min = min,fun.max = max,geom = "linerange")
# There is still one camera that has a deployment period earlier than all the
# rest, around march 2019 rather than 2020.
slice_min(data_raw, timestamp,n=10)
# The camera is MS#58, also called Ticho-Planes, from the Corcovado data. All
# the records are consistent, it seems this camera was just out a year before
# the rest

#### station names ####
unique(data_raw$deployment_id)
# There are 121 stations in the stations db, but only 118 in the records.
unique(data_raw$deployment_id)[which(!unique(data_raw$deployment_id) %in% unique(stn_data_raw$deployment_id))]
# all the stations in the records are in the stn database, written the same,
# with the exception of records with unknown station.
count(filter(data_raw, is.na(deployment_id)))
# there are 668 records with unknown station, could be worth determining which station they come from.
unique(stn_data_raw$deployment_id)[which(!unique(stn_data_raw$deployment_id) %in% unique(data_raw$deployment_id))]
# There are three stations that were not represented in the records, 111, 72,
# and 175. 
filter(data_raw,is.na(deployment_id)) %>% reframe(range(timestamp))
filter(stn_data_raw,deployment_id %in% c("MS#111","MS#175","MS#72"))
# The dates for NA cameras go from 25-1-20 to 6-7-20. These dates match the
# deployment period of MS#111, I will change the NAs to this station
data1 <- mutate(data_raw, deployment_id = replace(deployment_id, is.na(deployment_id),"MS#111"))
which(is.na(data1$deployment_id))

# remove # from codes
data1 <- mutate(data1, deployment_id = gsub("#","", deployment_id))

#### species names ####
count(data1, common_name) %>% view()
# there are a lot of species that I can filter out in a single step, but there
# are also a lot that need to be corrected or renamed
# there are 116 different common names
data1 <- mutate(data1, common_name = tolower(common_name)) %>% # first we change everything to lower case, this reduces to 112 names
  mutate(common_name = gsub("[[:space:]]", "_", common_name)) %>% # change spaces to underscore
  mutate(common_name = gsub("-", "_", common_name)) %>% # change dash to underscore
  mutate(common_name = case_when(common_name == "little_tinamou"~"tinamou_little", 
                                 common_name == "unidentified"~"uid",
                                 common_name == "crested_guan"~"guan_crested",
                                 common_name == "raccon_uid"~"raccoon_uid",
                                 common_name == "marbled_wood_quail"~"wood_quail_marbled",
                                 common_name == "opossum"~"common_opossum",
                                 common_name == "great_tinamou"~"tinamou_great",
                                 common_name %in% c("human","ranger","tourist","hunter")~"people",
                                 .default = common_name)) 

count(data1,common_name) %>% view()

# list of species to keep for analyses, first is with people and domestic, second without.
terr_verts <- c("baird's tapir","cat_domestic","cat_uid",
                "central_american_agouti","central_american_red_brocket",
                "collared_peccary","common_opossum", "cottontail_dices",
                "cow","coyote","crab_eating_raccoon","dog",
                "great_curassow","greater_grison",
                "guan_crested","guan_black","horse","people","jaguar","jaguarundi",
                "margay","nine_banded_armadillo","northern_raccon",
                "northern_tamandua","ocelot","oncilla","opossum_four_eyed",
                "opossum_uid","peccary_uid","puma","raccoon_uid", "spotted_paca",
                "striped_hog_nosed_skunk","tayra","tinamou_great","tinamou_little",
                "white_lipped_peccary","white_nosed_coati")
wild_verts <- c("baird's_tapir","cat_uid",
                "central_american_agouti","central_american_red_brocket",
                "collared_peccary","common_opossum", "cottontail_dices",
                "coyote","crab_eating_raccoon",
                "great_curassow","greater_grison",
                "guan_crested","guan_black","jaguar","jaguarundi",
                "margay","nine_banded_armadillo","northern_raccon",
                "northern_tamandua","ocelot","oncilla","opossum_four_eyed",
                "opossum_uid","peccary_uid","puma","raccoon_uid", "spotted_paca",
                "striped_hog_nosed_skunk","tayra","tinamou_great","tinamou_little",
                "white_lipped_peccary","white_nosed_coati")

# Filter data to keep only wild vertebrates
data_wild <- filter(data1, common_name %in% wild_verts)
data_wild
# This gives us a database with 17172 entries.

#### Independent events #### 
# We need to keep only independent events. We will set a five minute threshold to classify two detections as indep
data_wild <- arrange(data_wild,deployment_id,common_name,timestamp) %>% group_by(deployment_id, common_name) %>% 
  mutate(tdif = difftime(lead(timestamp), timestamp, units = "mins")) %>% 
  filter(is.na(tdif) | tdif>duration(5,"mins")) %>% 
  ungroup()

# this reduces to 8074 independent events, of 32 different taxa. This includes
# unidentified felines (Leopardus sp.), raccoons (Procyon sp.), opossums, and
# peccaries. We have then 28 different species, 5 birds and 23 mammals
count(data_wild, common_name) %>% arrange(desc(n)) %>% mutate(logn = log10(n)) %>% view()

#### rename and export ####
select(data_wild, deployment_id, common_name, timestamp) %>% 
  rename(site = deployment_id) %>% 
  write_csv("Data/records_wild.csv")

#### Variable names across sheets ####
names(read.csv("Data/ct_deployments.csv"))
names(read.csv("Data/site_temp.csv"))
names(read.csv("Data/soil_chemistry.csv"))

# in the chemistry and temp sheets sites are called sites, which is preferable
# to deployment_id. I will change all to site. Additionally, the codes don't
# have the '#' as in the ct_depoyments database, and the sites from Corcovado
# are written as PNCxxx, not MSxxx.

# site names
soilsites <- unique(read.csv("Data/soil_chemistry.csv")[,1])
tempsites <- unique(read.csv("Data/site_temp.csv")[,1])
ctsites <- gsub("#","",read.csv("Data/ct_deployments.csv")[,2])

# which ones are missing?
soilsites[!soilsites %in% ctsites]
tempsites[!tempsites %in% ctsites]

# There are 7 temperature sites that are not in the camera-trap sites. MS63, 71,
# 124, 130, 159, 169, and PBNP1_new. 63, 71, and 124 broke down and have no
# data.There are more soil chemistry sites that are not in the camera-trap data.
# For some it's because they have a different code (e.g. PNC instead of MS). I
# need to find out the identity of the PNC sites.

# The location of the temp loggers and chem sites is not always the same as the
# camera location. The data for temps and chem doesn't have coordinates, so I
# need to extract these from somewhere else.

