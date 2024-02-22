#' Juan S. Vargas
#' February 2024
#' 

# load libraries
library(tidyverse)

# Read in data
data_raw <- read_csv("Data/ct_records.csv",
                     col_types = cols_only(deployment_id = 'c', common_name = 'c', timestamp = col_datetime("%Y-%m-%d %H:%M") ))
stn_data_raw <- read_csv("Data/ct_deployments.csv", 
                         col_select = c(deployment_id, placename, longitude, latitude,start_date, end_date))
data_raw
stn_data_raw

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
                                 .default = common_name)) 

distinct(select(data1,common_name)) %>% view()



#### check times and dates ####
# check times
ggplot(data_raw, aes(hour(timestamp)))+geom_histogram()
# dates
ggplot(data_raw, aes(timestamp, deployment_id))+stat_summary(fun.min = min,fun.max = max,geom = "linerange")
# there is a wrong date, that has year 2035, let's find and correct.
data_raw %>% slice_max(timestamp,n=10)
# The problem is three records from station 129, possibly due to a malfunction.
# The records are of nothing so can be ignored
filter(data_raw, common_name!="Nothing") %>% 
  ggplot(aes(timestamp, deployment_id))+stat_summary(fun.min = min,fun.max = max,geom = "linerange")
# There is still one camera that has a deployment period earlier than all the rest, around march 2019 rather than 2020.
slice_min(data_raw, timestamp,n=10)
# The camera is MS#58, also called Ticho-Planes, from the Corcovado data. All
# the records are consistent, it seems this camera was just out a year before
# the rest

#### OLD CODE ####
#' This code reads in, cleans and puts all the camera-trap data in a single format.
#' THe original data came in MS EXCEL files in different formats from different partners.
# library(tidyverse)
# 
# PNC_data_raw <- readxl::read_excel("Data/Corcovado 2020.xlsx")
# PNC_data_raw
# 
# PILA_data_raw <- readxl::read_excel("Data/Data for Osa Conservation 10302020.xlsx",sheet = "DATA")
# PILA_sp_list <- readxl::read_excel("Data/Data for Osa Conservation 10302020.xlsx", sheet = "Species Names")[-c(68,70,71,72),c(2,3)]
# PILA_effort <- readxl::read_excel("Data/Data for Osa Conservation 10302020.xlsx", sheet = "EFFORT") %>% 
#   rename("Site" = site_name) %>% select(Site,long,lat,start_date,end_date) %>% 
#   rename("long" = lat, "lat" = long)
# 
# 
# PNC_effort <- PNC_data_raw %>% distinct(`Trap Station Name`, `Trap Station Latitude`, `Trap Station Longitude`) %>% 
#   rename("Site" = `Trap Station Name`, "lat" = `Trap Station Latitude`, "long" = `Trap Station Longitude`)
# PILA_data <- PILA_data_raw %>% mutate(datetime = as.POSIXct(paste(Date,Time),format = "%Y-%m-%d %H.%M.%S")) %>% 
#   rename("data_name" = Species) %>% 
#   left_join(select(PILA_sp_list,species_name,data_name)) %>% 
#   rename("Species" = species_name) %>% 
#   select(Site,Species,datetime)
# PNC_data <- PNC_data_raw %>% mutate(Species = paste(Genus, Species)) %>% 
#   select(`Trap Station Name`, `Species`, `Media Capture Timestamp`) %>% 
#   rename("Site"=`Trap Station Name`, "datetime" = `Media Capture Timestamp`)
# 
# all_stn_effort <- full_join(PILA_effort, PNC_effort)
# write_csv(all_stn_effort, "Data/station_info_all.csv")
# 
# records_all <- full_join(PILA_data, PNC_data) %>% filter(!grepl(x = Species, pattern = "unknown"), 
#                                                          !Species %in% c("Nothocercus bonapartei", "Chamaepetes unicolor", "Geotrygon chiriquensis", 
#                                                                          "Geotrygon montana", "Sciurus granatensis", "Catharus ustulatus",
#                                                                          "Parabuteo unicinctus","Buteogallus subtilis","Odontophorus guttatus" ,
#                                                                          "Pseudastur albicollis", "Columbina talpacoti", "Tigrisoma mexicanum", 
#                                                                          "Homo sapiens","Leptotila cassini","Patagioenas nigrirostris",
#                                                                          "Crypturellus soui","Formicarius analis", "Tinamus major",  
#                                                                          "Gymnocichla nudiceps","Mycteria americana", "Odontophorus gujanensis", 
#                                                                          "- -", "na", "NA NA"), 
#                                                          !is.na(Species)) 
# write_csv(records_all, "Data/combined_records.csv")
