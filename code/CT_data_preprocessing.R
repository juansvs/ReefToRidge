#' This code reads in, cleans and puts all the camera-trap data in a single format.
#' THe original data came in MS EXCEL files in different formats from different partners.

library(tidyverse)

PNC_data_raw <- readxl::read_excel("Data/Corcovado 2020.xlsx")
PNC_data_raw

PILA_data_raw <- readxl::read_excel("Data/Data for Osa Conservation 10302020.xlsx",sheet = "DATA")
PILA_sp_list <- readxl::read_excel("Data/Data for Osa Conservation 10302020.xlsx", sheet = "Species Names")[-c(68,70,71,72),c(2,3)]
PILA_effort <- readxl::read_excel("Data/Data for Osa Conservation 10302020.xlsx", sheet = "EFFORT") %>% 
  rename("Site" = site_name) %>% select(Site,long,lat,start_date,end_date) %>% 
  rename("long" = lat, "lat" = long)


PNC_effort <- PNC_data_raw %>% distinct(`Trap Station Name`, `Trap Station Latitude`, `Trap Station Longitude`) %>% 
  rename("Site" = `Trap Station Name`, "lat" = `Trap Station Latitude`, "long" = `Trap Station Longitude`)
PILA_data <- PILA_data_raw %>% mutate(datetime = as.POSIXct(paste(Date,Time),format = "%Y-%m-%d %H.%M.%S")) %>% 
  rename("data_name" = Species) %>% 
  left_join(select(PILA_sp_list,species_name,data_name)) %>% 
  rename("Species" = species_name) %>% 
  select(Site,Species,datetime)
PNC_data <- PNC_data_raw %>% mutate(Species = paste(Genus, Species)) %>% 
  select(`Trap Station Name`, `Species`, `Media Capture Timestamp`) %>% 
  rename("Site"=`Trap Station Name`, "datetime" = `Media Capture Timestamp`)

all_stn_effort <- full_join(PILA_effort, PNC_effort)
write_csv(all_stn_effort, "Data/station_info_all.csv")

records_all <- full_join(PILA_data, PNC_data) %>% filter(!grepl(x = Species, pattern = "unknown"), 
                                                         !Species %in% c("Nothocercus bonapartei", "Chamaepetes unicolor", "Geotrygon chiriquensis", 
                                                                         "Geotrygon montana", "Sciurus granatensis", "Catharus ustulatus",
                                                                         "Parabuteo unicinctus","Buteogallus subtilis","Odontophorus guttatus" ,
                                                                         "Pseudastur albicollis", "Columbina talpacoti", "Tigrisoma mexicanum", 
                                                                         "Homo sapiens","Leptotila cassini","Patagioenas nigrirostris",
                                                                         "Crypturellus soui","Formicarius analis", "Tinamus major",  
                                                                         "Gymnocichla nudiceps","Mycteria americana", "Odontophorus gujanensis", 
                                                                         "- -", "na", "NA NA"), 
                                                         !is.na(Species)) 
write_csv(records_all, "Data/combined_records.csv")
