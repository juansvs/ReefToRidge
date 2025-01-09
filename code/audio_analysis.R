library(tidyverse)
library(ggpubr)

# read in data
audiosdb_sub <- read.csv("Data/soundscape_file_subset.csv")
# features database, combining multiple csv outputs from Google colab code
audio_features_db <- do.call(rbind, lapply(
  list.files("Data", "acousticfeatures", full.names = T), 
  read.csv, header = F))
# give names to db
names(audio_features_db) <- c(
  'ZCR',  'MEANt',  'VARt', 'SKEWt',  'KURTt','LEQt', 'BGNt', 'SNRt', 'MED','Ht', 'ACTtFraction', 'ACTtCount',  'ACTtMean', 'EVNtFraction', 'EVNtMean', 'EVNtCount',## temporal indices
  'MEANf','VARf', 'SKEWf', 'KURTf', 'NBPEAKS', 'LEQf', 'ENRf', 'BGNf', 'SNRf','Hf', 'EAS','ECU','ECV','EPS','EPS_KURT','EPS_SKEW',
  'ACI','NDSI', 'rBA','AnthroEnergy', 'BioEnergy','BI', 'ROU','ADI','AEI','LFC','MFC','HFC','ACTspFract', 'ACTspCount', 'ACTspMean', 
  'EVNspFract', 'EVNspMean','EVNspCount', 'TFSD', 'H_Havrda', 'H_Renyi','H_pairedShannon',  'H_gamma',  'H_GiniSimpson','RAOQ', 'AGI',
  'ROItotal', 'ROIcover', ## spectral indices
  'File.Paths',# file name
  'start_s','end_s' # segment lims (s)
  )

audio_features_db <- mutate(audio_features_db, site = )

# see relationship between indices of interest
plot(select(audio_features_db, NBPEAKS, ACI, MFC, NDSI, BI, Ht, EPS, ECU))

# join with info about file
acoustic_jointdb <- right_join(audio_features_db, audiosdb_sub) %>% 
  inner_join(allsites_cov_db, by = join_by(stn == deployment_id)) %>% 
  mutate(dttm = ymd_hms(date))

# aggregate the index calculated into weekly values, and only taking sunrise and sunset values
summarised_MFC <- acoustic_jointdb %>% filter(hour(dttm)>5) %>% # keep only daytime recordings
  mutate(wk = week(dttm)) %>% # get the week of the recording
  summarise(median_MFC=median(MFC), .by = c(site, wk)) %>% 
  left_join(select(allsites_cov_db, site, alt1,lc1))
ggplot(summarised_MFC, aes(alt1, median_MFC, color = factor(lc1)))+geom_point()
# this plot looks very similar to taking a linerange with median and inter-quartile distance

#### plot ####
  # mutate(npk_scl = datawizard::standardise(NBPEAKS), mfc_scl = datawizard::standardise(MFC)) %>% 
  # mutate(index = (mfc_scl+npk_scl)/2) %>% 
  # slice_max(order_by = EPS, by = stn) %>% 
filter(acoustic_jointdb, hour(dttm)>5) %>% 
ggplot(aes(alt1,MFC))+
  stat_summary(aes(color = factor(lc1)), fun = median, fun.max = \(a) quantile(a, prob=0.75), fun.min = \(a) quantile(a, prob=0.25), show.legend = F) +
  labs(x = "Altitude (masl)", y = "Acoustic index (MFC)", color = "Land cover")+
  scale_color_manual(values = c("darkgreen", "goldenrod"), labels = c("Forest", "Pasture/Grassland"))+
  #theme_classic2(base_size = 16)+
  theme(aspect.ratio = 1)

