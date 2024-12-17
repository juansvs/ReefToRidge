library(tidyverse)

# read in database of audio files
audiosdb <- read_csv("Data/soundscape_file_db.csv")
names(audiosdb)

tibble(audiosdb,wk = week(audiosdb$date)) %>% distinct(stn,wk) %>% count(stn) %>% view()

# there are too many files to process, so we take only a subset of the data. We
# select only from Feb 2020 on, and the first 4 distinct weeks of sampling since
# that date. We also subset to select only every third day.
audiosdb_sub <- filter(audiosdb, date>=as.POSIXct("2020-02-01")) %>% 
  mutate(sw = as.integer(difftime(date,as.POSIXct("2020-02-01"),units = "weeks"))) %>% # sampling week
  mutate(gid = cur_group_id(), .by=c(stn,sw)) %>% # number of sampling week (combined)
  mutate(sw2 = gid-min(gid), .by=stn) %>% # diff wrt first sampling week by stn
  filter(sw2<4, as.integer(difftime(date,as.POSIXct("2020-02-01"), units = "days"))%%3==1) # every third day, first four (not necessarily consecutive) sampling weeks
# export
write.csv(audiosdb_sub, "Data/soundscape_file_subset.csv", quote = F, row.names = F)
