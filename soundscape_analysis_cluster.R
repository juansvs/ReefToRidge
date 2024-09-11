library(tuneR)
library(seewave)
library(audio)
library(phonTools)
library(tibble)
library(soundecology)
# library(kableExtra)
library(dplyr)
# library(pbapply)
# library(ggplot2)
library(stringr)
# library(foreach)
# library(doParallel)
library(progress)
library(doParallel)
library(googledrive)

# Variables

# Second to start analysis from recordings
start <- 5
# Second to end analysis from recordings 
end <- 290
# Threshold in decibels to ignore noise or background noise from sounds
threshdB <- -35
# Max frequency to be analysed
max_freq <- 12000
# Width of frequency bins used to make the spectrogram
freq_step <- 1000

# Lower freq for high pass filter
lowfreq <- 200
# Number of cores to make parallel processing
# numCores <- parallel::detectCores(logical = F)-1
numCores <- length(future::availableWorkers())

#--------------------Read files-------------------------------------
# site <- "MySite"
# audios <- list.files(paste0("E:/Data/Audios/", site),
#                      "*.wav",
#                      full.names = T,
#                      include.dirs = T,
#                      recursive = T)
filedb <- read.csv("Data/drive_audio_files.csv")

#-------------Google drive------------------
drive_auth("juans.vs@gmail.com")

# -------------------------Define functions---------------------------
filter_fun <- function(audio,
                       filtering = "none",
                       # frequency = 8000,
                       lower = 200,
                       higher = 9000){
  
  frequency <- audio@samp.rate
  
  # Filter low-pass
  if(filtering == "none"){
    audio2 <- audio
  }
  if(filtering == "low"){
    audio2 <- ffilter(audio,
                      f = frequency,
                      to = higher,
                      rescale = F)
  }
  if(filtering == "high"){
    audio2 <- ffilter(audio,
                      f = frequency,
                      from = lower,
                      rescale = F)
  }
  if(filtering == "band-pass"){
    audio2 <- ffilter(audio,
                      f = frequency,
                      from = lower,
                      to = higher,
                      rescale = F)
  }
  if(filtering == "band-stop"){
    audio2 <- ffilter(audio,
                      f = frequency,
                      from = lower,
                      to = higher,
                      rescale = F,
                      bandpass = FALSE)
  }
  if(filtering != "none"){
    audio <- Wave(audio2,
                  samp.rate = audio@samp.rate,
                  bit = audio@bit)
  }
  audio
}

alpha_ind <- function(audio){
  
  # Spectral entropy, calculated from the spectrogram, scaled between 0 and 1, also known as Pielou's eveness index,
  # SE = sh(spec(audio,
  #              f = 24000,
  #              wl = 512,
  #              flim = c(0,12)))
  # dividing the spectrogram into bins (default 10, each one of 1000 Hz) and taking the proportion of the signals in each bin above a threshold (default -50 dBFS). The ADI is the result of the Shannon index applied to these bins.
  AD = acoustic_diversity(audio,
                          max_freq = max_freq,
                          db_threshold = threshdB,
                          freq_step = freq_step)$adi_left
  
  tibble(AD)
}
  
calc_f <- function(x){
    audio1 <- readWave(x,
                       from = start,
                       to = end,
                       units = "seconds")
    audio_f1 <- filter_fun(audio1,
                           filtering = "high",
                           lower = lowfreq)#,
    #higher = 9000)
    resul1 <- alpha_ind(audio1)
    resul1 |>
      mutate(file = x) |>
      select(file, AD) |>
      as.data.frame()
}

# ------------------------Start parallel process--------------------------------

# registerDoParallel(numCores)

cl <- makeCluster(numCores)
registerDoParallel(cl)

clusterExport(cl, "filedb")

# df <- foreach(x = audios,
system.time(
df <- foreach(x = 1:32,
              .combine = rbind,
              .packages = c("seewave",
                            "dplyr",
                            "tuneR",
                            "soundecology",
                            "googledrive")) %dopar% {
                aud <- drive_download(as_id(filedb$id[x]))
                try(calc_f(aud$local_path))
                # file.remove(aud$local_path)
              }
)
stopCluster(cl)

write.csv(df,
          paste0("Results/All",site,"_thresh",threshdB,"maxFreq",max_freq,".csv"))