library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(eegkit)
library(readr)

source('preprocess_eeg.R')

eeg.files <- list.files(path="data", pattern=".bdf", full.names = T)

for(file in eeg.files){
  print(paste("Working on ", file))
  subject <- substr(file, 14, 15)
  preprocessed.data <- preprocess_eeg(file = file, subject_id = subject)
  write_csv(preprocessed.data, file=paste0("data/preprocessed/subject-", subject,"-epochs.csv"))
}

# compress one step further into RDS

all.eeg.epoch.files <- list.files(path="data/preprocessed", pattern=".csv", full.names = T)

all.eeg.epochs <- lapply(all.eeg.epoch.files, read_csv)


all.eeg.epochs.df <- bind_rows(all.eeg.epochs)

saveRDS(all.eeg.epochs.df, file="data/preprocessed/epochs.rds")




