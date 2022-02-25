preprocess_eeg <- function(file, subject_id){
  data <- read_eeg_file_tidy(file) %>%
    linked_ears_rereference() %>%
    bandpass_filter(low=0.1, high=40, sampling_rate=500) %>%
    segment(start = -100, end = 600, sampling_rate=500) %>%
    artifact_rejection(max_range=200) %>%
    baseline_correct() %>%
    mutate(subject = subject_id)
    
  return(data)
}

read_eeg_file_tidy <- function(file) {
  head <- edfReader::readEdfHeader(file)
  signals <- edfReader::readEdfSignals(head, signals="Ordinary")
  
  data <- map_df(signals, "signal")
  data$sample_id <- 1:nrow(data)
  
  events <- data %>% select(sample_id, TRIGGER) %>% dplyr::filter(TRIGGER < 256)
  
  data$`Packet Counter` <- NULL
  data$ACC21 <- NULL
  data$ACC22 <- NULL
  data$ACC23 <- NULL
  data$`ExG 1` <- NULL
  data$TRIGGER <- NULL
  
  tidy_data <- pivot_longer(data, 1:19, names_to = "electrode", values_to = "v") %>% select(sample_id, electrode, v, A2)
  
  return(list(signals=tidy_data, events=events))
}

linked_ears_rereference <- function(data){
  df <- data$signals
  df <- df %>% 
    mutate(v = v - (A2/2)) %>%
    select(-A2)
  data$signals <- df
  return(data)
}

bandpass_filter <- function(data, low, high, sampling_rate){
  order <- 3 * sampling_rate / low
  
  df <- data$signals
  df <- df %>% 
    group_by(electrode) %>% 
    mutate(v = eegfilter(v, sampling_rate, low, high, method="fir1", order=order)[,1]) %>%
    ungroup()
  data$signals <- df
  return(data)
}

segment <- function(data, start, end, sampling_rate){
  events <- data$events
  signals <- data$signals
  
  ms_per_sample <- 1000/sampling_rate
  min_sample <- start/ms_per_sample
  max_sample <- end/ms_per_sample
  
  epochs <- events %>% mutate(event_id = 1:n()) %>%
    group_by(event_id) %>%
    summarize(
      t=seq(start, end, ms_per_sample), 
      sample_id=seq(sample_id+min_sample, sample_id+max_sample, 1), 
      event_type=TRIGGER) %>%
    ungroup() %>%
    mutate(event_type=factor(event_type)) %>%
    left_join(signals, by="sample_id") %>%
    select(-sample_id)
  
  return(epochs)    
}

artifact_rejection <- function(epochs, max_range){
  d <- epochs %>% group_by(event_id, electrode) %>%
    summarize(r = max(v) - min(v)) %>%
    mutate(in_threshold = r <= max_range) %>%
    group_by(event_id) %>%
    summarize(good_segment = all(in_threshold))
  
  epochs.ar <- epochs %>%
    left_join(d, by=c("event_id"))
  
  return(epochs.ar)
}

baseline_correct <- function(epochs){
  baseline.means <- epochs %>%
    group_by(electrode, event_id) %>%
    dplyr::filter(t <= 0) %>%
    summarize(baseline.mean = mean(v))
  
  epoch.bc <- epochs %>%
    left_join(baseline.means, by=c("electrode", "event_id")) %>%
    mutate(v = v - baseline.mean) %>%
    select(-baseline.mean)
  
  return(epoch.bc)
}