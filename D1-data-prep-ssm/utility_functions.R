#' Split track at gaps
#' 
#' @param data Data frame with (at least) columns for "ID" and "time"
#' @param max_gap Longest allowed gap, in minutes (track will be split at longer gaps)
#' @param shortest_track Shortest track to keep after splitting, in minutes. Shorter
#' tracks will be removed from the output data set.
#' 
#' @return Data frame with identical structure as input, where ID column
#' has been replaced by new ID for split tracks. Old ID still accessible as
#' ID_old column
split_at_gap <- function(data, max_gap = 60, shortest_track = 0) {
  # Number of tracks
  n_tracks <- length(unique(data$ID))
  
  # Save old ID and reinitialise ID column
  data$ID_old <- data$ID
  data$ID <- character(nrow(data))
  
  # Loop over tracks (i.e., over IDs)
  for(i_track in 1:n_tracks) {
    # Indices for this track
    ind_this_track <- which(data$ID_old == unique(data$ID_old)[i_track])
    track_length <- length(ind_this_track)
    
    # Time intervals in min
    dtimes <- difftime(data$time[ind_this_track[-1]], 
                       data$time[ind_this_track[-track_length]],
                       units = "mins")
    
    # Indices of gaps longer than max_gap
    ind_gap <- c(0, which(dtimes > max_gap), track_length)
    
    # Create new ID based on split track
    subtrack_ID <- rep(1:(length(ind_gap) - 1), diff(ind_gap))
    data$ID[ind_this_track] <- paste0(data$ID_old[ind_this_track], "-", subtrack_ID)
  }
  
  # Only keep sub-tracks longer than some duration
  track_lengths <- sapply(unique(data$ID), function(id) {
    ind <- which(data$ID == id)
    difftime(data$time[ind[length(ind)]], data$time[ind[1]], units = "min")
  })
  ID_keep <- names(track_lengths)[which(track_lengths >= shortest_track)]
  data <- subset(data, ID %in% ID_keep)
  
  return(data)
}

#' Function to estimate proportion of missing locations 
#' 
#' @param time vector of times
#' @param res chosen resolution of data
#' 
#' @return proportion of total observations that would be missing (for chosen `res`)

p_na <- function(time, res) {
  if(res <= 60){
    time <- round_date(time, paste(res,"min")) # round to nearest resolution
  }else{
    time <- round_date(time, paste(res/60,"hour")) # round to nearest resolution
  }
  
  time <- na.omit(time[time != lead(time)]) # remove duplicate time
  # calculate maximum number of locations          
  max_n_loc <- length(seq(time[1], tail(time, 1) + res*60, 
                          by = as.difftime(res, units = "mins")))
  n_NA <- max_n_loc - (length(time)+1)
  return(n_NA / max_n_loc)
}



