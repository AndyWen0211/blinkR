#' @export
ExtendBlinks <- function (df, pupil = 'pupil_size', prior = 100,
                          post = 200, samp_rate = 500) {
  # add another column to monitor extend status
  df$is_blink_extended <- NA
  data <- df[[pupil]]
  # get the full indices of blinks
  index <- which(data == 0)
  if (length(index) == 0) {
    message("No blinks detected.")
    return(as_tibble(df))
  }
  # create vectors to save the StartBlink and EndBlink indices
  sblink <- numeric(0)
  eblink <- numeric(0)
  # assign the first index to StartBlink
  sblink[1] <- index[1]
  # detect each start and end of the blink
  for (i in 2:(length(index)-1)) {
    value1 <- index[i] - index[i-1]
    value2 <- index[i+1] - index[i]
    if (value1 != 1) sblink <- append(index[i], sblink, after = 0)
    if (value2 != 1) eblink <- append(index[i], eblink, after = 0)
  }
  # assign the last index to EndBlink
  eblink <- append(index[length(index)], eblink, after = 0)
  # The length of sblink and endblink should be the same
  if (length(sblink) != length(eblink))
    stop("The length of StartBlink and EndBlink doesn't match. Please tell Andy to check where goes wrong.")
  # Based on the sampling rate, get the number of sample points before and
  # after blink that might be affected by noise due to eyelid movements.
  interval <- 1000/samp_rate
  forward <- floor(prior/interval)
  backward <- floor(post/interval)
  # Reset the pupil size that is potentially distorted by eyelid movement
  n <- length(eblink)
  for (i in 1:n) {
    df[[pupil]][max(1, sblink[i]-forward):(sblink[i]-1)] <- 0
    df[[pupil]][(eblink[i]+1):min(eblink[i]+backward, nrow(df))] <- 0
    # Tag extended blinks
    # NOTE: is_blink_extended = NA -> True data, not a blink
    # is_blink_extended = FALSE -> True blinks
    # is_blink_extended = TRUE -> Data potentially being affected by a blink but not a blink originally
    df$is_blink_extended[max(1, sblink[i]-forward):(sblink[i]-1)] <- TRUE
    df$is_blink_extended[(eblink[i]+1):min(eblink[i]+backward, nrow(df))] <- TRUE
    df$is_blink_extended[sblink[i]:eblink[i]] <- FALSE
  }
  return(as_tibble(df))
}
