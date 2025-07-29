#' @export
DetectBlinkBlock <- function (df, pupil_col, blink_value = 0) {
  data <- df[[pupil_col]]
  # set two vectors to save the start and the end of a blink
  sblink <- numeric(0)
  eblink <- numeric(0)
  # get the full indices of blinks
  ifelse(is.na(blink_value), index <- which(is.na(data)),
         index <- which(data == blink_value))
  if (length(index) == 0) {
    message("No blinks detected. Check if the blink_value is set correctly.")
    return(list(startat = sblink, endat = eblink))
  }
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
    stop("The length of StartBlink and EndBlink doesn't match.")
  return(list(startat = sblink, endat = eblink))
}
