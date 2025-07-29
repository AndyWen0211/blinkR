#' @export
InterpolateLinearPupil <- function (df, pupil_col, sbj_col, trial_col,
                                    maxgap = Inf, samp_rate = 500) {
  # set a temporary feature managing delete status.
  df$delete <- FALSE
  # Detect blink blocks that are larger than the threshold and remove them
  data <- df[[pupil_col]]
  if (any(data == 0, na.rm = TRUE)) {
    warning("Detected 0s in data. Make sure they are properly converted to NA if needed.")
  }
  blink_boundary <- DetectBlinkBlock(df, pupil_col = pupil_col, blink_value = NA)
  sblink <- blink_boundary$startat
  if (length(sblink) == 0) {
    warning("The current data doesn't contain any blinks. Did you set all blinks to be NA?")
    return(as_tibble(df[,1:(ncol(df)-1)]))
  }
  eblink <- blink_boundary$endat
  # interval between two sample points
  interval <- 1000/samp_rate
  # set a counter checking the percentage of long period of blinks
  n <- 0
  for (i in 1:length(sblink)) {
    if ((eblink[i] - sblink[i])*interval > maxgap) {
      df$delete[sblink[i]:eblink[i]] <- TRUE
      n<-n+1
    }
  }
  prop_delete <- n/length(sblink)
  prop_retain <- 1-prop_delete
  df <- dplyr::filter(df, delete == FALSE)
  # Conduct linear interpolation grouped by trial
  # Removing potential single-sample blinks enables interpolation safer.
  subject_sym <- rlang::sym(sbj_col)
  trial_sym <- rlang::sym(trial_col)
  pupil_sym <- rlang::sym(pupil_col)
  df <- df %>% group_by(!!subject_sym,!!trial_sym) %>%
    mutate(!!pupil_sym := na.approx(!!pupil_sym, rule = 2))
  ungroup(df)
  message(paste0(round(prop_delete, 5)*100,"% blink block(s) deleted due to excessive long period."))
  message(paste0(round(prop_retain, 5)*100,"% blink block(s) linearly interpolated."))
  return(as_tibble(df[,1:(ncol(df)-1)]))
}
