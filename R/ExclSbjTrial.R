#' @export
ExclSbjTrial <- function(df, sbj_col, trial_col, prop_sbj = 0.2, pupil_col,
                         prop_trial = 0.2, export.path = NULL,
                         remove = TRUE, quiet = FALSE) {
  sbj_sym <- rlang::sym(sbj_col)
  trial_sym <- rlang::sym(trial_col)
  pupil_sym <- rlang::sym(pupil_col)
  if (!quiet) message("Starting trial and subject exclusion...")
  df <- df %>% dplyr::group_by(!!sbj_sym, !!trial_sym)
  trial_missing_summary <- df %>%
    summarise(n = n(), blink = sum(!!pupil_sym == 0), missing = sum(is.na(!!pupil_sym)),
              missing_prop = round((blink + missing)/n, 4),
              badtrial = (missing_prop > prop_trial), .groups = "keep")
  # export trial check summary
  if (!is.null(export.path)) {
    if (!quiet) message("Trial summary complete. Exporting results...")
    write.csv(trial_missing_summary, export.path)
    if (!quiet) message(paste0("Trial summary has been successfully saved as:\n", export.path))
  }
  # get info per group
  summary_list <- trial_missing_summary %>% dplyr::group_by(!!sbj_sym) %>% group_split()
  # save bad trials
  bad <- data.frame()
  # Process each subject, respectively
  for (i in 1:length(summary_list)) {
    data <- summary_list[[i]]
    current_sbj <- unique(data[[sbj_col]])
    if (!quiet) message(paste0("Trial Check Summary for "), current_sbj, ":")
    bad_trials <- data %>% dplyr::filter(badtrial)
    bad <- dplyr::bind_rows(bad, bad_trials)
    sbjlost_prep <- round(nrow(bad_trials)/nrow(data),4)
    if (!quiet) message(paste0("Bad trial percentage: ", sbjlost_prep*100, "%"))
    if (sbjlost_prep > prop_sbj) {
      stop(paste0("Subject ", current_sbj, " excluded, processing stops."))
    } else message("Subject kept. Keep preprocesing...")
  }
  trial_missing_summary <- tibble::as_tibble(trial_missing_summary) %>% dplyr::ungroup()
  bad <- tibble::as_tibble(bad) %>% dplyr::ungroup()
  if (remove) {
    df <- df %>% dplyr::filter(!(paste(!!sbj_sym, !!trial_sym) %in%
                                   paste(bad[[sbj_col]], bad[[trial_col]])))
    df <- tibble::as_tibble(dplyr::ungroup(df))
    if (!quiet) message("Bad trials have been successfully removed from the data.")
  }
  return(list(data = df, excl.summary = trial_missing_summary, bad_trials = bad))
}
