#' @export
BinifyPupil <- function(df, baseline, sbj_col, trial_col, time_col, stage_col,
                        onset, old_samp_rate = 500,
                        new_sample_rate = 50, method = c("mean", "median", "rms"),
                        stage.label = c("majority", "middle", "first.appear", "last", "mix"),
                        additional.onsets = NULL, name.bin = "bin",
                        basicInfoOnly = FALSE ,quiet.msg = FALSE, quiet.warning = FALSE,
                        quiet.debug.msg = FALSE, quiet.debug.warning = FALSE,...) {
  sbj_sym <- rlang::sym(sbj_col)
  trial_sym <- rlang::sym(trial_col)
  stage_sym <- rlang::sym(stage_col)
  time_sym <- rlang::sym(time_col)
  bin_sym <- rlang::sym(name.bin)
  base_sym <- rlang::sym(baseline)
  stage.label <- if (is.character(stage.label)) match.arg(stage.label) else stage.label
  if (is.null(stage.label) && !basicInfoOnly)
    stop("To skip stage tagging and timestamp assignment, please set basicInfoOnly = TRUE.")
  solve <- if (isTRUE(list(...)$append.missing) || is.null(list(...)$append.missing)) TRUE else FALSE
  df <- df %>% group_by(!!sbj_sym, !!trial_sym)
  # define a root mean square function
  rms_func <- function (x,...) return(sqrt(mean(x^2,...)))
  if (!quiet.msg) message("Begin adding bin info...")
  # calculate how many samples are in the same bin
  amount <- floor(old_samp_rate/new_sample_rate)
  df <- df %>% group_modify(~{
    n <- nrow(.x) - nrow(.x) %% amount
    # drop the tails
    dplyr::slice_head(.x, n = n)
  })
  # assigning bin#
  df <- df %>% group_modify(~{
    dplyr::mutate(.x, !!bin_sym := rep(1:(nrow(.x)/amount), each = amount))
  })
  df <- dplyr::relocate(df, !!bin_sym, .after = !!trial_sym)
  df <- df %>% dplyr::ungroup() %>% group_by(!!sbj_sym, !!trial_sym, !!bin_sym)
  if (!basicInfoOnly) {
    # save the processing outcome for later use
    df.binned <- df
    # Add column-related info
    info <- df %>% summarise(across(everything(), ~n_distinct(.x) == 1), .groups = "keep") %>%
      dplyr::ungroup() %>% dplyr::select(-!!sbj_sym, -!!trial_sym, -!!bin_sym)
    cols_to_add <- names(info %>%  dplyr::select(where(~all(.x, na.rm = TRUE))))
    df_to_add <- df %>% summarise(across(all_of(cols_to_add), ~unique(.x)), .groups = "keep")
  }
  # re-sampling
  if (is.character(method)) {
    method <- match.arg(method)
    if (method == "rms") func <- rms_func
    else func <- get(method)
  }
  if (is.function(method)) func <- method
  df <- df %>% summarise(!!base_sym := func(!!base_sym,...), .groups = "drop")
  if (!quiet.msg) {
    if (is.function(method)) message(paste0("Bins created via function: ", substitute(method), " across every ", floor(1000/new_sample_rate), " ms per subject and trial."))
    else {
      if (method == "mean" || method == "median") message(paste0("Bins created using the ", method, " value of every ", floor(1000/new_sample_rate), " ms per subject and trial."))
      if (method == "rms") message(paste0("Bins created using the root mean square of every ", floor(1000/new_sample_rate), " ms per subject and trial."))
    }
  }
  if (basicInfoOnly) {
    if (!quiet.warning)
      warning("Get core info only. Parameters below are ignored:\n'time_col', 'stage_col', and 'stage.label'.")
    if (!quiet.msg) message("Finish processing!")
    return(tibble::as_tibble(dplyr::ungroup(df)))
  }
  df <- dplyr::left_join(df, df_to_add, by = c(sbj_col, trial_col, name.bin))
  # add stage info
  if (is.character(stage.label)) {
    if (stage.label == "majority") {
      set_stage <- function(x) {
        # get the majority item of a vector
        return(unique(x)[which.max(tabulate(match(x,unique(x))))])
      }
      stage.info <- df.binned %>%
        summarise(!!stage_sym := set_stage(!!stage_sym),.groups = "keep")
      df <- dplyr::left_join(df, stage.info, by = c(sbj_col, trial_col, name.bin))
      df <- dplyr::relocate(df, !!stage_sym, .after = !!base_sym)
      if (!quiet.msg) message("Stage labels computed using the majority rule: the most frequent stage of each bin (per subject and trial).")
    }
    if (stage.label == "middle") {
      set_stage <- function(x) {
        # get the middle item of a vector
        return(x[ceiling(length(x)/2)])
      }
      stage.info <- df.binned %>%
        summarise(!!stage_sym := set_stage(!!stage_sym),.groups = "keep")
      df <- dplyr::left_join(df, stage.info, by = c(sbj_col, trial_col, name.bin))
      df <- dplyr::relocate(df, !!stage_sym, .after = !!base_sym)
      if (!quiet.msg) message("Stage labels computed by selecting the stage at the temporal midpoint within each bin (per subject and trial).")
    }
    if (stage.label == "first.appear") {
      set_stage <- function(x) {
        # get the first item of a vector
        return(x[1])
      }
      stage.info <- df.binned %>%
        summarise(!!stage_sym := set_stage(!!stage_sym),.groups = "keep")
      df <- dplyr::left_join(df, stage.info, by = c(sbj_col, trial_col, name.bin))
      df <- dplyr::relocate(df, !!stage_sym, .after = !!base_sym)
      if (!quiet.msg) message("Stage labels are assigned using the first occurring stage within each bin (per subject and trial).")
    }
    if (stage.label == "last") {
      set_stage <- function(x) {
        # get the last item of a vector
        return(x[length(x)])
      }
      stage.info <- df.binned %>%
        summarise(!!stage_sym := set_stage(!!stage_sym),.groups = "keep")
      df <- dplyr::left_join(df, stage.info, by = c(sbj_col, trial_col, name.bin))
      df <- dplyr::relocate(df, !!stage_sym, .after = !!base_sym)
      if (!quiet.msg) message("Stage labels are assigned using the last occurring stage within each bin (per subject and trial).")
    }
    if (stage.label == "mix") {
      set_stage <- function(x) {
        if (length(unique(x)) == 1) return(unique(x))
        else return("Mix")
      }
      stage.info <- df.binned %>%
        summarise(!!stage_sym := set_stage(!!stage_sym),.groups = "keep")
      df <- dplyr::left_join(df, stage.info, by = c(sbj_col, trial_col, name.bin))
      df <- dplyr::relocate(df, !!stage_sym, .after = !!base_sym)
      if (!quiet.msg) message("Stage labels are assigned as 'mix' when multiple stages occur within the same bin; otherwise, the stage is preserved.")
    }
  }
  if (is.function(stage.label)) {
    stage.info <- df.binned %>%
      summarise(!!stage_sym := stage.label(!!stage_sym),.groups = "keep")
    df <- dplyr::left_join(df, stage.info, by = c(sbj_col, trial_col, name.bin))
    df <- dplyr::relocate(df, !!stage_sym, .after = !!base_sym)
    if (!quiet.msg) message(paste0("Stage labels are assigned using the rule defined by function ", substitute(stage.label), "."))
  }
  # save a copy of df
  df.copy <- df
  # Assign a representative timestamp for each bin based on provided target stage
  # customize a function for assigning timestamps based on a given stage
  assign_timestamps <- function(df, onset, time_col, stage_col, samp_rate,
                                sbj_col, trial_col, bin_col) {
    time_sym <- rlang::sym(time_col)
    stage_sym <- rlang::sym(stage_col)
    sbj_sym <- rlang::sym(sbj_col)
    trial_sym <- rlang::sym(trial_col)
    bin_sym <- rlang::sym(bin_col)
    interval <- floor(1000/samp_rate)
    df <- df %>% dplyr::ungroup() %>% group_by(!!sbj_sym, !!trial_sym)
    # Check if there are groups that don't contain the target stage specified
    missingInfo <- df %>% dplyr::filter(!any(!!stage_sym == onset)) %>% summarise(.groups = "keep")
    # Filter out these missing trials
    df <- df %>% dplyr::filter((any(!!stage_sym == onset)))
    # Adding timestamps
    df <- df %>% group_modify(~{
      # add a timestamp column
      .x <- dplyr::mutate(.x, !!time_sym := NA_real_)
      # get the first appearing target stage of each group
      idx <- which(.x[[stage_col]] == onset)[1]
      # set the values of the timestamp
      .x[[time_col]][1:idx] <- rev(seq(from = 0, by = -interval, length.out = idx))
      .x[[time_col]][(idx+1):nrow(.x)] <- seq(from = interval, by = interval, length.out = nrow(.x)-idx)
      return(.x)
    })
    df <- df %>% relocate(!!time_sym, .after = !!bin_sym)
    return(list(data = df, missingInfo = missingInfo))
  }
  missingInfo <- assign_timestamps(df, onset = onset, time_col = time_col, stage_col = stage_col,
                                   sbj_col = sbj_col, trial_col = trial_col, bin_col = name.bin,
                                   samp_rate = new_sample_rate)$missingInfo
  df <- assign_timestamps(df, onset = onset, time_col = time_col, stage_col = stage_col,
                          sbj_col = sbj_col, trial_col = trial_col, bin_col = name.bin,
                          samp_rate = new_sample_rate)$data
  if (nrow(missingInfo) == 0 && !is.null(additional.onsets) && !quiet.warning)
    warning("No missing groups to be added. Ignored the additional.onsets argument.")
  if (nrow(missingInfo) > 0) {
    if (is.null(additional.onsets) && !quiet.msg && solve) {
      message(paste0("Detected group(s) that don't contain the target onset ", onset, " you put.\nBelow is a summary:"))
      print(missingInfo)
      message("These groups will be missing in your final output.")
      message("You should check these groups separately and provide\n**one or more** additional onsets for these groups, which will be operated sequentially.")
      message("Ignore these messages if you decide to drop these groups.")
    }
    if (is.null(additional.onsets) && !quiet.warning && solve)
      warning(paste0("Detected group(s) that don't contain the target onset ", onset, " you put.\nCheck the messages above for more details."))
    if (!solve) {
      if (!quiet.warning) warning("There are group(s) that don't include your specified onset stage which are unprocessed.\nDid you unintentionally set append.missing = FALSE?")
      if (!quiet.msg) message("Relevant column info has been added.")
      if (any(is.na(df[[baseline]])) && !quiet.warning)
        warning("Detect NA(s) in the final outcome of baseline pupil size column, which indicates that\nyour orginal baseline pupil size column contains NA(s). Considering set na.rm = TRUE as an argument.\nIf you feed 'method' with a customized function, please make sure you handle NA(s) properly when devising it.")
      if (!quiet.msg) message("Finish processing!")
      return(tibble::as_tibble(dplyr::ungroup(df)))
    }
    # Start appending missing groups base on the additional onsets
    if (!is.null(additional.onsets)){
      if(!quiet.debug.msg) message("Start recovering missing trials based on additional onset(s) you provide...")
      current_missing_info <- missingInfo
      df_recover <- data.frame()
      for (i in 1: length(additional.onsets)) {
        if (!quiet.debug.msg) message(paste0("Working on additional stage "), additional.onsets[i],"...")
        current_missing_df <- df.copy %>% group_by(!!sbj_sym, !!trial_sym) %>%
          dplyr::filter(!any(!!stage_sym == onset))
        df_to_add <- assign_timestamps(current_missing_df, onset = additional.onsets[i], time_col = time_col, stage_col = stage_col,
                                       sbj_col = sbj_col, trial_col = trial_col, bin_col = name.bin,
                                       samp_rate = new_sample_rate)$data
        current_missing_info <- assign_timestamps(current_missing_df, onset = additional.onsets[i],
                                                  time_col = time_col, stage_col = stage_col,
                                                  sbj_col = sbj_col, trial_col = trial_col,
                                                  bin_col = name.bin,samp_rate = new_sample_rate)$missingInfo
        if (nrow(current_missing_info) > 0 && !quiet.debug.msg)
          message(paste0("Detected group(s) that don't include ", additional.onsets[i],"."))
        df_recover <- dplyr::bind_rows(df_recover, df_to_add)
        if (!quiet.debug.msg) message(paste0("Finish additional onset stage "), additional.onsets[i],"...")
        if (nrow(current_missing_info) == 0) {
          if (i < length(additional.onsets) && !quiet.debug.warning)
            warning(paste0("No missing groups anymore using your first ", i," additional onset(s), Ignore redundant additional stage(s) you provided."))
          break
        }
      }
      if (nrow(current_missing_info) > 0) {
        if (!quiet.debug.msg) {
          if (!quiet.debug.msg) message("There are still group(s) that don't contain both your target onset stage and all of your additional onset stage(s).\nBelow is the info:")
          print(current_missing_info)
          if (!quiet.debug.msg) message("Check these groups and add more additional.onsets and run the function again.\nIgnore the message if you decide to drop these groups.")
        }
        if (!quiet.debug.warning) warning("Detected group(s) that are still missing stage(s). Check the message above for more info.\nIgnore this warning if you decide to drop these groups.")
      }
      if (nrow(current_missing_info) == 0 && !quiet.debug.msg)
        message("All groups in the original data are preserved.")
      df <- df %>% dplyr::bind_rows(df_recover) %>%
        dplyr::arrange(!!sbj_sym, !!trial_sym, !!bin_sym)
      if (!quiet.msg) message("Relevant column info has been added.")
      # Check for potential NAs
      if (any(is.na(df[[baseline]])) && !quiet.warning)
        warning("Detect NA(s) in the final outcome of baseline pupil size column, which indicates that\nyour orginal baseline pupil size column contains NA(s). Considering set na.rm = TRUE as an argument.\nIf you feed 'method' with a customized function, please make sure you handle NA(s) properly when devising it.")
      if (!quiet.msg) message("Finish processing!")
      return(tibble::as_tibble(dplyr::ungroup(df)))
    }
    if (!quiet.msg) message("Relevant column info has been added.")
    # Check for potential NAs
    if (any(is.na(df[[baseline]])) && !quiet.warning)
      warning("Detect NA(s) in the final outcome of baseline pupil size column, which indicates that\nyour orginal baseline pupil size column contains NA(s). Considering set na.rm = TRUE as an argument.\nIf you feed 'method' with a customized function, please make sure you handle NA(s) properly when devising it.")
    if (!quiet.msg) message("Finish processing!")
    return(tibble::as_tibble(dplyr::ungroup(df)))
  }
  if (!quiet.msg) message("Relevant column info has been added.")
  # Check for potential NAs
  if (any(is.na(df[[baseline]])) && !quiet.warning)
    warning("Detect NA(s) in the final outcome of baseline pupil size column, which indicates that\nyour orginal baseline pupil size column contains NA(s). Considering set na.rm = TRUE as an argument.\nIf you feed 'method' with a customized function, please make sure you handle NA(s) properly when devising it.")
  if (!quiet.msg) message("Finish processing!")
  return(tibble::as_tibble(dplyr::ungroup(df)))
}
