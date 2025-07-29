#' @export
BaselinePupil <- function (df, target.stage, howlong,sbj_col, pupil_col,
                           trial_col, stage_col, time_col, name = "baseline",
                           additional.stages = NULL,
                           direction = c("forward", "backward"),
                           value = c("mean", "median", "rms"),
                           samp_rate = 500, get_missing_groups = 0,
                           get_subtract_value = FALSE,
                           quiet.msg = FALSE, quiet.warning = FALSE,
                           quiet.debug.msg = FALSE, quiet.debug.warning = FALSE,...) {
  # save a copy of original df
  df.copy <- df
  # define a root mean square function
  rms <- function (x,...) return(sqrt(mean(x^2,...)))
  # save a copy of value
  method <- if (is.function(value) || is.numeric(value)) value
  samp_num <- round(howlong/(1000/samp_rate))
  direction <- match.arg(direction)
  sbj_sym <- rlang::sym(sbj_col)
  name_sym <- rlang::sym(name)
  stage_sym <- rlang::sym(stage_col)
  pupil_sym <- rlang::sym(pupil_col)
  trial_sym <- rlang::sym(trial_col)
  time_sym <- rlang::sym(time_col)
  group_list <- list(sbj_sym,trial_sym)
  group_vec <- c("subject", "trial")
  solve <- if (isTRUE(list(...)$append.missing) || is.null(list(...)$append.missing)) TRUE else FALSE
  df <- df %>% group_by(!!!group_list)
  if (get_missing_groups %% 1 > 0) stop("get_missing_groups must be an integar.")
  if (get_missing_groups < 0) stop("get_missing_groups must not be negative.")
  if (get_missing_groups > 0 && is.null(additional.stages)) {
    get_missing_groups <- 1
    if (!quiet.warning) warning("No aditional stages defined. get_missing_groups is forced to be 1.")
  }
  if (get_missing_groups > 1) {
    n <- min((get_missing_groups-1),length(additional.stages))
    check_list <- c(target.stage,c(additional.stages[1:n]))
    return.info <- df %>% dplyr::filter(!any(!!stage_sym %in% check_list)) %>% summarise(.groups = "keep")
    if (!quiet.warning)
      warning(paste0("You are getting the info of groups that don't contain stage(s): ", paste(check_list, collapse = ", "), ".\nAll parameters except for 'df', 'target.stage', 'additional.stages', 'sbj_col', 'stage_col' and 'trial_col' are ignored."))
    return(tibble::as_tibble(dplyr::ungroup(return.info)))
  }
  if (get_missing_groups == 1) {
    if (!quiet.warning)
      warning(paste0("You are getting the groups info that don't contain ", target.stage,".\nAll parameters except for 'df', 'target.stage', 'sbj_col', 'stage_col' and 'trial_col' are ignored."))
    return.info <- df %>% dplyr::filter(!any(!!stage_sym == target.stage)) %>% summarise(.groups = "keep")
    return(tibble::as_tibble(dplyr::ungroup(return.info)))
  }
  # Check non-numeric timestamp and pupil size
  if (!is.numeric(df[[pupil_col]])) {
    if (!quiet.msg) message("Converting pupil size column to numeric.")
    df <- dplyr::mutate(df, !!pupil_sym := as.numeric(df[[pupil_col]]))
  }
  if (!is.numeric(df[[time_col]])) {
    if (!quiet.msg) message("Converting timestamp column to numeric.")
    df <- dplyr::mutate(df, !!time_sym := as.numeric(df[[time_col]]))
  }
  if (get_subtract_value && !quiet.msg) message("You are getting the value used to be subtracted from the original pupil size per group.")
  if (is.function(value)) {
    if (!quiet.msg) {
      if (direction == "forward")
        message("Baseline via function: '", substitute(value),"' ", howlong, " ms before ", target.stage, " ends...")
      if (direction == "backward")
        message("Baseline via function: '", substitute(value),"' ", howlong, " ms after ", target.stage, " starts...")
    }
  }
  if (missing("value") || is.character(value)) {
    value <- match.arg(value)
    # save a copy of value
    method <- value
    if (!quiet.msg) {
      if (value == "mean" && direction == "forward")
        message("Baseline via mean pupil size of ", howlong, " ms before ", target.stage, " ends...")
      if (value == "mean" && direction == "backward")
        message("Baseline via mean pupil size of ", howlong, " ms after ", target.stage, " starts...")
      if (value == "median" && direction == "forward")
        message("Baseline via the median value of ", howlong, " ms before ", target.stage, " ends...")
      if (value == "median" && direction == "backward")
        message("Baseline via the median value of ", howlong, " ms after ", target.stage, " starts...")
      if (value == "rms" && direction == "forward")
        message("Baseline via the root mean square of ", howlong, " ms before ", target.stage, " ends...")
      if (value == "rms" && direction == "backward")
        message("Baseline via the root mean square of ", howlong, " ms after ", target.stage, " starts...")
    }
    value <- get(value)
  }
  # Detect if there is no data for a given stage in a group-wise manner
  missingInfo <- df %>% dplyr::filter(!any(!!stage_sym == target.stage)) %>% summarise(.groups = "keep")
  if (!is.null(list(...)$getmissingInfo) && list(...)$getmissingInfo)
    return(tibble::as_tibble(missingInfo))
  if (nrow(missingInfo) > 0) {
    if (!quiet.msg) message(paste0("Detected trial(s) that don't include ", target.stage,"."))
    if (!solve && !quiet.warning) warning("There are trials that don't include your specified stage which are unprocessed.\nDid you unintentionally set append.missing = FALSE?")
    # drop groups that don't have the target stage
    df <- df %>% dplyr::filter(any(!!stage_sym == target.stage))
  }
  if ((nrow(missingInfo) > 0 && is.null(additional.stages)) && solve) {
    if (!quiet.msg) {
      message(paste0("Some groups of data don't contain stage: ",target.stage, ". Below is a summary:"))
      print(missingInfo)
      message("These groups are hence DROPPED for baseline correction. Read carefully for what you should do next:\n")
      message("These groups may not contain your target stage because the target stage of these groups suffer extensive blinks\nthat greater than maxgap you defined during the interpolation part, and hence, they are dropped by the function automatically.\n")
      message("You should find another proper stage to conduct baseline correction for these groups.\nMaybe the stage before or after the current stage you put in the target.stage argument.\n")
      message("You should rerun the function but add additional.stages = 'newstage1'")
      message("If there are still groups that don't have stage 'newtarget1', you will get noticed.\nIf this happens, you can rerun the function again by adding another stage (e.g., newstage2) by using:\n additional.stages = c('newstage1','newstage2')")
      message("Repeat this until there is no missing groups.")
      message("Ignore all these messages if you intend to drop these groups.\n")
    }
    if (!quiet.warning)
      warning("Detected at least 1 group that doesn't contain your specified target.stage. Read the info above carefully.\nIgnore this warning if you intend to drop these groups.")
  }
  if (direction == "forward") {
    if (is.numeric(value)) {
      if (get_subtract_value && !solve)
        return(tibble::as_tibble(dplyr::ungroup(df %>% summarise(value = value, .groups = "keep"))))
      # add a temporary column to the original data frame
      if (!(get_subtract_value || quiet.warning))
        warning("The argument 'value' is set to be a fixed value, all pupil size is subtracted by this value.\nThe direction argument is ignored.")
      # add a temporary column to the original data frame
      df$value <- value
      subtract_value <- df %>% summarise(value = unique(value), .groups = "keep")
    } else {
      # Calculate the subtract value
      tmpdf <- df %>% dplyr::filter(!!stage_sym == target.stage) %>%
        group_modify(~dplyr::mutate(.x, "target.time" = .x[[time_col]][1]))
      vals <- unique(tmpdf$target.time)
      # Propagate the target.time across the full data
      tmpdf <- df  %>% dplyr::mutate("target.time" = vals[cur_group_id()])
      subtract_value <- tmpdf %>%
        dplyr::filter(!!time_sym < target.time) %>%
        group_modify(~ dplyr::slice_tail(.x, n = min(1+ samp_num, nrow(.x)))) %>%
        summarise(value = value(!!pupil_sym,...), .groups = "keep")
      if (get_subtract_value && !solve) {
        if (any(is.na(subtract_value)) && !quiet.warning)
          warning("Since your original dataset contains NA(s), please set rm.na == TRUE as an argument.\n If you feed a function to 'value', please make sure you handle NA(s) properly so that no groups have a baseline correction value of NA.")
        return(tibble::as_tibble(dplyr::ungroup(subtract_value)))
      }
      # merge into the original data frame
      df <- left_join(df, subtract_value, by = group_vec)
    }
    df <- df %>% mutate(!!name_sym := !!pupil_sym - value) %>%
      dplyr::relocate(!!name_sym, .after = !!pupil_sym)
    # drop the temporary column
    df <- df %>% dplyr::select(!value)
    if (!solve) {
      if (!quiet.msg) message("Complete Baseline Correction!")
      if (any(is.na(df[[name]])) && !quiet.warning)
        warning("Since your original dataset contains NA(s), please set rm.na == TRUE as an argument.\n If you feed a function to 'value', please make sure you handle NA(s) properly so that no groups have a baseline correction value of NA.")
      return(tibble::as_tibble(dplyr::ungroup(df)))
    }
  }
  if (direction == "backward") {
    if (is.numeric(value)) {
      if (get_subtract_value && !solve)
        return(tibble::as_tibble(dplyr::ungroup(df %>% summarise(value = value, .groups = "keep"))))
      # add a temporary column to the original data frame
      if (!(get_subtract_value || quiet.warning))
        warning("The argument 'value' is set to be a fixed value, all pupil size is subtracted by this value.\nThe direction argument is ignored.")
      df$value <- value
      subtract_value <- df %>% summarise(value = unique(value), .groups = "keep")
    } else {
      # Calculate the subtract value
      tmpdf <- df %>% dplyr::filter(!!stage_sym == target.stage) %>%
        group_modify(~dplyr::mutate(.x, "target.time" = .x[[time_col]][nrow(.x)]))
      vals <- unique(tmpdf$target.time)
      # Propagate the target.time across the full data
      tmpdf <- df  %>% dplyr::mutate("target.time" = vals[cur_group_id()])
      subtract_value <- tmpdf %>%
        dplyr::filter(!!time_sym > target.time) %>%
        group_modify(~ dplyr::slice_head(.x, n = min(1+ samp_num, nrow(.x)))) %>%
        summarise(value = value(!!pupil_sym,...), .groups = "keep")
      if (get_subtract_value && !solve) {
        if (any(is.na(subtract_value)) && !quiet.warning)
          warning("Since your original dataset contains NA(s), please set rm.na == TRUE as an argument.\n If you feed a function to 'value', please make sure you handle NA(s) properly so that no groups have a baseline correction value of NA.")
        return(tibble::as_tibble(dplyr::ungroup(subtract_value)))
      }
      # merge into the original data frame
      df <- left_join(df, subtract_value, by = group_vec)
    }
    df <- df %>% mutate(!!name_sym := !!pupil_sym - value) %>%
      dplyr::relocate(!!name_sym, .after = !!pupil_sym)
    # drop the temporary column
    df <- df %>% dplyr::select(!value)
    if (!solve) {
      if (!quiet.msg) message("Complete Baseline Correction!")
      if (any(is.na(df[[name]])) && !quiet.warning)
        warning("Since your original dataset contains NA(s), please set rm.na == TRUE as an argument.\n If you feed a function to 'value', please make sure you handle NA(s) properly so that no groups have a baseline correction value of NA.")
      return(tibble::as_tibble(dplyr::ungroup(df)))
    }
  }
  if (nrow(missingInfo) == 0 && !is.null(additional.stages) && !quiet.warning)
    warning("There is no missing group(s) to be added. Ignore the additional.stages argument.")
  # Add back the previous missing trials based on additional stages provided
  if (!is.null(additional.stages)) {
    if(!quiet.debug.msg) message("Start recovering missing trials based on additional target(s) you provide...")
    current_missing_info <- missingInfo
    df_recover <- data.frame()
    for (i in 1:length(additional.stages)) {
      if (nrow(current_missing_info) > 0) {
        if (!quiet.debug.msg) message(paste0("Working on additional stage "), additional.stages[i],"...")
        current_missing_df <- df.copy %>% dplyr::filter(paste(!!sbj_sym, !!trial_sym)
                                                        %in% paste(current_missing_info[[sbj_col]],current_missing_info[[trial_col]]))
        # get the missing trials after using each additional stage provided
        if (is.character(method)) {
          if (!quiet.debug.msg) {
            if (method == "mean" && direction == "forward")
              message("Baseline via mean pupil size of ", howlong, " ms before ", additional.stages[i], " ends...")
            if (method == "mean" && direction == "backward")
              message("Baseline via mean pupil size of ", howlong, " ms after ", additional.stages[i], " starts...")
            if (method == "median" && direction == "forward")
              message("Baseline via the median value of ", howlong, " ms before ", additional.stages[i], " ends...")
            if (method == "median" && direction == "backward")
              message("Baseline via the median value of ", howlong, " ms after ", additional.stages[i], " starts...")
            if (method == "rms" && direction == "forward")
              message("Baseline via the root mean square of ", howlong, " ms before ", additional.stages[i], " ends...")
            if (method == "rms" && direction == "backward")
              message("Baseline via the root mean square of ", howlong, " ms after ", additional.stages[i], " starts...")
          }
        } else {
          if (!quiet.debug.msg) {
            if (direction == "forward")
              message("Baseline via function: '", substitute(value),"' ", howlong, " ms before ", additional.stages[i], " ends...")
            if (direction == "backward")
              message("Baseline via function: '", substitute(value),"' ", howlong, " ms after ", additional.stages[i], " starts...")
          }
        }
        current_missing_info <- BaselinePupil(df = current_missing_df, target.stage = additional.stages[i],
                                              howlong = howlong, sbj_col = sbj_col,
                                              pupil_col = pupil_col, trial_col = trial_col,
                                              stage_col = stage_col,time_col = time_col,
                                              direction = direction, value = method, name = name,
                                              get_subtract_value = get_subtract_value, getmissingInfo = TRUE,
                                              quiet.msg = TRUE, quiet.warning = TRUE,...)
        current_df_recover <- BaselinePupil(df = current_missing_df, target.stage = additional.stages[i],
                                            howlong = howlong, sbj_col = sbj_col,
                                            pupil_col = pupil_col, trial_col = trial_col,
                                            stage_col = stage_col,time_col = time_col,
                                            direction = direction, value = method, name = name,
                                            get_subtract_value = get_subtract_value,
                                            append.missing = FALSE,
                                            quiet.msg = TRUE, quiet.warning = TRUE,...)
        df_recover <- df_recover %>% dplyr::bind_rows(current_df_recover)
        if (!quiet.debug.msg) message(paste0("Finish additional stage "), additional.stages[i],"...")
      }
      if (nrow(current_missing_info) > 0 && !quiet.debug.msg)
        message(paste0("Detected trial(s) that don't include ", additional.stages[i],"."))
      if (nrow(current_missing_info) == 0) {
        if (i < length(additional.stages) && !quiet.debug.warning)
          warning(paste0("No missing groups anymore using your first ", i," additional stage(s), Ignore redundant additional stage(s) you provided."))
        break
      }
    }
    if (nrow(current_missing_info) > 0) {
      if (!quiet.debug.msg) message("There are still group(s) that don't contain both your target stage and all of your additional stage(s).\nBelow is the info:")
      print(current_missing_info)
      if (!quiet.debug.msg) message("Check these groups and add more additional.stages and run the function again.\nIgnore the message if you decide to drop these groups.")
      if (!quiet.debug.warning) warning("Detected group(s) that are still missing stage(s). Check the message above for more info.\nIgnore this warning if you decide to drop these groups.")
    }
    if (nrow(current_missing_info) == 0 && !quiet.debug.msg) {
      if (get_subtract_value)
        message("The subtract values for ALL groups in the original data is calculated.") else
          message("All groups in the original data are preserved and baselined.")
    }
    if (get_subtract_value) {
      full_subtract_value <- subtract_value %>% dplyr::bind_rows(df_recover) %>%
        dplyr::arrange(!!sbj_sym, !!trial_sym) %>% dplyr::ungroup()
      if (any(is.na(full_subtract_value)) && !quiet.warning)
        warning("Since your original dataset contains NA(s), please set rm.na == TRUE as an argument.\n If you feed a function to 'value', please make sure you handle NA(s) properly so that no groups have a baseline correction value of NA.")
      return(tibble::as_tibble(full_subtract_value))
    }
    full_df <- df %>% dplyr::bind_rows(df_recover) %>%
      dplyr::arrange(!!sbj_sym, !!trial_sym) %>% dplyr::ungroup()
    if (!quiet.msg) message("Complete Baseline Correction!")
    if (any(is.na(full_df[[name]])) && !quiet.warning)
      warning("Since your original dataset contains NA(s), please set rm.na == TRUE as an argument.\n If you feed a function to 'value', please make sure you handle NA(s) properly so that no groups have a baseline correction value of NA.")
    return(tibble::as_tibble(full_df))
  }
  if (get_subtract_value) {
    if (any(is.na(subtract_value)) && !quiet.warning)
      warning("Since your original dataset contains NA(s), please set rm.na == TRUE as an argument.\n If you feed a function to 'value', please make sure you handle NA(s) properly so that no groups have a baseline correction value of NA.")
    return(tibble::as_tibble(dplyr::ungroup(subtract_value)))
  }
  if (!quiet.msg) message("Complete Baseline Correction!")
  if (any(is.na(df[[name]])) && !quiet.warning)
    warning("Since your original dataset contains NA(s), please set rm.na == TRUE as an argument.\n If you feed a function to 'value', please make sure you handle NA(s) properly so that no groups have a baseline correction value of NA.")
  return(tibble::as_tibble(dplyr::ungroup(df)))
}
