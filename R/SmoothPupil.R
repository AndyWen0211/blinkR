#' @export
SmoothPupil <- function(df, family = c("average","loess","Savitzky","Gaussian"),
                        grouping = c("both","subject","trial","none"), sbj_col, pupil_col,
                        trial_col, weight = FALSE,...) {
  family <- match.arg(family)
  grouping <- match.arg(grouping)
  sbj_sym <- rlang::sym(sbj_col)
  pupil_sym <- rlang::sym(pupil_col)
  trial_sym <- rlang::sym(trial_col)
  # grouping data frame
  group_list <- switch(grouping,
                       "none" = NULL,
                       "subject" = list(sbj_sym),
                       "trial" = list(trial_sym),
                       "both" = list(sbj_sym,trial_sym)
  )
  df <- df %>% group_by(!!!group_list)
  if (weight) {
    args <- list(...)
    sides <- if (!is.null(args$sides)) args$sides else 2
    window <- if (!is.null(args$window)) args$window else 11
    tofill <- if(!is.null(args$tofill)) args$tofill else TRUE
    if (!is.null(args$return_tibble)) {
      warning("Please don't set a 'return_tibble' argument within SmoothPupil.\nIn this case, your 'return_tibble' argument is ignored and forced to be FALSE.")
      args$return_tibble <- FALSE
    }
    warning("Ignoring the 'family' argument - weight is TRUE, hanning filter is automatically applied.")
    message("You are applying a hanning filter...\n")
    message(paste0("Apply sides = ", sides, ".\nInclude sides in the argument to change the parameter.\nRun ?stats::filter for more info on this parameter.\n"))
    message(paste0("Apply window size = ", window, ".\nInclude window in the argument to change the window size.\n"))
    if (tofill) {
      message("Edge values are automatically filled. Set 'tofill = FALSE' to change the behavior.\n")
    }
    df <- df %>% mutate(!!pupil_sym := HanningFilter(!!pupil_sym, !!!args))
    return(as_tibble(df))
  }
  if (family == "average") {
    args <- list(...)
    sides <- if (!is.null(args$sides)) args$sides else 2
    window <- if (!is.null(args$window)) args$window else 5
    tofill <- if(!is.null(args$tofill)) args$tofill else TRUE
    if (!is.null(args$return_tibble)) {
      warning("Please don't set a 'return_tibble' argument within SmoothPupil.\nIn this case, your 'return_tibble' argument is ignored and forced to be FALSE.")
      args$return_tibble <- FALSE
    }
    message("You are applying a moving average filter...\n")
    message(paste0("Apply sides = ", sides, ".\nInclude sides in the argument to change the parameter.\nRun ?stats::filter for more info on this parameter.\n"))
    message(paste0("Apply window size = ", window, ".\nInclude window in the argument to change the window size.\n"))
    if (tofill) {
      message("Edge values are automatically filled. Set 'tofill = FALSE' to change the behavior.\n")
    }
    df <- df %>% mutate(!!pupil_sym := MovingAverageFilter(!!pupil_sym,!!!args))
    return(as_tibble(df))
  }
  if (family == "Savitzky") {
    args <- list(...)
    p <- if (!is.null(args$p)) args$p else 3
    window <- if (!is.null(args$window)) args$window else p + 3 - p%%2
    if (!is.null(args$return_tibble)) {
      warning("Please don't set a 'return_tibble' argument within SmoothPupil.\nIn this case, your 'return_tibble' argument is ignored and forced to be FALSE.")
      args$return_tibble <- FALSE
    }
    message("You are applying a Savitzky Golay smoothing filter...\n")
    message(paste0("Apply polynomial term p = ", p, ".\nInclude p in the argument to change the polynomial term.\n"))
    message(paste0("Apply window size = ", window, ".\nInclude window in the argument to change the window size.\n"))
    message("Edge values are automatically filled.\n")
    df <- df %>% mutate(!!pupil_sym := SavitzkyGolayFilter(!!pupil_sym,!!!args))
    return(as_tibble(df))
  }
  if (family == "loess") {
    args <- list(...)
    if (is.null(args$time_col)) stop('You are applying a loess filter, please indicate the time column using "time_col = "')
    weights <- args$weights
    if (!is.null(args$return_tibble)) {
      warning("Please don't set a 'return_tibble' argument within SmoothPupil.\nIn this case, your 'return_tibble' argument is ignored and forced to be FALSE.")
      args$return_tibble <- FALSE
    }
    span <- if (!is.null(args$span)) args$span else 0.75
    if (is.null(weights)) message("You are applying an un-weighted Loess filter...\n")
    else message("You are applying a weighted Loess filter...\n")
    message(paste0("Apply span = ", span, "\nLarger span means a better global trend fit but worse local pattern fit.\n"))
    message("Edge values are automatically filled.\n")
    df <- df %>% mutate(!!pupil_sym := do.call("LoessFilter", args = c(list(df = pick(everything()),
                                                                            pupil_col = rlang::as_string(pupil_sym)),!!!args)))
    return(as_tibble(df))
  }
  if (family == "Gaussian") {
    args <- list(...)
    sides <- if (!is.null(args$sides)) args$sides else 2
    window <- if (!is.null(args$window)) args$window else 11
    tofill <- if(!is.null(args$tofill)) args$tofill else TRUE
    if (!is.null(args$return_tibble)) {
      warning("Please don't set a 'return_tibble' argument within SmoothPupil.\nIn this case, your 'return_tibble' argument is ignored and forced to be FALSE.")
      args$return_tibble <- FALSE
    }
    message("You are applying a Gaussian Kernel filter...\n")
    message(paste0("Apply sides = ", sides, ".\nInclude sides in the argument to change the parameter.\nRun ?stats::filter for more info on this parameter.\n"))
    message(paste0("Apply window size = ", window, ".\nInclude window in the argument to change the window size.\n"))
    if (tofill) {
      message("Edge values are automatically filled. Set 'tofill = FALSE' to change the behavior.\n")
    }
    df <- df %>% mutate(!!pupil_sym := GaussianKernelFilter(!!pupil_sym,!!!args))
    return(as_tibble(df))
  }
}
