#' @export
SavitzkyGolayFilter <- function (x, pupil_col = NULL, p = 3, window = p + 3 - p%%2,
                                 m = 0, ts = 1, return_tibble = FALSE) {
  if (!requireNamespace("signal", quietly = TRUE))
    stop('The "signal" package is required for Savitzky-Golay filtering. Please install it with:\ninstall.packages("signal")')
  if (window %% 2 == 0) stop("Window size must set to be an odd number.")
  if (is.numeric(x)) {
    if (!missing(return_tibble) && return_tibble == TRUE)
      message("Input is a vector, ignore 'return_tibble' argument.")
    x <- signal::sgolayfilt(x, p, n = window, m = m, ts = ts)
    return(as.numeric(x))
  }
  if (is.data.frame(x) || tibble::is.tibble(x)) {
    data <- x[[pupil_col]]
    if (!is.numeric(data)) {
      stop(paste0("Column \'", pupil_col, "\' must be numeric."))
    }
    data <- signal::sgolayfilt(data, p, n = window, m = m, ts = ts)
    x <- x %>% mutate(!!rlang::sym(pupil_col) := data)
    ifelse(return_tibble, return(tibble::as_tibble(x)), return(as.numeric(x[[pupil_col]])))
  }
}
