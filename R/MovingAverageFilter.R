#' @export
MovingAverageFilter <- function(x, pupil_col = NULL, window = 5, return_tibble = FALSE,
                                tofill = TRUE, align = c("center", "left", "right")) {
  if (window %% 2 == 0) stop("Window size must set to be an odd number.")
  align <- match.arg(align)
  filler <- if (tofill) "extend" else NA
  if (is.numeric(x)) {
    if (!missing(return_tibble) && return_tibble == TRUE)
      message("Input is a vector, ignore 'return_tibble' argument.")
    x <- zoo::rollmean(x, k = window, fill = filler, align = align)
    return(as.numeric(x))
  }
  if (is.data.frame(x) || tibble::is.tibble(x)) {
    data <- x[[pupil_col]]
    if (!is.numeric(data)) {
      stop(paste0("Column \'", pupil_col, "\' must be numeric."))
    }
    data <- zoo::rollmean(data, k = window, fill = filler, align = align)
    x <- x %>% mutate(!!rlang::sym(pupil_col) := data)
    ifelse(return_tibble, return(tibble::as_tibble(x)), return(as.numeric(x[[pupil_col]])))
  }
}
