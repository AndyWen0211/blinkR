#' @export
GaussianKernelFilter <- function(x, pupil_col = NULL, window = 11,
                                 sigma = window/6, tofill = TRUE,
                                 return_tibble = FALSE,
                                 method = c("convolution", "recursive"),
                                 sides = 2) {
  method <- match.arg(method)
  if (window %% 2 == 0) stop("Window size must set to be an odd number.")
  center <- (window-1)/2
  kernel <- seq(-center, center, 1)
  weights <- exp(-0.5 * (kernel / sigma)^2)
  weights <- weights / sum(weights)
  # Apply weights to the data
  if (is.numeric(x)) {
    if (!missing(return_tibble) && return_tibble == TRUE)
      message("Input is a vector, ignore 'return_tibble' argument.")
    x <- stats::filter(x, weights, method = method, sides = sides)
    # The filter automatically leave the two edges as NA, so we interpolate these values.
    if (tofill) x <- zoo::na.approx(x, rule = 2)
    return(as.numeric(x))
  }
  if (is.data.frame(x) || tibble::is.tibble(x)) {
    data <- x[[pupil_col]]
    if (!is.numeric(data)) {
      stop(paste0("Column \'", pupil_col, "\' must be numeric."))
    }
    data <- stats::filter(data, weights, method = method, sides = sides)
    if (tofill) data <- zoo::na.approx(data, rule = 2)
    x <- x %>% mutate(!!rlang::sym(pupil_col) := data)
    ifelse(return_tibble, return(tibble::as_tibble(x)), return(as.numeric(x[[pupil_col]])))
  }
}
