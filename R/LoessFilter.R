#' @export
LoessFilter <- function(df, time_col, pupil_col, weights = NULL,
                        na.action = c("na.exclude", "na.omit", "na.fail", "na.pass"),
                        span = 0.75, return_tibble = FALSE, degree = 2,
                        method = c("symmetric", "gaussian"),...) {
  na.action <- match.arg(na.action)
  method <- match.arg(method)
  weights <- if (missing(weights)) NULL else weights
  data <- predict(stats::loess(formula = as.formula(paste0(pupil_col,"~",time_col)),
                               data = df, weights = weights, na.action = get(na.action),
                               span = span, degree = degree, family = method))
  df <- df %>% dplyr::mutate(!!rlang::sym(pupil_col) := data)
  if (return_tibble) return(tibble::as_tibble(df))
  return(as.numeric(df[[pupil_col]]))
}
