#' @export
ExtractColumnInfo <- function (singleline, sep) {
  splitline <- strsplit(singleline, sep)[[1]]
  info <- splitline[length(splitline)]
  return(info)
}
