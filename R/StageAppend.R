#' @export
StageAppend <- function (tmpdf, start.index, end.index, rawlines, label) {
  df <- ExtractPupilData(start.index, end.index, rawlines)
  if (nrow(df)>0) {
    df$stage <- label
    tmpdf <- dplyr::bind_rows(tmpdf,df)
  }
  return(tmpdf)
}
