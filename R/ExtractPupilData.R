#' @export
ExtractPupilData <- function (start.index, end.index, rawlines) {
  # extract raw data lines
  data.lines <- rawlines[(start.index+1):(end.index-1)]
  # keep only lines that start with a digit
  data.lines <- data.lines[grepl("^\\d", data.lines)]
  # check if no data has been found
  if (length(data.lines) == 0) {
    warning("No sample data found between these indices. If you believe this is an error, tell Andy to check it.")
    tempdf <- data.frame("timestamp" = numeric(0), "pupil_size" = numeric(0))
    return(tempdf)
  }
  else {
    # create a temporary file in the working directory to read the data back in
    tmpfile.path <- file.path(getwd(),"cleaned_sample.txt")
    writeLines(data.lines, tmpfile.path)
    # transfer into a data frame, but we only keep the time stamp and pupil size.
    tempdf <- data.frame(fread(tmpfile.path))[,c(1,4)]
    colnames(tempdf) <- c("timestamp", "pupil_size")
    # delete the temporary file
    unlink(tmpfile.path)
    return(tempdf)
  }
}
