#' split_well_col
#' This function will take a well coordinate (formatted like A01) and return the col ie 01.
#' @param well this is a vector of well coordinates (Specified as LetterNumber)
#'
#' @return a character vector
#' @export
#'
#' @examples
#' wells = c("A1","A2")
#' cols = split_well_col(wells)
split_well_col= function(well){
  out = as.vector(strsplit(well, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=TRUE))
  out = sapply(out,"[[",2)
  return(out)
}
