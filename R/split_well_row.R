#' split_well_row
#' This function will take a well coordinate (formatted like A01) and return the row ie A.
#'
#' @param well this is a vector of well coordinates (Specified as LetterNumber)
#'
#' @return a character vector
#' @export
#'
#' @examples
#' wells = c("A1","A2")
#' rows = split_well_row(wells)
split_well_row= function(well){
  out = as.vector(strsplit(well, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=TRUE))
  out = sapply(out,"[[",1)
  return(out)
}
