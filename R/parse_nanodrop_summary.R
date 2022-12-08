#' parse_nanodrop_summary
#'
#' This function allows you to easily read in the summary csv file exported by the nanodrop one
#' For some reason, the Nanodrop outputs the summary csv file in a weird format.
#' It is actually tab separated (instead of comma) and is UTF-16 encoded.
#' This is all stuff I tend to forget, so here is a function to easily read this data in.
#'
#' @param filepath this is the file path to the csv that was exported from the nandrop.
#'
#' @return a tibble with the first 7 columns (only important ones) from the file.
#' @export
#'
#' @examples
#' data = parse_nanodrop_summary("tests/testdata/parse_nanodrop_summary/nanodrop_summary.csv")
parse_nanodrop_summary = function(filepath){
  out = readr::read_tsv(filepath,locale = locale(encoding = "UTF-16")) %>%
    #selecting only the important columns
    dplyr::select(1:7)
  return(out)
}
