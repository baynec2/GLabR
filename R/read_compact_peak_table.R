#'read_compact_peak_table
#'easily read in data from compact peak table export from Tapestation.
#' @param filepath
#'
#' @return a tibble that has been processed to remove
#' @export
#'
#' @examples
#' data = read_compact_peak_table("tests/read_compact_peak_table/compactPeakTable.csv")
read_compact_peak_table = function(filepath){
  data = readr::read_csv(filepath,skip = 1,
         col_names = c("filename",
                       "well",
                       "sample_description",
                       "size",
                       "calibrated_conc_ng_uL",
                       "assigned_conc_ng_uL",
                       "integrated_area",
                       "peak_molarity_nmol_l",
                       "From_bp_to_bp",
                       "peak_comment",
                       "observations")) %>%
  #removing the 100 bp reference peak and ladder
  dplyr::filter(size != 100,
                sample_description != "Ladder") %>%
  #just saying anything > 60kb is 60kb
  dplyr::mutate(size = as.numeric(gsub(">","",size)))

  return(data)
}
