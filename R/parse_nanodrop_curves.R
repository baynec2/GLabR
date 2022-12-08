#' parse_nanodrop_curves
#'
#' This function allows you to easily read in the wavelength/intensity tsv file exported by the Nanodrop One
#' For some reason, the Nanodrop outputs this file in a completely unusable format.
#' This function parses that file and extracts the usefull information.
#'
#' @param filepath
#'
#' @return a tibble with the following columns: sample_id, wavelength, and values.
#' @export
#'
#' @examples
#' data = parse_nanodrop_curves("tests/testdata/parse_nanodrop_cruves/nanodrop_curves.tsv")
parse_nanodrop_curves = function(filepath){

  #Reading in data
  col = readr::read_tsv(filepath,col_names = c("col")) %>%
    dplyr::select(col = 1) %>%
    dplyr::mutate(check =grepl(".*Absorbance.*",col),
                  start = grepl("220.0.*",col))

  #Determining how many samples there are
  nsample = sum(col$check)

  #Determining how many rows there are
  nrow = nrow(col)

  #Determining how many rows per sample there are
  row_per_sample = nrow/nsample

  #Determining the index of where the sample names are
  indexes = c(1,1:(nsample-1) * row_per_sample + 1)

  #Preparing a vector of sample names
  sample_names = col[indexes,1] %>%
    dplyr::pull()

  #Determining when the data actually starts, 13 rows down from sample name. shifts 2 down after every sample.
  data_start = which(col$start) - seq(1,-(nsample*2),by = -2)[1:nsample]

  data  = purrr::map_df(data_start,~readr::read_tsv(filepath,skip = .x, n_max = 261,col_names = c("wavelength","value")))

  #addding names
  sample_id = rep(sample_names,each = 261)

  out = tibble::tibble(sample_id,data)

  return(out)
}
