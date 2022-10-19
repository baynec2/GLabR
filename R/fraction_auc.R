#' fraction_auc
#' this function is designed to calculate the auc corresponding to each of the fractions from the . It is my hope that this will provide some sort of peptide quantification that can be useful when determining how much to shoot on the MS (unproven hypothesis)
#'
#' Note that this function has a few method specific details hard coded into it as these don't really change in the Gonzalez lab set up. They are as follows:
#'
#' 1. A + B rows are concatenated into a single fraction for each . ie A1 + B1 get combined, A2 + B2 get combined, and so on.
#' 2. The collector is set up to move column wise A1 -> A12 followed by B1 -> B12- then repeat.
#' 3. The time of the method isnt evenly divisible by 37 sec, this means there is an extra amount of time where it seems like A1 gets added to at the end of the method.
#'
#' @param parsed_hplc_chromatography_data this is the parsed hplc chromatography data, designed to be output of hplc_parser.
#' @param start this is the time in seconds when the method starts default 660 = 11 min
#' @param fraction_interval: this is the time between fraction collection intervals. Default = 37 sec.
#'
#' @return a tibble with auc calculated for each fraction.
#' @export
#'
#' @examples
#' parsed_hplc_chromatography_data = hplc_parser("tests/testdata/hplc_parser/CB004")
#' fraction_auc = fraction_auc(parsed_hplc_chromatography_data)
fraction_auc = function(parsed_hplc_chromatography_data, start = 660, fraction_interval = 37){
  #shorten lengthy name
  data = parsed_hplc_chromatography_data
  #Converting to min to match hplc export
  start_min = start/60
  fraction_interval_min = fraction_interval/60
  # filtering to only inclue time fractions are collected
  data = data %>%
    dplyr::filter(time_min >= start_min) %>%
    dplyr::mutate(cycle_number =  ceiling((time_min-11)/ (fraction_interval_min)))
  #Creating well names
  rows =rep(LETTERS[1:2],each = 12)
  cols = rep(1:12,2)
  well_order = paste0(rows,cols)

  #Determining how many
  n_iterations = data %>%
    dplyr::group_by(cycle_number) %>%
    dplyr::summarise(n = n())
  #Cnverting well_order to well id. note there seems to be an extra iteration, hardcoding it here as A1
  n_iterations$well = c(rep(well_order,nrow(n_iterations)/length(well_order)),"A1")
  # Joining with original data. Hardcoding to keep it explicit as to what samples we concationate
  data2 = dplyr::inner_join(data,n_iterations, by = "cycle_number") %>%
    dplyr::mutate(fraction = dplyr::case_when(well %in% c("A1","B1") ~ 1,
                                              well %in% c("A2","B2") ~ 2,
                                              well %in% c("A3","B3") ~ 3,
                                              well %in% c("A4","B4") ~ 4,
                                              well %in% c("A5","B5") ~ 5,
                                              well %in% c("A6","B6") ~ 6,
                                              well %in% c("A7","B7") ~ 7,
                                              well %in% c("A8","B8") ~ 8,
                                              well %in% c("A9","B9") ~ 9,
                                              well %in% c("A10","B10") ~ 10,
                                              well %in% c("A11","B11") ~ 11,
                                              well %in% c("A12","B12") ~ 12,
                                              TRUE ~ NA_real_))
  #calculating AUC
  auc = data2 %>%
    dplyr::group_by(fraction,plex) %>%
    dplyr::summarise(auc = MESS::auc(time_min,mau_value,type = "linear")) %>%
    dplyr::ungroup()
  return(auc)
}
