#' normalize_1plex : Normalize values when only a single TMT plex is ran. This function normalizes to the median of the sample/overall median.
#'
#' @param data : a tibble containing combined psms from an experiment with only one plex and no bridge chennel.
#' @param data_format : This is a character string specifying whether you would like the data in long or wide format. Note that if data_format ==  wide only data
#' corresponding to final_norm is returned. The intermediate normalization steps are disregarded in this case.
#' @return a tibble containing the normalized data.
#' @export
#'
#' @examples
#' data = read_delim("tests/testdata/combine_psm_fractions/PCB002_PSMs_Proteodiscover_output.txt") %>%
#' combine_psm_fractions() %>%
#' normalize_1plex()
normalize_1plex = function(data,data_format = "long"){
  #Getting average for each protein
  protein_avg = data %>%
    dplyr::group_by(ProteinID) %>%
    dplyr::summarise(protein_avg = mean(value,na.rm = T))

  #Getting the median of the averaged proteins
  median_protein_avg = protein_avg %>%
    dplyr::pull(protein_avg) %>%
    median(na.rm = T)

  #Getting the overall median for the entire dataset.
  median_all = data %>%
    dplyr::pull(value) %>%
    median(na.rm = T)

  #Appending median from each sample/plex combination, as well as the average or each protein and then performing normalization by dividing value by column (plex) median
  output = data %>%
    dplyr::inner_join(protein_avg,by = c("ProteinID")) %>%
    dplyr::mutate(intermediate_norm = value/(protein_avg/median_protein_avg))

  #Making the intermediate normalization

  medians = output %>%
    dplyr::group_by(Sample,TMT) %>%
    dplyr::summarise(median_of_sample_plex = median(intermediate_norm,na.rm = T))

  #finishing the normalization
  output2 = dplyr::inner_join(output,medians,by = c("Sample","TMT")) %>%
    dplyr::mutate(final_norm = intermediate_norm / (median_of_sample_plex/median_all))

  #Adding option to export data in long or wide format
  if(data_format == "long"){
    return(output2)
  }else if(data_format == "wide"){
    output3 = output2 %>%
      dplyr::select(Sample,TMT,ProteinID,final_norm) %>%
      tidyr::pivot_wider(names_from = c("Sample","TMT"), values_from = final_norm)
    return(output3)
  }else{
    print("data_format must be either long or wide")
  }
}

