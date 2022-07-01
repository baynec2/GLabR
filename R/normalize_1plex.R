#' normalize_1plex : Normalize values when only a single TMT plex is ran. This function normalizes to the median of the sample/overall median.
#'
#' @param data : a tibble containing combined psms from an experiment with only one plex and no bridge chennel.
#' @return a tibble containing the normalized data.
#' @export
#'
#' @examples
#' Here we are analyzing the results of the 1plex results.
#' dnorm = normalize_1plex("../")
normalize_1plex = function(data){
  #Getting average for each protein
  protein_avg = data %>%
    dplyr::group_by(ProteinID) %>%
    dplyr::summarise(protein_avg = mean(value))

  #Getting the median of the averaged proteins
  median_protein_avg = protein_avg %>%
    dplyr::pull(protein_avg) %>%
    median()

  #Getting the overall median for the entire dataset.
  median_all = data %>%
    pull(value) %>%
    median()

  #Appending median from each sample/plex combination, as well as the average or each protein and then performing normalization by dividing value by column (plex) median
  output = data %>%
    dplyr::inner_join(protein_avg,by = c("ProteinID")) %>%
    dplyr::mutate(intermediate_norm = value/(protein_avg/median_protein_avg))

  #Making the intermediate normalization

  medians = output %>%
    dplyr::group_by(Sample,TMT) %>%
    dplyr::summarise(median_of_sample_plex = median(intermediate_norm))

  #finishing the normalization
  output2 = dplyr::inner_join(output,medians,by = c("Sample","TMT")) %>%
    dplyr::mutate(final_norm = intermediate_norm / (median_of_sample_plex/median_all))

  return(output2)

  }
