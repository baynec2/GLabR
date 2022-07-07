#' normalize_to_bridge : Normalize PSM combine results to the bridge channel (batch correction).
#'
#'
#' @param data : This represents a data fame containing the PSMs. This is intended to be the output of combine_psm_fractions().
#' As such, the columns of this dataframe must follow certain naming conventions for the code to work. They are listed below.
#'
#' There must be a column named Sample, this contains the name of the plex that was run.
#' There must also be a column named TMT, this contains the isobaric tag that was used.
#' There must be a column named ProteinID that has an ID that corresponds to a protein identifier.
#' There must be a column named value that has all of the intensity values.
#'
#' @param bridge_channel_plex : This is the isobaric channel that contains the bridge channel. Defaults to 126 as that is often used in the lab.
#'
#' @param data_format : This is a character string specifying whether you would like the data in long or wide format. Note that if data_format ==  wide only data
#' corresponding to final_norm is returned. The intermediate normalization steps are disregarded in this case.
#'
#' @return: a tibble containing the normalized results.
#' There are a number of different normalization steps here.
#'
#' intermediate_norm = these values represent the intensity values for each sample where each value is divided by the corresponding bridge channel
#' value that has also been divided by the overall median of all the total bridge channel data.
#'
#' value = contains the original, unnormalized value reported in the exported thermodiscover file.
#' final_norm = these values represent the intermediate_norm values that have been divided by the median of each sample/plex combination which have
#' also been divided by the median of the overall data (discounting the bridge channel)
#'
#' These results can be further normalized using la_box_cox_norm(), which is leigh-ana's method for making the data even more normally distributed.
#' @export
#'
#' @examples
normalize_to_bridge = function(data = dataframe, bridge_channel_plex = "126",data_format = "long"){
#The normalization is done as follows:
  #First we need to normalize the bridge channel by dividing each value by the median of the overall values within the bridge
  #Then we need to normalize our data to the normalized bridge values. We do this by dividing the values by the normalized bridge values
  #Lastly, we need to normalize all of the plexes. We can do this by dividing the value that has been normalized to the bridge channel plex by the column medians (for each plex) and then the overall median.
  Bridge_median = data %>%
    dplyr::filter(TMT == bridge_channel_plex) %>%
    dplyr::pull() %>%
    median(na.rm=TRUE)

  #filtering for just the bridge values so we can append them later
  bridge_values = data %>%
    dplyr::filter(TMT == bridge_channel_plex) %>%
    dplyr::select(-TMT) %>%
    dplyr::rename(bridge_values = value)

  #appending bridge values. For each sample/TMT combination the
  norm_bridge = dplyr::inner_join(data,bridge_values,by = c("ProteinID","Sample")) %>%
    dplyr::mutate(intermediate_norm = value / (bridge_values /Bridge_median)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(TMT != bridge_channel_plex)

  #here data has been normalized to bridge. This matches the results from subset2 in original code

  #Here we are pulling out the data from all samples/plex combinations and taking the median for each.
  median_overall = norm_bridge %>%
    dplyr::pull(intermediate_norm) %>%
    median(na.rm=T)

  # medians per sample and plex
  sample_plex_medians = norm_bridge %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Sample,TMT) %>%
    dplyr::summarise(sample_plex_medians = median(intermediate_norm,na.rm=T))


  #combining medians with data frame
  norm_bridge = norm_bridge %>%
    dplyr::inner_join(sample_plex_medians,by=c("TMT","Sample"))

  #Doing the final Normalization
  output = norm_bridge %>%
    dplyr::mutate(final_norm = intermediate_norm / (sample_plex_medians / median_overall))

  #Adding option to export data in long or wide format
  if(data_format == "long"){
    return(output)
  }else if(data_format == "wide"){
    output2 = output %>%
      dplyr::select(Sample,TMT,ProteinID,final_norm) %>%
      tidyr::pivot_wider(names_from = c("Sample","TMT"), values_from = final_norm)
    return(output2)
  }else{
    print("data_format must be either long or wide")
  }
}
