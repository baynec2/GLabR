#' calc_phospho_ratio
#'
#' This function calculates the ratio of phosphosite to total proteins using box cox normalized values.
#'
#' @param proteomics_data this is a tibble containing proteomics data that has been read into R
#' @param phospho_data this is a tibble containing phosphoproteomics data that has been read into R (phospho enriched)
#' @param metadata this is the metadata that. Note the Plex/TMT combo needs to be concat
#' @param col_identifying_match this is the column in the data that identifies samples as being the same.
#' This is done using the metadata, and is important because this is how the ratio is actually calculated.
#' @param num_plexs this is the number of plexs that each dataset contains. If it is >1 things are normalized using the bridge channel
#' @param bridge_channel_plex this is the TMT label that contains the bridge channel to normalize to if num_plexes > 1.
#'
#' @return
#' @export
#'
#' @examples
calc_phospho_ratio = function(proteomics_data,phospho_data,metadata,col_identifying_match,num_plexs,bridge_channel_plex = "126"){

# Different code is used for 1 plex vs bridge channel plexes
if(num_plexs == 1){
proteomics_data_2 = GLabR::combine_psm_fractions(proteomics_data) %>%
  GLabR::normalize_1plex() %>%
  GLabR::la_box_cox_norm() %>%
  dplyr::rename(proteomics_box_cox_scaled_values = box_cox_scaled_values) %>%
  dplyr::select(Sample,TMT,ProteinID,proteomics_box_cox_scaled_values) %>%
  dplyr::mutate(SampleID = paste0(Sample,".",TMT))


phospho_data_2 = GLabR::psm_phospho_mod(phospho_data) %>%
  GLabR::combine_psm_fractions(analysis_type = "phospho") %>%
  GLabR::normalize_1plex() %>%
  GLabR::la_box_cox_norm() %>%
  tidyr::separate(ProteinID,into = c("ProteinID","Annotated_Sequence","ptmRS") ,sep = "_") %>%
  dplyr::rename(phospho_box_cox_scaled_values = box_cox_scaled_values) %>%
  dplyr::mutate(SampleID = paste0(Sample,".",TMT))

# Combining the different data streams. Calculating ratio.
# The function needs to know what sample IDs are the same between the two studies. There is no way to know this without the sampleIDs. These are contained in the metadata.

proteomics_data_2 = dplyr::inner_join(proteomics_data_2,metadata,by = "SampleID") %>%
  dplyr::select(!!col_identifying_match,ProteinID,proteomics_box_cox_scaled_values)
phospho_data_2 = dplyr::inner_join(phospho_data_2,metadata,by = "SampleID") %>%
  dplyr::select(!!col_identifying_match,ProteinID,Annotated_Sequence,ptmRS,phospho_box_cox_scaled_values)

combined = dplyr::inner_join(phospho_data_2,proteomics_data_2,by = c(col_identifying_match,"ProteinID")) %>%
  dplyr::mutate(Phospho_Prot_ratio = phospho_box_cox_scaled_values / proteomics_box_cox_scaled_values)

return(combined)


}else{
  proteomics_data_2 = GLabR::combine_psm_fractions(proteomics_data) %>%
    GLabR::normalize_to_bridge(bridge_channel_plex = bridge_channel_plex) %>%
    GLabR::la_box_cox_norm() %>%
    dplyr::rename(proteomics_box_cox_scaled_values = box_cox_scaled_values) %>%
    dplyr::select(Sample,TMT,ProteinID,proteomics_box_cox_scaled_values) %>%
    dplyr::mutate(SampleID = paste0(Sample,".",TMT))

  phospho_data_2 = GLabR::psm_phospho_mod(phospho_data) %>%
    GLabR::combine_psm_fractions(analysis_type = "phospho") %>%
    GLabR::normalize_to_bridge(bridge_channel_plex = bridge_channel_plex) %>%
    GLabR::la_box_cox_norm() %>%
    tidyr::separate(ProteinID,into = c("ProteinID","Annotated_Sequence","ptmRS") ,sep = "_") %>%
    dplyr::rename(phospho_box_cox_scaled_values = box_cox_scaled_values) %>%
    dplyr::mutate(SampleID = paste0(Sample,".",TMT))

  proteomics_data_2 = dplyr::inner_join(proteomics_data_2,metadata,by = "SampleID") %>%
    dplyr::select(!!col_identifying_match,ProteinID,proteomics_box_cox_scaled_values)
  phospho_data_2 = dplyr::inner_join(phospho_data_2,metadata,by = "SampleID") %>%
    dplyr::select(!!col_identifying_match,ProteinID,Annotated_Sequence,ptmRS,phospho_box_cox_scaled_values)

  combined = dplyr::inner_join(phospho_data_2,proteomics_data_2,by = c(col_identifying_match,"ProteinID")) %>%
    dplyr::mutate(Phospho_Prot_ratio = phospho_box_cox_scaled_values / proteomics_box_cox_scaled_values)

return(combined)
}

}





