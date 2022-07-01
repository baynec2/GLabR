#' combine_psm_fractions_replica : This function exactly replicates the old script that was used to combine and normalize the PSM data.
#' Should generally not be used except for comparison purposes as this method was found to be flawed.
#'
#' @param filepath
#'
#' @return
#' @export
#'
#' @examples
combine_psm_fractions_replica = function(data){

  #Filtering to only keep high quality PSMs based on the criteria from Jacob's script
  fdata = data %>%
    #Replacing the same columns as the previous code with 1s. I am pretty sure that this is a bug and shouldn't be done.
    mutate(`Abundance: 132N` = tidyr::replace_na(`Abundance: 132N`,replace = 1),
           `Abundance: 132C` = tidyr::replace_na(`Abundance: 132C`,replace = 1),
           `Abundance: 133N` = tidyr::replace_na(`Abundance: 133N`,replace = 1),
           `Abundance: 133C` = tidyr::replace_na(`Abundance: 133C`,replace = 1),
           `Abundance: 134N` = tidyr::replace_na(`Abundance: 134N`,replace = 1)) %>%
    #Taking the average across all of the TMT channels
    #is this correct? Not sure why this makes sense to do. Is this really signal to noise?
    mutate(n_NA = rowSums(is.na(select(., contains("Abundance"))))) %>%
    mutate(AvgSN = rowMeans(select(., contains("Abundance")), na.rm = TRUE)) %>%
    filter(`PSM Ambiguity` != "Rejected",
           `Isolation Interference [%]` < 30,
           AvgSN >= 10) %>%
    #Spliting the file name into plexs
    mutate(Sample = stringr::str_split(`Spectrum File`,"_",simplify = T)[,1],
           Fraction = stringr::str_extract(`Spectrum File`, "F[:digit:]")) %>%
    #Taking only the first Protein Accession listed in the Master.Protein.Accession Column
    mutate(ProteinID = stringr::str_split(`Master Protein Accessions`,";",simplify = T)[,1]) %>%
    #This really makes the output more similar. Looks like Previous script just obmits all NAs.
    filter(n_NA == 0)

  #Combining all of the PSMs from each sample
  output = fdata %>%
    select(Sample,Fraction,ProteinID,contains("Abundance")) %>%
    tidyr::pivot_longer(4:length(.)) %>%
    dplyr::group_by(Sample,ProteinID,name) %>%
    dplyr::summarise(value = sum(value,na.rm=T)) %>%
    dplyr::mutate(TMT = stringr::str_extract(name,"\\d\\d\\d[n|c|N|C]?")) %>%
    dplyr::select(Sample,TMT,ProteinID,value) %>%
    dplyr::ungroup()


  return(output)
}
