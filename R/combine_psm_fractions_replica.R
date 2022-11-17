#' combine_psm_fractions_replica : This function exactly replicates the old script that was used to combine and normalize the PSM data.
#' Should generally not be used except for comparison purposes as this method was found to be flawed.
#'
#' @param data This is a tibble containing the contents of the file that was exported from proteomediscoverer as a .txt file.
#' Note that the file names must be named following the convention of "Plex_Fraction" for this to work properly.
#' ie if the plex was PCB001 and the fraction number was F1, it would need to named PCB001_F1.raw.
#' This is because the code uses the file name to determine what plexes should be grouped together.
#'
#' @return a tibble
#' @export
#'
#' @examples
#'data = read_delim("tests/testdata/normalize_to_bridge/PSM_output.txt") %>%
#'combine_psm_fractions_replica()
combine_psm_fractions_replica = function(data){

  #Filtering to only keep high quality PSMs based on the criteria from Jacob's script
  fdata = data %>%
    #Replacing the same columns as the previous code with 1s. I am pretty sure that this is a bug and shouldn't be done.
    dplyr::mutate(`Abundance: 132N` = tidyr::replace_na(`Abundance: 132N`,replace = 1),
           `Abundance: 132C` = tidyr::replace_na(`Abundance: 132C`,replace = 1),
           `Abundance: 133N` = tidyr::replace_na(`Abundance: 133N`,replace = 1),
           `Abundance: 133C` = tidyr::replace_na(`Abundance: 133C`,replace = 1),
           `Abundance: 134N` = tidyr::replace_na(`Abundance: 134N`,replace = 1)) %>%
    #Taking the average across all of the TMT channels
    #is this correct? Not sure why this makes sense to do. Is this really signal to noise?
    dplyr::mutate(n_NA = rowSums(is.na(dplyr::select(., dplyr::contains("Abundance"))))) %>%
    dplyr::mutate(AvgSN = rowMeans(dplyr::select(., dplyr::contains("Abundance")), na.rm = TRUE)) %>%
    dplyr::filter(`PSM Ambiguity` != "Rejected",
           `Isolation Interference [%]` < 30,
           AvgSN >= 10) %>%
    #Spliting the file name into plexs
    dplyr::mutate(Sample = stringr::str_split(`Spectrum File`,"_",simplify = T)[,1],
           Fraction = stringr::str_extract(`Spectrum File`, "F[:digit:]")) %>%
    #Taking only the first Protein Accession listed in the Master.Protein.Accession Column
    dplyr::mutate(ProteinID = stringr::str_split(`Master Protein Accessions`,";",simplify = T)[,1]) %>%
    #This really makes the output more similar. Looks like Previous script just obmits all NAs.
    dplyr::filter(n_NA == 0)


  #Combining all of the PSMs from each sample
  output = fdata %>%
    dplyr::select(Sample,Fraction,ProteinID,dplyr::contains("Abundance")) %>%
    tidyr::pivot_longer(4:length(.)) %>%
    dplyr::group_by(Sample,ProteinID,name) %>%
    dplyr::summarise(value = sum(value,na.rm=T)) %>%
    dplyr::mutate(TMT = stringr::str_extract(name,"\\d\\d\\d[n|c|N|C]?")) %>%
    dplyr::select(Sample,TMT,ProteinID,value) %>%
    dplyr::ungroup()


  return(output)
}
