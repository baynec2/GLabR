#' combine_psm_fractions: Combine PSM matches from different fractions in a proteomics run.
#'
#' This is an analogous function to the part of the TMT16/10_norm script that outputs nonnormalizeddataall.txt.
#' Briefly, it filters out the PSMs if their average is less than 10 across all of the plexes, if the isolation interference is > 30, or if the PSM Ambiguity == Rejected.
#'
#' @param data : This is a tibble containing the contents of the file that was exported from proteomediscoverer as a .txt file.
#' Note that the file names must be named following the convention of "Plex_Fraction" for this to work properly.
#' ie if the plex was PCB001 and the fraction number was F1, it would need to named PCB001_F1.raw.
#' This is because the code uses the file name to determine what plexes should be grouped together.
#'
#' @param analysis_type : This is a string representing the type of analysis that is being done. Currently, proteomics and phospho are supported.
#'
#'
#' @return a tibble containing the summed abundance for each plex.
#' @export
#'
#' @examples
#' data = read_delim("tests/testdata/normalize_to_bridge/PSM_output.txt") %>%
#'combine_psm_fractions()
combine_psm_fractions = function(data,analysis_type = "proteomics"){
    #Filtering to only keep high quality PSMs based on the criteria from Jacob's script
  fdata = data %>%
    #Note that original script replaced nas with 1s, but a better approach is to leave NAs as NA.
    #Taking the average across all of the TMT channels
    dplyr::mutate(n_NA = rowSums(is.na(dplyr::select(., dplyr::contains("Abundance"))))) %>%
    #is this correct? Not sure why this makes sense to do. Is the sum across a row this really signal to noise?
    #Doesn't seem like the best to me but including because that is what has been traditionally done.
    dplyr::mutate(AvgSN = rowMeans(dplyr::select(., dplyr::contains("Abundance")), na.rm = TRUE)) %>%
    #Filtering out the rows that contain the following
    dplyr::filter(`PSM Ambiguity` != "Rejected",
           `Isolation Interference [%]` < 30,
           AvgSN >= 10) %>%
    #Spiting the file name into plexs
    dplyr::mutate(Sample = stringr::str_split(`Spectrum File`,"_",simplify = T)[,1],
           Fraction = stringr::str_extract(`Spectrum File`, "F[:digit:]"))

  #Adapting code to deal with both proteomics and phospho via psm_phospho_mod
  if(analysis_type == "proteomics"){
  #Taking only the first Protein Accession listed in the Master.Protein.Accession Column
    fdata = fdata %>%
      dplyr::mutate(ProteinID = stringr::str_split(`Master Protein Accessions`,";",simplify = T)[,1])
  }else if(analysis_type == "phospho"){
  # We don't need to change the columns if we are doing a phospho analysis
    fdata = fdata

  }else{
    print("analysis_type must be either proteomics or phospho")
  }
  #Combining all of the PSMs from each sample
  output = fdata %>%
    dplyr::select(Sample,Fraction,ProteinID,dplyr::contains("Abundance")) %>%
    tidyr::pivot_longer(4:length(.)) %>%
    dplyr::group_by(Sample,ProteinID,name) %>%
    # INote this was modified to prevent coercion to 0s in the case that all are NA.
    dplyr::summarise(value = sumna(value)) %>%
    dplyr::mutate(TMT = stringr::str_extract(name,"\\d\\d\\d[n|c|N|C]?")) %>%
    dplyr::select(Sample,TMT,ProteinID,value) %>%
    dplyr::ungroup()

    return(output)
  }
