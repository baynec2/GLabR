#' psm_phospho_mod
#' This code modifies the PSM file that is output from proteome discoverer to facilitate phospho PTM analysis.
#' @param data This is the=PSM proteome discoverer output file that has been into R as a tibble
#'
#' @return a tibble with a Protein ID column that is compatible with downstream processing through the combine_psm_fractions() and other functions.
#'
#' Note that the concatination of Master Accession, Annotated Sequence and ptmRS are concatinated into the ProteinID column.
#' Even though the name ProteinID isn't named exactly correctly, it was much easier to do it this way since that was the convention susequent code uses.
#'
#' @export
#'
#' @examples
psm_phospho_mod = function(data){

mod_data = data %>%
  dplyr::mutate(ProteinID = stringr::str_split(`Master Protein Accessions`,";",simplify = T)[,1],
                `Annotated Sequence` = toupper(`Annotated Sequence`)) %>%
  #Even though we are naming the column ProteinID, it isn't really a ProteinID - instead it is an ID. Done to be compatible with other functions easily.
  dplyr::mutate(ProteinID = paste0(ProteinID,"_",`Annotated Sequence`,"_",`ptmRS: Best Site Probabilities`)) %>%
  dplyr::select(ProteinID,dplyr::everything()) %>%
  dplyr::mutate(GLabR_step = "psm_phospho_mod")

return(mod_data)
}


