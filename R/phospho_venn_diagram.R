#' phospho_venn_diagram
#'
#' This function takes the results from a non-phospho enriched proteomics experiment (Sequest with dynamic phospho mode and ptmRS node enabled), looks at uniwue Protein-Sequence-PTM combinations
#' and compares them to the results from a phospho enriched proteomics experiment (Sequest with dynamic phospho mode and ptmRS node enabled).
#'
#' @param proteomics_data output from a proteome discoverer exported text file acquured from a non phospho enriched experiment. Intended to have been read into R using read_delim().
#' Note that the workflow used in proteome discoverer must include the ptmRS node and dynamic phospho mods in Sequest.
#' @param phospho_enriched_data from a paired phospho-enriched experient (same samples as proteomics data). Intended to have beenread into R using read_delim().
#' Note that the workflow used must include the ptmRS node and dynamic phospho mods in Sequest.
#' @param bridge_channel this is the bridge channel used if applicable. Default is NULL meaning that data will be processed as a single plex unless otherwise indicated.
#'
#' @return a ggplot2 object
#' @export
#' @examples
#'
#' phospho_enriched = readr::read_delim("tests/testdata/psm_phospho_mod/PCB001_PSM.txt")
#' proteomics = readr::read_delim("tests/testdata/phospho_venn_diagram/PCB002_proteomics_ptmRS_data.txt")
#'phospho_venn_diagram(proteomics,phospho_enriched)
phospho_venn_diagram = function(proteomics_data, phospho_enriched_data,bridge_channel = NULL){

  if(is.null(bridge_channel)){

    proteomics_data_2 = GLabR::psm_phospho_mod(proteomics_data) %>%
      GLabR::combine_psm_fractions(analysis_type = "phospho") %>%
      GLabR::normalize_1plex() %>%
      GLabR::la_box_cox_norm() %>%
      tidyr::separate(ProteinID,into = c("ProteinID","Annotated_Sequence","ptmRS") ,sep = "_") %>%
      dplyr::mutate(PTMID = paste0(ProteinID,"_",Annotated_Sequence,"_",ptmRS)) %>%
      dplyr::mutate(SampleID = paste0(Sample,".",TMT)) %>%
      dplyr::filter(ptmRS != "NA") %>%
      dplyr::pull(PTMID) %>%
      unique()


    phospho_data_2 = GLabR::psm_phospho_mod(phospho_enriched_data) %>%
      GLabR::combine_psm_fractions(analysis_type = "phospho") %>%
      GLabR::normalize_1plex() %>%
      GLabR::la_box_cox_norm() %>%
      tidyr::separate(ProteinID,into = c("ProteinID","Annotated_Sequence","ptmRS") ,sep = "_") %>%
      dplyr::mutate(PTMID = paste0(ProteinID,"_",Annotated_Sequence,"_",ptmRS)) %>%
      dplyr::mutate(SampleID = paste0(Sample,".",TMT)) %>%
      dplyr::filter(ptmRS != "NA") %>%
      dplyr::pull(PTMID) %>%
      unique()

    list = list(proteomics = proteomics_data_2,phosho_enriched = phospho_data_2)


    p1 = ggVennDiagram::ggVennDiagram(list)

    return(p1)
  }else{
    proteomics_data_2 = GLabR::psm_phospho_mod(phospho_enriched_data) %>%
      GLabR::combine_psm_fractions(analysis_type = "phospho") %>%
      GLabR::normalize_bridge(bridge_channel) %>%
      GLabR::la_box_cox_norm() %>%
      tidyr::separate(ProteinID,into = c("ProteinID","Annotated_Sequence","ptmRS") ,sep = "_") %>%
      dplyr::mutate(PTMID = paste0(ProteinID,"_",Annotated_Sequence,"_",ptmRS)) %>%
      dplyr::mutate(SampleID = paste0(Sample,".",TMT)) %>%
      dplyr::filter(ptmRS != "NA") %>%
      dplyr::pull(PTMID) %>%
      unique()

    phospho_data_2 = GLabR::psm_phospho_mod(phospho_enriched_data) %>%
      GLabR::combine_psm_fractions(analysis_type = "phospho") %>%
      GLabR::normalize_bridge(bridge_channel) %>%
      GLabR::la_box_cox_norm() %>%
      tidyr::separate(ProteinID,into = c("ProteinID","Annotated_Sequence","ptmRS") ,sep = "_") %>%
      dplyr::mutate(PTMID = paste0(ProteinID,"_",Annotated_Sequence,"_",ptmRS)) %>%
      dplyr::mutate(SampleID = paste0(Sample,".",TMT)) %>%
      dplyr::filter(ptmRS != "NA") %>%
      dplyr::pull(PTMID) %>%
      unique()

    list = list(proteomics = proteomics_data_2,phosho_enriched = phospho_data_2)


    p1 = ggVennDiagram::ggVennDiagram(list)

    return(p1)

  }
}
