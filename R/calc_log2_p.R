#' calc_log2_p
#'
#' This function can be used to give the log2 fold change and the pvalues from a t test of proteins between groups.
#'
#' @param data This is a data frame containing Protein ID
#' @param column_split_by This is the factor you would like to split data and compare by. Must only have 2 levels.
#' @return a tibble
#' @export
#'
#' @examples
#'# Reading and processing data
#'data = readr::read_delim("tests/testdata/combine_psm_fractions/PCB002_PSMs_Proteodiscover_output.txt") %>%
#'combine_psm_fractions() %>%
#' normalize_1plex() %>%
#'la_box_cox_norm() %>%
#'dplyr::mutate(Sample_ID = paste0(Sample,".",TMT))
#'
#'#Reading in metadata
#'md = readr::read_csv("tests/testdata/metadata.csv") %>%
#'#Removed redundant columns to prevent .x and .y columns in data_md
#'dplyr::select(-Sample,-TMT)
#'
#'# Appending md to data
#'data_md = dplyr::inner_join(data,md,by = "Sample_ID")
#'
#'#Filtering as appropriately
#'f_data_md = data_md %>%
#'dplyr::filter(`Mayo_Endoscopic_Sub_Score` %in% c("Healthy_control","3: Severe disease (spontaneous bleeding, ulceration)")) %>%
#'dplyr::mutate_if(is.numeric,~tidyr::replace_na(.,0))
#'
#'#generating stats
#'stats = calc_log2_p(f_data_md,"Mayo_Endoscopic_Sub_Score")
calc_log2_p = function(data,column_split_by) {
  # For Phospho Data
  if ("Phospho_Prot_ratio" %in% names(data)) {
    stat = data %>% dplyr::group_by(ProteinID, Annotated_Sequence,
                                    ptmRS) %>% rstatix::t_test(as.formula(paste(
                                      "Phospho_Prot_ratio",
                                      "~", column_split_by
                                    ))) %>% rstatix::adjust_pvalue(p.col = "p",
                                                                   output.col = "p.adj_fdr",
                                                                   method = "fdr")
    log2 = data %>% dplyr::group_by_at(dplyr::vars("ProteinID",
                                                   column_split_by)) %>% dplyr::summarise(mean = mean(Phospho_Prot_ratio,
                                                                                                      na.rm = T)) %>% tidyr::pivot_wider(names_from = 2,
                                                                                                                                         values_from = 3)
    first_col_name = names(log2)[2]
    second_col_name = names(log2)[3]
    total_header = paste0("Log2FC(", first_col_name, "/",
                          second_col_name, ")")
    log2 = log2 %>% dplyr::mutate(`:=`(!!total_header, log2(
      !!rlang::sym(first_col_name) / !!rlang::sym(second_col_name)
    )))
    d1 = dplyr::inner_join(stat, log2, by = "ProteinID") %>%
      dplyr::mutate(pi_score = GLabR::calc_pi_score(p.adj_fdr, !!rlang::sym(total_header)))
    return(d1)
  }
  else {
    stat = data %>% dplyr::group_by(ProteinID) %>% rstatix::t_test(as.formula(paste(
      "box_cox_scaled_values",
      "~", column_split_by
    ))) %>% rstatix::adjust_pvalue(p.col = "p",
                                   output.col = "p.adj_fdr",
                                   method = "fdr")
    log2 = data %>% dplyr::group_by_at(dplyr::vars("ProteinID",
                                                   column_split_by)) %>% dplyr::summarise(mean = mean(box_cox_scaled_values,
                                                                                                      na.rm = T)) %>% tidyr::pivot_wider(names_from = 2,
                                                                                                                                         values_from = 3)
    first_col_name = names(log2)[2]
    second_col_name = names(log2)[3]
    total_header = paste0("Log2FC(", first_col_name, "/",
                          second_col_name, ")")
    log2 = log2 %>% dplyr::mutate(`:=`(!!total_header, log2(
      !!rlang::sym(first_col_name) / !!rlang::sym(second_col_name)
    )))
    d1 = dplyr::inner_join(stat, log2, by = "ProteinID") %>%
      dplyr::mutate(pi_score = GLabR::calc_pi_score(p.adj_fdr, !!rlang::sym(total_header)))
    return(d1)
  }
}
