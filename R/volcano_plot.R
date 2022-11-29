#' volcano_plot
#' this plots a volcano plot with proteins that meet the specified criteria highlighted.
#' @param data = this is intended to be the output of calc_log2_p.
#' @param column_split_by = this is the name of the column you would like to split data by in order to compare groups
#' Note that it must contain only 2 groups, and the column name cannot have any spaces.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' # Reading and processing data
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
#'
#'#plotting volcano plot
#'volcano plot(stats,p_theshold = 0.05, log2fc_threshold = 1)
volcano_plot = function(calc_log2_p_output,p_threshold = 0.05, log2fc_threshold = 1){

  d1 = calc_log2_p_output
  total_header = names(d1)[grepl("Log2FC.*",names(d1))]
  labels = d1 %>% dplyr::filter(p.adj_fdr <= p_threshold &
                                  (!!rlang::sym(total_header) >= fc_threshold | !!rlang::sym(total_header) <=
                                     -fc_threshold))
  p1 = d1 %>%
    ggplot2::ggplot(ggplot2::aes(!!rlang::sym(total_header),
                                 -log10(p.adj_fdr))) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = -log10(0.05),
                        linetype = "dashed", color = "red") +
    ggplot2::geom_vline(xintercept = -fc_threshold,
                        linetype = "dashed", color = "red") +
    ggplot2::geom_vline(xintercept = fc_threshold,
                        linetype = "dashed", color = "red") +
    ggrepel::geom_text_repel(data = labels,
                                                                                                                                                                           ggplot2::aes(!!rlang::sym(total_header), -log10(p.adj_fdr),
                                                                                                                                                                                                                                                                                                             label = ProteinID))
  return(p1)
}
