#' extract_sig_proteins
#'
#' This function extracts the protein IDs for all proteins that have significant differences between the groups being compared.
#'
#' @param data this is the data
#' @param column_split_by this is the column to split by
#' @param p_threshold this is the p threshold to use
#' @param fc_threshold this is the fc threshold to use
#' @param pi_score_theshold this is the pi score threshold to use if applicable
#'
#' @return a tibble
#' @export
#'
#' @examples
#' #' # Reading and processing data
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
#'#extracting the sig proteins
#'extract_sig_proteins(f_data_md,"Mayo_Endoscopic_Sub_Score",p_threshold = 0.05,fc_threshold = 1)
extract_sig_proteins = function(data,column_split_by,p_threshold = 0.05,fc_threshold = 1,pi_score_theshold = NULL){

  #Generating the stats
  stat = data %>%
    dplyr::group_by(ProteinID) %>%
    rstatix::t_test(as.formula(paste("box_cox_scaled_values",'~',column_split_by))) %>%
    rstatix::adjust_pvalue(p.col = "p",output.col = "p.adj_fdr", method = "fdr")

  #Calculating log2FC
  log2 = data %>%
    dplyr::group_by_at(dplyr::vars("ProteinID",column_split_by)) %>%
    dplyr::summarise(mean = mean(box_cox_scaled_values,na.rm = T)) %>%
    tidyr::pivot_wider(names_from = 2,values_from = 3)

  #Pulling out name of first and second column
  first_col_name = names(log2)[2]
  second_col_name = names(log2)[3]

  #Defining new mutate column name
  total_header=paste0("Log2FC(",first_col_name,"/", second_col_name,")")

  log2 = log2 %>%
    dplyr::mutate(!!total_header := log2(!!rlang::sym(first_col_name)/!!rlang::sym(second_col_name)))

  d1 = dplyr::inner_join(stat,log2,by = "ProteinID") %>%
    dplyr::mutate(pi_score = GLabR::calc_pi_score(p.adj_fdr,!!rlang::sym(total_header)))

  #Labeling what
  labels = d1 %>%
    dplyr::filter(p.adj_fdr <= p_threshold & (!!rlang::sym(total_header) >= fc_threshold | !!rlang::sym(total_header) <= -fc_threshold))


  return(labels$ProteinID)
}
