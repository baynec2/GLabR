#' t_test_log2
#'
#' This function reports the log2 change, t test p values/ p.adj values when comparing two groups.
#'
#' @param data : This is the inbut data. Intended to be the output of la_box_cox_norm()
#' @param column_split_by :  this is the name of the column you would like to split data by in order to compare groups. Should only contain 2 levels.
#' Note that it must contain only 2 groups, and the column name cannot have any spaces
#'
#' @return a tibble
#' @export
#'
#' @examples
#' data = readr::read_csv("tests/testdata/t_test_log2/t_test_data.csv")
#' t_test = t_test_log2(PCB002_norm,"Donor_type")
t_test_log2 = function(data,column_split_by){
  if("Phospho_Prot_ratio" %in% names(data)){

    stat = data %>%
      dplyr::group_by(ProteinID,Annotated_Sequence,ptmRS) %>%
      rstatix::t_test(as.formula(paste("Phospho_Prot_ratio",'~',column_split_by))) %>%
      rstatix::adjust_pvalue(p.col = "p",output.col = "p.adj_fdr", method = "fdr")

    #Calculating log2FC
    log2 = data %>%
      dplyr::group_by_at(dplyr::vars("ProteinID",column_split_by)) %>%
      dplyr::summarise(mean = mean(Phospho_Prot_ratio,na.rm = T)) %>%
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

    return(d1)
  }else{
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

    return(d1)
  }
}
