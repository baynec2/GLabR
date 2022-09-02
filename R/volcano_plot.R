#' volcano_plot
#' this plots a volcano plot with proteins that meet the specified criteria highlighted.
#' @param data = this is a tibble that has been loaded into R.
#' @param column_split_by = this is the name of the column you would like to split data by in order to compare groups
#' Note that it must contain only 2 groups, and the column name cannot have any spaces.
#' @param p_threshold = this is the adjusted p value threshold that is used to generate the dashed line on the volcano plot
#' @param fc_threshold = this is the fc_threshold that is used to generate the dashed line on the volcano plot.
#' @param pi_score_theshold = if used, the indicated pi score threshold is used instead of the p_threshold and fc_threshold.
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' This is where I will put an example.
volcano_plot = function(data,column_split_by,p_threshold = 0.05,fc_threshold = 1){

  # Handling phospho data appropriately (detecting this data based on column name)
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

    #Labeling what
    labels = d1 %>%
      dplyr::filter(p.adj_fdr <= p_threshold & (!!rlang::sym(total_header) >= fc_threshold | !!rlang::sym(total_header) <= -fc_threshold))

    # Doing the plotting and labeling the significant proteins.
    p1 = d1 %>%
      ggplot2::ggplot(ggplot2::aes(!!rlang::sym(total_header),-log10(p.adj_fdr)))+
      ggplot2::geom_point()+
      ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")+
      ggplot2::geom_vline(xintercept = -fc_threshold, linetype = "dashed", color = "red")+
      ggplot2::geom_vline(xintercept = fc_threshold, linetype = "dashed", color = "red")+
      ggrepel::geom_text_repel(data = labels,ggplot2::aes(!!rlang::sym(total_header),-log10(p.adj_fdr),label = ProteinID))

    return(p1)
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

    #Labeling what
    labels = d1 %>%
      dplyr::filter(p.adj_fdr <= p_threshold & (!!rlang::sym(total_header) >= fc_threshold | !!rlang::sym(total_header) <= -fc_threshold))

    # Doing the plotting and labeling the significant proteins.
    p1 = d1 %>%
      ggplot2::ggplot(ggplot2::aes(!!rlang::sym(total_header),-log10(p.adj_fdr)))+
      ggplot2::geom_point()+
      ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")+
      ggplot2::geom_vline(xintercept = -fc_threshold, linetype = "dashed", color = "red")+
      ggplot2::geom_vline(xintercept = fc_threshold, linetype = "dashed", color = "red")+
      ggrepel::geom_text_repel(data = labels,ggplot2::aes(!!rlang::sym(total_header),-log10(p.adj_fdr),label = ProteinID))

    return(p1)
  }
}
