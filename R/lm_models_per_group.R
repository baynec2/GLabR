#' lm_models_per_group
#'
#' This function groups by a user specified variable and then conducts linear models using user specified x and y values for each group.
#' The p values of each linear model are then adjusted using a fdr correction.
#' I have found this to be useful for proteomics, and metagenomics data in the past.
#' In the case of proteomics, it is useful to group by protein, and then see whether the intensity of each protein correlates with an outcome such as tumor volume.
#' In the case of metagenomics, it can be useful to group by some level of taxonomy and then see how relative abundance correlates with an outcome such as tumor volume.
#'
#' @param data this is the data you will be operating on.
#' @param group this is the group you would like to do linear correlations for each member of.
#' @param yvalue this is the y value.
#' @param xvalue this is the x value.
#'
#' @return a tibble with the r squared values, p value, fdr adjusted pvalue, and direction of the correlation.
#' @export
#'
#' @examples
#'
#' data = data.frame(group = rep(1:10,each = 50), yvalue = rnorm(500,10), xvalue = rnorm(500,10))
#' mod = lm_models_per_group(data,"group","yvalue","xvalue")
lm_models_per_group = function(data,group,yvalue,xvalue){

  # This will give the summary terms from each model
  fitted_models = data %>%
    dplyr::group_by(!!!syms(group)) %>%
    do(broom::glance(lm(!!!syms(yvalue) ~ !!!syms(xvalue), data = .))) %>%
    dplyr::ungroup() %>%
    rstatix::adjust_pvalue(method = "fdr")

  # Above doesn't give information about the direction of correlation, getting that below
  dir = data %>%
    dplyr::group_by(!!!syms(group)) %>%
    do(broom::tidy(lm(!!!syms(yvalue) ~ !!!syms(xvalue), data = .))) %>%
    dplyr::filter(!grepl(".*Intercept.*",term)) %>%
    dplyr::mutate(direction = case_when(estimate <= 0 ~ "negative",
                                        TRUE ~ "postive")) %>%
    dplyr::select(!!!syms(group),direction)

  # Combining the model stats and directionality.
  fitted_models = dplyr::inner_join(fitted_models,dir,by = group) %>%
    dplyr::arrange(p.value.adj) %>%
    dplyr::select(!!!syms(group),r.squared,p.value,p.value.adj,direction)

}
