#' calc_pi_score
#'
#' This function returns the pi score
#'
#' @param pval : p value
#' @param Log2FC : log 2 fold change
#'
#' @return The pi score
#' @export
#'
#' @examples
#' pval = 0.001
#' Log2FC = 2
#' pi_score = calc_pi_score(pval,Log2FC)
calc_pi_score = function(pval,Log2FC){
  output = -log10(pval) * abs(Log2FC)
  return(output)
}
