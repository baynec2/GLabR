#' adjust_nanopore_barcodes
#'
#' this function is designed to take the proportion of nanopore barcodes observed in a previous run, and calculate how the user should readjust the barcodes to achieve the same number of basepairs across samples.
#' this is probably most usefull in the case where a flongle has been used to observe the barcode distribution of pooled sample prior to a minion run.
#' After seeing this distribution, the barcodes could be re-pooled to achieve an even distribution.
#'
#' @param barcodes this is the names of the barcodes.
#' @param proportion_of_bases this is the proportion of base pairs that belong to each barcode. Note this should be in the same order as the names of the barcodes.
#' @param volume_of_each_barcode_added this is the volume that the protocol indicates each barcode to be added at. Note this should be in the same order as the names of the barcodes.
#' @param same_total_volume do you want the same total volume? Or do you want a different total volume? Different total volume will have the same amount of total DNA while same total volume will be even but scaled down.
#' @return a tibble with intermediate calculations and final volume to be added.
#' @export
#'
#' @examples
#'values = runif(90)
#'x = values/sum(values)
#'barcodes = test$barcode
#'proportion_of_bases = test$proportion_of_bases
#'volume_of_each_barcode_added = 1.25
#'adjust_nanopore_barcodes(barcodes,proportion_of_bases,volume_of_each_barcode_added = 1.25,same_total_volume == TRUE)
adjust_nanopore_barcodes = function(barcodes,proportion_of_bases,volume_of_each_barcode_added = 1.25,same_total_volume = TRUE){

# What proportion would each barcode be if everything worked such that each barcode was evenly distributed.
expected_proportion = 1/ length(barcodes)

# what would the adjustment factor have to be to match the expected proportions given each proportion that was actually observed?
adjustments = expected_proportion / proportion_of_bases

# how would this translate into the amount of volume that needs to be added?


if(same_total_volume == TRUE){
  # This gives the same volume at the end as what the original protocol calls for.
  # You could also argue that you don't need to do that, you could just use a greater volume since there is an ampure cleanup.
  # factor to multiply by adjustment = total_volume / sum of proportions
  # adjusted volume = factor to multiply by adjustment * adjustment

  total_volume = volume_of_each_barcode_added * length(barcodes)
  sum_of_proportions = sum(adjustments)

  x = (total_volume / sum_of_proportions)

adjusted_volume = x * adjustments
}else{
  adjusted_volume = adjustments * volume_of_each_barcode_added
}
df = tibble::tibble(barcodes = barcodes,
                    proportion_of_bases = proportion_of_bases,
                    expected_proportion = expected_proportion,
                    adjustment = adjustments,
                    adjusted_volume_uL= adjusted_volume)
return(df)
}

