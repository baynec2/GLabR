#' randomize_wells
#'
#' This function will randomize well location given a data frame of metadata.
#'
#' @param metadata this is a dataframe containing the metadata.
#' Note that it must have a column titled plate which contains the plate number of each sample.
#'
#' @return
#' @export
#'
#' @examples
#' metadata = data.frame(sample_id = 1:108,plate = c(rep(1,each = 96),rep(2,each = 12)))
#' metadata_with_randomly_assigned_wells= randomize_wells(metadata)
randomize_wells = function(metadata){

  set.seed(2)
  rows = rep(LETTERS[1:8],each = 12)
  columns = rep(1:12,8)
  wells = paste0(rows,columns)

  gdf = metadata %>%
    dplyr::group_by(plate) %>%
    dplyr::mutate(well = sample(wells,size = n()))

  return(gdf)
}

metadata = data.frame(sample_id = 1:108,plate = c(rep(1,each = 96),rep(2,each = 12)))
metadata_with_randomly_assigned_wells= randomize_wells(metadata)
