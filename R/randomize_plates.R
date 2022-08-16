#' randomize_plates
#'
#' this function will randomize 96 well plate location given a dataframe containing metadata.
#' Note that this will fill as many plates as can be completely and then add samples to the last (n < 96) remaining plate randomly.
#'
#' @param metadata this is the metadata for your experiment.
#'
#' @return metadata dataframe with randomized plate assignment per sample.
#' @export
#'
#' @examples
#'
#' metadata = data.frame(sample_id = 1:108)
#' metadata_with_randomly_assigned_plates= randomize_plates(metadata)
#'
randomize_plates = function(metadata){
  n_samples = nrow(metadata)
  plates= c(1,ceiling(n_samples/96))
  nplate = length(plates)
  last_plate = plates[nplate]
  full_plates = plates[1:nplate-1]

  remainder = n_samples %% 96

  plate_vector = c(rep(full_plates,each = 96),rep(last_plate, each = remainder))

  random = sample(plate_vector,n_samples)

  out = metadata %>%
    dplyr::bind_cols(plate = random)

  return(out)
}
