#' randomize_tmt
#'
#' This function will provide randomized tmt labels for given plex groups.
#'
#' @param metadata this is a data frame containing at the minimum a column named plex_number.
#' This column is
#' Note that this should include the bridge channel for each plex.
#'
#' @param tmt_plex this specifies whether a 10 or 16 plex is being used.
#'
#' @return a data frame with randomly assigned plex numbers
#' @export
#'
#' @examples
#'
#' metadata = data.frame(sample_id = 1:30,plex_num = rep(c(1,2), each = 15))
#' metadata_with_randomly_assigned_tmts = randomize_tmt(metadata)
#'
randomize_tmt = function(metadata, tmt_plex = 16){
  set.seed(2)
  if(tmt_plex == 10){
    labels = c("126","127N","127C","128N",
               "128C","129N","129C","130N",
               "130C","131")
    gmetadata = metadata %>%
      dplyr::group_by(plex_num) %>%
      dplyr::mutate(tmt_label = sample(labels,size = n()))

    return(gmetadata)
    }else if(tmt_plex == 16){
    labels = c("126","127N","127C","128N",
               "128C","129N","129C","130N",
               "130C","131N","131C","132N",
               "132C","133N","133C","134N")

    gmetadata = metadata %>%
      dplyr::group_by(plex_num) %>%
      dplyr::mutate(tmt_label = sample(labels,size = n()))

    return(gmetadata)
    }else{
    print("tmt_plex must be either 10 or 16")
  }
}


