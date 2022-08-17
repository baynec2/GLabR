#' assign_plex_num
#'
#' this function assigns plex number for your metadata. This function also adds a bridge channel.
#' This is necessary for experiments that have more samples than can fit in a single plex.
#` You can specify that members of a given group (plex_grouping_factor) all have the same plex.
#` This is particularly useful for paired samples. `
#' @param metadata this is a data frame containing the metadata.
#' @param plex_grouping_factor all samples with the same value will be assigned the same plex.
#' This is useful for ensuring that pared samples are in the same plex.
#' @param tmt_plex the n of TMT plex that you are using, either 10 or 16.
#'
#' @return a data frame with asigned plex numbers
#' @export
#'
#' @examples
#'
#' metadata = data.frame(sample_id = as.character(1:108),
#'                       plate = c(rep(1,each = 96),rep(2,each = 12)),
#'                       mouse_id = rep(paste0("Mouse",1:27),each = 4)
#'                       )
#'
#'metadata_with_assigned_plex_n = assign_plex_number(metadata,"mouse_id")

assign_plex_num = function(metadata,plex_grouping_factor,tmt_plex = 16){

  #making sure metadata is arranged properly (necessary for adding plexes)
    o = stringr::str_order(metadata %>% dplyr::pull(plex_grouping_factor),numeric = T)
    metadata = metadata[o,]
  #determining how many samples there are
  n_samples = nrow(metadata)
  #determining how many samples are in each group
  n_per_group = n_samples / nrow(unique(metadata[plex_grouping_factor]))
  #How many groups can fit evenly in a plex
  n_sample_per_plex = floor((tmt_plex-1)/n_per_group) * n_per_group
  #determining how many bridge channels we will need
  n_bridge = ceiling(n_samples/(n_sample_per_plex))
  #calculating the total number of samples
  totaln = n_samples + n_bridge
  #determining how many extra samples there are (within a not complete plex)
  extra = totaln %% (n_sample_per_plex+1)
  #determining how many plexes we will need
  n_complete_plexes = ceiling(totaln/(n_sample_per_plex+1))

  #assigning a plex number to all of the samples
  plex_num = c(rep(1:n_complete_plexes,each = n_sample_per_plex),
               rep(n_bridge,each = extra))

  #assembling metadata with plex number for each sample
  md_plex = metadata %>%
    dplyr::bind_cols(plex_num = plex_num)

  #adding info about the bridge channel
  bridge_channels_df = data.frame(sample_id = paste0("bridge_channel_",1:n_bridge),
                                  plex_num =1:n_bridge)

  #Assembling the final data frame
  all_md = md_plex %>%
    dplyr::bind_rows(bridge_channels_df)


  return(all_md)
}


