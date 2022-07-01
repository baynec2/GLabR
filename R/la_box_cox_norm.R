#' la_box_cox_norm
#' here is leigh-ana's method for conducting box_cox_normalization on proteomics data.
#'
#' @param data : This is a tibble containing at least the columns named Sample , TMT, and final_norm. This is intended to be the output of normalize_to_bridge() or noramlize_1plex()
#'
#' @return a normalized csv file in the same directory.
#' @export
#'
#' @examples
la_box_cox_norm = function(data){

 #Had to modify data to make compatible with leighana's script. Removed all nas, infinte values, and 0s.
 # Transformed to wide format to do lm in column format (didn't feel like figuring out how to do this within the tidyverse)
  mod_data = data %>%
    dplyr::select(Sample,TMT,ProteinID,final_norm) %>%
    dplyr::filter(is.finite(final_norm),
                  final_norm != 0 ) %>%
    tidyr::pivot_wider(names_from = c("Sample","TMT"),values_from = "final_norm") %>%
    na.omit() %>%
    as.data.frame()

  # Leigh-ana's script
  transformed_data <- matrix(data=NA,nrow=length(mod_data$ProteinID),ncol=length(mod_data)-1)
  rownames(transformed_data) <- mod_data$ProteinID
  colnames(transformed_data) <- colnames(mod_data)[2:length(mod_data)]

  for (i in 2:(length(mod_data))) {
    temp_data = mod_data[,i]
    b <-MASS::boxcox(lm(temp_data ~ 1),plotit = FALSE,interp = TRUE)
    lambda <- b$x[which.max(b$y)]
    new_data <- (temp_data^lambda - 1)/lambda
    scaled_data <- reshape::rescaler(new_data,type="range")
    scaled_data <- scaled_data/mean(scaled_data,na.rm = TRUE)
    transformed_data[,i-1] <- scaled_data
  }


  #tranforming final data into long data format
  output = transformed_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "ProteinID") %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(2:length(.)) %>%
    tidyr::separate(name, into = c("Sample", "TMT"),sep = "_") %>%
    dplyr::select(Sample,TMT,ProteinID,box_cox_scaled_values = value) %>%
    dplyr::inner_join(data,by = c("Sample", "TMT", "ProteinID")) %>%
    dplyr::select(Sample,TMT,ProteinID,final_norm,box_cox_scaled_values)

  return(output)

}

