#' la_box_cox_norm
#' here is leigh-ana's method for conducting box_cox_normalization on proteomics data.
#'
#' @param data : This is a tibble containing at least the columns named Sample , TMT, and final_norm. This is intended to be the output of normalize_to_bridge() or noramlize_1plex().
#' Note that this expects data to be in the long format
#'
#'
#' @param data_format :This is a character string specifying whether you would like the data in long or wide format. Note that if data_format =  wide only data
#' corresponding to box_cox_scaled_values is returned. The intermediate normalization steps are disregarded in this case.
#'
#' @return a tibble containing the normalized data.
#' @export
#'
#' @examples
#' data = readr::read_delim("tests/testdata/combine_psm_fractions/PCB002_PSMs_Proteodiscover_output.txt") %>%
#' combine_psm_fractions() %>%
#' normalize_to_bridge(bridge_channel_plex = 126) %>%
#' la_box_cox_norm()
la_box_cox_norm = function(data,data_format = "long"){
  # Had to modify data to make compatible with leighana's script. Removed all nas, infinte values, and 0s.
  #Transformed to wide format to do lm in column format (didn't feel like figuring out how to do this within the tidyverse)
  mod_data = data %>%
    dplyr::select(Sample,TMT,ProteinID,final_norm) %>%
    dplyr::mutate(final_norm = dplyr::case_when(is.infinite(final_norm) ~ NA_real_,
                                                final_norm == 0 ~ NA_real_,
                                                TRUE ~ final_norm)) %>%
    tidyr::pivot_wider(names_from = c("Sample","TMT"),values_from = "final_norm") %>%
    as.data.frame()

  # Leigh-ana's script
  transformed_data <- matrix(data=NA,nrow=length(mod_data$ProteinID),ncol=length(mod_data)-1)
  rownames(transformed_data) <- mod_data$ProteinID
  colnames(transformed_data) <- colnames(mod_data)[2:length(mod_data)]

  for(i in 2:length(mod_data)) {
    temporary_data <-mod_data[,i]
    #Seems like the MASS function has a problem with the environment. I added y = TRUE and qr = TRUE based on this:
    #https://stackoverflow.com/questions/39728374/r-passing-linear-model-to-another-function-inside-a-function.
    b <- MASS::boxcox(lm(temporary_data ~ 1,y=TRUE, qr=TRUE),plotit = FALSE,interp = TRUE)
    lambda <- b$x[which.max(b$y)]
    new_data <- (temporary_data^lambda - 1)/lambda
    scaled_data <- reshape::rescaler(new_data,type="range")
    scaled_data <- scaled_data/mean(scaled_data,na.rm = TRUE)
    transformed_data[,i-1] <- scaled_data
  }

  #transforming final data into long data format
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

  # Adding option to export data in long or wide format
  if(data_format == "long"){
    return(output)
  }else if(data_format == "wide"){
    output2 = output %>%
      dplyr::select(Sample,TMT,ProteinID,box_cox_scaled_values) %>%
      tidyr::pivot_wider(names_from = c("Sample","TMT"), values_from = box_cox_scaled_values)
    return(output2)
  }else{
    print("format must be either long or wide")
  }

}

