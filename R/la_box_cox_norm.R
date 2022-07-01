#' la_box_cox_norm
#' here is leigh-ana's method for conducting box_cox_normalization on proteomics data.
#'
#' @param filepath : This is the filepath to a file containing the protein names in the first column and other columns containing values for each sample.
#'
#' @return a normalized csv file in the same directory.
#' @export
#'
#' @examples
la_box_cox_norm = function(filepath){

  data <- read.csv(filepath,header=TRUE)

  # first column is protein names, other columns are values for each sample
  transformed_data <- matrix(data=NA,nrow=length(data$ProteinID),ncol=length(data)-1)
  rownames(transformed_data) <- data$ProteinID
  colnames(transformed_data) <- colnames(data)[2:length(data)]

  for (i in 2:(length(data))) {
    print(i)
    temp_data = data[,i]
    b <-MASS::boxcox(lm(temp_data ~ 1))
    lambda <- b$x[which.max(b$y)]
    new_data <- (temp_data^lambda - 1)/lambda
    scaled_data <- reshape::rescaler(new_data,type="range")
    scaled_data <- scaled_data/mean(scaled_data,na.rm = TRUE)
    transformed_data[,i-1] <- scaled_data
    car::qqPlot(scaled_data)
  }

  output_file = gsub(".csv","",filepath)

  write.table(transformed_data, paste0(output_file,"_normalized.csv"),sep=",")

}

