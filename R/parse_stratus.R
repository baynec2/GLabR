#' parse_stratus
#'
#' This function parses the output from the stratus plate reader by cerillo (https://cerillo.bio/stratus/) into a usable form.
#' It takes the  csv file exported by the plate reader and transforms it into a tidy data format.
#'
#' @param input this is the .csv file exported from the stratus plate reader.
#'
#' @return
#' a tibble with the following columns.
#' time_hr -> the amount of time in hrs that have elapsed since the start of the read
#' datetime -> the date time of the measurement.
#' temperature -> the temperature that was recorded at the time of the measurement in C
#' well -> the well corresponding to a measurement.
#' paste0("od",wavelength) -> values of the measured wavelength. Our machine measured 600nm wavelength so default is od600.
#' @export
#'
#' @examples
#' parsed_OD600_data = parse_stratus("tests/testdata/parse_stratus/stratus_data_export.csv")
parse_stratus = function(input,wavelength_nm =  600){
  #Necessary to specify this encoding because of the degree sign special character
  #Otherwise, reading it in will error.
  data = read.csv(input,skip = 9,fileEncoding = "ISO-8859-1") %>%
    dplyr::rename(unix_time =1, temperature = 2) %>%
    dplyr::mutate(datetime = lubridate::as_datetime(unix_time)) %>%
    select(datetime, temperature, dplyr::everything(), -unix_time)

  #Defining the start time as first value in datetime column/
  start_time = data$datetime[[1]]

  #Calculating amount of time past and tidying the data.
  output = data %>%
    dplyr::mutate(time_difference = datetime - start_time,
                  time_difference = as.numeric(time_difference,units = "hours")) %>%
    dplyr::select(time_hr = time_difference, datetime,temperature, dplyr::everything()) %>%
    tidyr::pivot_longer(4:length(.),names_to = "well",values_to = paste0("od",wavelength_nm))

  return(output)
}
