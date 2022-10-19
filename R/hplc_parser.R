#' hplc_parser
#'
#' This function is designed to take the output of a .txt file that is exported from chromeleon and format it into a data frame.
#' To get the most out of this function, the files should be named with the following conventions:
#' 1. an index representing the order should be placed at the beginning of the file ie 00_plex1
#' 2. blank samples should be run between each experimental sample and should be labeled index_BLANK.
#' 3. plexes should be names in a way that is easy to map with metadata. I use the convention PCB001 (Plex, Charlie Bayne, plex # 0001)
#'
#' @param file_name: this is the file name
#'
#' @return a formated tibble with chromatography data
#' @export
#'
#' @examples
#' chromatography_data = hplc_parser("tests/testdata/hplc_parser/CB004")
#'
hplc_parser = function(directory){
  #Listing files in the directory
  files = list.files(directory,full.names = TRUE)
  #Defining parsing function
  parse = function(file_path){
  #Extracting file_name from file path
  file_name = basename(file_path) %>%
    gsub(".txt","",.)
  #Reading files and tidying up
  d1 = readr::read_delim(file_path,skip = 42) %>%
    dplyr::mutate(order = gsub("_.*", "" ,file_name),
                  plex = gsub("^.*_","",file_name),
                  type = dplyr::case_when(grepl("*BLANK.*",file_name)~ "BLANK",
                                   TRUE ~ "experimental")) %>%
    #fixing columnn names
    dplyr::rename(time_min = `Time (min)`,
           step = `Step (s)` ,
           mau_value = `Value (mAU)`)
  return(d1)
  }
  # iterating through all of the files
  d1 = purrr::map_df(files,parse)
  return(d1)
}
