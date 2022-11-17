#' annotate_pubchem
#'
#' this function allows you to search a vector containing terms against pubchem's API.
#'
#' @param input this is a list of the things to be searched
#'
#' @return a data frame with the input specified by user and output
#' @export
#'
#' @examples
#'#'pub_chem_smiles = "CN1C2CCC1C(C(C2)OC(=O)C3=CC=CC=C3)C(=O)OC"
#'annotate_pubchem(pub_chem_smiles)
annotate_pubchem = function(input){
  t = purrr::map(input,~tryCatch({
    n_count = 1
    pco = readr::read_lines(pubchem_url_request(.x),n_max = n_count)
    while(!grepl("[A-Z]",pco)){
      n_count = n_count + 1
      pco = readr::read_lines(pubchem_url_request(.x),n_max = n_count)[[n_count]]
    }
    data = dplyr::tibble(input = .x, output = pco,)
    return(data)
  },
  # upon error do this...
  error = function(.x){
    message(paste("request does not seem to exist:", .x))
    return(NA)
  }))
  #collapsing list into dataframe
  df=do.call("rbind", t)

  return(df)
}

