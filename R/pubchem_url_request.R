#' pubchem_url_request
#'
#' This function generates a pubchem url request for user specified info. designed to be used in combination with the
#' @param input this is the input
#' @param search1 this is the first search string
#' @param search2 this is the second search
#' @param operation this is the operation that you would like to perform on the data
#' @return a url that hits pubchem's api and can be used to search for user specified data.
#' @export
#'
#' @examples
#'
#'
pubchem_url_request = function(input,
                            search1 = "compound",
                            search2 = "smiles",
                            operation = "synonyms"){
  base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
  request = paste0(base_url,search1,"/",search2,"/",input,"/",operation,"/TXT")
  return(request)
}





