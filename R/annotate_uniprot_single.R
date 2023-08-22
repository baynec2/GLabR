#' annotate_uniprot_single
#' This function uses Uniprot's API to annotate proteins one protein at a time. This was made because of problems with the batch version.
#' @param protein_list This is a list of the protein accessions that you would like to annotate
#' @param columns This is a string  of the column names that you would like to include followed by commas.
#' Column names accepted by the API can be found here: https://www.uniprot.org/help/return_fields.
#' If not specified, the default columns are returned.
#'
#' @return A tibble containing the protein information for the specified accessions.
#' @export
#'
#' @examples
#' protein_list = c("P00642")
#' annotate_uniprot_single(protein_list)
annotate_uniprot_single = function(protein_list,columns = NULL){
  # BaseURL for request
  baseUrl <- "https://rest.uniprot.org/uniprotkb/search?query=accession:"
  # Adding OR term to do multiple proteins per request

  #Initializing a data frame to hold results
  output = data.frame()
  for( i in 1:length(protein_list)){
    #Iterating through the batches
    temp_list = protein_list[[i]]
    temp_string_list =  paste(temp_list, sep="", collapse="")
    # Defining columns that we want to return. Defaults to default columns.
    if(is.null(columns)){
      request = paste0(baseUrl,temp_string_list,"&format=tsv","&size=",length(temp_list))
    }else{
      request = paste0(baseUrl,temp_string_list,"&format=tsv","&fields=",columns,"&size=",length(temp_list))
    }
    #Note this error handling will mean that the entire batch gets assigned as NA, not only any entries that may not exist.
    #This is a problem, but it only seems to affect a very small amount of entries so I have decided to leave it as is.
    returned_data = tryCatch({
      # Try
      read.csv(request, header = TRUE, sep = '\t')
    },
    # upon error do this...
    error = function(cond){
      message(paste("request does not seem to exist:", request))
      message(paste("n_iteration: ", i))
      return(NA)
    }
    )
    output = rbind(output,returned_data)
  }

  output = tibble::as_tibble(output)

  return(output)

}
