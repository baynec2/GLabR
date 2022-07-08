#' annotate_proteins
#' This function uses Uniprot's API to annotate proteins.
#' @param protein_list This is a list of the protein accessions that you would like to annotate
#' @param columns This is a string  of the column names that you would like to include followed by commas.
#' Column names accepted by the API can be found here: https://www.uniprot.org/help/return_fields.
#' If not specified, the default columns are returned.
#'
#' @return A tibble containing the protein information for the specified accessions.
#' @export
#'
#' @examples
annotate_proteins = function(protein_list,columns = NULL){
  # BaseURL for request
  baseUrl <- "https://rest.uniprot.org/uniprotkb/search?query=accession:"
  # Adding OR term to do multiple proteins per request

  if(length(protein_list) > 300){
  master_list = paste0(protein_list,"+OR+")
  #Determining how many times we will have to iterate to get all accessions on initial list. 300 seems to be a good number per batch
  num_iterations = floor(length(master_list)/300)
  #calculating the end of each iteration
  end = c((1:num_iterations) * 300,length(master_list))
  #calculating the start of each iteration
  start = c((end - 299)[1:length(end) - 1],end[length(end)-1] +1)
  }else{
    master_list = paste0(protein_list,"+OR+")
    start = 1
    end = length(protein_list)
}
  #Initializing a data frame to hold results
  output = data.frame()
    for( i in 1:length(start)){
      #Iterating through the batches
      temp_list = master_list[start[i]:end[i]]
      #Getting rid of last OR term
      temp_list[length(temp_list)] = gsub("\\+OR\\+","",temp_list[length(temp_list)])
      temp_string_list =  paste(temp_list, sep="", collapse="")
      # Defining columns that we want to return. Defaults to default columns.
      if(is.null(columns)){
        request = paste0(baseUrl,temp_string_list,"&format=tsv","&size=",length(temp_list))
      }else{
        request = paste0(baseUrl,temp_string_list,"&format=tsv","&fields=",columns,"&size=",length(temp_list))
      }
      returned_data = read.csv(request, header = TRUE, sep = '\t')
      output = rbind(output,returned_data)
    }

  output = tibble::as_tibble(output)

  return(output)

}


