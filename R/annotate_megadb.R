#' annotate_megadb
#'
#' annotate proteins against the megadb (hosted as a postgres sql database on AWS).
#'
#' @param proteinids
#'
#' @return a tibble with proteinid, best_taxa_level match, description, and the database from which the annotation was derived.
#' note food and human are derived from uniprot.
#' @export
#'
#' @examples
annotate_megadb = function(proteinids){
#Connecting to megadb AWS postgresql database
connec = DBI::dbConnect(RPostgres::Postgres(),
               host = "megadb.c86w0dpwhuey.us-west-1.rds.amazonaws.com",
               port = 5432,
               user = "postgres",
               password = "megadbpass")

#Defining Negate function for later use.
`%!in%` = Negate(`%in%`)

#Determining which proteins are in which database based on naming conventions.
phage_ids = proteinids[grepl("^uvig.*",proteinids)]
bacteria_ids = proteinids[grepl("^MGY.*",proteinids)]

#Defining uniprot ids as being those that are not in the phage or bacteria database
uniprot_ids = proteinids[proteinids %!in% c(phage_ids,bacteria_ids)]

#Searching against phage database
phage_results = dplyr::tbl(connec,"phage") %>%
  dplyr::select(proteinid,best_tax_level,description = eggnog_free_text_description) %>%
  dplyr::filter(proteinid %in% phage_ids) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(database = "phage")

#Searching against bacteria database
bacteria_results = dplyr::tbl(connec,"bacteria") %>%
  dplyr::select(proteinid,best_tax_level = max_annot_lvl,description) %>%
  dplyr::filter(proteinid %in% bacteria_ids) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(database = "bacteria")

#Searching against uniprot API
uniprot_results = GLabR::annotate_uniprot(uniprot_ids) %>%
  select(proteinid = Entry,best_tax_level = Organism,description = Protein.names) %>%
  #If it isn't human, it is food
  dplyr::mutate(database = dplyr::case_when(grepl("Homo sapiens.*",best_tax_level) ~ "human",
                                            TRUE ~ "food"))
#Combining data into cohesive table
#Contains proteinid, taxa level, and description.
output = dplyr::bind_rows(uniprot_results,bacteria_results,phage_results) %>%
  dplyr::mutate(common_name = stringr::str_extract(string = best_tax_level,pattern = "\\(.*?\\)"),
                common_name = gsub("\\(|\\)","",common_name)) %>%
  dplyr::select(proteinid,database,best_tax_level,common_name,description)

#Output
return(output)
}





