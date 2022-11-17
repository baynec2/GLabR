#' annotate_food
#' This function can used to annotate foods into higher groupings.
#' ie a pig can be classified as red meat, animal, etc depending on the level of granularity.
#'
#' @param common_names a list of the common names in a database
#'
#' @return a tibble  containing
#' @export
#'
#' @examples
#' #list of protein IDs
#' proteinids = read_csv("tests/testdata/megadb/proteinids.csv") %>%
#'  pull(ProteinID)
#'
#' megadb_results = annotate_megadb(proteinids)
#'
#' #pulling common names
#' common_names = megadb_results %>%
#' dplyr::pull(common_name) %>%
#' dplyr::unique()
#'
#' #Now we can figure out what these foods are
#' foods = annotate_food(common_names)
annotate_food = function(common_names){

  connec = DBI::dbConnect(RPostgres::Postgres(),
                          host = "megadb.c86w0dpwhuey.us-west-1.rds.amazonaws.com",
                          port = 5432,
                          user = "postgres",
                          password = "megadbpass")


  output = dplyr::tbl(connec,"food_variety") %>%
    dplyr::filter(common_name %in% common_names) %>%
    dplyr::as_tibble()

  return(output)
}
