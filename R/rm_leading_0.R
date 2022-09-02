#' rm_leading_0
#'
#' useful function to remove leading 0s from well IDS.
#' This can be usefull when trying to map data from sample.IDs with leading 0s to other data without leading 0s.
#' Created it as a function so I don't have to look up the regular expression every time I use it.
#'
#' @param Wells: vector of well IDs
#'
#' @return vector of well IDs modified to remove the leading 0.
#' @export
#'
#' @examples
#' Wells = c("A01","A02", "B04")
#'
#' rm_leading_0(Wells)
#' ##Output###
#' "A1","A2","B4"

rm_leading_0 = function(Wells){

    output = stringr::str_remove(Wells,"(?<![0-9])0+")
    return(output)

}
