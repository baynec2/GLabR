#' gibson_calculator
#' This function is useful for calculating how to set up the Gibson reaction given user input lists of parts,concentrations, and lengths.
#'
#' @param names a vector of the names of each of the parts. Make sure the vector is first, followed by the inserts.
#' @param concentrations a vector of concentrations. Make sure that they are in the same order as names
#' @param lengths a vector of lengths. Make sure these are in the same order as the names
#'
#' @return
#' a tibble containing all the ingredients at the max molarity according to the NEB protocol that is possible given the constraints imposed by the concentration of the parts. .
#'
#' @examples
#' concentrations = c(100,10,10,10)
#' lengths = c(5000,1000,1000,1000)
#' names = c("pUC19","insert_1","insert_2","insert_3")
#'
#' gibson_calculator(concentrations,lengths,names)
#'
#'
concentrations = c(100,10,10,10)
lengths = c(5000,1000,1000,1000)
names = c("pUC19","CB001","CB002","CB003")
gibson_calculator = function(concentrations,lengths,names){

 num_of_pieces = length(lengths)
  #determining what ratio to use
  if(num_of_pieces <= 3){
    ratio = 0.5
    max_vector_pmol = (0.2 / num_of_pieces) * ratio
    max_insert_pmol =   (0.2 - max_vector_pmol)/ (num_of_pieces-1)

  }else{
    ratio = 1
    max_vector_pmol = 0.5/num_of_pieces
    max_insert_pmol  = 0.5/num_of_pieces
  }

 #determining mass of DNA needed to get the max molarity
 #moles of dsDNA (mol) x ((length of dsDNA (bp) x 617.96 g/mol/bp) + 36.04 g/mol)
  vector_mass_for_max_molarity = max_vector_pmol * 1E-12 * (lengths[1] * 617.96 + 36.04) /1E-9
  insert_masses_for_max_molarity = max_insert_pmol * 1E-12 * (lengths[2:num_of_pieces] * 617.96 + 36.04) / 1E-9

  #can we achieve the max given volume constraints?
  uL_vector = vector_mass_for_max_molarity / concentrations[1]
  uL_insert = insert_masses_for_max_molarity / concentrations[2:num_of_pieces]

  total_uL = uL_vector + sum(uL_insert)

  #Adjusting the total volume to be 10 if it is not.
  if(total_uL > 10){
    adjustment = total_uL/10
    uL_vector = uL_vector / adjustment
    uL_insert = uL_insert / adjustment
    vector_pmol = max_vector_pmol / adjustment
    insert_pmol = max_insert_pmol / adjustment
  }else{
    vector_pmol = max_vector_pmol
    insert_pmol = max_insert_pmol
  }

  pmol = c(vector_pmol,rep(insert_pmol,num_of_pieces-1))
  uL = c(uL_vector,uL_insert)

  df = tibble::tibble(names,pmol,uL)

  #adding other components
  water_uL = 10-  sum(df$uL)
  water = list(names = "water",pmol = NA_real_,uL = water_uL)
  MM = list(names = "2x MM", pmol = NA_real_, uL = 10)

  df2 = dplyr::bind_rows(df,water,MM)

  return(df2)
}
