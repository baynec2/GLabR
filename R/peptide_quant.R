#' peptide_quant
#' this function takes the export from the plate reader and formats data into a tidy format while calculating predicted concentrations of unknowns using a linear model of the standard curve.
#'
#' @param filepath this is the filepath to an exported .txt file from either the peptide_quant_od480_single_plate.spr or peptide_quant_od480_single_plate.spr protocols.
#' Note that when doing the export for this file you need to only check the standards and unknowns boxes.
#'
#' @return a tibble containing processed data. Note that predicted values are calculated for unknown sample concentrations using a linear model.
#' @export
#'
#' @examples
peptide_quant = function(filepath){
  #Reading in the file from plate reader. Note that it has a weird encoding.
  data = read_tsv(filepath,
                  locale = locale(encoding ="utf-16le"),
                  skip =2,
                  col_types = "ccncnccc" )
  #number of  rows for plate # detectopn
  nrow = nrow(data)
  if(nrow <= 106){
    #If there is only one plate, code below is executed (one plate will only have 106 rows max)
    #Subseting standards section
    standards = data[1:24,] %>%
      dplyr::mutate(Sample_Type = "standards",
                    Conc = as.numeric(Conc)) %>%
      dplyr::select(Sample_Type,Well = Wells, Conc_mg_mL = Conc, OD480 = Value)

    #Pulling blank value
    blank = standards %>%
      dplyr::filter(Conc_mg_mL == 0) %>%
      dplyr::pull()

    # Taking unknowns
    unknowns = data[31:102,]
    #fixing column names
  }else{
    #If there is data from more than one plate in the file, code below is executed
    standards = data[1:24,] %>%
      dplyr::mutate(Sample_Type = "standards",
                    Conc = as.numeric(Conc)) %>%
      dplyr::select(Sample_Type,Well = Wells, Conc_mg_mL = Conc, OD480 = Value)

    unknowns = data[35:nrow-4,]

  }

  colnames(unknowns) =c("Sample","Well","Value","MeanValue","R","Result","SD","CV","..9")

  unknowns = unknowns %>%
    dplyr::mutate(Sample_Type = "unknown") %>%
    dplyr::select(Plate = Sample,Well,OD480=Value,Sample_Type) %>%
    tidyr::fill(Plate) %>%
    dplyr::mutate(Blank_OD480 = OD480 - blank,
                  Well = as.character(Well))
  #Grouping the standards together.
  standards = standards %>%
    dplyr::mutate(Blank_OD480 = OD480 - blank)

  #Constructing the linear model
  lm = lm(Conc_mg_mL ~ Blank_OD480,data = std_sum)

  predicted_conc = predict(lm,newdata = unknowns)

  unknowns= unknowns %>%
    dplyr::bind_cols(Conc_mg_mL = predicted_conc)

  all = dplyr::bind_rows(standards,unknowns) %>%
    dplyr::mutate(Plate = dtringr::str_remove(Plate,"^0")) %>%
    dplyr::select(Plate,Well,Conc_mg_mL,OD480,Blank_OD480,Sample_Type)

  return(all)
}




