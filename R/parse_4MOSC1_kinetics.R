#' parse_4MOSC1_kinetics
#'
#' This function is designed to parse the excel files the Gutkind lab uses to
#' keep track of their mouse kinetic data into tidy data.
#'
#' Will output a tibble with the following columns
#' 1. experiment (this is the name of the tab)
#' 2. mouse_id this is the identifier of the mice
#' 3. cage_id: this is the experiment/cage_number combination
#' 4. cage_number: this is the per experiment cage number (cage number is calc
#' by assigning first five mice to cage 1 and so on. If the total number of mice
#' is not divisable by 5, the last cage gets less mice).
#' 5. measurement: the measurement taken
#' 6. value: the value of the measurement.
#' 7. responder_status: the assigned responder status based on final tumor volume
#       * 0-5 = responder
#       * 5-115 = stable
#       * 115 = non responder
#'
#' @param filepath this is the filepath to the completed excel sheet template
#' containing the data to parse
#' @param sheets this is a description of the sheets that you want to parse
#' can either be a vector of numbers or names.
#'
#' @return a tidy tibble
#' @export
#'
#' @examples
#' parsed_data = parse_4MOSC1_kinetic("tests/testdata/parse_4MOSC1_kinetics/LiBIO_kinetics.xlsx",
#' sheets = c(2,3,6))
parse_4MOSC1_kinetic = function(filepath,sheets){
# Defining the parse function
parse = function(sheet){
#Reading in the data
data = readxl::read_excel(filepath,sheet)
#Defining the mouse_ids
mouse_id = data[2:nrow(data),2]
colnames(mouse_id) = "mouse_id"
#Finding out how many samples there are
is_not_na = !is.na(mouse_id)
n_samples = sum(is_not_na)

mouse_id = mouse_id[1:n_samples,] %>%
  as.data.frame()

#Finding out how many days there are
colnames = colnames(data)

days = stringr::str_extract(colnames,"day.*")
days = days[!is.na(days)]
days = gsub("day ","",days)

days = data.frame(day = days)

#Adding cage_number
cage_number = rep(1:(floor(nrow(mouse_id)/5)),5)
remainder = nrow(mouse_id) - length(cage_number)
if(remainder > 0){
add = 1:remainder
cage_number = c(cage_number,add) %>% sort()
} else{
cage_number = cage_number %>% sort()
}
#Adding response status


#Extracting the sheet name
current_sheet = readxl::excel_sheets(filepath)[sheet]

#Binding the sheet name to to make unique cage numbers across experiments
experiment = current_sheet

#Parsing the measurements
measurements = as.character(data[1,])

#getting the colums for neck circumference
NC_index = which(grepl("Neck Circumf",measurements))
#getting the column indexes we need
L_index = NC_index + 1
W_index = L_index + 1
V_index = W_index + 1

#extracting the data
NC = data[2:(nrow(mouse_id)+1),NC_index]
L = data[2:(nrow(mouse_id)+1),L_index]
W = data[2:(nrow(mouse_id)+1),W_index]
V = data[2:(nrow(mouse_id)+1),V_index]

# Flipping data structure
data_w = data.frame(neck_circumference = NC, l_mm = L, W_mm = W, vol_mm3 = V) %>%
  tidyr::pivot_longer(1:length(.)) %>%
  mutate(name = gsub("\\..*","",name))

# Adding missing data
data_w$days = as.numeric(rep(as.character(days[,1]),nrow(mouse_id)* 4))

data_w$mouse_id = rep(as.character(mouse_id[,1]),each = nrow(days) * 4)

data_w$cage_number = as.numeric(rep(cage_number,each = nrow(days) * 4))

data_w$value = as.numeric(data_w$value)
#assigning responder status based on final volume
# 0-5 = responder
# 5-115 = stable
# 115 = non responder
max_days = max(as.numeric(as.character(days[,1])))
assignment = data_w %>%
  dplyr::filter(days == max_days,
                name == "vol_mm3") %>%
  dplyr::mutate(responder_status = case_when(value < 5 ~ "responder",
                                             value >= 5 & value <= 115 ~ "stable",
                                             TRUE ~ "non responder")) %>%
  dplyr::select(mouse_id,responder_status)

data_f = data_w %>%
  dplyr::mutate(experiment = experiment,
                cage_id = paste0(experiment,"_",cage_number)) %>%
  dplyr::select(experiment,mouse_id,cage_id,cage_number,day = days,measurement = name, value)

#adding the responder status
data_f = dplyr::inner_join(data_f,assignment,by = "mouse_id")
return(data_f)
}

#iterating through all of the sheets
out = purrr::map_df(sheets,parse)

return(out)
}
