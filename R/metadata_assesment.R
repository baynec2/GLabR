#' metadata_assesment(): this is a R function designed to replicate the analysis found in the paper:
#' Mortality Risk Profiling of Staphylococcus aureus Bacteremia by Multi-omic Serum Analysis Reveals Early Predictive
#' and Pathogenic Signatures
#'
#' Brief description per the method sections below.
#' First, the metadata was split into two groups, categorical metadata and continuous metadata.
#' Categorical metadata associations were determined using MWU test (2 categories) or Kruskall-Wallis test (> 2 categories).
#' Continuous metadata associations were determined using Pearson correlation.
#'
#' @param data = a tibble contining data and metadata from a. filter it on your own to only contain the things you are interested in assessing.
#' @param sample_id = name of the column containing a sample identifier (unquoted)
#' @param metric  = name of the column contining the metric of interest (gene, row_id etc- unquoted)
#' @param features = a vector of features of the metric of interest contained in the data of interest you are interested in
#' @param outcome = the outcome metric you are interested in measuring (unquoted)
#' @param outcome_c = the coutcome metric you are interested in mearuring (quoted)
#'
#' @return a tibble conting the unadjusted p values from the test.
#' @export
#'
#' @examples
metadata_assesment = function(data,
                              sample_id = sample_id,
                              metric = row_id,
                              features,
                              outcome = norm_value,
                              outcome_c = "norm_value"){


  # Splitting the inital data frame into the data columns and the
  input = data %>%
    dplyr::select({{sample_id}},{{metric}},{{outcome}})

  md_in_all = data %>%
    dplyr::select({{sample_id}},{{metric}})

  md = data %>%
    dplyr::select(-{{sample_id}},-{{metric}},-{{outcome}})

  # What columns are categorical
  categorical_cols = purrr::map_lgl(md,~!is.numeric(.))

  categorical_md = md %>%
    dplyr::select(which(categorical_cols))

  # number of distinct things in each categorical category
  col_numbers_in_cat = categorical_md %>%
    purrr::map_dbl(~length(unique(na.omit(.))))

  ##Of the categories that are categorical, which have 2 categories
  # Sample ID should be in the first column. Adding it here.
  cols_2 = categorical_md[,which(col_numbers_in_cat == 2)]

  ##Of the categories that are categorical, which have more than 2 categories.
  cols_g_2_cat = categorical_md[,which(col_numbers_in_cat > 2)]

  # What columns are continuous
  continuous_cols = purrr::map_lgl(md,is.numeric)

  continuous_md = md[,which(continuous_cols)]


  input_md_continuous = dplyr::bind_cols(input,continuous_md)
  input_md_2_col = dplyr::bind_cols(input,cols_2)
  input_md_g_2_col = dplyr::bind_cols(input,cols_g_2_cat)


  # Lets do the continuous stats first #
  stats = data.frame()

  # Outer loop with the different variables to test
  for(x in names(which(continuous_cols))){

    for(i in features){
      loop = dplyr::filter(input_md_continuous,{{metric}} == i)
      #Let's look at the continuous stuff first
      c_mod = loop %>%
        rstatix::cor_test({{outcome}},x,method = "pearson",use="complete.obs") %>%
        dplyr::select(variable = var2,p) %>%
        dplyr::mutate({{metric}} := i)

      stats = bind_rows(stats,c_mod)
    }
  }

  # Outer loop
  for(x in names(cols_2)){
    # Categorical two categories
    for(i in features){
      loop = dplyr::filter(input_md_2_col,{{metric}} == i)

      #Let's do the stats for the things with only two values first.
      inner_loop = rstatix::wilcox_test(data = loop,as.formula(paste0(outcome_c," ~ ",x))) %>%
        dplyr::mutate(variable = x,
                      {{metric}} := i) %>%
        select(variable,p,{{metric}})

      stats = bind_rows(stats,inner_loop)

    }

  }

  # Outer loop
  for(x in names(cols_g_2_cat)){
    # Categorical two categories
    for(i in features){
      loop = dplyr::filter(input_md_g_2_col,{{metric}} == i)

      #Let's do the stats for the things with more than two values now.
      inner_loop = rstatix::kruskal_test(data = loop,as.formula(paste0(outcome_c," ~ ",x))) %>%
        dplyr::mutate(variable = x,{{metric}} := i) %>%
        dplyr::select(variable,p,{{metric}})


      stats = dplyr::bind_rows(stats,inner_loop)

    }

  }
  return(stats)
}

