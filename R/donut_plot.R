#' donut_plot
#'
#' this plots a donut plot
#' @param data this is a tibble containing your data.
#' @param col this is a string containing the name of the column that you will be grouping by
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples
#' # Plotting a donut plot
#' data = readr::read_csv("tests/testdata/donut_plot/donut_plot_test_data.csv")
#' donut_plot(data)
donut_plot = function(data,col = "sample_type_group_3"){

  #counting n of observations in each.
  fdata = data %>%
    dplyr::group_by(.data[[col]]) %>%
    dplyr::summarise(count = n())
  #Total n
  total = sum(fdata$count)
  #Fractions
  fdata = fdata %>%
    dplyr::mutate(fraction = (count/total))
  # Compute the cumulative percentages (top of each rectangle)
  fdata$ymax=cumsum(fdata$fraction)
  # Compute the bottom of each rectangle
  fdata$ymin=c(0, head(fdata$ymax, n=-1))
  # Compute label position
  fdata$labelPosition = (fdata$ymax + fdata$ymin) / 2
  # Make the plot
  p1 = ggplot2::ggplot(fdata,ggplot2::aes_string(ymax="ymax", ymin="ymin", xmax=4, xmin=3,fill=col)) +
    ggplot2::geom_rect() +
    viridis::scale_fill_viridis(discrete = TRUE)+
    ggrepel::geom_label_repel(x=3.5,ggplot2::aes_string(label = col,y = "labelPosition"),size=4)+
    ggplot2::coord_polar(theta="y") +
    ggplot2::xlim(c(2, 4)) +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "none")

  return(p1)

}
