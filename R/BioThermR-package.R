#' @keywords internal
"_PACKAGE"

#' @keywords internal
"_PACKAGE"

#' @import ggplot2
#' @importFrom stats median sd quantile IQR density runif setNames
#' @importFrom grDevices colorRampPalette
#' @importFrom ggrepel geom_text_repel geom_label_repel
#' @importFrom plotly ggplotly plot_ly layout
#' @importFrom ggsci scale_fill_npg scale_color_npg
NULL

utils::globalVariables(c(
  ".data",
  "final_x", "final_y", "val", "x", "y", "label",
  "temp", "x_shifted", "y_shifted","auto","manual"
))
