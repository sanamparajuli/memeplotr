#' A clean theme for motif figures
#'
#' A light `ggplot2` theme with no panel grid in the motif direction, left
#' aligned sequence labels and a horizontal legend. Every element can be
#' overridden afterwards with the usual `+ theme(...)`.
#'
#' @param base_size Base font size in points.
#' @param base_family Base font family.
#' @param grid One of `"x"` (vertical guides only, the default), `"none"`,
#'   `"y"`, `"both"`.
#' @param legend Legend position, passed to [ggplot2::theme()].
#' @param axis_text_y_face Font face for sequence labels, e.g. `"italic"`.
#'
#' @return A `ggplot2` theme object.
#' @examples
#' res <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
#' gg_motif_map(res) + theme_motif(base_size = 13, grid = "none")
#' @export
theme_motif <- function(base_size = 11, base_family = "",
                        grid = c("x", "none", "y", "both"),
                        legend = "bottom", axis_text_y_face = "plain") {
  grid <- match.arg(grid)
  gx <- grid %in% c("x", "both")
  gy <- grid %in% c("y", "both")
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major.x = if (gx) ggplot2::element_line(linewidth = 0.25, colour = "grey88") else ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = if (gy) ggplot2::element_line(linewidth = 0.25, colour = "grey88") else ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(hjust = 1, face = axis_text_y_face,
                                          size = ggplot2::rel(0.95)),
      axis.ticks.x = ggplot2::element_line(linewidth = 0.3, colour = "grey50"),
      axis.title.y = ggplot2::element_blank(),
      legend.position = legend,
      legend.key.height = grid::unit(0.8, "lines"),
      legend.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.9)),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.title.position = "plot",
      strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.9))
    )
}

#' @rdname theme_motif
#' @export
theme_tree_motif <- function(base_size = 11, base_family = "") {
  ggplot2::theme_void(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(5, 2, 5, 5),
      legend.position = "none"
    )
}
