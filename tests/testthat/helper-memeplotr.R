## The motif blocks are drawn as polygons, so one block spans several rows of
## built data; identify the layer by its fill aesthetic and count groups.
block_data <- function(p) {
  bl <- ggplot2::ggplot_build(p)$data
  cand <- Filter(function(d) "fill" %in% names(d) && "group" %in% names(d), bl)
  if (!length(cand)) return(bl[[1]])
  cand[[which.max(vapply(cand, function(d) length(unique(d$group)), integer(1)))]]
}

n_blocks <- function(p) length(unique(block_data(p)$group))

block_fills <- function(p) {
  d <- block_data(p)
  unique(d$fill[!is.na(d$fill)])
}

fill_labels <- function(p) {
  sc <- ggplot2::ggplot_build(p)$plot$scales$get_scales("fill")
  if (is.null(sc)) return(character())
  as.character(sc$get_labels())
}

y_labels <- function(p) {
  pp <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]
  as.character(pp$y$get_labels())
}
