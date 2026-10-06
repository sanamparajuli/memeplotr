#' Tip order of a phylogeny
#'
#' The sequence names in the order [gg_phylo()] draws them, top to bottom.
#' Use it to put your own plot on the same vertical order as a tree panel
#' before composing the two with [gg_tree_panel()].
#'
#' @param tree An `ape::phylo` object, Newick file, or Newick string.
#' @param ladderize,right Passed to the internal layout; must match the values
#'   used for the tree panel.
#'
#' @return A character vector of tip labels, first element at the top.
#' @examples
#' tree_tip_order(system.file("extdata", "example_tree.nwk", package = "memeplotr"))
#' @export
tree_tip_order <- function(tree, ladderize = TRUE, right = TRUE) {
  .tip_order(read_tree(tree), ladderize = ladderize, right = right)
}

#' Compose a tree panel with any tip-ordered plot
#'
#' Places a [gg_phylo()] panel to the left of another plot and aligns the two
#' vertically. The right-hand plot must put tip *i* of
#' [tree_tip_order()] at y = `n - i + 1` and use the same `expand_y`, which is
#' what [gg_motif_map()], [gg_motif_heatmap()] and [gg_motif_counts()] do when
#' given `order = "tree"`.
#'
#' @param tree_plot A `ggplot` from [gg_phylo()].
#' @param panel A `ggplot` to place on the right, or a list of them.
#' @param widths Relative panel widths, recycled over
#'   `1 + length(panel)` panels.
#' @param guides Passed to [patchwork::plot_layout()]; `"collect"` merges the
#'   legends of all panels.
#' @param ... Further arguments for [patchwork::plot_layout()].
#'
#' @return A `patchwork` object.
#' @examples
#' res <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
#' nwk <- system.file("extdata", "example_tree.nwk", package = "memeplotr")
#' gg_tree_panel(
#'   gg_phylo(nwk),
#'   gg_motif_map(res, order = "tree", tree = nwk, keep_empty = TRUE) +
#'     ggplot2::theme(axis.text.y = ggplot2::element_blank())
#' )
#' @seealso [gg_tree_motifs()]
#' @export
gg_tree_panel <- function(tree_plot, panel, widths = c(1, 2.4),
                          guides = "collect", ...) {
  .need("patchwork")
  panels <- if (inherits(panel, "ggplot")) list(panel) else panel
  ps <- c(list(tree_plot), panels)
  widths <- rep_len(widths, length(ps))
  patchwork::wrap_plots(ps, nrow = 1) +
    patchwork::plot_layout(widths = widths, guides = guides, ...)
}

#' Motif map aligned to a phylogeny
#'
#' The figure MEME cannot draw: a motif composition diagram with the sequences
#' ordered and aligned to a phylogenetic tree read from a Newick file. Both
#' halves are ordinary `ggplot2` objects composed with `patchwork`, so you can
#' restyle either one, or pull them apart and rebuild the layout yourself.
#'
#' Tips present in the tree but absent from the motif output keep an empty
#' track, so the rows always line up with the branches. Sequences in the other
#' direction - present in the MEME output but not in the tree - are dropped with
#' a warning naming them, because the tree defines the rows; a row the tree does
#' not have would shift every label off the track it names. A warning here
#' almost always means a tip label and a sequence name differ in spelling.
#'
#' @param x A [meme_result], `motif_table`, or data frame of hits.
#' @param tree An `ape::phylo` object, path to a Newick file, or a Newick
#'   string.
#' @param widths Relative widths of the tree and map panels.
#' @param show_labels Where to draw sequence names. `"axis"` (default) puts
#'   them in the tree panel's right-hand axis, where grid reserves their full
#'   width, so long names can never overlap the motif blocks. `"tree"` draws
#'   them inside the tree panel next to each tip - compact, but names longer
#'   than the room given by `gg_phylo(xlim_mult =)` will run over the map.
#'   `"map"` uses the map panel's y axis; `"none"` omits them. Named
#'   `show_labels` rather than `labels` so that `label = TRUE` for block labels
#'   passes cleanly through `...` to [gg_motif_map()].
#' @param gap Extra space in points between the tree panel (labels included)
#'   and the motif map.
#' @param expand_y Vertical padding, in track units, applied to both panels.
#'   Raise it if blocks are clipped at the top or bottom.
#' @param ladderize,right,cladogram,align_tips,branch_colour,branch_linewidth,scale_bar
#'   Passed to [gg_phylo()].
#' @param tree_args Further arguments for [gg_phylo()], as a named list.
#' @param title,subtitle Title and subtitle for the composed figure.
#' @param guides Legend handling, passed to [patchwork::plot_layout()].
#' @param ... Passed to [gg_motif_map()].
#'
#' @return A `patchwork` object. Component plots are available as
#'   `attr(p, "tree")` and `attr(p, "map")`.
#' @examples
#' res <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
#' nwk <- system.file("extdata", "example_tree.nwk", package = "memeplotr")
#' gg_tree_motifs(res, nwk)
#' gg_tree_motifs(res, nwk, shape = "hexagon", palette = "muted",
#'                align_tips = TRUE, widths = c(1, 3))
#' @seealso [gg_phylo()], [gg_motif_map()], [gg_tree_panel()]
#' @export
gg_tree_motifs <- function(x, tree, widths = c(1, 2.4),
                           show_labels = c("axis", "tree", "map", "none"),
                           gap = 6,
                           expand_y = 0.6, ladderize = TRUE, right = TRUE,
                           cladogram = FALSE, align_tips = FALSE,
                           branch_colour = "grey25", branch_linewidth = 0.5,
                           scale_bar = FALSE, tree_args = list(),
                           title = NULL, subtitle = NULL, guides = "collect",
                           ...) {
  .need("patchwork")
  show_labels <- match.arg(show_labels)
  phy <- read_tree(tree)

  tp <- do.call(gg_phylo, utils::modifyList(
    list(tree = phy, ladderize = ladderize, right = right, cladogram = cladogram,
         align_tips = align_tips, branch_colour = branch_colour,
         branch_linewidth = branch_linewidth, scale_bar = scale_bar,
         tip_labels = show_labels %in% c("tree", "axis"),
         tip_label_position = if (show_labels == "axis") "axis" else "panel",
         expand_y = expand_y),
    tree_args))
  ## keep the attributes gg_phylo set; `+` on a ggplot drops them
  if (gap > 0) {
    tat <- attributes(tp)[c("tip_order", "n_tip")]
    tp <- tp + ggplot2::theme(plot.margin = ggplot2::margin(5, 2 + gap, 5, 5))
    attr(tp, "tip_order") <- tat$tip_order
    attr(tp, "n_tip") <- tat$n_tip
  }

  mp <- gg_motif_map(x, order = "tree", tree = phy, keep_empty = TRUE,
                     drop_untreed = TRUE, expand_y = expand_y, ...)
  ## both panels must span the same rows, or patchwork aligns two different
  ## y scales and the labels drift off the tracks they name
  if (length(attr(mp, "seq_levels")) != length(attr(tp, "tip_order"))) {
    .abort("Tree and map panels disagree on row count; this is a bug - please report it.")
  }
  if (show_labels != "map") {
    mp <- mp + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                              axis.ticks.y = ggplot2::element_blank())
  }

  tips <- attr(tp, "tip_order")
  lv <- attr(mp, "seq_levels")
  if (!identical(as.character(lv)[seq_along(tips)], as.character(tips))) {
    rlang::warn("Tree tip order and map track order disagree; check tip labels against sequence names.")
  }

  out <- gg_tree_panel(tp, mp, widths = widths, guides = guides)
  if (!is.null(title) || !is.null(subtitle)) {
    out <- out + patchwork::plot_annotation(
      title = title, subtitle = subtitle,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        plot.title.position = "plot"))
  }
  attr(out, "tree") <- tp
  attr(out, "map") <- mp
  out
}
