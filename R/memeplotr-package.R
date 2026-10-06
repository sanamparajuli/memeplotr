#' memeplotr: flexible ggplot2 visualisation of MEME Suite output
#'
#' @description
#' `memeplotr` turns MEME Suite XML into tidy tables and gives you ordinary
#' `ggplot2` objects back, so every colour, shape, scale, label and theme is
#' modifiable with the usual `+` grammar. The main entry points are:
#'
#' * [read_meme()] - parse `meme.xml` / `streme.xml` into a [meme_result] object
#' * [gg_motif_map()] - motif composition (block) diagram
#' * [gg_motif_logo()] / [gg_motif_logos()] - sequence logos from the PWMs
#' * [gg_phylo()] / [gg_tree_motifs()] - Newick tree, and tree-aligned motif map
#' * [gg_motif_heatmap()], [gg_motif_counts()], [gg_motif_cooccurrence()]
#' * [run_memeplotr()] - the interactive Shiny front end
#'
#' @keywords internal
"_PACKAGE"

#' @import ggplot2
#' @importFrom rlang .data %||%
#' @importFrom stats setNames reorder median
#' @importFrom utils head modifyList
NULL

utils::globalVariables(c("."))
