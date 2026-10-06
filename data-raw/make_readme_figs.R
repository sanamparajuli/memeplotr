## Figures for README.md and the delivered demo set.
library(memeplotr); library(ggplot2); library(patchwork)
xml <- system.file("extdata", "example_meme.xml", package = "memeplotr")
nwk <- system.file("extdata", "example_tree.nwk", package = "memeplotr")
res <- relabel_motifs(read_meme(xml),
         c("R2R3-a", "R2R3-b", "SANT", "LHEQLE", "C-term", "Gly-rich"))
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path("man/figures", f)

ggsave(out("README-map.png"),
  gg_motif_map(res, label = TRUE, label_size = 2.1, length_label = TRUE,
               title = "Motif composition of a MYB-like family",
               subtitle = "arrow blocks oriented by match strand; MEME XML read directly"),
  width = 9.2, height = 4.4, dpi = 150, bg = "white")

ggsave(out("README-tree.png"),
  gg_tree_motifs(res, nwk, align_tips = TRUE, widths = c(1, 2.6), shape = "hexagon",
                 palette = "bright",
                 title = "The same diagram aligned to a Newick phylogeny",
                 subtitle = "rows follow the tree's tip order; tips without hits keep their track"),
  width = 10, height = 4.6, dpi = 150, bg = "white")

ggsave(out("README-logos.png"),
  gg_motif_logos(res, motifs = c("R2R3-a", "R2R3-b", "SANT"), ncol = 1,
                 show_consensus = TRUE),
  width = 8.5, height = 5.6, dpi = 150, bg = "white")

ggsave(out("README-summary.png"),
  (gg_motif_heatmap(res, order = "tree", tree = nwk, keep_empty = TRUE,
                    value = "count", label = TRUE, title = "occurrences per sequence") |
   gg_motif_cooccurrence(res, measure = "jaccard", label = TRUE,
                         title = "motif co-occurrence")) +
    plot_annotation(tag_levels = "A"),
  width = 11, height = 4.6, dpi = 150, bg = "white")

ggsave(out("README-custom.png"),
  (gg_motif_map(res, sequences = c("Ath_MYB1", "Osa_MYB1", "Zma_MYB1", "Sly_MYB6"),
                colours = c("LHEQLE" = "#D55E00"), palette = "grey",
                title = "one motif pinned, rest in grey") |
   (gg_motif_map(res, sequences = c("Ath_MYB1", "Osa_MYB1", "Zma_MYB1", "Sly_MYB6"),
                 shape = "capsule", backbone_colour = "grey80") +
      scale_fill_viridis_d(option = "plasma", end = 0.9) +
      labs(title = "your own scale and theme") +
      theme_minimal(base_size = 10) +
      theme(panel.grid.major.y = element_blank(), legend.position = "none"))),
  width = 11, height = 3.2, dpi = 150, bg = "white")

ggsave(out("README-shapes.png"),
  wrap_plots(lapply(motif_shapes(), function(s)
    gg_motif_map(res, sequences = c("Ath_MYB1", "Zma_MYB1", "Sly_MYB6"),
                 shape = s, title = s) +
      theme(legend.position = "none", plot.title = element_text(size = 9),
            axis.title.x = element_blank())), ncol = 3),
  width = 11, height = 5, dpi = 150, bg = "white")

cat("wrote:", paste(list.files("man/figures"), collapse = ", "), "\n")
