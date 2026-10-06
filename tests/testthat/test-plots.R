r    <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
nwk  <- system.file("extdata", "example_tree.nwk", package = "memeplotr")
built <- function(p) ggplot2::ggplot_build(p)

test_that("gg_motif_map() builds for every shape", {
  for (s in motif_shapes()) {
    p <- gg_motif_map(r, shape = s)
    expect_s3_class(p, "ggplot")
    expect_no_error(built(p))
  }
})

test_that("gg_motif_map() draws one block per site", {
  expect_equal(n_blocks(gg_motif_map(r)), nrow(r$sites))
  expect_equal(n_blocks(gg_motif_map(r, source = "contributing")),
               nrow(r$contributing))
})

test_that("filters subset the data", {
  p <- gg_motif_map(r, motifs = c("MEME-1", "MEME-2"), backbone = FALSE)
  expect_setequal(fill_labels(p), c("MEME-1", "MEME-2"))
  expect_equal(n_blocks(p), sum(r$sites$motif %in% c("MEME-1", "MEME-2")))

  seqs <- utils::head(unique(r$sequences$sequence), 3)
  p2 <- gg_motif_map(r, sequences = seqs)
  expect_setequal(y_labels(p2), seqs)
  expect_setequal(attr(p2, "seq_levels"), seqs)

  cut <- stats::median(r$sites$pvalue, na.rm = TRUE)
  p3 <- gg_motif_map(r, pvalue_max = cut)
  expect_lt(n_blocks(p3), nrow(r$sites))
  expect_equal(n_blocks(p3), sum(r$sites$pvalue <= cut, na.rm = TRUE))
})

test_that("sequence ordering options are honoured", {
  ## "seq_levels" is the order the tracks are drawn in, first element on top
  u <- unique(r$sequences$sequence)
  expect_equal(attr(gg_motif_map(r, order = "input"), "seq_levels"), u)
  expect_equal(attr(gg_motif_map(r, order = "name"), "seq_levels"), sort(u))
  expect_equal(attr(gg_motif_map(r, seq_order = rev(u)), "seq_levels"), rev(u))
  expect_equal(attr(gg_motif_map(r, order = "tree", tree = nwk), "seq_levels"),
               tree_tip_order(nwk))

  len <- stats::setNames(r$sequences$length, r$sequences$sequence)
  expect_true(all(diff(len[attr(gg_motif_map(r, order = "length"), "seq_levels")]) <= 0))

  ## the first track must sit at the largest y, so that row 1 is at the top
  p <- gg_motif_map(r, order = "name")
  pp <- ggplot2::ggplot_build(p)$layout$panel_params[[1]]
  expect_equal(pp$y$get_breaks()[1], length(u))
  expect_equal(y_labels(p)[1], sort(u)[1])
})

test_that("colours can be pinned per motif", {
  f <- block_fills(gg_motif_map(r, colours = c("MEME-1" = "#111111")))
  expect_true("#111111" %in% f)
  expect_equal(length(f), nrow(r$motifs))
})

test_that("unknown palettes and motifs error informatively", {
  expect_error(gg_motif_map(r, palette = "not_a_palette"), "palette")
  expect_error(gg_motif_logo(r, "nope"), "Unknown motif|not found")
})

test_that("gg_phylo() and tree composites align", {
  tp <- gg_phylo(nwk)
  expect_s3_class(tp, "ggplot")
  expect_no_error(built(tp))
  expect_equal(tree_tip_order(nwk), attr(tp, "tip_order"))

  cp <- gg_tree_motifs(r, nwk)
  expect_s3_class(cp, "patchwork")
  expect_equal(attr(attr(cp, "tree"), "tip_order"),
               attr(attr(cp, "map"), "seq_levels"))
  expect_no_error(built(gg_phylo(nwk, cladogram = TRUE, node_labels = TRUE,
                                tip_points = TRUE, scale_bar = TRUE)))
})

test_that("keep_empty retains tracks for tips without hits", {
  r2 <- r
  drop <- r$sequences$sequence[1]
  r2$sites <- r$sites[r$sites$sequence != drop, ]
  lv <- attr(gg_motif_map(r2, order = "tree", tree = nwk, keep_empty = TRUE),
             "seq_levels")
  expect_length(lv, length(read_tree(nwk)$tip.label))
  expect_true(drop %in% lv)
})

test_that("summary plots build", {
  ps <- list(gg_motif_heatmap(r), gg_motif_heatmap(r, value = "count", label = TRUE),
             gg_motif_heatmap(r, motif_colours = TRUE),
             gg_motif_counts(r), gg_motif_counts(r, position = "dodge"),
             gg_motif_cooccurrence(r), gg_motif_cooccurrence(r, measure = "jaccard",
                                                             label = TRUE, upper = TRUE),
             gg_motif_positions(r), gg_motif_positions(r, geom = "histogram"),
             gg_motif_positions(r, geom = "jitter", scale = "relative"))
  for (p in ps) {
    expect_s3_class(p, "ggplot")
    expect_no_error(built(p))
  }
})

test_that("plots remain modifiable with the usual grammar", {
  p <- gg_motif_map(r) +
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    ggplot2::labs(title = "mine") +
    ggplot2::theme_minimal()
  expect_s3_class(p, "ggplot")
  expect_no_error(built(p))
  expect_equal(p$labels$title, "mine")
})

test_that("geom_motif() works as a standalone layer", {
  d <- data.frame(sequence = c("a", "a"), motif = c("M1", "M2"),
                  start = c(1, 20), stop = c(10, 30))
  p <- ggplot2::ggplot(d, ggplot2::aes(xmin = start, xmax = stop, y = sequence,
                                       fill = motif)) +
    geom_motif()
  expect_no_error(built(p))
})

test_that("axis tip labels sit outside the panel and keep their attributes", {
  xrange <- function(p) ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x.range
  geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], character(1))
  tips <- tree_tip_order(nwk)

  ## on the axis: names become y-scale labels, grid reserves their width, and
  ## nothing is drawn inside the panel that could run into a neighbour
  ta <- gg_phylo(nwk, tip_label_position = "axis")
  expect_setequal(y_labels(ta), tips)
  expect_false("GeomText" %in% geoms(ta))

  ## in the panel: a text layer, no tip names on the axis, and the x range is
  ## padded by xlim_mult to make room for them
  tp <- gg_phylo(nwk, tip_label_position = "panel")
  expect_true("GeomText" %in% geoms(tp))
  expect_false(any(tips %in% y_labels(tp)))
  expect_lt(diff(xrange(ta)), diff(xrange(tp)))

  ## gg_tree_motifs() defaults to the axis and the gap margin must not drop
  ## the tip_order attribute the panels are aligned on
  cp <- gg_tree_motifs(r, nwk)
  expect_setequal(y_labels(attr(cp, "tree")), tips)
  expect_equal(attr(attr(cp, "tree"), "tip_order"), tips)
  expect_equal(attr(attr(cp, "tree"), "tip_order"),
               attr(attr(cp, "map"), "seq_levels"))
  expect_no_error(built(gg_tree_motifs(r, nwk, show_labels = "tree", gap = 0)))
})

test_that("sequences missing from the tree cannot knock rows out of register", {
  ## rename one sequence so it no longer matches any tip, as a misspelled or
  ## differently-prefixed label does in real data
  r2 <- r
  bad <- as.character(r2$sites$sequence[1])
  r2$sites$sequence[r2$sites$sequence == bad] <- "Not_In_Tree_xyz"
  r2$sequences$sequence[r2$sequences$sequence == bad] <- "Not_In_Tree_xyz"
  tips <- tree_tip_order(nwk)

  ## standalone map keeps it, appended at the bottom
  expect_warning(m <- gg_motif_map(r2, order = "tree", tree = nwk), "absent from the tree")
  expect_true("Not_In_Tree_xyz" %in% attr(m, "seq_levels"))

  ## the composite drops it, so both panels span exactly the tree's tips
  expect_warning(cp <- gg_tree_motifs(r2, nwk), "absent from the tree")
  expect_equal(attr(attr(cp, "map"), "seq_levels"), tips)
  expect_equal(attr(attr(cp, "tree"), "tip_order"), tips)
  yr <- function(p) ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y.range
  expect_equal(yr(attr(cp, "tree")), yr(attr(cp, "map")))
  expect_no_error(built(cp))
})
