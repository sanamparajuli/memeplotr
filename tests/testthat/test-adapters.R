test_that("read_fimo() produces a motif_table", {
  f <- system.file("extdata", "example_fimo.tsv", package = "memeplotr")
  skip_if(!nzchar(f))
  d <- read_fimo(f)
  expect_s3_class(d, "motif_table")
  expect_true(all(c("sequence", "motif", "start", "stop", "strand") %in% names(d)))
  expect_true(all(d$stop >= d$start))
})

test_that("as_motif_table() adapts an arbitrary data frame", {
  df <- data.frame(seq_id = c("a", "a", "b"), mot = c("M1", "M2", "M1"),
                   from = c(1, 20, 5), to = c(10, 29, 14),
                   sense = c("+", "-", "+"), p = c(1e-5, 1e-3, 1e-4))
  d <- as_motif_table(df, sequence = "seq_id", motif = "mot", start = "from",
                      stop = "to", strand = "sense", pvalue = "p")
  expect_s3_class(d, "motif_table")
  expect_equal(nrow(d), 3L)
  expect_equal(d$width, c(10, 10, 10))
  expect_equal(levels(factor(d$motif)), c("M1", "M2"))
})

test_that("a motif_table can be plotted like a meme_result", {
  df <- data.frame(sequence = rep(c("a", "b"), each = 2),
                   motif = c("M1", "M2", "M1", "M2"),
                   start = c(1, 20, 3, 25), stop = c(10, 30, 12, 35),
                   strand = "+", seq_length = 60)
  d <- as_motif_table(df)
  p <- gg_motif_map(d)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("missing required columns are reported", {
  expect_error(as_motif_table(data.frame(a = 1)), "sequence|motif|start")
})
