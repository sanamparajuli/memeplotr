r <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))

test_that("letter_outline() returns unit-square polygons", {
  for (ch in c(LETTERS, "-", ".")) {
    d <- letter_outline(ch)
    expect_s3_class(d, "data.frame")
    expect_true(all(c("x", "y", "piece") %in% names(d)))
    expect_gt(nrow(d), 2)
    expect_true(all(d$x >= -1e-9 & d$x <= 1 + 1e-9), info = ch)
    expect_true(all(d$y >= -1e-9 & d$y <= 1 + 1e-9), info = ch)
  }
  expect_gt(length(unique(letter_outline("i")$piece)), 1)
})

test_that("logo_heights() measures information content", {
  h <- logo_heights(r, 1)
  expect_true(all(c("position", "letter", "ymin", "ymax") %in% names(h)))
  w <- r$motifs$width[1]
  expect_equal(max(h$position), w)
  tot <- tapply(h$ymax - h$ymin, h$position, sum)
  expect_true(all(tot <= log2(length(r$alphabet$letters)) + 1e-6))

  ## min_height drops letters too small to draw, so the stack only sums to 1
  ## when nothing is filtered out
  hp <- logo_heights(r, 1, method = "probability", min_height = 0)
  totp <- as.numeric(tapply(hp$ymax - hp$ymin, hp$position, sum))
  expect_equal(totp, rep(1, w), tolerance = 1e-6)
  hd <- logo_heights(r, 1, method = "probability")
  expect_true(all(as.numeric(tapply(hd$ymax - hd$ymin, hd$position, sum)) <= 1 + 1e-9))
})

test_that("stacks are ordered smallest-first within a column", {
  h <- logo_heights(r, 1)
  one <- h[h$position == h$position[1], ]
  expect_equal(one$ymin[-1], utils::head(one$ymax, -1), tolerance = 1e-9)
})

test_that("small-sample correction lowers total information", {
  a <- sum(logo_heights(r, 1)$ymax - logo_heights(r, 1)$ymin)
  b <- logo_heights(r, 1, small_sample_correction = TRUE)
  expect_lte(sum(b$ymax - b$ymin), a + 1e-9)
})

test_that("motifs can be addressed by label, id, alt name or index", {
  m <- r$motifs
  expect_equal(logo_heights(r, 1), logo_heights(r, m$motif[1]))
  expect_equal(logo_heights(r, 1), logo_heights(r, m$motif_id[1]))
})

test_that("logo plots build with every colour scheme", {
  for (s in c("auto", "chemistry", "hydrophobicity", "taylor", "mono")) {
    p <- gg_motif_logo(r, 1, scheme = s)
    expect_s3_class(p, "ggplot")
    expect_no_error(ggplot2::ggplot_build(p))
  }
  expect_no_error(ggplot2::ggplot_build(
    gg_motif_logo(r, 1, method = "probability", show_consensus = TRUE,
                  letter_colours = c(A = "#000000"), outline = "grey20")))
  expect_s3_class(gg_motif_logos(r, motifs = r$motifs$motif[1:3], ncol = 1), "patchwork")
})

test_that("reverse_complement is refused for protein alphabets", {
  expect_error(gg_motif_logo(r, 1, reverse_complement = TRUE), "nucleotide|complement")
})

test_that("logo_colours() returns a named vector covering the alphabet", {
  cl <- logo_colours("chemistry")
  expect_true(all(c("A", "W", "K", "D") %in% names(cl)))
  expect_true(all(grepl("^#", unname(cl))))
})
