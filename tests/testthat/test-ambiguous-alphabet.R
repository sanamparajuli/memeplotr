test_that("ambiguous letters are excluded from the alphabet and PWM columns", {
  dna <- read_meme(system.file("extdata", "example_meme_dna.xml", package = "memeplotr"))

  expect_identical(dna$meta$alphabet, "dna")
  expect_identical(dna$alphabet$letters, c("A", "C", "G", "T"))
  for (m in dna$pwm) expect_identical(colnames(m), c("A", "C", "G", "T"))
  for (m in dna$pssm) expect_identical(colnames(m), c("A", "C", "G", "T"))
  expect_true(all(abs(rowSums(dna$pwm[[1]]) - 1) < 1e-6))

  prot <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
  expect_length(prot$alphabet$letters, 20L)
  expect_false(any(c("B", "Z", "J", "X") %in% prot$alphabet$letters))
  expect_identical(ncol(prot$pwm[[1]]), 20L)
})

test_that("ambiguous letter_ids still resolve in contributing sites", {
  dna <- read_meme(system.file("extdata", "example_meme_dna.xml", package = "memeplotr"))
  expect_true(nrow(dna$contributing) > 0)
  expect_true(all(grepl("^[ACGT]+$", dna$contributing$site_seq)))
  # the site string must reproduce the motif consensus
  expect_identical(
    sort(unique(dna$contributing$site_seq)),
    sort(unique(dna$motifs$name))
  )
})

test_that("a ragged alphabet_matrix is reported, not silently reshaped", {
  f <- tempfile(fileext = ".xml")
  x <- readLines(system.file("extdata", "example_meme_dna.xml", package = "memeplotr"))
  i <- grep("letter_id=\"A\">0.03", x)[1]
  x[i] <- sub("<value letter_id=\"T\">[0-9.]+</value>", "", x[i])
  writeLines(x, f)
  expect_error(read_meme(f), "ragged")
})

test_that("every plot entry point works on DNA input", {
  dna <- read_meme(system.file("extdata", "example_meme_dna.xml", package = "memeplotr"))
  nwk <- system.file("extdata", "example_tree_dna.nwk", package = "memeplotr")
  expect_s3_class(gg_motif_map(dna), "ggplot")
  expect_s3_class(gg_motif_logo(dna, 1), "ggplot")
  expect_s3_class(gg_motif_heatmap(dna), "ggplot")
  expect_s3_class(gg_motif_counts(dna), "ggplot")
  expect_s3_class(gg_motif_cooccurrence(dna), "ggplot")
  expect_s3_class(gg_motif_positions(dna), "ggplot")
  expect_no_error(ggplot2::ggplot_build(gg_tree_motifs(dna, nwk)))
})
