ex  <- function() system.file("extdata", "example_meme.xml", package = "memeplotr")
ex4 <- function() system.file("extdata", "example_meme_v4.xml", package = "memeplotr")

test_that("read_meme() returns a complete meme_result", {
  r <- read_meme(ex())
  expect_s3_class(r, "meme_result")
  expect_true(all(c("motifs", "sites", "contributing", "sequences", "pwm",
                    "pssm", "background", "alphabet", "meta") %in% names(r)))
  expect_gt(nrow(r$motifs), 0)
  expect_gt(nrow(r$sites), 0)
  expect_true(all(c("motif", "motif_id", "width", "consensus") %in% names(r$motifs)))
  expect_true(all(c("sequence", "motif", "start", "stop", "strand") %in% names(r$sites)))
})

test_that("sites are consistent with the motif and sequence tables", {
  r <- read_meme(ex())
  expect_true(all(r$sites$sequence %in% r$sequences$sequence))
  expect_true(all(r$sites$motif %in% r$motifs$motif))
  expect_true(all(r$sites$start >= 1))
  len <- stats::setNames(r$sequences$length, r$sequences$sequence)
  expect_true(all(r$sites$stop <= len[r$sites$sequence]))
  expect_true(all(r$sites$stop >= r$sites$start))
  expect_true(all(r$sites$strand %in% c("+", "-", "*")))
})

test_that("PWMs are proper probability matrices", {
  r <- read_meme(ex())
  expect_length(r$pwm, nrow(r$motifs))
  for (m in r$pwm) {
    expect_true(is.matrix(m))
    expect_equal(unname(rowSums(m)), rep(1, nrow(m)), tolerance = 1e-6)
    expect_true(all(m >= 0))
    expect_identical(colnames(m), r$alphabet$letters)
  }
  expect_equal(vapply(r$pwm, nrow, integer(1)), r$motifs$width,
               ignore_attr = TRUE)
})

test_that("the MEME 4.x layout parses to the same content", {
  a <- read_meme(ex())
  b <- read_meme(ex4())
  expect_equal(nrow(a$motifs), nrow(b$motifs))
  expect_equal(nrow(a$sites), nrow(b$sites))
  expect_equal(a$motifs$consensus, b$motifs$consensus)
  expect_equal(a$alphabet$letters, b$alphabet$letters)
})

test_that("contributing sites are a separate table carrying the site sequence", {
  r <- read_meme(ex())
  expect_gt(nrow(r$contributing), 0)
  expect_true(all(c("site_seq", "left_flank", "right_flank") %in% names(r$contributing)))
  expect_true(all(r$sites$source == "scanned"))
  expect_true(all(r$contributing$source == "contributing"))
  expect_true(all(nchar(r$contributing$site_seq) ==
                    r$motifs$width[match(r$contributing$motif, r$motifs$motif)]))
})

test_that("print and summary methods work", {
  r <- read_meme(ex())
  expect_output(print(r), "meme_result")
  s <- summary(r)
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), nrow(r$motifs))
  expect_s3_class(as.data.frame(r), "data.frame")
})

test_that("relabel_motifs() renames consistently across tables", {
  r <- read_meme(ex())
  r2 <- relabel_motifs(r, "index")
  expect_equal(r2$motifs$motif, paste("Motif", seq_len(nrow(r$motifs))))
  expect_true(all(r2$sites$motif %in% r2$motifs$motif))
  expect_true(all(r2$contributing$motif %in% r2$motifs$motif))
  lab <- paste0("box", seq_len(nrow(r$motifs)))
  r3 <- relabel_motifs(r, lab)
  expect_equal(r3$motifs$motif, lab)
  expect_setequal(unique(r3$sites$motif), lab)
  expect_error(relabel_motifs(r, c("alpha", "beta")), "one label per motif")
})

test_that("bad input is rejected with a useful message", {
  expect_error(read_meme(tempfile()), "exist")
  f <- tempfile(fileext = ".xml")
  writeLines("<not_meme><a/></not_meme>", f)
  expect_error(read_meme(f), "MEME")
})
