test_that("the Shiny app builds every panel and emits runnable code", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  app_dir <- system.file("shiny", package = "memeplotr")
  skip_if(!nzchar(app_dir))

  shiny::testServer(app_dir, {
    session$setInputs(
      src = "example", tsrc = "example", order = "input", source = "sites",
      motifs = character(), seqs = character(), pmax = 0, shape = "arrow",
      palette = "okabe_ito", custom_cols = FALSE, height = 0.6, alpha = 1,
      strand_shape = FALSE, label = FALSE, length_label = FALSE, backbone = TRUE,
      cladogram = FALSE, align_tips = TRUE, tree_w = 1, title = "", subtitle = "",
      w = 11, h = 5.5, dpi = 300, tab = "Motif map", tbl = "sites",
      logo_motifs = character(), logo_method = "bits", logo_scheme = "auto",
      logo_cons = TRUE, logo_rc = FALSE, logo_ssc = FALSE, logo_ncol = 1,
      logo_outline = 0, sum_kind = "heatmap", sum_label = FALSE)

    expect_s3_class(res(), "meme_result")
    expect_s3_class(the_map(), "ggplot")
    expect_s3_class(the_tree_map(), "patchwork")
    expect_s3_class(the_summary(), "ggplot")
    expect_gt(nrow(tbl()), 0)

    session$setInputs(order = "tree", shape = "hexagon", label = TRUE)
    expect_s3_class(the_map(), "ggplot")
    expect_s3_class(the_tree_map(), "patchwork")

    session$setInputs(tab = "Summary", sum_kind = "counts", sum_position = "dodge")
    expect_s3_class(the_summary(), "ggplot")
    session$setInputs(sum_kind = "co-occurrence", sum_measure = "jaccard")
    expect_s3_class(the_summary(), "ggplot")
    session$setInputs(sum_kind = "positions", sum_geom = "histogram")
    expect_s3_class(the_summary(), "ggplot")

    ## the code tab has to parse, and the call it names has to exist
    session$setInputs(tab = "Motif map")
    code <- output$code
    expect_no_error(parse(text = code))
    expect_match(code, "read_meme\\(")
    expect_match(code, "gg_motif_map\\(")
  })
})

test_that("an uploaded XML is never paired with the bundled example tree", {
  skip_if_not_installed("shiny")
  app <- system.file("shiny", "app.R", package = "memeplotr")
  skip_if(!nzchar(app))
  up <- file.path(tempdir(), "uploaded.xml")
  file.copy(system.file("extdata", "example_meme_dna.xml", package = "memeplotr"),
            up, overwrite = TRUE)

  shiny::testServer(shiny::shinyAppFile(app), {
    ## example data + example tree: the tree must load
    session$setInputs(src = "example", tsrc = "example", exset = "protein")
    expect_false(is.null(phy()))

    ## uploaded data + example tree: the tree must be dropped, not matched to
    ## zero tips (the bundled tree describes the bundled sequences only)
    session$setInputs(src = "upload", tsrc = "example")
    session$setInputs(xml = list(name = "uploaded.xml", size = file.size(up),
                                 type = "text/xml", datapath = up))
    expect_s3_class(res(), "meme_result")
    expect_null(tree_path())
    expect_null(phy())
  })
})

test_that("both bundled example sets load in the app", {
  skip_if_not_installed("shiny")
  app <- shiny::shinyAppFile(system.file("shiny", "app.R", package = "memeplotr"))
  for (set in c("protein", "dna")) {
    shiny::testServer(app, {
      session$setInputs(src = "example", exset = set, tsrc = "example")
      r <- res()
      expect_s3_class(r, "meme_result")
      expect_gt(nrow(r$motifs), 0L)
      expect_false(is.null(tree_path()))
    })
  }
})

test_that("every file named in the app's example table is installed", {
  for (f in c("example_meme.xml", "example_tree.nwk",
              "example_meme_dna.xml", "example_tree_dna.nwk")) {
    expect_true(nzchar(system.file("extdata", f, package = "memeplotr")), info = f)
  }
})

test_that("ambiguity codes in the alphabet do not break the PWM column names", {
  r <- read_meme(system.file("extdata", "example_meme_dna.xml", package = "memeplotr"))
  expect_gt(length(r$pwm), 0L)
  expect_true(all(vapply(r$pwm, function(m) ncol(m) == 4L, logical(1))))
  expect_setequal(colnames(r$pwm[[1]]), c("A", "C", "G", "T"))

  p <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
  expect_gt(length(p$pwm), 0L)
  expect_true(all(vapply(p$pwm, function(m) ncol(m) == 20L, logical(1))))
  expect_setequal(colnames(p$pwm[[1]]), strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]])
})

test_that("bundle_shiny_app() writes a loadable, self-contained app directory", {
  src <- test_path("..", "..")           # the package source under test
  skip_if_not(file.exists(file.path(src, "inst", "shiny", "app.R")),
              "not running from a source tree")
  skip_if_not_installed("pkgload")
  out <- file.path(tempdir(), "memeplotr-bundle-test")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  expect_message(bundle_shiny_app(out, pkg_source = src), "Shiny bundle")
  expect_true(all(file.exists(file.path(out,
    c("app.R", "app_body.R", "memeplotr/DESCRIPTION", "memeplotr/NAMESPACE",
      "memeplotr/R", "memeplotr/inst/extdata/example_meme.xml")))))

  ## the whole point: no library(memeplotr), since nothing is installed there
  body <- readLines(file.path(out, "app_body.R"), warn = FALSE)
  expect_false(any(grepl("^\\s*library\\(memeplotr\\)", body)))
  ## and tests/vignettes stay out, so the server does not install testthat
  expect_false(dir.exists(file.path(out, "memeplotr", "tests")))
  expect_false(dir.exists(file.path(out, "memeplotr", "vignettes")))

  expect_error(bundle_shiny_app(out, pkg_source = src, overwrite = FALSE),
               "already exists")
  expect_error(bundle_shiny_app(out, pkg_source = tempfile()), "does not exist")
})
