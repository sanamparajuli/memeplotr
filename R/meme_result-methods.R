#' Methods for `meme_result` objects
#'
#' @param x A `meme_result` from [read_meme()].
#' @param object A `meme_result` from [read_meme()].
#' @param row.names,optional Ignored; present for S3 consistency.
#' @param ... Passed on to methods.
#'
#' @return `print()` returns `x` invisibly; `summary()` returns a tibble with one
#'   row per motif plus occupancy statistics; `as.data.frame()` returns the site
#'   table as a plain data frame.
#'
#' @examples
#' res <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
#' print(res)
#' summary(res)
#' head(as.data.frame(res))
#'
#' @name meme_result
NULL

#' @rdname meme_result
#' @export
print.meme_result <- function(x, ...) {
  m <- x$meta
  cat(sprintf("<meme_result>  %s %s   alphabet: %s\n", m$tool, m$version %||% "?", m$alphabet))
  cat(sprintf("  motifs     : %d  (widths %s)\n", nrow(x$motifs),
              if (nrow(x$motifs)) paste(range(x$motifs$width, na.rm = TRUE), collapse = "-") else "-"))
  cat(sprintf("  sequences  : %d\n", nrow(x$sequences)))
  cat(sprintf("  sites      : %d (%s)\n", nrow(x$sites),
              paste(unique(x$sites$source), collapse = "/")))
  if (nrow(x$motifs)) {
    tb <- x$motifs[, c("motif", "consensus", "width", "nsites", "evalue")]
    tb$consensus <- ifelse(nchar(tb$consensus %||% "") > 24,
                           paste0(substr(tb$consensus, 1, 21), "..."), tb$consensus)
    cat("\n")
    print(as.data.frame(utils::head(tb, 10)), row.names = FALSE)
    if (nrow(tb) > 10) cat(sprintf("  ... %d more motifs\n", nrow(tb) - 10))
  }
  invisible(x)
}

#' @rdname meme_result
#' @export
summary.meme_result <- function(object, ...) {
  s <- object$sites
  nseq <- max(nrow(object$sequences), 1L)
  occ <- if (nrow(s)) {
    dplyr::summarise(
      dplyr::group_by(s, .data$motif_id),
      n_sites = dplyr::n(),
      n_seqs = dplyr::n_distinct(.data$sequence),
      min_pvalue = suppressWarnings(min(.data$pvalue, na.rm = TRUE)),
      median_start = stats::median(.data$start, na.rm = TRUE),
      .groups = "drop"
    )
  } else {
    tibble::tibble(motif_id = character(0), n_sites = integer(0), n_seqs = integer(0),
                   min_pvalue = numeric(0), median_start = numeric(0))
  }
  out <- dplyr::left_join(object$motifs, occ, by = "motif_id")
  out$n_sites <- ifelse(is.na(out$n_sites), 0L, out$n_sites)
  out$n_seqs <- ifelse(is.na(out$n_seqs), 0L, out$n_seqs)
  out$pct_seqs <- round(100 * out$n_seqs / nseq, 1)
  out[, c("motif", "consensus", "width", "evalue", "n_sites", "n_seqs",
          "pct_seqs", "min_pvalue", "median_start")]
}

#' @rdname meme_result
#' @export
as.data.frame.meme_result <- function(x, row.names = NULL, optional = FALSE, ...) {
  as.data.frame(x$sites, row.names = row.names, optional = optional, ...)
}

#' Relabel motifs
#'
#' Replaces the `motif` display label used by every plotting function.
#'
#' @param x A [meme_result].
#' @param by One of `"alt"` (e.g. `MEME-1`), `"name"`, `"consensus"`, `"regex"`,
#'   `"id"`, `"index"` (e.g. `Motif 1`), or a character vector of new labels the
#'   same length as `x$motifs`, optionally named by `motif_id`.
#' @param prefix Prefix used when `by = "index"`. Default `"Motif "`.
#'
#' @return The `meme_result` with updated `motif` labels in `$motifs`, `$sites`
#'   and `$contributing`.
#'
#' @examples
#' res <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
#' relabel_motifs(res, "consensus")$motifs$motif
#' relabel_motifs(res, "index")$motifs$motif
#' @export
relabel_motifs <- function(x, by = c("alt", "name", "consensus", "regex", "id", "index"),
                           prefix = "Motif ") {
  stopifnot(inherits(x, "meme_result"))
  keywords <- c("alt", "name", "consensus", "regex", "id", "index")
  if (is.character(by) && length(by) > 1L && !all(by %in% keywords) &&
      length(by) != nrow(x$motifs)) {
    .abort(c(sprintf("`by` must be one of %s, or one label per motif.",
                     paste0('"', keywords, '"', collapse = ", ")),
             i = sprintf("Got %d labels for %d motifs.", length(by), nrow(x$motifs))))
  }
  if (is.character(by) && length(by) == nrow(x$motifs) && !identical(length(by), 1L)) {
    new <- if (!is.null(names(by))) unname(by[x$motifs$motif_id]) else by
  } else {
    by <- match.arg(by)
    new <- switch(by,
      alt = ifelse(is.na(x$motifs$alt), x$motifs$motif, x$motifs$alt),
      name = ifelse(is.na(x$motifs$name), x$motifs$motif, x$motifs$name),
      consensus = x$motifs$consensus,
      regex = x$motifs$regex,
      id = x$motifs$motif_id,
      index = paste0(prefix, x$motifs$index)
    )
  }
  new <- ifelse(is.na(new) | !nzchar(new), x$motifs$motif_id, new)
  if (anyDuplicated(new)) new <- make.unique(new, sep = "_")
  map <- stats::setNames(new, x$motifs$motif_id)
  x$motifs$motif <- unname(map[x$motifs$motif_id])
  for (nm in c("sites", "contributing")) {
    if (nrow(x[[nm]])) x[[nm]]$motif <- unname(map[x[[nm]]$motif_id])
  }
  x
}
