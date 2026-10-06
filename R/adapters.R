#' Coerce motif hits to the standard memeplotr site table
#'
#' Every plotting function in `memeplotr` consumes the same tidy "motif table".
#' `as_motif_table()` builds one from a [meme_result], a FIMO/MAST result, or
#' any data frame of motif hits - which is how output from non-MEME tools
#' (HMMER, InterProScan, custom scanners) enters the package.
#'
#' @param x A [meme_result], a data frame, or a tibble of hits.
#' @param sequence,motif,start,stop,strand,pvalue Column names in `x` (given as
#'   strings) holding each field. Defaults match the memeplotr convention.
#' @param width Optional column name holding motif width; if absent it is
#'   derived as `stop - start + 1`.
#' @param seq_length Optional: either a column name in `x`, or a named numeric
#'   vector of sequence lengths, used to draw the sequence backbone. If `NULL`
#'   the backbone is drawn to the furthest motif end of each sequence.
#' @param source Which site table to use for a [meme_result]: `"sites"` (all
#'   scanned hits, the default) or `"contributing"` (only the sites MEME used to
#'   build each motif).
#' @param ... Unused.
#'
#' @return A tibble of class `motif_table` with columns `sequence`, `motif`,
#'   `start`, `stop`, `width`, `strand`, `pvalue`, `seq_length`.
#'
#' @examples
#' res <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
#' head(as_motif_table(res))
#'
#' # any scanner's output, adapted
#' hits <- data.frame(seq = c("a", "a", "b"), mot = c("M1", "M2", "M1"),
#'                    from = c(3, 40, 11), to = c(12, 49, 20))
#' as_motif_table(hits, sequence = "seq", motif = "mot",
#'                start = "from", stop = "to")
#' @export
as_motif_table <- function(x, ...) UseMethod("as_motif_table")

#' @rdname as_motif_table
#' @export
as_motif_table.meme_result <- function(x, source = c("sites", "contributing"), ...) {
  source <- match.arg(source)
  s <- x[[source]]
  if (!nrow(s)) {
    rlang::warn(sprintf("`%s` table is empty for this %s result.", source, x$meta$tool))
  }
  len <- stats::setNames(x$sequences$length, x$sequences$sequence)
  out <- tibble::tibble(
    sequence = s$sequence, motif = s$motif, motif_id = s$motif_id,
    start = as.numeric(s$start), stop = as.numeric(s$stop),
    width = as.numeric(s$width), strand = s$strand, pvalue = s$pvalue,
    seq_length = unname(len[s$sequence])
  )
  .finish_motif_table(out, len)
}

#' @rdname as_motif_table
#' @export
as_motif_table.data.frame <- function(x, sequence = "sequence", motif = "motif",
                                      start = "start", stop = "stop",
                                      strand = "strand", pvalue = "pvalue",
                                      width = "width", seq_length = NULL, ...) {
  need <- c(sequence, motif, start)
  miss <- need[!need %in% names(x)]
  if (length(miss)) {
    .abort(paste0("Columns not found in `x`: ", paste(miss, collapse = ", "), "."),
           i = paste0("Available: ", paste(names(x), collapse = ", ")))
  }
  pick <- function(nm, default = NA) if (!is.null(nm) && nm %in% names(x)) x[[nm]] else default
  w <- pick(width)
  st <- as.numeric(x[[start]])
  sp <- pick(stop)
  if (all(is.na(sp))) {
    if (all(is.na(w))) .abort("Supply either a `stop` or a `width` column.")
    sp <- st + as.numeric(w) - 1
  }
  sp <- as.numeric(sp)
  len <- NULL
  if (is.character(seq_length) && length(seq_length) == 1L && seq_length %in% names(x)) {
    len <- stats::setNames(as.numeric(x[[seq_length]]), x[[sequence]])
    len <- len[!duplicated(names(len))]
  } else if (is.numeric(seq_length)) {
    len <- seq_length
  }
  ## values are resolved before tibble() is called: inside tibble() the growing
  ## data mask would shadow the `motif` / `strand` argument names.
  v_seq <- as.character(x[[sequence]])
  v_mot <- as.character(x[[motif]])
  v_str <- .norm_strand(pick(strand, "*"))
  v_pv  <- as.numeric(pick(pvalue, NA_real_))
  v_w   <- if (all(is.na(w))) sp - st + 1 else as.numeric(w)
  v_len <- if (is.null(len)) NA_real_ else unname(len[v_seq])
  out <- tibble::tibble(
    sequence = v_seq, motif = v_mot, motif_id = v_mot,
    start = st, stop = sp, width = v_w,
    strand = v_str, pvalue = v_pv, seq_length = v_len
  )
  .finish_motif_table(out, len)
}

#' @keywords internal
#' @noRd
.finish_motif_table <- function(out, len = NULL) {
  if (nrow(out) && anyNA(out$seq_length)) {
    fallback <- stats::setNames(
      as.numeric(tapply(out$stop, out$sequence, max, na.rm = TRUE)),
      names(tapply(out$stop, out$sequence, max, na.rm = TRUE)))
    out$seq_length[is.na(out$seq_length)] <- unname(fallback[out$sequence[is.na(out$seq_length)]])
  }
  out$motif <- factor(out$motif, levels = unique(out$motif))
  structure(out, class = c("motif_table", class(tibble::tibble())))
}

#' Read FIMO output
#'
#' Reads the TSV produced by `fimo --text` or `fimo.tsv` into the standard
#' memeplotr site table, so that motifs discovered by MEME and then scanned
#' across a larger sequence set can be plotted with the same grammar.
#'
#' @param file Path to `fimo.tsv`.
#' @param q_max Optional maximum q-value (FDR) to retain.
#' @param p_max Optional maximum p-value to retain.
#'
#' @return A `motif_table` (see [as_motif_table()]).
#' @examples
#' f <- system.file("extdata", "example_fimo.tsv", package = "memeplotr")
#' if (nzchar(f)) head(read_fimo(f))
#' @export
read_fimo <- function(file, q_max = NULL, p_max = NULL) {
  x <- utils::read.delim(file, comment.char = "#", stringsAsFactors = FALSE)
  x <- x[nzchar(trimws(x[[1]])), , drop = FALSE]
  nm <- names(x)
  seq_col <- intersect(c("sequence_name", "sequence.name", "sequence"), nm)[1]
  mot_col <- intersect(c("motif_alt_id", "motif_id", "motif.alt.id", "motif"), nm)[1]
  if (is.na(seq_col) || is.na(mot_col)) {
    .abort("Could not find sequence/motif columns in the FIMO file.",
           i = paste0("Columns present: ", paste(nm, collapse = ", ")))
  }
  if (!is.null(p_max) && "p.value" %in% nm) x <- x[x$p.value <= p_max, , drop = FALSE]
  if (!is.null(q_max) && "q.value" %in% nm) x <- x[x$q.value <= q_max, , drop = FALSE]
  as_motif_table(x, sequence = seq_col, motif = mot_col, start = "start",
                 stop = "stop", strand = "strand", pvalue = "p.value")
}
