#' Read MEME Suite XML output
#'
#' Parses the XML produced by `meme`, `streme` or `dreme` into a tidy
#' [meme_result] object. Both the MEME 4.x and 5.x XML layouts are supported;
#' the differences in how letters, sequences and motif names are encoded are
#' resolved internally.
#'
#' @param file Path (or URL, or [xml2::read_xml()]-compatible connection) to a
#'   `meme.xml` / `streme.xml` / `dreme.xml` file.
#' @param sites_sequence Logical; extract the matched sub-sequence (and flanks)
#'   for each contributing site. Set to `FALSE` for a faster parse of very large
#'   files. Default `TRUE`.
#' @param p_thresh Optional numeric. If supplied, scanned sites with a p-value
#'   above this threshold are dropped at read time. MEME has usually already
#'   applied its own threshold (reported in `x$meta$site_p_thresh`).
#'
#' @return An object of class `meme_result`: a list with elements
#'   \describe{
#'     \item{meta}{one-row tibble of run metadata (tool, version, alphabet,
#'       command line, model type, thresholds).}
#'     \item{motifs}{tibble, one row per motif: `motif_id`, `motif`, `name`,
#'       `alt`, `consensus`, `regex`, `width`, `nsites`, `evalue`, `pvalue`,
#'       `ic`, `re`, `llr`, `index`.}
#'     \item{sites}{tibble of motif hits, one row per site:
#'       `sequence`, `sequence_id`, `motif`, `motif_id`, `start`, `stop`,
#'       `width`, `strand`, `pvalue`, `source`.}
#'     \item{contributing}{tibble of the sites MEME used to build each motif,
#'       with `site_seq`, `left_flank`, `right_flank` when available.}
#'     \item{sequences}{tibble of input sequences: `sequence_id`, `sequence`,
#'       `length`, `weight`, `seq_pvalue`, `n_sites`.}
#'     \item{pwm}{named list of position probability matrices
#'       (rows = position, columns = letters).}
#'     \item{pssm}{named list of log-odds score matrices, when present.}
#'     \item{background}{named numeric vector of background letter frequencies.}
#'     \item{alphabet}{list with `type`, `letters` and `colours`.}
#'   }
#'
#' @details
#' MEME reports site positions 0-based; `memeplotr` converts them to 1-based
#' inclusive `start`/`stop` coordinates so that `stop - start + 1 == width`.
#'
#' `streme.xml` and `dreme.xml` contain motifs and PWMs but no per-sequence site
#' table. For those files `sites` is returned empty with a warning; combine them
#' with [read_fimo()] or [as_motif_table()] to draw motif maps.
#'
#' @examples
#' xml <- system.file("extdata", "example_meme.xml", package = "memeplotr")
#' res <- read_meme(xml)
#' res
#' head(res$sites)
#'
#' @seealso [read_fimo()], [as_motif_table()], [gg_motif_map()]
#' @export
read_meme <- function(file, sites_sequence = TRUE, p_thresh = NULL) {
  doc <- xml2::read_xml(file)
  root <- xml2::xml_root(doc)
  tool <- xml2::xml_name(root)
  res <- switch(
    toupper(tool),
    "MEME"   = .read_meme_doc(root, sites_sequence = sites_sequence),
    "STREME" = .read_streme_doc(root),
    "DREME"  = .read_streme_doc(root),
    .abort(sprintf("Unrecognised root element <%s>. Expected <MEME>, <STREME> or <DREME>.", tool))
  )
  if (!is.null(p_thresh) && nrow(res$sites)) {
    res$sites <- res$sites[is.na(res$sites$pvalue) | res$sites$pvalue <= p_thresh, , drop = FALSE]
  }
  res
}

# ---- MEME ------------------------------------------------------------------

#' @keywords internal
#' @noRd
.read_meme_doc <- function(root, sites_sequence = TRUE) {
  ver <- xml2::xml_attr(root, "version")

  ## ---- alphabet (4.x: id="letter_A"; 5.x: id="A") -------------------------
  alph_node <- xml2::xml_find_first(root, ".//training_set/alphabet")
  lnodes <- xml2::xml_find_all(root, ".//training_set/alphabet/letter")
  lid <- xml2::xml_attr(lnodes, "id")
  lsym <- xml2::xml_attr(lnodes, "symbol")
  lsym[is.na(lsym)] <- sub("^letter_", "", lid[is.na(lsym)])
  lmap <- stats::setNames(lsym, lid)
  lcol <- xml2::xml_attr(lnodes, "colour")
  if (all(is.na(lcol))) lcol <- xml2::xml_attr(lnodes, "color")
  lcol <- stats::setNames(ifelse(is.na(lcol), NA_character_, paste0("#", lcol)), lsym)

  # MEME lists ambiguous letters (DNA: N, V, H, ...; protein: B, Z, J, X) in the
  # same <alphabet> block as the core ones, marked by an `equals` attribute, but
  # the PWM/PSSM matrices only carry the core letters. Keep the full id->symbol
  # map for <letter_ref> lookups; everything alphabet-shaped uses core only.
  core <- is.na(xml2::xml_attr(lnodes, "equals"))
  if (!any(core)) core <- rep(TRUE, length(lnodes))
  lsym <- lsym[core]
  lcol <- lcol[core]

  a_like <- xml2::xml_attr(alph_node, "like")
  a_name <- xml2::xml_attr(alph_node, "name") %||% xml2::xml_attr(alph_node, "id")
  atype <- .alphabet_type(a_like, a_name, lsym)

  ## ---- sequences -----------------------------------------------------------
  snodes <- xml2::xml_find_all(root, ".//training_set/sequence")
  sequences <- tibble::tibble(
    sequence_id = xml2::xml_attr(snodes, "id"),
    sequence    = xml2::xml_attr(snodes, "name"),
    length      = .as_int(xml2::xml_attr(snodes, "length")),
    weight      = .as_num(xml2::xml_attr(snodes, "weight"))
  )
  smap <- stats::setNames(sequences$sequence, sequences$sequence_id)

  ## ---- model / background --------------------------------------------------
  gettxt <- function(xp) {
    n <- xml2::xml_find_first(root, xp)
    if (inherits(n, "xml_missing")) NA_character_ else trimws(xml2::xml_text(n))
  }
  bgn <- xml2::xml_find_all(root, ".//model/background_frequencies/alphabet_array/value")
  background <- if (length(bgn)) {
    stats::setNames(.as_num(xml2::xml_text(bgn)), unname(lmap[xml2::xml_attr(bgn, "letter_id")]))
  } else stats::setNames(numeric(0), character(0))

  ss_sum <- xml2::xml_find_first(root, ".//scanned_sites_summary")
  meta <- tibble::tibble(
    tool          = "MEME",
    version       = ver,
    alphabet      = atype,
    nmotifs       = .as_int(gettxt(".//model/nmotifs")),
    model_type    = gettxt(".//model/type"),
    nsequences    = nrow(sequences),
    command_line  = gettxt(".//model/command_line"),
    datafile      = xml2::xml_attr(xml2::xml_find_first(root, ".//training_set"), "primary_sequences") %||%
                    xml2::xml_attr(xml2::xml_find_first(root, ".//training_set"), "datafile"),
    site_p_thresh = if (inherits(ss_sum, "xml_missing")) NA_real_ else .as_num(xml2::xml_attr(ss_sum, "p_thresh"))
  )

  ## ---- motifs --------------------------------------------------------------
  mnodes <- xml2::xml_find_all(root, ".//motifs/motif")
  nL <- length(lsym)

  read_matrix <- function(m, which) {
    arrs <- xml2::xml_find_all(m, sprintf("./%s/alphabet_matrix/alphabet_array", which))
    if (!length(arrs)) return(NULL)
    ids <- xml2::xml_attr(xml2::xml_find_all(arrs[[1]], "./value"), "letter_id")
    cols <- unname(lmap[ids])
    cols[is.na(cols)] <- ids[is.na(cols)]
    if (all(is.na(cols))) cols <- lsym[seq_along(cols)]
    vals <- xml2::xml_find_all(arrs, "./value")
    v <- .as_num(xml2::xml_text(vals))
    if (length(v) != length(arrs) * length(cols)) {
      rlang::abort(sprintf(
        paste0("Motif '%s' has a ragged <%s> matrix: %d values over %d positions ",
               "x %d letters. The XML may be truncated."),
        xml2::xml_attr(m, "id") %||% "?", which, length(v), length(arrs), length(cols)))
    }
    matrix(v, nrow = length(arrs), byrow = TRUE, dimnames = list(NULL, cols))
  }

  motif_id <- xml2::xml_attr(mnodes, "id")
  m_name <- .blank_to_na(xml2::xml_attr(mnodes, "name"))
  m_alt  <- .blank_to_na(xml2::xml_attr(mnodes, "alt"))
  width  <- .as_int(xml2::xml_attr(mnodes, "width"))
  regex  <- vapply(mnodes, function(m) {
    n <- xml2::xml_find_first(m, "./regular_expression")
    if (inherits(n, "xml_missing")) NA_character_ else gsub("\\s+", "", xml2::xml_text(n))
  }, character(1))

  pwm  <- lapply(mnodes, read_matrix, which = "probabilities")
  pssm <- lapply(mnodes, read_matrix, which = "scores")
  names(pwm) <- names(pssm) <- motif_id

  consensus <- vapply(seq_along(mnodes), function(i) {
    if (!is.na(m_name[i]) && grepl("^[A-Za-z.]+$", m_name[i]) &&
        !is.na(width[i]) && nchar(m_name[i]) == width[i]) m_name[i]
    else .consensus_from_pwm(pwm[[i]])
  }, character(1))

  ## display label: MEME 5 puts "MEME-1" in alt and the consensus in name;
  ## MEME 4 puts the motif number in name and has no alt.
  label <- ifelse(!is.na(m_alt), m_alt,
                  ifelse(!is.na(m_name), paste0("MEME-", m_name), motif_id))

  motifs <- tibble::tibble(
    motif_id = motif_id,
    motif    = label,
    name     = m_name,
    alt      = m_alt,
    consensus = consensus,
    regex    = regex,
    width    = width,
    nsites   = .as_int(xml2::xml_attr(mnodes, "sites")),
    evalue   = .as_num(xml2::xml_attr(mnodes, "e_value")),
    pvalue   = .as_num(xml2::xml_attr(mnodes, "p_value")),
    ic       = .as_num(xml2::xml_attr(mnodes, "ic")),
    re       = .as_num(xml2::xml_attr(mnodes, "re")),
    llr      = .as_num(xml2::xml_attr(mnodes, "llr")),
    index    = seq_along(mnodes)
  )
  mmap <- stats::setNames(motifs$motif, motifs$motif_id)
  wmap <- stats::setNames(motifs$width, motifs$motif_id)

  ## ---- contributing sites ---------------------------------------------------
  ## NB: xml2 de-duplicates nodesets returned by xml_parent(), so the parent of
  ## each site is resolved by walking motif / scanned_sites nodes explicitly.
  contributing <- dplyr::bind_rows(lapply(seq_along(mnodes), function(i) {
    cs <- xml2::xml_find_all(mnodes[[i]], "./contributing_sites/contributing_site")
    if (!length(cs)) return(NULL)
    cs_seqid <- xml2::xml_attr(cs, "sequence_id")
    pos <- .as_int(xml2::xml_attr(cs, "position"))
    w <- rep(width[i], length(cs))
    site_seq <- left_f <- right_f <- rep(NA_character_, length(cs))
    if (isTRUE(sites_sequence)) {
      site_seq <- vapply(cs, function(s) {
        lr <- xml2::xml_find_all(s, "./site/letter_ref")
        if (!length(lr)) return(NA_character_)
        paste(unname(lmap[xml2::xml_attr(lr, "letter_id")]), collapse = "")
      }, character(1))
      gettx <- function(s, tag) {
        n <- xml2::xml_find_first(s, paste0("./", tag))
        if (inherits(n, "xml_missing")) NA_character_ else trimws(xml2::xml_text(n))
      }
      left_f  <- vapply(cs, gettx, character(1), tag = "left_flank")
      right_f <- vapply(cs, gettx, character(1), tag = "right_flank")
    }
    tibble::tibble(
      sequence_id = cs_seqid,
      sequence    = unname(smap[cs_seqid]),
      motif_id    = motif_id[i],
      motif       = label[i],
      start       = pos + 1L,
      stop        = pos + w,
      width       = w,
      strand      = .norm_strand(xml2::xml_attr(cs, "strand")),
      pvalue      = .as_num(xml2::xml_attr(cs, "pvalue")),
      site_seq    = site_seq,
      left_flank  = left_f,
      right_flank = right_f,
      source      = "contributing"
    )
  }))
  if (!nrow(contributing)) {
    contributing <- tibble::tibble(
      sequence_id = character(0), sequence = character(0), motif_id = character(0),
      motif = character(0), start = integer(0), stop = integer(0), width = integer(0),
      strand = character(0), pvalue = numeric(0), site_seq = character(0),
      left_flank = character(0), right_flank = character(0), source = character(0))
  }

  ## ---- scanned sites --------------------------------------------------------
  ssn <- xml2::xml_find_all(root, ".//scanned_sites_summary/scanned_sites")
  if (length(ssn)) {
    sites <- dplyr::bind_rows(lapply(seq_along(ssn), function(i) {
      sn <- xml2::xml_find_all(ssn[[i]], "./scanned_site")
      if (!length(sn)) return(NULL)
      sc_motif <- xml2::xml_attr(sn, "motif_id")
      ## a few MEME builds write the bare motif number instead of "motif_1"
      if (!all(sc_motif %in% motifs$motif_id)) {
        alt_key <- stats::setNames(motifs$motif_id, motifs$name)
        hit <- !sc_motif %in% motifs$motif_id & sc_motif %in% names(alt_key)
        sc_motif[hit] <- unname(alt_key[sc_motif[hit]])
      }
      pos <- .as_int(xml2::xml_attr(sn, "position"))
      w <- unname(wmap[sc_motif])
      sid <- xml2::xml_attr(ssn[[i]], "sequence_id")
      tibble::tibble(
        sequence_id = rep(sid, length(sn)),
        sequence    = unname(smap[rep(sid, length(sn))]),
        motif_id    = sc_motif,
        motif       = unname(mmap[sc_motif]),
        start       = pos + 1L,
        stop        = pos + w,
        width       = w,
        strand      = .norm_strand(xml2::xml_attr(sn, "strand")),
        pvalue      = .as_num(xml2::xml_attr(sn, "pvalue")),
        source      = "scanned"
      )
    }))
    seq_stats <- tibble::tibble(
      sequence_id = xml2::xml_attr(ssn, "sequence_id"),
      seq_pvalue  = .as_num(xml2::xml_attr(ssn, "pvalue")),
      n_sites     = .as_int(xml2::xml_attr(ssn, "num_sites"))
    )
  } else {
    sites <- contributing[, c("sequence_id", "sequence", "motif_id", "motif",
                              "start", "stop", "width", "strand", "pvalue", "source")]
    if (nrow(sites)) sites$source <- "contributing"
    seq_stats <- tibble::tibble(sequence_id = character(0), seq_pvalue = numeric(0),
                                n_sites = integer(0))
  }

  sequences <- dplyr::left_join(sequences, seq_stats, by = "sequence_id")
  if (!"n_sites" %in% names(sequences)) sequences$n_sites <- NA_integer_
  if (nrow(sites)) {
    cnt <- dplyr::count(sites, .data$sequence_id, name = ".n")
    sequences <- dplyr::left_join(sequences, cnt, by = "sequence_id")
    sequences$n_sites <- ifelse(is.na(sequences$n_sites), sequences$.n, sequences$n_sites)
    sequences$.n <- NULL
  }
  sequences$n_sites[is.na(sequences$n_sites)] <- 0L

  .new_meme_result(meta, motifs, sites, contributing, sequences, pwm, pssm,
                   background, list(type = atype, letters = lsym, colours = lcol))
}

# ---- STREME / DREME ---------------------------------------------------------

#' @keywords internal
#' @noRd
.read_streme_doc <- function(root) {
  tool <- toupper(xml2::xml_name(root))
  ver <- xml2::xml_attr(root, "version")
  lnodes <- xml2::xml_find_all(root, ".//alphabet/letter")
  lsym <- xml2::xml_attr(lnodes, "symbol")
  lsym[is.na(lsym)] <- sub("^letter_", "", xml2::xml_attr(lnodes, "id")[is.na(lsym)])
  core <- is.na(xml2::xml_attr(lnodes, "equals"))
  if (!any(core)) core <- rep(TRUE, length(lnodes))
  lsym <- lsym[core]
  atype <- .alphabet_type(xml2::xml_attr(xml2::xml_find_first(root, ".//alphabet"), "like"),
                          xml2::xml_attr(xml2::xml_find_first(root, ".//alphabet"), "name"), lsym)

  mnodes <- xml2::xml_find_all(root, ".//motifs/motif")
  pwm <- lapply(mnodes, function(m) {
    pos <- xml2::xml_find_all(m, "./pos")
    if (!length(pos)) return(NULL)
    mat <- vapply(lsym, function(l) .as_num(xml2::xml_attr(pos, l)), numeric(length(pos)))
    matrix(mat, nrow = length(pos), dimnames = list(NULL, lsym))
  })
  mid <- xml2::xml_attr(mnodes, "id")
  alt <- .blank_to_na(xml2::xml_attr(mnodes, "alt"))
  names(pwm) <- mid
  width <- .as_int(xml2::xml_attr(mnodes, "width"))
  width[is.na(width)] <- vapply(pwm, function(p) if (is.null(p)) NA_integer_ else nrow(p), integer(1))[is.na(width)]

  motifs <- tibble::tibble(
    motif_id = mid,
    motif    = ifelse(!is.na(alt), alt, mid),
    name     = .blank_to_na(xml2::xml_attr(mnodes, "name")),
    alt      = alt,
    consensus = vapply(pwm, .consensus_from_pwm, character(1)),
    regex    = NA_character_,
    width    = width,
    nsites   = .as_int(xml2::xml_attr(mnodes, "train_pos_count")),
    evalue   = .as_num(xml2::xml_attr(mnodes, "test_evalue")),
    pvalue   = .as_num(xml2::xml_attr(mnodes, "test_pvalue")),
    ic       = NA_real_, re = NA_real_, llr = NA_real_,
    index    = seq_along(mnodes)
  )
  rlang::warn(sprintf(
    "%s XML contains motifs but no per-sequence site table; `sites` is empty. Supply hits via read_fimo() or as_motif_table() to draw motif maps.", tool))

  empty_sites <- tibble::tibble(
    sequence_id = character(0), sequence = character(0), motif_id = character(0),
    motif = character(0), start = integer(0), stop = integer(0), width = integer(0),
    strand = character(0), pvalue = numeric(0), source = character(0))

  .new_meme_result(
    tibble::tibble(tool = tool, version = ver, alphabet = atype,
                   nmotifs = nrow(motifs), model_type = NA_character_,
                   nsequences = NA_integer_,
                   command_line = trimws(xml2::xml_text(xml2::xml_find_first(root, ".//command_line"))),
                   datafile = NA_character_, site_p_thresh = NA_real_),
    motifs, empty_sites, empty_sites,
    tibble::tibble(sequence_id = character(0), sequence = character(0),
                   length = integer(0), weight = numeric(0),
                   seq_pvalue = numeric(0), n_sites = integer(0)),
    pwm, stats::setNames(list(), character(0)),
    stats::setNames(numeric(0), character(0)),
    list(type = atype, letters = lsym, colours = stats::setNames(rep(NA_character_, length(lsym)), lsym))
  )
}

#' @keywords internal
#' @noRd
.alphabet_type <- function(like, name, letters) {
  x <- tolower(paste(c(like, name), collapse = " "))
  if (grepl("protein|amino", x)) return("protein")
  if (grepl("rna", x)) return("rna")
  if (grepl("dna|nucle", x)) return("dna")
  if (length(letters) > 8) return("protein")
  if ("U" %in% letters) return("rna")
  "dna"
}

#' @keywords internal
#' @noRd
.new_meme_result <- function(meta, motifs, sites, contributing, sequences,
                             pwm, pssm, background, alphabet) {
  structure(
    list(meta = meta, motifs = motifs, sites = sites, contributing = contributing,
         sequences = sequences, pwm = pwm, pssm = pssm, background = background,
         alphabet = alphabet),
    class = "meme_result"
  )
}
